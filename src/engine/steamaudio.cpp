#include "steamaudio.h"

#include "audioscene.h"
#include "dynlib.h"
#include "utils/logger.h"

#include <phonon.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <thread>

// The functions used, looked up by name in the library (their declarations in phonon.h give the types).
#define SA_FUNCS(X) \
    X(iplContextCreate) X(iplContextRelease) X(iplHRTFCreate) X(iplHRTFRelease) \
    X(iplBinauralEffectCreate) X(iplBinauralEffectRelease) X(iplBinauralEffectReset) X(iplBinauralEffectApply) \
    X(iplDirectEffectCreate) X(iplDirectEffectRelease) X(iplDirectEffectReset) X(iplDirectEffectApply) \
    X(iplSceneCreate) X(iplSceneRelease) X(iplSceneCommit) \
    X(iplSimulatorCreate) X(iplSimulatorRelease) X(iplSimulatorSetScene) X(iplSimulatorCommit) \
    X(iplSimulatorSetSharedInputs) X(iplSimulatorRunDirect) \
    X(iplSourceCreate) X(iplSourceRelease) X(iplSourceAdd) X(iplSourceRemove) X(iplSourceSetInputs) X(iplSourceGetOutputs) \
    X(iplReflectionEffectCreate) X(iplReflectionEffectRelease) X(iplReflectionEffectApply)

namespace SteamAudio {
namespace {
    struct Fns {
#define SA_DECL(n) decltype(&::n) n = nullptr;
        SA_FUNCS(SA_DECL)
#undef SA_DECL
    } fn;
    void* lib = nullptr;
    IPLContext context = nullptr;
    IPLHRTF hrtf = nullptr;
    bool active = false;
    std::string status = "not loaded";

    void IPLCALL Log(IPLLogLevel level, const char* message) {
        if (level == IPL_LOGLEVEL_ERROR || level == IPL_LOGLEVEL_WARNING) LOG_WARN("Audio", "steam audio: %s\n", message);
    }
    IPLAudioSettings AudioSettings() { IPLAudioSettings s{}; s.samplingRate = kSampleRate; s.frameSize = kBlock; return s; }

    // ---- Occlusion simulation state (#1293) ----
    // Results per voice, written by the simulation thread, read by the audio thread.
    std::atomic<float> occTarget[kMaxVoices];
    std::atomic<float> transTarget[kMaxVoices][3];
    std::atomic<uint32_t> resultSeq[kMaxVoices];
    std::atomic<bool> simOn{ false };

    IPLScene simScene = nullptr;
    IPLSimulator simulator = nullptr;
    IPLSource sources[kMaxVoices] = {};
    std::thread simThread;
    std::mutex inMutex;
    std::condition_variable inCv;
    bool stopSim = false;
    struct Inputs {
        std::shared_ptr<const AudioScene> scene;
        DVec3 listener{};
        float ahead[3] = { 0, 0, -1 }, up[3] = { 0, 1, 0 };
        SimVoice voices[kMaxVoices];
    } inputs;
    std::mutex statsMutex;
    SimStats stats;

    // What the ray callbacks trace against during one RunDirect (set by the simulation thread before it).
    struct RayContext { std::shared_ptr<const AudioScene> scene; Vector3 listener{}; } rays;

    // The acoustic materials by AudioMat (their values: AudioMatPropsOf, audioscene.cpp).
    IPLMaterial kMaterials[(int)AudioMat::Count];
    void InitMaterials() {
        for (int m = 0; m < (int)AudioMat::Count; ++m) {
            const AudioMatProps& p = AudioMatPropsOf((AudioMat)m);
            for (int b = 0; b < 3; ++b) { kMaterials[m].absorption[b] = p.absorption[b]; kMaterials[m].transmission[b] = p.transmission[b]; }
            kMaterials[m].scattering = p.scattering;
        }
    }

    // ---- Reverb (#1294): one listener-centric parametric reverb on a send bus ----
    // How enclosed the listener is (AudioScene::enclosure, 64 rays about 3 times a second on the simulation thread)
    // sets its decay times and wet level; the audio thread runs the effect on the voices' summed sends.
    IPLReflectionEffect reverbFx = nullptr;
    std::atomic<float> reverbRt[3] = { 0.1f, 0.1f, 0.1f };
    std::atomic<float> reverbWet{ 0.0f };
    std::atomic<bool> reverbOn{ false };
    float reverbIn[kBlock], reverbL[kBlock], reverbR[kBlock];

    void IPLCALL ClosestHit(const IPLRay* ray, IPLfloat32 minD, IPLfloat32 maxD, IPLHit* hit, void*) {
        hit->distance = INFINITY;
        hit->triangleIndex = hit->objectIndex = hit->materialIndex = -1;
        if (maxD <= minD || !rays.scene) return;
        AudioMat m = AudioMat::Generic;
        const Vector3 o = { ray->origin.x, ray->origin.y, ray->origin.z }, d = { ray->direction.x, ray->direction.y, ray->direction.z };
        const float t = rays.scene->closestHit(o, d, minD, maxD, &m, &rays.listener);
        if (t < 0.0f) return;
        hit->distance = t;
        hit->normal = { -d.x, -d.y, -d.z }; // (direct simulation does not use it)
        hit->material = &kMaterials[(int)m];
        hit->materialIndex = (int)m;
    }
    void IPLCALL AnyHit(const IPLRay* ray, IPLfloat32 minD, IPLfloat32 maxD, IPLuint8* occluded, void*) {
        if (maxD <= minD || !rays.scene) { *occluded = maxD <= minD ? 1 : 0; return; }
        const Vector3 o = { ray->origin.x, ray->origin.y, ray->origin.z }, d = { ray->direction.x, ray->direction.y, ray->direction.z };
        *occluded = rays.scene->anyHit(o, d, minD, maxD, &rays.listener) ? 1 : 0;
    }

    void SimLoop() {
        bool added[kMaxVoices] = {};
        int periodMs = 100;
        long long iteration = 0;
        for (;;) {
            Inputs in;
            {
                std::unique_lock<std::mutex> lock(inMutex);
                inCv.wait_for(lock, std::chrono::milliseconds(periodMs), [] { return stopSim; });
                if (stopSim) break;
                in = inputs;
            }
            if (!in.scene) continue;
            const auto t0 = std::chrono::steady_clock::now();
            const DVec3 lrel = in.listener - in.scene->origin;
            rays.scene = in.scene;
            rays.listener = { (float)lrel.x, (float)lrel.y, (float)lrel.z };
            bool commit = false;
            int n = 0;
            for (int v = 0; v < kMaxVoices; ++v) {
                const bool want = in.voices[v].active;
                n += want ? 1 : 0;
                if (want == added[v]) continue;
                if (want) fn.iplSourceAdd(sources[v], simulator); else fn.iplSourceRemove(sources[v], simulator);
                added[v] = want;
                commit = true;
            }
            if (commit) fn.iplSimulatorCommit(simulator);
            // (#1294) The room around the listener, every third run (~3 times a second), whether or not anything plays.
            if (reverbFx && iteration++ % 3 == 0) {
                const AudioEnclosure enc = in.scene->enclosure(rays.listener, 64);
                float rt[3], wet = 0.0f;
                AudioReverbFromEnclosure(enc, rt, wet);
                for (int b = 0; b < 3; ++b) reverbRt[b].store(rt[b], std::memory_order_relaxed);
                reverbWet.store(wet, std::memory_order_relaxed);
                std::lock_guard<std::mutex> lock(statsMutex);
                stats.enclosure = enc.hitFraction;
                stats.meanDistance = enc.meanDistance;
                for (int b = 0; b < 3; ++b) stats.rt60[b] = rt[b];
                stats.wet = wet;
            }
            if (n == 0) { rays.scene.reset(); continue; }
            IPLSimulationSharedInputs sh{};
            const IPLVector3 ahead = { in.ahead[0], in.ahead[1], in.ahead[2] }, up = { in.up[0], in.up[1], in.up[2] };
            sh.listener.ahead = ahead;
            sh.listener.up = up;
            sh.listener.right = { ahead.y * up.z - ahead.z * up.y, ahead.z * up.x - ahead.x * up.z, ahead.x * up.y - ahead.y * up.x };
            sh.listener.origin = { rays.listener.x, rays.listener.y, rays.listener.z };
            sh.numRays = 1;
            sh.numBounces = 1;
            sh.duration = 0.1f;
            sh.order = 0;
            sh.irradianceMinDistance = 1.0f;
            fn.iplSimulatorSetSharedInputs(simulator, IPL_SIMULATIONFLAGS_DIRECT, &sh);
            for (int v = 0; v < kMaxVoices; ++v) {
                if (!added[v]) continue;
                const DVec3 r = in.voices[v].pos - in.scene->origin;
                IPLSimulationInputs si{};
                si.flags = IPL_SIMULATIONFLAGS_DIRECT;
                si.directFlags = (IPLDirectSimulationFlags)(IPL_DIRECTSIMULATIONFLAGS_OCCLUSION | IPL_DIRECTSIMULATIONFLAGS_TRANSMISSION);
                si.source.origin = { (float)r.x, (float)r.y, (float)r.z };
                si.source.ahead = { 0, 0, -1 };
                si.source.up = { 0, 1, 0 };
                si.source.right = { 1, 0, 0 };
                si.occlusionType = IPL_OCCLUSIONTYPE_RAYCAST;
                si.numTransmissionRays = 2;
                fn.iplSourceSetInputs(sources[v], IPL_SIMULATIONFLAGS_DIRECT, &si);
            }
            fn.iplSimulatorRunDirect(simulator);
            float occ[kMaxVoices];
            for (int v = 0; v < kMaxVoices; ++v) {
                occ[v] = 1.0f;
                if (!added[v]) continue;
                IPLSimulationOutputs out{};
                fn.iplSourceGetOutputs(sources[v], IPL_SIMULATIONFLAGS_DIRECT, &out);
                occTarget[v].store(std::clamp(out.direct.occlusion, 0.0f, 1.0f), std::memory_order_relaxed);
                for (int b = 0; b < 3; ++b) transTarget[v][b].store(std::clamp(out.direct.transmission[b], 0.0f, 1.0f), std::memory_order_relaxed);
                resultSeq[v].fetch_add(1, std::memory_order_release);
                occ[v] = out.direct.occlusion;
            }
            rays.scene.reset();
            const float ms = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
            periodMs = ms > 50.0f ? (int)(ms * 2.0f) : 100; // never more than about half a core
            std::lock_guard<std::mutex> lock(statsMutex);
            ++stats.ticks;
            stats.lastMs = ms;
            stats.maxMs = std::max(stats.maxMs, ms);
            stats.sources = n;
            std::copy(occ, occ + kMaxVoices, stats.occlusion);
        }
        for (int v = 0; v < kMaxVoices; ++v) if (added[v]) fn.iplSourceRemove(sources[v], simulator);
        fn.iplSimulatorCommit(simulator);
    }
}

struct VoiceFx {
    int index = 0;
    IPLBinauralEffect binaural = nullptr;
    IPLDirectEffect direct = nullptr;
    float in[kBlock], mid[kBlock], left[kBlock], right[kBlock];
    // (#1293) the occlusion the audio thread applies, gliding to the simulation's
    bool fresh = true;              // no result yet for this sound: the first one is taken at once
    uint32_t seenSeq = 0;
    float occ = 1.0f, trans[3] = { 1, 1, 1 };
};

bool Load() {
    if (active) return true;
#if RAYLIBGAME_HEADLESS_PLATFORM
    status = "off (server build)";
    if ((true)) return false;
#endif
    if (const char* e = getenv("RAYLIBGAME_STEAM_AUDIO"); e && e[0] == '0') { status = "off (RAYLIBGAME_STEAM_AUDIO=0)"; return false; }
#if defined(_WIN32)
    lib = DynLib_Open("phonon.dll");
    if (!lib) { status = "not found (phonon.dll: tools/fetch_steamaudio.sh)"; return false; }
#else
    lib = DynLib_Open("libphonon.so");
    if (!lib) { status = "not found (libphonon.so: tools/fetch_steamaudio.sh)"; return false; }
#endif
    std::string missing;
#define SA_LOAD(n) fn.n = reinterpret_cast<decltype(fn.n)>(DynLib_Symbol(lib, #n)); if (!fn.n) missing += " " #n;
    SA_FUNCS(SA_LOAD)
#undef SA_LOAD
    if (!missing.empty()) { status = "library lacks:" + missing; LOG_WARN("Audio", "steam audio: %s\n", status.c_str()); DynLib_Close(lib); lib = nullptr; return false; }

    IPLContextSettings cs{};
    cs.version = STEAMAUDIO_VERSION;
    cs.logCallback = Log;
    cs.simdLevel = IPL_SIMDLEVEL_AVX2; // AVX-512 can throttle the clock (the SDK's own advice)
    if (fn.iplContextCreate(&cs, &context) != IPL_STATUS_SUCCESS) { status = "context failed (version mismatch?)"; DynLib_Close(lib); lib = nullptr; return false; }
    IPLAudioSettings as = AudioSettings();
    IPLHRTFSettings hs{};
    hs.type = IPL_HRTFTYPE_DEFAULT;
    hs.volume = 1.0f;
    hs.normType = IPL_HRTFNORMTYPE_NONE;
    if (fn.iplHRTFCreate(context, &as, &hs, &hrtf) != IPL_STATUS_SUCCESS) {
        status = "HRTF failed";
        fn.iplContextRelease(&context);
        DynLib_Close(lib);
        lib = nullptr;
        return false;
    }
    char buf[96];
    snprintf(buf, sizeof buf, "%d.%d.%d, HRTF", STEAMAUDIO_VERSION_MAJOR, STEAMAUDIO_VERSION_MINOR, STEAMAUDIO_VERSION_PATCH);
    status = buf;
    active = true;
    LOG_INFO("Audio", "steam audio %s loaded\n", buf);
    return true;
}

void Shutdown() {
    if (!active) return;
    StopSimulation();
    fn.iplHRTFRelease(&hrtf);
    fn.iplContextRelease(&context);
    DynLib_Close(lib);
    lib = nullptr;
    active = false;
    status = "shut down";
}

bool Active() { return active; }
const char* Status() { return status.c_str(); }

VoiceFx* CreateVoice(int index) {
    if (!active || index < 0 || index >= kMaxVoices) return nullptr;
    VoiceFx* v = new VoiceFx();
    v->index = index;
    IPLAudioSettings as = AudioSettings();
    IPLBinauralEffectSettings bs{};
    bs.hrtf = hrtf;
    IPLDirectEffectSettings ds{};
    ds.numChannels = 1;
    if (fn.iplBinauralEffectCreate(context, &as, &bs, &v->binaural) != IPL_STATUS_SUCCESS ||
        fn.iplDirectEffectCreate(context, &as, &ds, &v->direct) != IPL_STATUS_SUCCESS) {
        DestroyVoice(v);
        return nullptr;
    }
    return v;
}

void DestroyVoice(VoiceFx* v) {
    if (!v) return;
    if (active) {
        if (v->binaural) fn.iplBinauralEffectRelease(&v->binaural);
        if (v->direct) fn.iplDirectEffectRelease(&v->direct);
    }
    delete v;
}

void ResetVoice(VoiceFx* v) {
    if (!v || !active) return;
    fn.iplBinauralEffectReset(v->binaural);
    fn.iplDirectEffectReset(v->direct);
    v->fresh = true;
    v->seenSeq = resultSeq[v->index].load(std::memory_order_acquire);
    v->occ = 1.0f;
    for (float& t : v->trans) t = 1.0f;
}

void Process(VoiceFx* v, const float* mono, float* left, float* right, const Params& p) {
    memcpy(v->in, mono, sizeof v->in);
    float* inCh[1] = { v->in };
    float* midCh[1] = { v->mid };
    float* outCh[2] = { v->left, v->right };
    IPLAudioBuffer in{ 1, kBlock, inCh }, mid{ 1, kBlock, midCh }, out{ 2, kBlock, outCh };

    IPLDirectEffectParams dp{};
    dp.flags = IPL_DIRECTEFFECTFLAGS_APPLYAIRABSORPTION;
    dp.transmissionType = IPL_TRANSMISSIONTYPE_FREQINDEPENDENT;
    dp.distanceAttenuation = 1.0f; // (ours: the voice's volume)
    for (int b = 0; b < 3; ++b) { dp.airAbsorption[b] = p.air[b]; dp.transmission[b] = 1.0f; }
    dp.directivity = 1.0f;
    dp.occlusion = 1.0f;
    if (simOn.load(std::memory_order_relaxed)) { // (#1293)
        const int i = v->index;
        const uint32_t seq = resultSeq[i].load(std::memory_order_acquire);
        const bool haveNew = seq != v->seenSeq;
        if (!v->fresh || haveNew) {
            const float occT = occTarget[i].load(std::memory_order_relaxed);
            float trT[3];
            for (int b = 0; b < 3; ++b) trT[b] = transTarget[i][b].load(std::memory_order_relaxed);
            const float k = v->fresh ? 1.0f : 0.2f; // ~25 ms glide at 256-frame blocks
            v->occ += (occT - v->occ) * k;
            for (int b = 0; b < 3; ++b) v->trans[b] += (trT[b] - v->trans[b]) * k;
            v->fresh = false;
            v->seenSeq = seq;
        }
        dp.flags = (IPLDirectEffectFlags)(dp.flags | IPL_DIRECTEFFECTFLAGS_APPLYOCCLUSION | IPL_DIRECTEFFECTFLAGS_APPLYTRANSMISSION);
        dp.transmissionType = IPL_TRANSMISSIONTYPE_FREQDEPENDENT;
        dp.occlusion = v->occ;
        for (int b = 0; b < 3; ++b) dp.transmission[b] = v->trans[b];
    }
    fn.iplDirectEffectApply(v->direct, &dp, &in, &mid);

    IPLBinauralEffectParams bp{};
    bp.direction = { p.dir[0], p.dir[1], p.dir[2] };
    bp.interpolation = IPL_HRTFINTERPOLATION_BILINEAR;
    bp.spatialBlend = p.spatialBlend;
    bp.hrtf = hrtf;
    bp.peakDelays = nullptr;
    fn.iplBinauralEffectApply(v->binaural, &bp, &mid, &out);
    memcpy(left, v->left, sizeof v->left);
    memcpy(right, v->right, sizeof v->right);
}

// ---- Occlusion and transmission (#1293) ----

bool StartSimulation() {
    if (!active) return false;
    if (simulator) return true;
    if (const char* e = getenv("RAYLIBGAME_AUDIO_OCCLUSION"); e && e[0] == '0') return false;
    InitMaterials();
    IPLSceneSettings ss{};
    ss.type = IPL_SCENETYPE_CUSTOM;
    ss.closestHitCallback = ClosestHit;
    ss.anyHitCallback = AnyHit;
    if (fn.iplSceneCreate(context, &ss, &simScene) != IPL_STATUS_SUCCESS) { LOG_WARN("Audio", "steam audio: scene failed: no occlusion\n"); return false; }
    fn.iplSceneCommit(simScene);
    IPLSimulationSettings sim{};
    sim.flags = IPL_SIMULATIONFLAGS_DIRECT;
    sim.sceneType = IPL_SCENETYPE_CUSTOM;
    sim.reflectionType = IPL_REFLECTIONEFFECTTYPE_CONVOLUTION;
    sim.maxNumOcclusionSamples = 16;
    sim.maxNumRays = 1;
    sim.numDiffuseSamples = 1;
    sim.maxDuration = 0.1f;
    sim.maxOrder = 0;
    sim.maxNumSources = kMaxVoices;
    sim.numThreads = 1;
    sim.rayBatchSize = 1;
    sim.samplingRate = kSampleRate;
    sim.frameSize = kBlock;
    if (fn.iplSimulatorCreate(context, &sim, &simulator) != IPL_STATUS_SUCCESS) {
        fn.iplSceneRelease(&simScene);
        LOG_WARN("Audio", "steam audio: simulator failed: no occlusion\n");
        return false;
    }
    fn.iplSimulatorSetScene(simulator, simScene);
    fn.iplSimulatorCommit(simulator);
    IPLSourceSettings srcs{};
    srcs.flags = IPL_SIMULATIONFLAGS_DIRECT;
    for (int v = 0; v < kMaxVoices; ++v) {
        if (fn.iplSourceCreate(simulator, &srcs, &sources[v]) == IPL_STATUS_SUCCESS) continue;
        for (int k = 0; k < v; ++k) fn.iplSourceRelease(&sources[k]);
        fn.iplSimulatorRelease(&simulator);
        fn.iplSceneRelease(&simScene);
        LOG_WARN("Audio", "steam audio: sources failed: no occlusion\n");
        return false;
    }
    for (int v = 0; v < kMaxVoices; ++v) {
        occTarget[v].store(1.0f);
        for (int b = 0; b < 3; ++b) transTarget[v][b].store(1.0f);
    }
    // (#1294) The reverb: parametric (a feedback delay network driven by three decay times), cheap at any size.
    if (const char* e = getenv("RAYLIBGAME_AUDIO_REVERB"); !(e && e[0] == '0')) {
        IPLAudioSettings as = AudioSettings();
        IPLReflectionEffectSettings rs{};
        rs.type = IPL_REFLECTIONEFFECTTYPE_PARAMETRIC;
        rs.irSize = kSampleRate; // (unused by the parametric type)
        rs.numChannels = 1;
        if (fn.iplReflectionEffectCreate(context, &as, &rs, &reverbFx) == IPL_STATUS_SUCCESS) reverbOn.store(true);
        else { reverbFx = nullptr; LOG_WARN("Audio", "steam audio: reverb effect failed: dry\n"); }
    }
    stopSim = false;
    simThread = std::thread(SimLoop);
    simOn.store(true);
    status += reverbFx ? ", occlusion, reverb" : ", occlusion";
    LOG_INFO("Audio", "steam audio: occlusion and transmission on (custom scene: the density field and part cells)%s\n", reverbFx ? ", reverb" : "");
    return true;
}

bool ReverbOn() { return reverbOn.load(std::memory_order_relaxed); }
float ReverbWet() { return reverbOn.load(std::memory_order_relaxed) ? reverbWet.load(std::memory_order_relaxed) : 0.0f; }

void ReverbProcess(const float* mono, float* left, float* right) {
    if (!reverbFx || !reverbOn.load(std::memory_order_relaxed)) {
        memset(left, 0, sizeof(float) * kBlock);
        memset(right, 0, sizeof(float) * kBlock);
        return;
    }
    memcpy(reverbIn, mono, sizeof reverbIn);
    float* inCh[1] = { reverbIn };
    float* outCh[2] = { reverbL, reverbR };
    IPLAudioBuffer in{ 1, kBlock, inCh }, out{ 2, kBlock, outCh };
    IPLReflectionEffectParams rp{};
    rp.type = IPL_REFLECTIONEFFECTTYPE_PARAMETRIC;
    for (int b = 0; b < 3; ++b) { rp.reverbTimes[b] = reverbRt[b].load(std::memory_order_relaxed); rp.eq[b] = 1.0f; }
    rp.numChannels = 1;
    rp.irSize = kSampleRate;
    fn.iplReflectionEffectApply(reverbFx, &rp, &in, &out, nullptr);
    // (If the parametric reverb fills only its first channel, both ears get it.)
    float eR = 0.0f;
    for (int i = 0; i < kBlock; ++i) eR += reverbR[i] * reverbR[i];
    memcpy(left, reverbL, sizeof reverbL);
    memcpy(right, eR > 0.0f ? reverbR : reverbL, sizeof reverbR);
}

void StopSimulation() {
    if (!simulator) return;
    simOn.store(false);
    { std::lock_guard<std::mutex> lock(inMutex); stopSim = true; inputs = Inputs{}; }
    inCv.notify_all();
    if (simThread.joinable()) simThread.join();
    for (IPLSource& s : sources) if (s) fn.iplSourceRelease(&s);
    reverbOn.store(false); // (the reverb stream is gone by now: audio.cpp unloads its streams first)
    if (reverbFx) fn.iplReflectionEffectRelease(&reverbFx);
    reverbFx = nullptr;
    fn.iplSimulatorRelease(&simulator);
    fn.iplSceneRelease(&simScene);
    simulator = nullptr;
    simScene = nullptr;
}

bool Simulating() { return simOn.load(std::memory_order_relaxed); }

void SetSimInputs(std::shared_ptr<const AudioScene> scene, const DVec3& listener, const float ahead[3], const float up[3],
                  const SimVoice (&voices)[kMaxVoices]) {
    if (!simulator) return;
    std::lock_guard<std::mutex> lock(inMutex);
    inputs.scene = std::move(scene);
    inputs.listener = listener;
    for (int k = 0; k < 3; ++k) { inputs.ahead[k] = ahead[k]; inputs.up[k] = up[k]; }
    for (int v = 0; v < kMaxVoices; ++v) inputs.voices[v] = voices[v];
}

SimStats GetSimStats() {
    std::lock_guard<std::mutex> lock(statsMutex);
    return stats;
}

} // namespace SteamAudio
