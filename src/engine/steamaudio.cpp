#include "steamaudio.h"

#include "utils/logger.h"

#include <phonon.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif

// The functions used, looked up by name in the library (their declarations in phonon.h give the types).
#define SA_FUNCS(X) \
    X(iplContextCreate) X(iplContextRelease) X(iplHRTFCreate) X(iplHRTFRelease) \
    X(iplBinauralEffectCreate) X(iplBinauralEffectRelease) X(iplBinauralEffectReset) X(iplBinauralEffectApply) \
    X(iplDirectEffectCreate) X(iplDirectEffectRelease) X(iplDirectEffectReset) X(iplDirectEffectApply)

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

    void* OpenLib() {
#if defined(_WIN32)
        return (void*)LoadLibraryA("phonon.dll");
#else
        return dlopen("libphonon.so", RTLD_NOW | RTLD_LOCAL);
#endif
    }
    void* Sym(const char* name) {
#if defined(_WIN32)
        return reinterpret_cast<void*>(GetProcAddress((HMODULE)lib, name));
#else
        return dlsym(lib, name);
#endif
    }
    void CloseLib() {
        if (!lib) return;
#if defined(_WIN32)
        FreeLibrary((HMODULE)lib);
#else
        dlclose(lib);
#endif
        lib = nullptr;
    }
    IPLAudioSettings AudioSettings() { IPLAudioSettings s{}; s.samplingRate = kSampleRate; s.frameSize = kBlock; return s; }
}

struct VoiceFx {
    IPLBinauralEffect binaural = nullptr;
    IPLDirectEffect direct = nullptr;
    float in[kBlock], mid[kBlock], left[kBlock], right[kBlock];
};

bool Load() {
    if (active) return true;
#if RAYLIBGAME_HEADLESS_PLATFORM
    status = "off (server build)";
    if ((true)) return false;
#endif
    if (const char* e = getenv("RAYLIBGAME_STEAM_AUDIO"); e && e[0] == '0') { status = "off (RAYLIBGAME_STEAM_AUDIO=0)"; return false; }
    lib = OpenLib();
#if defined(_WIN32)
    if (!lib) { status = "not found (phonon.dll: tools/fetch_steamaudio.sh)"; return false; }
#else
    if (!lib) { status = "not found (libphonon.so: tools/fetch_steamaudio.sh)"; return false; }
#endif
    std::string missing;
#define SA_LOAD(n) fn.n = reinterpret_cast<decltype(fn.n)>(Sym(#n)); if (!fn.n) missing += " " #n;
    SA_FUNCS(SA_LOAD)
#undef SA_LOAD
    if (!missing.empty()) { status = "library lacks:" + missing; LOG_WARN("Audio", "steam audio: %s\n", status.c_str()); CloseLib(); return false; }

    IPLContextSettings cs{};
    cs.version = STEAMAUDIO_VERSION;
    cs.logCallback = Log;
    cs.simdLevel = IPL_SIMDLEVEL_AVX2; // AVX-512 can throttle the clock (the SDK's own advice)
    if (fn.iplContextCreate(&cs, &context) != IPL_STATUS_SUCCESS) { status = "context failed (version mismatch?)"; CloseLib(); return false; }
    IPLAudioSettings as = AudioSettings();
    IPLHRTFSettings hs{};
    hs.type = IPL_HRTFTYPE_DEFAULT;
    hs.volume = 1.0f;
    hs.normType = IPL_HRTFNORMTYPE_NONE;
    if (fn.iplHRTFCreate(context, &as, &hs, &hrtf) != IPL_STATUS_SUCCESS) {
        status = "HRTF failed";
        fn.iplContextRelease(&context);
        CloseLib();
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
    fn.iplHRTFRelease(&hrtf);
    fn.iplContextRelease(&context);
    CloseLib();
    active = false;
    status = "shut down";
}

bool Active() { return active; }
const char* Status() { return status.c_str(); }

VoiceFx* CreateVoice() {
    if (!active) return nullptr;
    VoiceFx* v = new VoiceFx();
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

} // namespace SteamAudio
