// audio/directions (#1724): positional sounds come from where they are, also while the planet turns under them. On the
// planet surface the player stands, level, facing local north (body-local tangent frame: north = toward the body's +Y pole
// projected onto the tangent plane; east = north x up, the player's right when facing north). Four emitters sit 10 m north,
// east, south and west at ear height, each set ONCE at the start of its case (as a creature call, a footstep or a machine
// does: a global position) and left alone for 2 s of sim time, while the ground turns about 100 m/s under it. Each plays
// alone; then the player turns east and it repeats. Per case it prints
//   (a) the direction Steam Audio is given (AudioSteamDirection: x right, y up, -z ahead), its angle off the true axis at
//       the start and after 2 s, the distance, and what a position kept as set (the old behaviour) would have given,
//   (b) the listener's stored forward and up against the true ones (the camera the player really looks through),
//   (c) the RMS of the left and right output of the voice, when Steam Audio runs a device (Audio_DebugVoiceLR), or else
//       of Steam Audio's own HRTF driven offline with white noise and that same direction (no device needed: Load(true)).
// Checks: the direction holds within 5 degrees and does not drift (3 degrees, 0.5 m over 2 s); a side source is louder in
// its ear by 3 dB, one ahead or behind is within 2 dB (when (c) is there). Needs RAYLIBGAME_AUDIO_V2=1 (the scenario sets
// it) for the live part; with no audio device (the headless build) the anchoring itself is tested against the true
// listener (Audio_AnchorOf / Audio_AnchorPos, which the emitters and one-shots use), plus the offline HRTF run.
#include "scenario.h"
#include "audio.h"
#include "audiospatial.h"
#include "player.h"
#include "steamaudio.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>

namespace {
struct Expect { char axis; int sign; }; // x: across the view (+ right), z: along it (-1 ahead, +1 behind)
constexpr const char* kName[4] = { "N", "E", "S", "W" };
constexpr Expect kExpect[2][4] = {
    { { 'z', -1 }, { 'x', +1 }, { 'z', +1 }, { 'x', -1 } },  // facing north: N ahead, E right, S behind, W left
    { { 'x', -1 }, { 'z', -1 }, { 'x', +1 }, { 'z', +1 } },  // facing east: N left, E ahead, S right, W behind
};
constexpr float kRange = 10.0f;     // m from the listener
constexpr double kSettle = 0.6;     // s of wall time before the window: the fade-in and the HRTF's start
constexpr double kWindow = 1.0;     // s of wall time the RMS covers

struct World4 {
    Vector3 up{}, facing{}, dirs[4]{}, ear{};
    DVec3 pos[4];
};

struct St {
    int k = 0;                      // the case, 0..7: phase = k / 4 (0 north, 1 east), emitter = k % 4
    double t0 = -1.0, tReset = -1.0, T0 = 0.0, err0 = 0.0, dist0 = 0.0;
    bool init = false, sampled = false;
    long frames = 0;
    int em[4] = { -1, -1, -1, -1 };
    bool live = false, lr = false, offline = false, loadedHere = false;
    SteamAudio::VoiceFx* fx = nullptr;
    DVec3 start[4];                 // where each emitter was set, once, at the case's start
    AudioAnchor anc[4];
} S;

World4 Compute(ScenarioCtx& c, int phase) {
    World4 w;
    const Vector3 feet = c.body->toLocalD(player.position).toFloat();
    w.up = Vector3Normalize(feet);
    const Vector3 pole = { 0, 1, 0 };
    const Vector3 north = Vector3Normalize(Vector3Subtract(pole, Vector3Scale(w.up, Vector3DotProduct(pole, w.up))));
    const Vector3 east = Vector3CrossProduct(north, w.up);
    w.dirs[0] = north; w.dirs[1] = east; w.dirs[2] = Vector3Negate(north); w.dirs[3] = Vector3Negate(east);
    w.facing = phase == 0 ? north : east;
    w.ear = Vector3Add(feet, Vector3Scale(w.up, player.eyeHeight));
    for (int i = 0; i < 4; ++i) w.pos[i] = c.body->toWorldD(DVec3(Vector3Add(w.ear, Vector3Scale(w.dirs[i], kRange))));
    return w;
}

// The listener as the player's pose is right now (no game state between): for the no-device run.
AudioListener TruthListener(ScenarioCtx& c, const World4& w) {
    AudioListener l;
    l.pos = c.body->toWorldD(DVec3(w.ear));
    l.forward = c.body->dirToWorld(w.facing);
    l.up = c.body->dirToWorld(w.up);
    return l;
}

// Steam Audio's HRTF on white noise from a direction (the library's own, no device): the RMS of each ear.
bool OfflineLR(const Vector3& dir, float& l, float& r) {
    SteamAudio::ResetVoice(S.fx);
    SteamAudio::Params p;
    p.dir[0] = dir.x; p.dir[1] = dir.y; p.dir[2] = dir.z;
    p.spatialBlend = 1.0f;
    uint32_t rng = 0x9e3779b9u;
    float in[SteamAudio::kBlock], L[SteamAudio::kBlock], R[SteamAudio::kBlock];
    double el = 0.0, er = 0.0;
    long n = 0;
    for (int blk = 0; blk < 60; ++blk) {
        for (float& x : in) { rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5; x = ((rng & 0xffff) / 32768.0f - 1.0f) * 0.3f; }
        SteamAudio::Process(S.fx, in, L, R, p);
        if (blk < 20) continue; // the HRTF's start
        for (int i = 0; i < SteamAudio::kBlock; ++i) { el += L[i] * L[i]; er += R[i] * R[i]; }
        n += SteamAudio::kBlock;
    }
    l = (float)sqrt(el / n); r = (float)sqrt(er / n);
    return n > 0;
}

void Setup(ScenarioCtx& c) {
    Audio_DebugLR(true);
    const int snd = Audio_Sound("test_noise");
    S.live = snd >= 0;
    if (S.live) {
        for (int i = 0; i < 4; ++i) {
            S.em[i] = Audio_EmitterCreate(snd, AudioBus::Effects, AudioPriority::Normal, {}, true);
            if (S.em[i] < 0) S.live = false;
            else Audio_EmitterSetPlaying(S.em[i], false);
        }
    }
    S.lr = S.live && SteamAudio::Active();
    if (!SteamAudio::Active()) { // no device: Steam Audio's HRTF on its own, when the library is there
        S.loadedHere = SteamAudio::Load(true);
        if (S.loadedHere) { S.fx = SteamAudio::CreateVoice(0); S.offline = S.fx != nullptr; }
    }
    c.note(std::string("emitters: ") + (S.live ? "live (audio device)" : "none (no audio device: the anchoring is tested against the true listener)") + "; left/right: "
           + (S.lr ? "the voice's output" : S.offline ? "Steam Audio's HRTF offline" : "none (no steam audio)")
           + "; backend: " + Audio_SpatialBackend() + "; steam audio: " + SteamAudio::Status());
}

void Finish() {
    for (int i = 0; i < 4; ++i) if (S.em[i] >= 0) Audio_EmitterDestroy(S.em[i]);
    if (S.fx) SteamAudio::DestroyVoice(S.fx);
    if (S.loadedHere) SteamAudio::Shutdown();
    Audio_DebugLR(false);
}
}

SCENARIO("audio/directions") {
    s.tags = { "audio", "planet" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 300.0;
    s.record = false;
    s.env = { { "RAYLIBGAME_AUDIO_V2", "1" } };

    s.step("settle", 150.0, [](ScenarioCtx& c) {
        const bool idle = WorldIdle();
        c.waitingFor(std::string("2 s of play, an idle world (") + (idle ? "yes" : "no") + ") and the player on the ground (" + (player.grounded ? "yes" : "no") + ")");
        static bool stood = false;
        if (idle && !stood) { stood = true; player.flyMode = false; player.spawnOn(c.world, *c.body, c.spawnDirWorld()); } // stands on the ground, as the player does (the camera is theirs)
        return stood && c.t >= 3.0 && player.grounded;
    });

    s.step("four emitters, facing north then east, the ground turning", 200.0, [](ScenarioCtx& c) {
        if (!S.init) { S.init = true; Setup(c); }
        const int phase = S.k / 4, i = S.k % 4;
        const World4 w = Compute(c, phase);
        player.forward = c.body->dirToWorld(w.facing); // level, as the player looks
        if (S.t0 < 0.0) { // a new case: the emitters are set ONCE (as a creature call or a machine does), then left to the planet's spin
            S.t0 = GetTime();
            S.T0 = c.T;
            S.tReset = -1.0;
            S.frames = 0;
            S.sampled = false;
            for (int j = 0; j < 4; ++j) {
                S.start[j] = w.pos[j];
                S.anc[j] = Audio_AnchorOf(w.pos[j]);
                if (S.em[j] >= 0) { Audio_EmitterSetPosition(S.em[j], w.pos[j]); Audio_EmitterSetPlaying(S.em[j], j == i); }
            }
        }
        ++S.frames;
        const Expect ex = kExpect[phase][i];
        // What the game hears now: the listener it stored and the emitter's resolved place (live), or the true listener now
        // and the anchor resolved now (no device: the anchoring itself).
        const AudioListener truth = TruthListener(c, w);
        const AudioListener& L = S.live ? Audio_GetListener() : truth;
        const DVec3 now = S.live ? Audio_DebugEmitterPos(S.em[i]) : Audio_AnchorPos(S.anc[i], S.start[i]);
        auto angleOf = [&](const DVec3& pos, double& dist) {
            const Vector3 d = AudioSteamDirection(L, pos);
            const Vector3 ideal = ex.axis == 'x' ? Vector3{ (float)ex.sign, 0, 0 } : Vector3{ 0, 0, (float)ex.sign };
            dist = (pos - L.pos).length();
            return (double)acosf(std::clamp(Vector3DotProduct(d, ideal), -1.0f, 1.0f)) * 57.29578;
        };
        double dist = 0.0;
        const double err = angleOf(now, dist);
        if (!S.sampled && S.frames >= 4) { S.sampled = true; S.err0 = err; S.dist0 = dist; }

        bool done = S.sampled && S.frames >= 8 && c.T - S.T0 >= 2.0; // two seconds of the planet turning
        const double el = GetTime() - S.t0;
        if (S.live && S.lr) {
            if (S.tReset < 0.0 && el >= kSettle) {
                float l, r;
                Audio_DebugVoiceLR(Audio_DebugEmitterVoice(S.em[i]), &l, &r, true); // (zero the sums)
                S.tReset = el;
            }
            done = done && S.tReset >= 0.0 && el >= S.tReset + kWindow;
        }
        if (!done) return false;

        // ---- the case is over: read everything ----
        const AudioListener& LG = Audio_GetListener();
        const Vector3 dir = AudioSteamDirection(L, now);
        const Vector3 fwdW = c.body->dirToWorld(w.facing), upW = c.body->dirToWorld(w.up);
        double staleDist = 0.0;
        const double staleErr = angleOf(S.start[i], staleDist); // the old behaviour: the position kept as it was set
        const std::string tag = std::string(phase == 0 ? "facing N, " : "facing E, ") + kName[i] + ": ";
        char buf[600];
        snprintf(buf, sizeof buf, "%sdir (right %.3f, up %.3f, back %.3f), %.2f m away, off the true axis by %.2f deg (at the start %.2f deg, %.2f m); the position kept as set would be %.2f deg off, %.1f m away; listener forward (%.3f %.3f %.3f) vs true (%.3f %.3f %.3f), up (%.3f %.3f %.3f) vs true (%.3f %.3f %.3f)",
                 tag.c_str(), dir.x, dir.y, dir.z, dist, err, S.err0, S.dist0, staleErr, staleDist, LG.forward.x, LG.forward.y, LG.forward.z, fwdW.x, fwdW.y, fwdW.z, LG.up.x, LG.up.y, LG.up.z, upW.x, upW.y, upW.z);
        c.note(buf);
        c.metric(tag + "the position kept as set: angle off after 2 s (deg)", staleErr, "deg");
        c.gt(tag + "listener forward matches the player's view (dot)", Vector3DotProduct(Vector3Normalize(LG.forward), fwdW), 0.99);
        c.gt(tag + "listener up matches the ground's up (dot)", Vector3DotProduct(Vector3Normalize(LG.up), upW), 0.99);
        c.between(tag + "distance to the emitter (m)", dist, kRange - 0.5, kRange + 0.5);
        c.lt(tag + "distance drift over 2 s (m)", fabs(dist - S.dist0), 0.5);
        c.lt(tag + "angle drift over 2 s (deg)", fabs(err - S.err0), 3.0);
        c.lt(tag + "angle off the true direction after 2 s (deg)", err, 5.0);
        c.gt(tag + (ex.axis == 'x' ? (ex.sign > 0 ? "direction to the right" : "direction to the left") : (ex.sign < 0 ? "direction ahead" : "direction behind")),
             ex.sign * (ex.axis == 'x' ? dir.x : dir.z), 0.7);

        float l = 0.0f, r = 0.0f;
        bool haveLR = false;
        if (S.live && S.lr) haveLR = Audio_DebugVoiceLR(Audio_DebugEmitterVoice(S.em[i]), &l, &r, true);
        else if (S.offline) haveLR = OfflineLR(dir, l, r);
        if (haveLR) {
            const double dB = 20.0 * log10(std::max(1e-9f, r) / std::max(1e-9f, l));
            snprintf(buf, sizeof buf, "%sleft %.5f, right %.5f (RMS), right - left %.2f dB%s", tag.c_str(), l, r, dB, S.live && S.lr ? " (the voice's output)" : " (offline HRTF)");
            c.note(buf);
            if (ex.axis == 'x') c.gt(tag + "louder in the near ear (dB)", ex.sign * dB, 3.0);
            else c.lt(tag + "balanced between the ears (|dB|)", fabs(dB), 2.0);
        } else c.note(tag + "no left/right (no steam audio)");

        S.t0 = -1.0;
        if (++S.k >= 8) { Finish(); return true; }
        return false;
    });
}
