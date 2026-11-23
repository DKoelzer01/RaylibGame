// audio/directions (#1724): positional sounds come from where they are. On the planet surface the player stands, level,
// facing local north (body-local tangent frame: north = toward the body's +Y pole projected onto the tangent plane; east =
// north x up, the player's right when facing north). Four emitters sit 10 m north, east, south and west at ear height. Each
// plays alone for about a second; then the player turns east and it repeats. Per case it prints
//   (a) the direction Steam Audio is given (AudioSteamDirection of the listener the game stored: x right, y up, -z ahead),
//   (b) the listener's stored forward and up against the true ones (the camera the player really looks through),
//   (c) the RMS of the left and right output of the voice, when Steam Audio runs a device (Audio_DebugVoiceLR), or else
//       of Steam Audio's own HRTF driven offline with white noise and that same direction (no device needed: Load(true)).
// Checks: the direction (always), left against right (when (c) is there): a side source is louder in its ear by 3 dB, one
// ahead or behind is within 2 dB. Needs RAYLIBGAME_AUDIO_V2=1 (the scenario sets it) for the live part; without an audio
// device (the headless build) only (a), (b) and the offline HRTF run.
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
    double t0 = -1.0, tReset = -1.0;
    long frames = 0;
    int em[4] = { -1, -1, -1, -1 };
    bool live = false, lr = false, offline = false, loadedHere = false;
    SteamAudio::VoiceFx* fx = nullptr;
    DVec3 prev[4];                  // the positions of the last frame: the listener the game stored is from its end
    bool havePrev = false;
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
    c.note(std::string("emitters: ") + (S.live ? "live (audio device)" : "none (no audio device)") + "; left/right: "
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
    s.timeout = 120.0;
    s.record = false;
    s.env = { { "RAYLIBGAME_AUDIO_V2", "1" } };

    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play, an idle world and the player on the ground"); return c.T >= 2.0 && WorldIdle() && player.grounded; });

    s.step("four emitters, facing north then east", 100.0, [](ScenarioCtx& c) {
        if (S.k == 0 && S.t0 < 0.0 && !S.havePrev) Setup(c);
        const int phase = S.k / 4, i = S.k % 4;
        const World4 w = Compute(c, phase);
        player.forward = c.body->dirToWorld(w.facing); // level, as the player looks
        if (S.t0 < 0.0) {
            S.t0 = GetTime();
            S.tReset = -1.0;
            S.frames = 0;
            for (int j = 0; j < 4; ++j) if (S.em[j] >= 0) Audio_EmitterSetPlaying(S.em[j], j == i);
        }
        ++S.frames;
        for (int j = 0; j < 4; ++j) if (S.em[j] >= 0) Audio_EmitterSetPosition(S.em[j], w.pos[j]);
        const double el = GetTime() - S.t0;
        bool done = false;
        if (S.live && S.lr) {
            if (S.tReset < 0.0 && el >= kSettle) {
                float l, r;
                Audio_DebugVoiceLR(Audio_DebugEmitterVoice(S.em[i]), &l, &r, true); // (zero the sums)
                S.tReset = el;
            }
            done = S.tReset >= 0.0 && el >= S.tReset + kWindow;
        } else {
            done = S.frames >= 8 && S.havePrev; // the listener the game stored is from the end of the last frame
        }
        const DVec3 posNow[4] = { w.pos[0], w.pos[1], w.pos[2], w.pos[3] };
        if (!done || !S.havePrev) { for (int j = 0; j < 4; ++j) S.prev[j] = posNow[j]; S.havePrev = true; return false; }

        // ---- the case is over: read everything ----
        const AudioListener& L = Audio_GetListener();
        const Vector3 dir = AudioSteamDirection(L, S.prev[i]);
        const Vector3 fwdW = c.body->dirToWorld(w.facing), upW = c.body->dirToWorld(w.up);
        const double dist = (S.prev[i] - L.pos).length();
        const std::string tag = std::string(phase == 0 ? "facing N, " : "facing E, ") + kName[i] + ": ";
        char buf[400];
        snprintf(buf, sizeof buf, "%sdir (right %.3f, up %.3f, back %.3f), %.2f m away; listener forward (%.3f %.3f %.3f) vs true (%.3f %.3f %.3f), up (%.3f %.3f %.3f) vs true (%.3f %.3f %.3f)",
                 tag.c_str(), dir.x, dir.y, dir.z, dist, L.forward.x, L.forward.y, L.forward.z, fwdW.x, fwdW.y, fwdW.z, L.up.x, L.up.y, L.up.z, upW.x, upW.y, upW.z);
        c.note(buf);
        c.gt(tag + "listener forward matches the player's view (dot)", Vector3DotProduct(Vector3Normalize(L.forward), fwdW), 0.99);
        c.gt(tag + "listener up matches the ground's up (dot)", Vector3DotProduct(Vector3Normalize(L.up), upW), 0.99);
        c.between(tag + "distance to the emitter (m)", dist, kRange - 0.5, kRange + 0.5);
        const Expect ex = kExpect[phase][i];
        const float along = ex.axis == 'x' ? dir.x : dir.z;
        c.gt(tag + (ex.axis == 'x' ? (ex.sign > 0 ? "direction to the right" : "direction to the left") : (ex.sign < 0 ? "direction ahead" : "direction behind")), ex.sign * along, 0.7);

        float l = 0.0f, r = 0.0f;
        bool haveLR = false;
        if (S.lr) haveLR = Audio_DebugVoiceLR(Audio_DebugEmitterVoice(S.em[i]), &l, &r, true);
        else if (S.offline) haveLR = OfflineLR(dir, l, r);
        if (haveLR) {
            const double dB = 20.0 * log10(std::max(1e-9f, r) / std::max(1e-9f, l));
            snprintf(buf, sizeof buf, "%sleft %.5f, right %.5f (RMS), right - left %.2f dB%s", tag.c_str(), l, r, dB, S.lr ? " (the voice's output)" : " (offline HRTF)");
            c.note(buf);
            if (ex.axis == 'x') c.gt(tag + "louder in the near ear (dB)", ex.sign * dB, 3.0);
            else c.lt(tag + "balanced between the ears (|dB|)", fabs(dB), 2.0);
        } else if (S.live) c.note(tag + "no left/right (no steam audio: own panning)");
        else c.note(tag + "no left/right (no steam audio)");

        S.t0 = -1.0;
        S.havePrev = false;
        if (++S.k >= 8) { Finish(); return true; }
        return false;
    });
}
