#include "cliptelemetry.h"
#include "player.h"
#include "plots.h"
#include "profiler.h"
#include "recorder.h"
#include "scenario.h"
#include "scene.h"
#include "voxel.h"
#include "worldtypes.h"

#include <raylib.h>
#include <raymath.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

namespace {
ClipStats S;
bool camIn = false, plIn = false;     // this frame
float camDepth = 0.0f;
double camSince = 0.0, plSince = 0.0; // game seconds when the current stay inside began

int Mode() { // -1 unset, 0 off, 1 forced on
    static const int m = [] { const char* e = getenv("RAYLIBGAME_CLIP_TELEMETRY"); return e && *e ? (e[0] == '0' ? 0 : 1) : -1; }();
    return m;
}
bool Wanted() {
    const int m = Mode();
    if (m >= 0) return m == 1;
    return benchMode || Scenario_Active() || Recorder_Active() || Profiler_On();
}

// The terrain body for a body-local test: the one whose surface bound the point is nearest.
const Planetoid* BodyFor(const Scene& world, Vector3 renderPos) {
    const Planetoid* best = nullptr;
    float bestD = 1e30f;
    for (const auto& o : world.objects) {
        const Planetoid* p = dynamic_cast<const Planetoid*>(o.get());
        if (!p || !p->field) continue;
        const float d = Vector3Length(p->toLocal(renderPos)) - p->field->maxRadius();
        if (d < bestD) { bestD = d; best = p; }
    }
    return best && bestD < 50.0f ? best : nullptr; // above every surface bound: nothing to be inside
}

bool CubeAt(const Planetoid& p, Vector3 local) {
    if (!p.plots) return false;
    Plot* pl = p.plots->plotFor(local, 0.0f, false);
    if (!pl) return false;
    const int id = pl->partAt(pl->cellAt(pl->toPlot(local)));
    return id >= 0 && id < (int)pl->parts.size() && pl->parts[id].type == PartType::Cube;
}

bool Solid(const Planetoid& p, Vector3 local) {
    return p.field->density(local, p.edits.get()) >= ISO_LEVEL || CubeAt(p, local); // (density >= ISO_LEVEL is rock: > 0 put the surface 2 m high)
}

// Metres under the open surface: up the local vertical in 0.5 m steps to the first open sample (60 m at most).
float DepthAt(const Planetoid& p, Vector3 local) {
    const Vector3 up = Vector3Normalize(local);
    for (float d = 0.5f; d <= 60.0f; d += 0.5f)
        if (!Solid(p, Vector3Add(local, Vector3Scale(up, d)))) return d;
    return 60.0f;
}

void Track(bool in, bool& was, double& since, long& frames, double& seconds, double& first, double& last, float& depthMax,
           float depth, double t, double dt, const char* who) {
    if (in) {
        ++frames; seconds += dt; last = t;
        if (first < 0.0) first = t;
        depthMax = std::max(depthMax, depth);
        if (!was) { since = t; printf("bench: %s inside solid (d=%.1f m) at t=%.2f s\n", who, depth, t); }
    } else if (was) {
        printf("bench: %s out of solid at t=%.2f s (%.2f s inside)\n", who, t, t - since);
    }
    was = in;
}
}

void ClipTelemetry_Frame(Scene& world, double t, double dt) {
    if (!Wanted()) { camIn = plIn = false; camDepth = 0.0f; return; }
    // The camera: its eye and six points 0.15 m around it (the near plane reaches that far: a wall that close shows its back).
    bool cam = false; float cd = 0.0f;
    const Vector3 eyeR = world.camera.position;
    if (const Planetoid* p = BodyFor(world, eyeR)) {
        const Vector3 eye = p->toLocal(eyeR);
        static const Vector3 kRing[] = { { 0.15f, 0, 0 }, { -0.15f, 0, 0 }, { 0, 0.15f, 0 }, { 0, -0.15f, 0 }, { 0, 0, 0.15f }, { 0, 0, -0.15f } };
        cam = Solid(*p, eye);
        for (int k = 0; k < 6 && !cam; ++k) cam = Solid(*p, Vector3Add(eye, p->dirToLocal(kRing[k])));
        if (cam) cd = DepthAt(*p, eye);
    }
    Track(cam, camIn, camSince, S.cameraFrames, S.cameraSeconds, S.cameraFirst, S.cameraLast, S.cameraDepthMax, cd, t, dt, "camera");
    camDepth = cam ? cd : 0.0f;
    // The player's capsule (not while flying: then it is the camera): its sphere centres at 0.3, 0.85 and 1.4 m over the feet.
    bool pin = false; float pd = 0.0f;
    if (!player.flyMode) {
        const Vector3 feetR = (player.position - world.renderOrigin).toFloat();
        if (const Planetoid* p = BodyFor(world, feetR)) {
            const Vector3 feet = p->toLocal(feetR), up = Vector3Normalize(feet);
            for (float h : { 0.3f, 0.85f, 1.4f }) {
                const Vector3 c = Vector3Add(feet, Vector3Scale(up, h));
                if (Solid(*p, c)) { pin = true; pd = std::max(pd, DepthAt(*p, c)); }
            }
        }
    }
    Track(pin, plIn, plSince, S.playerFrames, S.playerSeconds, S.playerFirst, S.playerLast, S.playerDepthMax, pd, t, dt, "player");
}

bool ClipTelemetry_CameraInside() { return camIn; }
float ClipTelemetry_CameraDepth() { return camDepth; }
const ClipStats& ClipTelemetry_Stats() { return S; }
void ClipTelemetry_Reset() { S = ClipStats{}; }

void ClipTelemetry_DrawTag() {
    if (!camIn) return;
    const char* text = TextFormat("CLIPPED d=%.1f m", camDepth);
    const int fs = std::max(20, GetScreenHeight() / 40), w = MeasureText(text, fs);
    DrawRectangle(GetScreenWidth() - w - 24, 12, w + 16, fs + 10, Fade(MAROON, 0.85f));
    DrawText(text, GetScreenWidth() - w - 16, 17, fs, RAYWHITE);
}
