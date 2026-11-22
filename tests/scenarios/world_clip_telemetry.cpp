// world/clip-telemetry (#1442): the clip telemetry sees a camera inside the ground and a camera above it. In the sandbox:
// 2 s with the eye 3 m under the plane (frames counted, a depth of ~3 m), then 2 s with it 5 m above (no new frames).
#include "scenario.h"
#include "cliptelemetry.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <cmath>

namespace {
struct State { long inside = 0; } S;
Vector3 Ground(ScenarioCtx& c) {
    const Vector3 up = c.spawnDirLocal();
    return Vector3Scale(up, c.field->radius + c.field->elevation(up));
}
Vector3 Side(Vector3 up) { return Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 })); }
}

SCENARIO("world/clip-telemetry") {
    s.tags = { "world" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.step("settle", 30.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T >= 1.0 && WorldIdle(); });
    s.step("the eye 3 m under the ground", 4.0, [](ScenarioCtx& c) {
        const Vector3 g = Ground(c), up = Vector3Normalize(g);
        const Vector3 eye = Vector3Subtract(g, Vector3Scale(up, 3.0f));
        c.camera(eye, Vector3Add(eye, Vector3Scale(Side(up), 10.0f)), up);
        if (c.t < 2.0) return false;
        const ClipStats& k = ClipTelemetry_Stats();
        c.gt("camera frames inside solid", (double)k.cameraFrames, 30.0);
        c.check("the camera is inside now", ClipTelemetry_CameraInside());
        c.near("its depth (m)", ClipTelemetry_CameraDepth(), 3.0, 0.75);
        S.inside = k.cameraFrames;
        return true;
    });
    s.step("the eye 5 m above it", 4.0, [](ScenarioCtx& c) {
        const Vector3 g = Ground(c), up = Vector3Normalize(g);
        const Vector3 eye = Vector3Add(g, Vector3Scale(up, 5.0f));
        c.camera(eye, Vector3Add(eye, Vector3Scale(Side(up), 10.0f)), up);
        if (c.t < 2.0) return false;
        c.check("the camera is out now", !ClipTelemetry_CameraInside());
        c.lt("camera frames inside solid since (one frame of lag at most)", (double)(ClipTelemetry_Stats().cameraFrames - S.inside), 2.0);
        return true;
    });
}
