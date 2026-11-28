// fluid/surface-probe: Fluid::SurfaceAbove (what swimming, pumps and buckets read) finds the water's real surface while the
// water moves. It used to stop at the first cell short of full going up from the point, and moving water - a pour landing
// in a pool, a column the pressure pass (#1168) is lifting - has such cells under the surface, so the probe read a false
// surface metres down. A 2 x 2 m pit 3 m deep is filled by a steady pour; every frame, SurfaceAbove from 0.3 m over the
// floor is compared with the top of the highest wet cell above that point (scanned down from over the rim). Counted: the
// frames where the probe reads more than 0.2 m under that top. (A guard for a still pool filling: with the pressure pass
// off this never showed the bug, before or after the fix. The test that fails without it is fluid/u-bend.)
#include "scenario.h"
#include "fluid.h"
#include "waterfield.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <string>

namespace {
constexpr float kDepth = 3.0f;
struct State {
    Vector3 home{}, ax{}, ay{};
    float r0 = 0.0f, poured = 0.0f;
    double nextPour = 0.0;
    int frames = 0, low = 0;
    float worst = 0.0f;
    Vector3 pt(float x, float y, float h) const {
        const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(home, r0), Vector3Add(Vector3Scale(ax, x), Vector3Scale(ay, y))));
        return Vector3Scale(d, r0 + h);
    }
    // The top of the highest wet cell over the pit's middle (m against the ground), scanning down from over the rim.
    float top(const Planetoid& p) const {
        const float cs = WaterField::CellSize();
        for (float h = 0.5f - 0.5f * cs; h > -kDepth; h -= cs) {
            const float w = Fluid::WaterAt(p, pt(0.0f, 0.0f, h), 0.0f) / WaterField::CellVolume();
            if (w > 0.02f) return h - 0.5f * cs + std::min(w, 1.0f) * cs;
        }
        return -kDepth;
    }
} S;
}

SCENARIO("fluid/surface-probe") {
    s.tags = { "fluid", "water" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 180.0;
    s.step("settle", 30.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("dig the pit", 1.0, [](ScenarioCtx& c) {
        S.home = c.field->findSpawnDir(c.body->seed);
        S.ax = Vector3Normalize(Vector3CrossProduct(S.home, { 0, 1, 0 })); S.ay = Vector3CrossProduct(S.ax, S.home);
        S.r0 = c.field->radius + c.field->elevation(S.home);
        TerrainEdit e{};
        e.center = S.pt(0.0f, 0.0f, -0.5f * kDepth + 0.5f);
        e.carve = true; e.box = true; e.half = { 1.0f, 1.0f, 0.5f * kDepth + 0.5f };
        e.ax = S.ax; e.ay = S.ay; e.az = Vector3Normalize(e.center);
        e.radius = Vector3Length(e.half);
        ApplyTerrainEdits(*c.body, &e, 1);
        c.camera(S.pt(0.0f, -6.0f, 3.0f), S.pt(0.0f, 0.0f, -1.5f), S.home);
        return true;
    });
    s.step("pour and probe every frame", 90.0, [](ScenarioCtx& c) {
        c.camera(S.pt(0.0f, -6.0f, 3.0f), S.pt(0.0f, 0.0f, -1.5f), S.home);
        c.waitingFor("the pit filled to 1 m under the rim");
        if (c.t <= 1.0) return false;
        if (c.t >= S.nextPour) { // a steady stream at one corner: water keeps moving under the surface
            S.nextPour = c.t + 0.1;
            S.poured += Fluid::Pour(*c.body, S.pt(0.6f, 0.6f, -0.25f), 0.08f);
        }
        const float top = S.top(*c.body);
        if (top > -kDepth + 0.6f) { // once there is a pool over the probe point
            const float s = Fluid::SurfaceAbove(*c.body, S.pt(0.0f, 0.0f, -kDepth + 0.3f), 6.0f); // (in the water: the floor's own cell is part ground)
            const float probe = s > 0.0f ? -kDepth + 0.3f + s : -kDepth;
            ++S.frames;
            const float under = top - probe;
            S.worst = std::max(S.worst, under);
            if (under > 0.2f) ++S.low;
        }
        if (top < -1.0f) return false;
        c.note(TextFormat("%d frames probed while it filled (%.1f m^3 poured): SurfaceAbove read more than 0.2 m under the water's top in %d, at worst %.2f m",
                          S.frames, S.poured, S.low, S.worst));
        c.metric("frames with a false surface", S.low, "frames");
        c.metric("worst false surface depth", S.worst, "m");
        c.gt("frames probed", S.frames, 50);
        c.lt("frames where SurfaceAbove read a false surface (> 0.2 m under the top)", S.low, 1);
        return true;
    });
}
