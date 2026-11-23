// water/river-bank-dig (#1520): digging near a river floods the hole even if there is no clear path for water to flow.
// The river nearest the start (seed 42, as water/river-trench finds it) gets a hole on its higher bank, 1.5 m beyond the
// water's edge, dug to 1.5 m under the river's level. Undug ground (>= ~0.9 m) lies between the hole and the water, so
// the hole must stay dry. After 20 s: is the hole's bottom wet, and how much water is within 1 m of it.
#include "scenario.h"
#include "fluid.h"
#include "rivers.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace {
struct BState {
    bool dug = false;
    double t0 = 0.0;
    Vector3 floorPt{}, up{};
    float ground = 0.0f, level = 0.0f, half = 0.0f;
    Vector3 hole{}, eye{}, target{};
};
BState B;
}

SCENARIO("water/river-bank-dig") {
    s.tags = { "water" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day" };
    s.timeout = 240.0;
    s.record = false;
    s.step("an idle world", 120.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T > 1.0 && WorldIdle(); });
    s.step("dig a hole 1.5 m off the river's edge", 5.0, [](ScenarioCtx& c) {
        B = BState{};
        const RiverNetwork* net = c.field->rivers();
        if (!net) { c.fail("no rivers"); return true; }
        const Vector3 home = c.spawnDirLocal();
        const Vector3 ax = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), ay = Vector3CrossProduct(ax, home);
        const float R = c.field->radius;
        Vector3 hit{}; bool found = false;
        for (float r = 0.0f; r < 3000.0f && !found; r += 10.0f)
            for (int k = 0; k < 64 && !found; ++k) {
                const float a = 2.0f * PI * k / 64.0f;
                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(ax, r * cosf(a) / R), Vector3Scale(ay, r * sinf(a) / R))));
                if (c.field->waterSourceAt(d).kind == 1 && net->nearest(d, 60.0f).depth >= 1.5f) { hit = d; found = true; }
            }
        if (!found) { c.fail("no river within 3 km"); return true; }
        const RiverNetwork::Hit rh = net->nearest(hit, 60.0f);
        const PlanetField::WaterSource ws = c.field->waterSourceAt(hit);
        const Vector3 across = Vector3Normalize(Vector3CrossProduct(hit, rh.down)); // cross(up, down)
        const float off = rh.halfWidth + 1.5f;
        auto at = [&](float sgn) { return Vector3Normalize(Vector3Add(hit, Vector3Scale(across, sgn * off / R))); };
        const Vector3 dP = at(1.0f), dN = at(-1.0f);
        const float eP = c.field->elevation(dP), eN = c.field->elevation(dN);
        const float sgn = eP >= eN ? 1.0f : -1.0f;
        const Vector3 dB = sgn > 0 ? dP : dN;
        const float ground = std::max(eP, eN);
        B.ground = ground; B.level = ws.level; B.half = rh.halfWidth; B.up = dB;
        // Stacked 0.6 m spheres from the ground down to 1.5 m under the river's level (centres 0.6 m over the floor at the lowest).
        const float floorE = ws.level - 1.5f;
        std::vector<TerrainEdit> cut;
        for (float e = ground; e >= floorE + 0.6f - 0.01f; e -= 0.6f)
            cut.push_back({ Vector3Scale(dB, R + e), 0.6f, -1.0f, true });
        cut.push_back({ Vector3Scale(dB, R + floorE + 0.6f), 0.6f, -1.0f, true });
        ApplyTerrainEdits(*c.body, cut.data(), cut.size());
        B.hole = Vector3Scale(dB, R + floorE + 0.3f);
        B.floorPt = B.hole;
        const Vector3 acrossB = Vector3Scale(across, sgn);
        // 6 m from the hole on the far side from the river, 3 m up, looking at the hole and the river.
        B.eye = Vector3Add(Vector3Scale(dB, R + ground + 3.0f), Vector3Scale(acrossB, 6.0f));
        B.target = Vector3Add(Vector3Scale(dB, R + ground - 0.5f), Vector3Scale(acrossB, -2.0f));
        c.note(TextFormat("river lat %.5f lon %.5f: level %.2f m, half width %.1f m; bank ground %.2f m (%s side, %.1f m over the level), hole floor %.2f m",
                          asinf(hit.y) * RAD2DEG, atan2f(hit.z, hit.x) * RAD2DEG, ws.level, rh.halfWidth, ground, sgn > 0 ? "+" : "-", ground - ws.level, floorE));
        c.camera(B.eye, B.target, dB);
        B.dug = true; B.t0 = c.T;
        return true;
    });
    s.step("20 s of water", 60.0, [](ScenarioCtx& c) {
        if (B.dug) c.camera(B.eye, B.target, B.up);
        c.waitingFor("20 s of water after the dig");
        if (c.T - B.t0 < 20.0) return false;
        const bool wet = Fluid::Wet(*c.body, B.floorPt);
        const float vol = Fluid::WaterAt(*c.body, B.hole, 1.0f);
        c.note(TextFormat("river level %.2f m, bank ground %.2f m, half width %.1f m; hole bottom wet %s, water within 1 m of the hole %.3f m^3",
                          B.level, B.ground, B.half, wet ? "YES" : "no", vol));
        c.metric("water within 1 m of the hole", vol, "m^3");
        c.check("dug the hole", B.dug);
        c.check("the hole by the river stays dry (no open path)", !wet);
        return true;
    });
}
