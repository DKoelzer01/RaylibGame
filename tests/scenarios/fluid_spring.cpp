// fluid/spring (#870, moved from --test-spring, Plan 5): the nearest water spring to the start (Fluid::AddSource, its head
// over the ground) runs off as a stream: after a while most of its water is well away from it. PASS as the flag: a source
// exists, the moving water exceeds max(3, 5 x cell size) m^3, and under 60% of it is within 3 m of the spring.
#include "scenario.h"
#include "ambience.h"
#include "fluid.h"
#include "player.h"
#include "waterfield.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>

namespace {
struct State { Vector3 spring{}; } S;
void SpringCamera(ScenarioCtx& c) {
    const Vector3 up = Vector3Normalize(S.spring);
    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 }));
    c.camera(Vector3Add(Vector3Add(S.spring, Vector3Scale(t1, 18.0f)), Vector3Scale(up, 26.0f)), S.spring, up); // well above it
}
}

SCENARIO("fluid/spring") {
    s.tags = { "fluid", "water", "planet" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 150.0;
    s.legacyFlag = "--seed 42 --test-spring";

    s.step("find", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (!WorldIdle()) return false; // (the flag's branch starts on an idle frame)
        const Vector3 home = Vector3Scale(c.field->findSpawnDir(c.body->seed), c.field->radius);
        std::vector<Vector3> found;
        float best = 1e30f;
        for (float range : { 3000.0f, 10000.0f }) {
            Ambience_SpringsNear(*c.body, home, range, found);
            for (const Vector3& sp : found) { const float d = Vector3Distance(sp, home); if (d < best) { best = d; S.spring = sp; } }
            if (best < 1e29f) break;
        }
        if (best > 1e29f) { c.fail("no spring within 10 km"); return false; }
        const Vector3 up = Vector3Normalize(S.spring);
        c.note("nearest spring " + std::to_string(best / 1000.0f) + " km from the start (" + BiomeName(c.field->biomeAt(up, c.field->elevation(up))) + ")");
        player.spawnOn(c.world, *c.body, up);
        player.flyMode = true;
        return true;
    });
    s.step("arrive", 10.0, [](ScenarioCtx& c) { SpringCamera(c); c.waitingFor("5 s at the spring"); return c.t > 5.0; });
    s.step("run", 40.0, [](ScenarioCtx& c) {
        SpringCamera(c);
        if (c.t < 0.05) Fluid::SetSpeed(6); // fast-forward, as the flag
        c.waitingFor("30 s of six-fold fluid speed");
        if (c.t <= 30.0) return false;
        Fluid::SetSpeed(1);
        const float total = (float)Fluid::GetStats().volume, nearSrc = Fluid::FluidAt(*c.body, S.spring, 3.0f, 0);
        const float want = std::max(3.0f, 5.0f * WaterField::CellSize());
        c.metric("moving water", total, "m^3", false);
        c.metric("within 3 m of the spring", nearSrc, "m^3");
        c.gt("sources", Fluid::SourceCount(), 0.5);
        c.gt("moving water (m^3)", total, want);
        c.lt("share within 3 m of the spring", total > 0.0f ? nearSrc / total : 1.0, 0.6);
        return true;
    });
}
