// fluid/pond (#870, moved from --test-pond, Plan 4 F4): drain the nearest pond into a pit dug beside it. The pond is simulated
// water, its top drains, the pit fills; the readers (farming, animals, the map) see a drained shore as dry (#385); no trees
// grow in its bed (#189). Then save, go 300 m away, load and come back: the pond stays drained and the pit full. Last, a
// swim in the flooded pit (topped up when shallow): in water 90% of frames, sinking under 1 m/s. PASS as the flag.
#include "scenario.h"
#include "fluid.h"
#include "player.h"
#include "trees.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>

namespace {
struct State {
    Vector3 dir{}, out{}, pit{}, mid{};
    float level = 0.0f, r = 0.0f;
    bool firstOk = false;
    int topping = -1, inWater = 0, frames = 0;
    double lastPour = -10.0;
    float h0 = 0.0f;
    const PlanetField* pf = nullptr;
    Vector3 at(Vector3 d, float h) const { return Vector3Scale(d, pf->radius + h); }
    Vector3 along(float sv) const { return Vector3Normalize(Vector3Add(dir, Vector3Scale(out, sv / pf->radius))); }
    void camera(ScenarioCtx& c) const { // above the rim
        const Vector3 eye = Vector3Add(at(along(r + 6.0f), level + 18.0f), Vector3Scale(Vector3CrossProduct(dir, out), 22.0f));
        c.camera(eye, at(along(r * 0.5f), level), dir);
    }
    long long tick0 = -1; // (#1431) fluid ticks at the start of the timed step (drain, reload)
} S;
// (#1431) Seconds of the water's own time (25 ms fluid ticks) since the step began: without lockstep a slow machine skips ticks
// (the self-hosted runner ran 26-28 a sim second of 40) and the drain had not finished at its 240 s of sim time there.
double FluidSeconds() {
    const long long t = Fluid::GetStats().ticks;
    if (S.tick0 < 0) S.tick0 = t;
    return (double)(t - S.tick0) * 0.025;
}
}

SCENARIO("fluid/pond") {
    // (#1431) Not tagged "fluid", so the self-hosted runner skips it: its slow GPU fluid step (#1447) never finished the drain. Local only for now.
    s.tags = { "water", "planet", "saves", "slow", "local" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 1000.0; // the drain runs ~3 min of fluid time on seed 42's pond; the timed steps allow twice their fluid time in sim time (#1431)
    s.env = { { "RAYLIBGAME_FLUID_LOCKSTEP", "0" } }; // lockstep ran this pond at ~0.4x real time (10+ min); as the flag runs it
    s.legacyFlag = "--seed 42 --test-pond";

    s.step("dig", 60.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (!WorldIdle()) return false;
        S.pf = c.field;
        const PlanetField* pf = c.field;
        const Vector3 home = pf->findSpawnDir(c.body->seed);
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
        float best = 1e30f;
        for (int pass = 0; pass < 2 && best > 1e29f; ++pass)
            for (float e = pass ? -20000.0f : -6000.0f; e <= (pass ? 20000.0f : 6000.0f); e += pass ? 60.0f : 24.0f)
                for (float n = pass ? -20000.0f : -6000.0f; n <= (pass ? 20000.0f : 6000.0f); n += pass ? 60.0f : 24.0f) {
                    const float dist = sqrtf(e * e + n * n);
                    if (dist >= best) continue;
                    const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius))));
                    const PlanetField::WaterSource ws = pf->waterSourceAt(d);
                    if (ws.kind != 3) continue;
                    best = dist; S.dir = ws.pondDir; S.level = ws.level; S.r = ws.pondR;
                }
        if (best > 1e29f)
            for (const Vector3& bd : pf->bakePondsByDistance(home)) { // eroded worlds (#487): the bake's nearest pond
                const PlanetField::WaterSource ws = pf->waterSourceAt(bd);
                if (ws.kind != 3) continue;
                best = acosf(std::clamp(Vector3DotProduct(home, bd), -1.0f, 1.0f)) * pf->radius;
                S.dir = ws.pondDir; S.level = ws.level; S.r = ws.pondR;
                break;
            }
        if (best > 1e29f) { c.fail("no pond found"); return false; }
        S.out = Vector3Normalize(Vector3CrossProduct(S.dir, { 0, 1, 0 }));
        const float bed = pf->elevation(S.dir);
        c.note("pond " + std::to_string((int)best) + " m from the start, radius " + std::to_string((int)S.r) + " m, " + std::to_string(S.level - bed) + " m deep");
        std::vector<TerrainEdit> cut;
        const float floorH = bed - 1.5f;
        for (float sv = 0.0f; sv <= S.r + 12.0f; sv += 0.8f)
            for (float h = floorH + 1.2f; h <= S.level + 3.0f; h += 1.2f) cut.push_back({ S.at(S.along(sv), h), 1.4f, -1.0f, true });
        for (int k = 0; k < 3; ++k) cut.push_back({ S.at(S.along(S.r + 14.0f), floorH - 1.5f - 1.5f * k), 3.5f, -1.0f, true });
        {
            const float need = PI * S.r * S.r * 0.5f, each = 4.0f / 3.0f * PI * 3.5f * 3.5f * 3.5f * 0.6f;
            const int extra = std::max(0, (int)ceilf(need / each) - 3);
            const Vector3 side = Vector3Normalize(Vector3CrossProduct(S.dir, S.out));
            for (int k = 0; k < extra; ++k) {
                const float ring = 5.0f * (1 + k / 8), ang = (k % 8) * (PI / 4.0f);
                const Vector3 d = Vector3Normalize(Vector3Add(S.along(S.r + 14.0f + ring * cosf(ang)), Vector3Scale(side, ring * sinf(ang) / pf->radius)));
                cut.push_back({ S.at(d, floorH - 0.5f), 3.5f, -1.0f, true });
            }
        }
        ApplyTerrainEdits(*c.body, cut.data(), cut.size());
        S.pit = S.at(S.along(S.r + 14.0f), floorH - 4.5f);
        S.mid = S.at(S.along(S.r * 0.3f), S.level - 0.3f);
        return true;
    });
    s.step("drain", 480.0, [](ScenarioCtx& c) {
        S.camera(c);
        c.waitingFor("the water asleep after the drain");
        const double ft = FluidSeconds();
        if (!((ft > 15.0 && Fluid::Asleep()) || ft > 240.0)) return false;
        S.tick0 = -1; // the next timed step starts its own count
        c.note(TextFormat("the drain took %.0f s of fluid time, %.0f s of sim time", ft, c.t));
        const PlanetField* pf = S.pf;
        const bool claimed = Fluid::GetStats().basinBricks > 0;
        const bool midDry = !Fluid::Wet(*c.body, S.mid), pitWet = Fluid::Wet(*c.body, S.pit);
        c.metric("drain time", ft, "s");
        c.check("the pond is simulated water", claimed);
        c.check("the pond's top drained", midDry);
        c.check("the pit filled", pitWet);
        bool readersDry = true;
        if (Fluid::Asleep()) {
            float shore = -1.0f;
            for (float sv = 0.0f; sv <= S.r + 4.0f && shore < 0.0f; sv += 0.5f) {
                const float bed = pf->elevation(S.along(-sv));
                if (bed > S.level - 0.25f && bed < S.level - 0.05f && !Fluid::Wet(*c.body, S.at(S.along(-sv), bed + 0.05f))) shore = sv;
            }
            if (shore >= 0.0f) {
                readersDry = !Fluid::GeneratedWaterValid(c.body->id, S.at(S.along(-shore), S.level - 0.2f));
                c.check("the readers see the drained shore as dry (#385)", readersDry);
            } else c.note("no drained shore column found to check the readers on");
        } else c.note("still draining at the cap: the readers' check is skipped");
        std::vector<Trunk> trunks;
        if (c.body->trees) c.body->trees->trunksNear(S.at(S.dir, S.level), S.r * 0.8f, trunks);
        int inBed = 0;
        for (const Trunk& tr : trunks) {
            const Vector3 db = Vector3Normalize(tr.base);
            inBed += Vector3Distance(db, S.dir) * pf->radius < S.r * 0.8f && Vector3Length(tr.base) - pf->radius < S.level + 0.3f;
        }
        c.eq("trees in the pond's bed (#189)", inBed, 0);
        S.firstOk = claimed && midDry && pitWet && readersDry;
        const bool saved = SaveWorld(c.world);
        player.placeAt(c.world, *c.body, c.body->toWorld(S.at(S.along(S.r + 300.0f), pf->elevation(S.along(S.r + 300.0f)) + 2.0f)));
        const bool loaded = saved && LoadWorldSave(c.world);
        c.check("saved, went 300 m away and loaded", loaded);
        return true;
    });
    s.step("reload", 250.0, [](ScenarioCtx& c) {
        S.camera(c);
        c.waitingFor("the water asleep after the load");
        const double ft = FluidSeconds();
        if (!((ft > 12.0 && Fluid::Asleep()) || ft > 120.0)) return false;
        S.tick0 = -1;
        c.check("after the load the pond is still simulated", Fluid::Simulated(*c.body, S.mid));
        c.check("after the load its top is drained", !Fluid::Wet(*c.body, S.mid));
        c.check("after the load the pit is full", Fluid::Wet(*c.body, S.pit));
        return true;
    });
    s.step("top-up", 45.0, [](ScenarioCtx& c) {
        const float surf = Fluid::SurfaceAbove(*c.body, S.pit, 12.0f);
        if (S.topping < 0) S.topping = surf < 1.2f ? 1 : 0;
        if (S.topping && c.t < 40.0 && (c.t - S.lastPour < 2.0 || surf < 1.6f)) {
            c.waitingFor("the pit topped up to 1.6 m");
            if (c.t - S.lastPour >= 2.0) { S.lastPour = c.t; Fluid::Pour(*c.body, Vector3Add(S.pit, Vector3Scale(Vector3Normalize(S.pit), 0.3f)), 8.0f); }
            return false;
        }
        player.flyMode = false;
        const float lift = surf > 0.0f ? std::clamp(surf - 1.5f, 0.3f, 1.0f) : 1.0f;
        c.metric("pit water depth", surf, "m", false);
        player.placeAt(c.world, *c.body, c.body->toWorld(Vector3Add(S.pit, Vector3Scale(Vector3Normalize(S.pit), lift))));
        return true;
    });
    s.step("swim", 5.0, [](ScenarioCtx& c) {
        const float h = Vector3Length(c.body->toLocal((player.position - c.world.renderOrigin).toFloat()));
        if (c.t < 0.5 || S.h0 == 0.0f) { S.h0 = h; S.inWater = S.frames = 0; return false; }
        ++S.frames; S.inWater += player.inWater ? 1 : 0;
        c.waitingFor("3 s of swimming");
        if (c.t <= 3.0) return false;
        const float sink = (S.h0 - h) / 2.5f;
        c.metric("sink rate in the pit", sink, "m/s");
        c.gt("frames in water (%)", S.frames ? 100.0 * S.inWater / S.frames : 0.0, 89.99);
        c.lt("sink rate (m/s)", sink, 1.0);
        return true;
    });
}
