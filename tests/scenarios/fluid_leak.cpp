// fluid/leak (#870, moved from --test-leak, #680): a 20 m reservoir 4 m wide, 2 m deep behind a 2 m wall of ground, an empty
// channel 6 m deep beyond; a 1 x 1 m hole is cut through the wall 0.5 m over the floor. The early outflow (the volume past
// the wall, 1-3 s) against an orifice's sqrt(2 g h) x area, and after 60 s the pool against an orifice-then-weir model of
// the drain to the hole's bottom. No water on the grass (a bad fill or a leak around the edges), none past the wall before
// the hole. PASS as the flag: h0 > 1.5 m, outflow 0.4-1.6 of ideal, final depth within 0.12 m of the model, grass < 0.05 m^3,
// past the wall before the hole < 0.2 m^3.
#include "scenario.h"
#include "fluid.h"
#include "waterfield.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>

namespace {
constexpr float kDrain = 60.0f; // s judged after the hole opens
struct State {
    Vector3 home{}, ax{}, ay{};
    float r0 = 0.0f, h0 = 0.0f, vol1 = -1.0f, vol3 = -1.0f, pastBefore = 0.0f, poured = 0.0f;
    double nextPour = 0.0;
    int k = 0;
    Vector3 pt(float x, float y, float h) const {
        const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(home, r0), Vector3Add(Vector3Scale(ax, x), Vector3Scale(ay, y))));
        return Vector3Scale(d, r0 + h);
    }
    TerrainEdit box(float x0, float x1, float y0, float y1, float h0b, float h1b) const {
        TerrainEdit e{};
        e.center = pt(0.5f * (x0 + x1), 0.5f * (y0 + y1), 0.5f * (h0b + h1b));
        e.carve = true; e.box = true; e.half = { 0.5f * (x1 - x0), 0.5f * (y1 - y0), 0.5f * (h1b - h0b) };
        e.ax = ax; e.ay = ay; e.az = Vector3Normalize(e.center);
        e.radius = Vector3Length(e.half);
        return e;
    }
    float depthAt(const Planetoid& p, float x) const { const float s = Fluid::SurfaceAbove(p, pt(x, 0.0f, -2.5f + 0.55f), 4.0f); return s > 0.0f ? s + 0.55f : 0.0f; }
    float past(const Planetoid& p) const { float v = 0.0f; for (float x = 4.0f; x < 40.0f; x += 4.0f) for (float h : { -5.0f, -1.0f }) v += Fluid::WaterAt(p, pt(x, 0.0f, h), 2.0f); return v; }
    float onGrass(const Planetoid& p) const { float v = 0.0f; for (float x = -19.0f; x < 40.0f; x += 3.0f) for (float y : { -3.5f, 3.5f }) v += Fluid::WaterAt(p, pt(x, y, 0.25f), 0.5f); return v; }
    // 45 degrees off the channel on its downstream side, 45 degrees up: the trench, the hole and the dam's edges (#691).
    void camera(ScenarioCtx& c) const { c.camera(pt(6.0f, -6.0f, 6.0f), pt(0.5f, 0.0f, -1.5f), home); }
} S;
}

SCENARIO("fluid/leak") {
    s.tags = { "fluid", "water", "physics-reference" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 200.0;
    s.legacyFlag = "--sandbox --day --test-leak";

    s.step("settle", 30.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("dig", 1.0, [](ScenarioCtx& c) {
        S.home = c.field->findSpawnDir(c.body->seed);
        S.ax = Vector3Normalize(Vector3CrossProduct(S.home, { 0, 1, 0 })); S.ay = Vector3CrossProduct(S.ax, S.home);
        S.r0 = c.field->radius + c.field->elevation(S.home);
        std::vector<TerrainEdit> cut = { S.box(-20.0f, 0.0f, -2.0f, 2.0f, -2.5f, 3.0f), S.box(2.0f, 40.0f, -2.0f, 2.0f, -6.0f, 3.0f) };
        ApplyTerrainEdits(*c.body, cut.data(), cut.size());
        c.note("a 20 m reservoir behind a 2 m wall, a channel beyond (cells " + std::to_string(WaterField::CellSize()) + " m)");
        S.camera(c);
        return true;
    });
    s.step("fill", 90.0, [](ScenarioCtx& c) {
        S.camera(c);
        c.waitingFor("a 2 m mean depth");
        if (c.t <= 1.0) return false;
        if (c.t >= S.nextPour) {
            S.nextPour = c.t + 0.1;
            const float x = -19.0f + 1.4f * (float)(S.k % 10), y = (S.k / 10) % 2 ? 1.0f : -1.0f;
            ++S.k;
            S.poured += Fluid::Pour(*c.body, S.pt(x, y, -2.3f), 0.25f);
        }
        const float mean = 0.2f * (S.depthAt(*c.body, -18.0f) + S.depthAt(*c.body, -14.0f) + S.depthAt(*c.body, -10.0f) + S.depthAt(*c.body, -6.0f) + S.depthAt(*c.body, -2.0f));
        return S.poured >= 200.0f || mean >= 2.0f;
    });
    s.step("calm", 45.0, [](ScenarioCtx& c) { S.camera(c); c.waitingFor("the water asleep"); return (c.t > 8.0 && Fluid::Asleep()) || c.t > 40.0; });
    s.step("hole", 1.0, [](ScenarioCtx& c) {
        S.camera(c);
        c.lt("water on the grass after the fill (m^3)", S.onGrass(*c.body), 0.05);
        S.h0 = S.depthAt(*c.body, -10.0f);
        S.pastBefore = S.past(*c.body);
        std::vector<TerrainEdit> cut = { S.box(-0.2f, 2.2f, -0.5f, 0.5f, -2.0f, -1.0f) }; // 1 x 1 m through the wall, 0.5-1.5 m over the floor
        ApplyTerrainEdits(*c.body, cut.data(), cut.size());
        c.cue();
        c.note("pool " + std::to_string(S.h0) + " m deep; the hole opened");
        return true;
    });
    s.step("drain", kDrain + 5.0, [](ScenarioCtx& c) {
        S.camera(c);
        if (c.t > 1.0 && S.vol1 < 0.0f) S.vol1 = S.past(*c.body);
        if (c.t > 3.0 && S.vol3 < 0.0f) S.vol3 = S.past(*c.body);
        c.waitingFor("60 s of draining");
        if (c.t <= kDrain) return false;
        const float lip = 0.5f, g = 9.81f, head = std::max(0.0f, S.h0 - (lip + 0.5f));
        const float qIdeal = sqrtf(2.0f * g * head), q = (S.vol3 - S.vol1) / 2.0f, final = S.depthAt(*c.body, -10.0f);
        const float poolArea = S.h0 > 0.1f ? (S.poured - Fluid::PourDropped()) / S.h0 : 80.0f;
        float d = S.h0;
        for (int k = 0; k < (int)(kDrain / 0.01f); ++k) {
            const float qd = d >= lip + 1.0f ? 0.6f * sqrtf(2.0f * g * (d - (lip + 0.5f))) : d > lip ? 1.7f * powf(d - lip, 1.5f) : 0.0f;
            d -= qd * 0.01f / poolArea;
        }
        c.metric("outflow / orifice ideal", q / std::max(qIdeal, 0.01f), "", false);
        c.metric("final depth", final, "m");
        c.gt("pool depth before the hole (m)", S.h0, 1.5);
        c.between("early outflow (m^3/s)", q, 0.4 * qIdeal, 1.6 * qIdeal);
        c.near("depth after 60 s against the orifice-weir model (m)", final, d, 0.12);
        c.lt("water on the grass at the end (m^3)", S.onGrass(*c.body), 0.05);
        c.lt("past the wall before the hole (m^3)", S.pastBefore, 0.2);
        return true;
    });
}
