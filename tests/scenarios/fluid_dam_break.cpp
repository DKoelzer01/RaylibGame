// fluid/dam-break (#870, moved from --test-dam-break, #680): a physics reference. A channel 4 m wide and 2.5 m deep; its first
// 20 m hold water 2 m deep behind a 2 m wall of ground, which is cut away. Ritter's solution for a dam break on a dry,
// frictionless bed: the front runs at 2 sqrt(g h0), the depth at the dam is 4/9 h0 and the discharge (8/27) sqrt(g) h0^1.5
// per metre of width. Measured in the first 3 s (the wave reflected off the reservoir's back reaches the dam at ~4.5 s).
// PASS (as the flag): the reservoir filled over 1.5 m and the mean front speed to 10 m is 0.4-1.3 of theory.
#include "scenario.h"
#include "fluid.h"
#include "waterfield.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>

namespace {
struct State {
    Vector3 home{}, ax{}, ay{};
    float r0 = 0.0f, h0 = 0.0f, floor = -2.5f, poured = 0.0f, depthDam = -1.0f, qDam = 0.0f, vol1 = -1.0f, vol2 = -1.0f;
    double front[3] = { -1.0, -1.0, -1.0 }, nextPour = 0.0;
    int k = 0;
    Vector3 pt(float x, float y, float h) const { // x along the channel, y across, h over the ground (m), body-local
        const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(home, r0), Vector3Add(Vector3Scale(ax, x), Vector3Scale(ay, y))));
        return Vector3Scale(d, r0 + h);
    }
    TerrainEdit box(float x0, float x1) const { // the channel's cut from x0 to x1, up through the air (#693)
        TerrainEdit e{};
        e.center = pt(0.5f * (x0 + x1), 0.0f, 0.25f);
        e.carve = true; e.box = true; e.half = { 0.5f * (x1 - x0), 2.0f, 2.75f };
        e.ax = ax; e.ay = ay; e.az = Vector3Normalize(e.center);
        e.radius = Vector3Length(e.half);
        return e;
    }
    float depthAt(const Planetoid& p, float x) const { const float s = Fluid::SurfaceAbove(p, pt(x, 0.0f, floor + 0.55f), 4.0f); return s > 0.0f ? s + 0.55f : 0.0f; }
    void camera(ScenarioCtx& c) const { c.camera(pt(-4.0f, -9.0f, 9.0f), pt(-4.0f, 0.0f, -2.0f), home); } // over the wall
} S;
}

SCENARIO("fluid/dam-break") {
    s.tags = { "fluid", "water", "physics-reference" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 150.0;
    s.legacyFlag = "--sandbox --day --test-dam-break";

    // The legacy test starts at 2 s of play once the world has been idle (#814: thread timing must not pick the phase).
    s.step("settle", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("2 s of play and an idle world");
        return c.T >= 2.0 && WorldIdle();
    });
    s.step("dig", 1.0, [](ScenarioCtx& c) {
        S.home = c.field->findSpawnDir(c.body->seed);
        S.ax = Vector3Normalize(Vector3CrossProduct(S.home, { 0, 1, 0 })); S.ay = Vector3CrossProduct(S.ax, S.home);
        S.r0 = c.field->radius + c.field->elevation(S.home);
        std::vector<TerrainEdit> cut = { S.box(-20.0f, 0.0f), S.box(2.0f, 50.0f) }; // a 2 m wall (#691)
        ApplyTerrainEdits(*c.body, cut.data(), cut.size());
        c.note("a 20 m reservoir behind a 2 m wall, a 50 m channel, 4 m wide, 2.5 m deep (cells " + std::to_string(WaterField::CellSize()) + " m)");
        S.camera(c);
        return true;
    });
    // Fill to 2 m: 0.25 m^3 every 0.1 s over the floor, kept 6.4 m from the wall (#691, #680); done at a 2 m mean depth or 200 m^3.
    s.step("fill", 90.0, [](ScenarioCtx& c) {
        S.camera(c);
        c.waitingFor("a 2 m mean depth (" + std::to_string((int)S.poured) + " m^3 poured)");
        if (c.t < 2.0) return false;
        if (c.t >= S.nextPour) {
            S.nextPour = c.t + 0.1;
            const float x = -19.0f + 1.4f * (float)(S.k % 10), y = (S.k / 10) % 2 ? 1.0f : -1.0f;
            ++S.k;
            S.poured += Fluid::Pour(*c.body, S.pt(x, y, -2.3f), 0.25f);
        }
        float sum = 0.0f;
        for (float x : { -18.0f, -14.0f, -10.0f, -6.0f, -2.0f }) { const float sf = Fluid::SurfaceAbove(*c.body, S.pt(x, 0.0f, -2.5f + 0.55f), 4.0f); sum += sf > 0.0f ? sf + 0.55f : 0.0f; }
        return sum / 5.0f >= 2.0f || S.poured >= 200.0f;
    });
    s.step("calm", 45.0, [](ScenarioCtx& c) {
        S.camera(c);
        c.waitingFor("the water asleep");
        return (c.t > 8.0 && Fluid::Asleep()) || c.t > 40.0;
    });
    s.step("break", 1.0, [](ScenarioCtx& c) {
        S.camera(c);
        for (float h = -1.0f; h > -3.0f; h -= 0.02f) if (c.body->densityAt(S.pt(-10.0f, 0.0f, h)) >= ISO_LEVEL) { S.floor = h; break; }
        S.h0 = S.depthAt(*c.body, -10.0f);
        const float beyond = Fluid::WaterAt(*c.body, S.pt(10.0f, 0.0f, S.floor + 0.5f), 4.0f);
        std::vector<TerrainEdit> cut = { S.box(-0.2f, 2.2f) }; // the wall goes
        ApplyTerrainEdits(*c.body, cut.data(), cut.size());
        c.cue();
        c.metric("reservoir depth", S.h0, "m", false);
        c.metric("poured", S.poured, "m^3", false);
        c.metric("past the wall before the break", beyond, "m^3");
        return true;
    });
    s.step("flow", 10.0, [](ScenarioCtx& c) {
        S.camera(c);
        const double t = c.t;
        const float dist[3] = { 5.0f, 10.0f, 20.0f };
        for (int i = 0; i < 3; ++i)
            if (S.front[i] < 0.0 && Fluid::WaterAt(*c.body, S.pt(2.0f + dist[i], 0.0f, S.floor + 0.5f), 0.5f) > 0.02f) S.front[i] = t;
        if (t > 1.4 && S.depthDam < 0.0f) S.depthDam = S.depthAt(*c.body, 1.0f);
        const float qNow = Vector3DotProduct(Fluid::FlowAt(*c.body, S.pt(2.5f, 0.0f, -1.8f), 2.0f), S.ax);
        if (t > 1.0 && t < 2.0) S.qDam = std::max(S.qDam, qNow);
        auto downstream = [&]() { float v = 0.0f; for (float x = 4.5f; x < 50.0f; x += 4.0f) v += Fluid::WaterAt(*c.body, S.pt(x, 0.0f, S.floor + 2.0f), 2.0f); return v; };
        if (t > 1.0 && S.vol1 < 0.0f) S.vol1 = downstream();
        if (t > 2.0 && S.vol2 < 0.0f) S.vol2 = downstream();
        c.waitingFor("3 s after the break");
        if (t <= 3.0) return false;
        const float g = 9.81f, ct = 2.0f * sqrtf(g * std::max(S.h0, 0.01f));
        const float dTheory = 4.0f / 9.0f * S.h0, qTheory = 8.0f / 27.0f * sqrtf(g) * powf(S.h0, 1.5f) * 4.0f;
        const float speed = S.front[1] > 0.0 ? 10.0f / (float)S.front[1] : 0.0f;
        c.metric("front at 5 m", S.front[0], "s"); c.metric("front at 10 m", S.front[1], "s"); c.metric("front at 20 m", S.front[2], "s");
        c.metric("front speed / theory", speed / ct, "", false);
        c.metric("depth at the dam / theory", S.depthDam / std::max(dTheory, 0.01f), "", false);
        c.metric("discharge / theory", S.qDam / std::max(qTheory, 0.01f), "", false);
        c.metric("discharge by volume / theory", (S.vol2 - S.vol1) / std::max(qTheory, 0.01f), "", false);
        // The legacy verdict: h0 > 1.5, the front reached 10 m, its mean speed 0.4-1.3 of theory.
        c.gt("reservoir depth (m)", S.h0, 1.5);
        c.check("front reached 10 m", S.front[1] > 0.0);
        c.between("front speed / theory", speed / ct, 0.4, 1.3);
        return true;
    });
}
