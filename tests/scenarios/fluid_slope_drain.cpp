// fluid/slope-drain (#454, playtest 2026-11-24: "Water clings to slopes as a thin film after all water has flowed away"): on the
// bench planet, the steepest open slope (9-35 degrees each 2 m, still falling 8 m down) within 80 m of the spawn; 0.5 m^3 poured at its
// top runs down. A film under the spread limit (~1 cm) never moves, and the bricks fall asleep after 2-3 s, so their films
// stayed for good. Checked: the water left on the slope (the simulation's cells in a box over its upper 6 m, 2 m wide) 150 s
// later is under 2 litres. The amount at 20 s (runoff done, films left) is reported for comparison.
#include "scenario.h"
#include "fluid.h"
#include "player.h"
#include "world.h"
#include <raymath.h>
#include <cmath>

namespace {
struct State {
    Vector3 top{}, down{}, up{}, side{};
    double nextPour = 0.0, at20 = -1.0;
    float poured = 0.0f, grade = 0.0f;
} S;
double SlopeWater(ScenarioCtx& c) { // m^3 on the slope: a box from the top 6 m down it, 2 m wide, 1.5 m over and under the ground's mean
    const Vector3 mid = Vector3Add(S.top, Vector3Scale(S.down, 3.0f));
    return Fluid::VolumeIn(*c.body, mid, S.down, S.side, S.up, { 3.0f, 1.0f, 1.5f + 3.0f * S.grade });
}
}

SCENARIO("fluid/slope-drain") {
    s.tags = { "fluid", "water" };
    s.timeout = 300.0;
    s.step("find a slope", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle() || !c.field) return false;
        const Vector3 me = Vector3Normalize(c.body->toLocalD(player.position).toFloat());
        const float R = c.field->radius;
        const Vector3 e = Vector3Normalize(Vector3CrossProduct(me, fabsf(me.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        const Vector3 n = Vector3CrossProduct(me, e);
        auto dirAt = [&](float x, float y) { return Vector3Normalize(Vector3Add(me, Vector3Add(Vector3Scale(e, x / R), Vector3Scale(n, y / R)))); };
        float best = 0.0f;
        for (float x = -80.0f; x <= 80.0f; x += 4.0f)
            for (float y = -80.0f; y <= 80.0f; y += 4.0f)
                for (int a = 0; a < 8; ++a) {
                    const float ang = a * PI / 4.0f, dx = cosf(ang), dy = sinf(ang);
                    const float h0 = c.field->elevation(dirAt(x, y));
                    if (h0 < 1.0f) continue; // over the sea's level
                    // Falling all the way down 8 m, each 2 m by a 9-35 degree grade.
                    bool ok = true;
                    float prev = h0;
                    for (int k = 1; k <= 4 && ok; ++k) {
                        const float h = c.field->elevation(dirAt(x + dx * 2.0f * k, y + dy * 2.0f * k));
                        const float g = (prev - h) / 2.0f;
                        ok = g > 0.15f && g < 0.7f;
                        prev = h;
                    }
                    const float g = (h0 - prev) / 8.0f;
                    if (!ok || g <= best || Fluid::AnyWaterNear(*c.body, Vector3Scale(dirAt(x, y), R + h0), 10.0f)) continue;
                    best = g;
                    S.top = Vector3Scale(dirAt(x, y), R + h0);
                    const Vector3 low = Vector3Scale(dirAt(x + dx * 8.0f, y + dy * 8.0f), R + prev);
                    S.up = dirAt(x, y);
                    S.down = Vector3Normalize(Vector3Subtract(Vector3Subtract(low, S.top), Vector3Scale(S.up, Vector3DotProduct(Vector3Subtract(low, S.top), S.up))));
                    S.side = Vector3CrossProduct(S.up, S.down);
                    S.grade = g;
                }
        c.gt("a 9-35 degree slope found within 80 m (grade)", best, 0.0f);
        if (best <= 0.0f) return true;
        c.note(TextFormat("slope grade %.2f (%.0f degrees)", S.grade, atanf(S.grade) * RAD2DEG));
        return true;
    });
    s.step("pour at its top and wait (150 s)", 200.0, [](ScenarioCtx& c) {
        if (S.grade <= 0.0f) return true;
        // The camera across the slope, a few metres off its side.
        const Vector3 mid = Vector3Add(S.top, Vector3Scale(S.down, 3.0f));
        c.camera(Vector3Add(Vector3Add(mid, Vector3Scale(S.side, 6.0f)), Vector3Scale(S.up, 3.0f)), mid, S.up);
        if (c.t >= S.nextPour && c.t < 2.0) { S.nextPour = c.t + 0.1; S.poured += Fluid::Pour(*c.body, Vector3Add(S.top, Vector3Scale(S.up, 0.3f)), 0.025f); }
        if (S.at20 < 0.0 && c.t >= 20.0) S.at20 = SlopeWater(c);
        return c.t >= 150.0;
    });
    s.step("report", 1.0, [](ScenarioCtx& c) {
        if (S.grade <= 0.0f) return true;
        const double left = SlopeWater(c);
        c.note(TextFormat("%.2f m^3 poured; on the slope %.4f m^3 at 20 s, %.4f m^3 at 150 s", S.poured, S.at20, left));
        c.metric("water left on the slope at 150 s", (float)(left * 1000.0), "litres");
        c.gt("water poured", S.poured, 0.4f);
        c.lt("the films on the slope dried (litres left at 150 s)", (float)(left * 1000.0), 2.0f);
        return true;
    });
}
