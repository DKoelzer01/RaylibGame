// visual/small-flora (#1995 phase 3): the small flora's size classes. Near the start, the nearest bush and the nearest
// loose stone of each class (small bush .. thicket; pebbles, stones, a rock, a boulder), each from a few metres at its
// height, with whatever stands in the view taken away; then the grass of an open savanna (tall drifts) and of the chalk
// downs (short turf). The check: every class grows near the start and each class's plants are bigger than the last's;
// the shots are for a person to judge (the sizes are GUESSes).
#include "scenario.h"
#include "biomevariants.h"
#include "foliage.h"
#include "trees.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace {
constexpr int kShots = 2 * kSmallFloraClasses; // bushes 0-3, then stones 0-3
const char* kNames[kShots] = { "small bush", "bush", "large bush", "thicket", "pebbles", "stones", "rock", "boulder" };
struct Visit { bool skip = false; Plant plant{}; Vector3 eye{}, target{}, up{}; bool cleared = false; double clearedAt = 0.0; };
Visit V[kShots];
Vector3 gHome{};
Vector3 Tangent(Vector3 up) { return Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 })); }
// Whether q stands between eye and p (within a metre and a half of the sight line, or by the eye), or hides it.
bool InTheWay(const Plant& p, Vector3 eye, const Plant& q) {
    if (q.id == p.id) return false;
    const Vector3 seg = Vector3Subtract(p.base, eye);
    const float t = std::clamp(Vector3DotProduct(Vector3Subtract(q.base, eye), seg) / std::max(Vector3LengthSqr(seg), 1e-3f), 0.0f, 1.2f);
    Vector3 off = Vector3Subtract(q.base, Vector3Add(eye, Vector3Scale(seg, t)));
    off = Vector3Subtract(off, Vector3Scale(p.up, Vector3DotProduct(off, p.up)));
    return Vector3Length(off) < 1.5f + q.radius || Vector3Distance(q.base, eye) < 2.0f + q.radius;
}
struct GrassShot { const char* variant; const char* name; bool skip = false; Vector3 eye{}, target{}, up{}; };
GrassShot G[2] = { { "savanna.grassland", "savanna grass" }, { "meadow.chalk_downs", "chalk downs turf" } };
}

SCENARIO("visual/small-flora") {
    s.tags = { "visual", "foliage", "planet" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day", "--no-hud" };
    s.env = { { "RAYLIBGAME_WEATHER_FORCE", "0" } }; // clear: rain streaks hide small things
    s.timeout = 600.0;
    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("pick the subjects", 90.0, [](ScenarioCtx& c) {
        const Vector3 home = c.field->findSpawnDir(c.body->seed);
        gHome = Vector3Scale(home, c.field->radius + c.field->elevation(home));
        const Vector3 up = home;
        c.camera(Vector3Add(gHome, Vector3Scale(up, 40.0f)), Vector3Add(gHome, Vector3Scale(Tangent(up), 30.0f)), up);
        c.waitingFor("the tree tiles around the start, an idle world");
        if (!c.body->trees) { c.fail("no trees on the planet"); return true; }
        if (c.t < 2.0 || ((c.body->trees->readyTilesNear(gHome, 300.0f) < 40 || !WorldIdle()) && c.t < 60.0)) return false;
        std::vector<Plant> plants;
        c.body->trees->plantsNear(gHome, 400.0f, plants);
        int counts[kShots] = {};
        float meanTop[kShots] = {};
        for (const Plant& p : plants) {
            if (p.tree()) continue;
            const int k = (p.stone() ? kSmallFloraClasses : 0) + p.smallClass();
            ++counts[k]; meanTop[k] += p.height;
        }
        std::string line;
        for (int k = 0; k < kShots; ++k) { meanTop[k] /= std::max(1, counts[k]); line += TextFormat(" %s %d (%.2f m)", kNames[k], counts[k], meanTop[k]); }
        c.note("small flora within 400 m (count, mean height):" + line);
        for (int k = 0; k < kShots; ++k) {
            c.check(std::string(kNames[k]) + " grows near the start", counts[k] > 0);
            if (k % kSmallFloraClasses) c.check(std::string(kNames[k]) + " is bigger than the class below", counts[k] == 0 || counts[k - 1] == 0 || meanTop[k] > meanTop[k - 1]);
        }
        // Each class's subject: the nearest to the start on gentle open ground, framed from a few metres at its height.
        for (int k = 0; k < kShots; ++k) {
            std::vector<Plant> cands;
            for (const Plant& p : plants) if (!p.tree() && p.stone() == (k >= kSmallFloraClasses) && p.smallClass() == k % kSmallFloraClasses) cands.push_back(p);
            std::sort(cands.begin(), cands.end(), [](const Plant& a, const Plant& b) { return Vector3Distance(a.base, gHome) < Vector3Distance(b.base, gHome); });
            if (cands.size() > 8) cands.resize(8);
            int bestCost = 1 << 30;
            for (const Plant& p : cands) {
                // From above the grass (an eye at least 1.8 m up), looking down at it.
                const float h = std::max(p.height, 0.25f), back = std::clamp(1.5f + 2.2f * std::max(h, p.radius), 2.5f, 9.0f);
                const Vector3 e1 = Tangent(p.up), e2 = Vector3CrossProduct(p.up, e1);
                for (int d = 0; d < 8; ++d) {
                    const float a = d * 0.785398f;
                    const Vector3 side = Vector3Add(Vector3Scale(e1, cosf(a)), Vector3Scale(e2, sinf(a)));
                    Vector3 eye = Vector3Add(p.base, Vector3Add(Vector3Scale(side, back), Vector3Scale(p.up, std::max(1.8f, 0.6f * h + 0.5f))));
                    const Vector3 ed = Vector3Normalize(eye);
                    const float ground = c.field->radius + c.field->elevation(ed) + 1.8f;
                    if (Vector3Length(eye) < ground) eye = Vector3Scale(ed, ground);
                    const Vector3 target = Vector3Add(p.base, Vector3Scale(p.up, 0.4f * h));
                    auto rock = [&](Vector3 q) { return c.field->density(q, nullptr) >= ISO_LEVEL; };
                    bool open = true;
                    for (int i = 1; i < 16 && open; ++i) open = !rock(Vector3Lerp(eye, target, i / 16.0f));
                    for (const Vector3& dv : { e1, Vector3Negate(e1), e2, Vector3Negate(e2), p.up }) open = open && !rock(Vector3Add(eye, Vector3Scale(dv, 1.0f)));
                    int cost = open ? 0 : 1000;
                    cost += (int)((Vector3Length(eye) - ground) * 4.0f);
                    for (const Plant& q : plants) cost += InTheWay(p, eye, q) ? (q.stone() ? 30 : 1) : 0;
                    if (cost < bestCost) { bestCost = cost; V[k].plant = p; V[k].eye = eye; V[k].target = target; V[k].up = p.up; }
                }
            }
            if (cands.empty()) { V[k].skip = true; continue; }
            int cleared = 0;
            for (const Plant& q : plants) if (InTheWay(V[k].plant, V[k].eye, q) && !q.stone()) { c.body->trees->remove(q.id); ++cleared; }
            c.note(TextFormat("%s: %.2f m tall, %.2f m reach, %.0f m from the start; %d plants cleared%s", kNames[k], V[k].plant.height, V[k].plant.radius,
                              Vector3Distance(V[k].plant.base, gHome), cleared, bestCost >= 1000 ? " (terrain in the sight line)" : ""));
        }
        return true;
    });
    for (int k = 0; k < kShots; ++k) {
        s.step(std::string("look at the ") + kNames[k], 60.0, [k](ScenarioCtx& c) {
            if (V[k].skip) return true;
            // The grass hid pebbles and seedlings (as in play): for judging sizes it fades out within a metre of the eye
            // (1.8 m up, so none shows); the grass shots below restore it.
            kGrassDrawDistance = 1.0f;
            c.camera(V[k].eye, V[k].target, V[k].up);
            c.waitingFor("the view settled (an idle world)");
            if (c.t < 2.0 || (!WorldIdle() && c.t < 30.0)) return false;
            if (!V[k].cleared) { // again from the eye: tiles that loaded only now hold plants the first pass never saw
                V[k].cleared = true;
                std::vector<Plant> near;
                c.body->trees->plantsNear(V[k].eye, 30.0f, near);
                int n = 0;
                for (const Plant& q : near) if (InTheWay(V[k].plant, V[k].eye, q) && !q.stone()) { c.body->trees->remove(q.id); ++n; }
                V[k].clearedAt = n ? c.t : 0.0;
                return false;
            }
            if (c.t < V[k].clearedAt + 2.0 || (!WorldIdle() && c.t < V[k].clearedAt + 30.0)) return false;
            c.shot(kNames[k], 100.0);
            return true;
        });
    }
    for (int g = 0; g < 2; ++g) {
        s.step(std::string("to the ") + G[g].name, 120.0, [g](ScenarioCtx& c) {
            kGrassDrawDistance = 100.0f; // (foliage.h's default)
            if (!G[g].eye.x && !G[g].eye.y && !G[g].eye.z) {
                const int v = BiomeVariants().find(G[g].variant);
                if (v < 0) { c.note(std::string("no variant ") + G[g].variant + " (skipped)"); G[g].skip = true; return true; }
                const Vector3 home = Vector3Normalize(gHome);
                // A flat spot in the patch: of points around its nearest bit, the one whose ground varies least.
                const Vector3 base = c.field->nearestVariantDir(home, v), b1 = Tangent(base), b2 = Vector3CrossProduct(base, b1);
                Vector3 best = base;
                float bestSlope = 1e9f;
                for (int j = 0; j < 17; ++j) {
                    const float a = (j - 1) * 0.785398f, r = j == 0 ? 0.0f : (j <= 8 ? 150.0f : 300.0f) / c.field->radius;
                    const Vector3 d = Vector3Normalize(Vector3Add(base, Vector3Add(Vector3Scale(b1, cosf(a) * r), Vector3Scale(b2, sinf(a) * r))));
                    if (c.field->variantAt(d, c.field->elevation(d)) != v || c.field->elevation(d) < 6.0f) continue;
                    const float e0 = c.field->elevation(d), k = 12.0f / c.field->radius;
                    const float slope = fabsf(c.field->elevation(Vector3Normalize(Vector3Add(d, Vector3Scale(b1, k)))) - e0) + fabsf(c.field->elevation(Vector3Normalize(Vector3Add(d, Vector3Scale(b2, k)))) - e0);
                    if (slope < bestSlope) { bestSlope = slope; best = d; }
                }
                const Vector3 at = Vector3Scale(best, c.field->radius + c.field->elevation(best)), e1 = Tangent(best);
                G[g].up = best;
                G[g].eye = Vector3Add(at, Vector3Scale(best, 1.6f));
                const Vector3 ahead = Vector3Normalize(Vector3Add(best, Vector3Scale(e1, 25.0f / c.field->radius)));
                G[g].target = Vector3Scale(ahead, c.field->radius + c.field->elevation(ahead) + 0.3f);
                c.note(TextFormat("%s: %.0f m from the start, slope %.1f m over 12 m", G[g].name, Vector3Distance(at, gHome), bestSlope));
            }
            if (G[g].skip) return true;
            c.camera(G[g].eye, G[g].target, G[g].up);
            c.waitingFor("the ground and its grass around the eye (an idle world)");
            if (c.t < 4.0 || (!WorldIdle() && c.t < 90.0)) return false;
            c.shot(G[g].name, 100.0);
            return true;
        });
    }
    s.step("done", 1.0, [](ScenarioCtx& c) { c.keepCapture("the small flora's size classes (#1995 phase 3): for a person to judge"); return true; });
}
