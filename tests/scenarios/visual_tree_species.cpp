// visual/tree-species (#1995 phase 2): each new tree species where it grows. For acacia, baobab, mangrove, willow,
// krummholz, snag, palm, kapok and tree fern: the nearest deep patch of a biome variant that lists it (biomes.txt
// `species`, the --bench-biome search), a mature tree of that species there (the largest class found, mature first), and
// the camera beside it at its mid-height, in daylight, with the biome's own tints (#2032). One shot a species, for a
// person to judge (the look is meant to change while the values are GUESSes); the check is that each species grows.
#include "scenario.h"
#include "biomevariants.h"
#include "trees.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <cmath>
#include <string>
#include <vector>

namespace {
struct Subject { const char* variant; WoodSpecies species; };
const Subject kSubjects[] = {
    { "savanna.acacia", WoodSpecies::Acacia },     { "savanna.kopje", WoodSpecies::Baobab },
    { "swamp.mangrove", WoodSpecies::Mangrove },   { "swamp.reed_marsh", WoodSpecies::Willow },
    { "alpine.krummholz", WoodSpecies::Krummholz }, { "taiga.burnt", WoodSpecies::Snag },
    { "jungle.karst_towers", WoodSpecies::Palm },  { "jungle.lowland", WoodSpecies::Kapok },
    { "jungle.cloud_forest", WoodSpecies::TreeFern },
};
constexpr int kN = sizeof(kSubjects) / sizeof(kSubjects[0]);
struct Visit { bool skip = false; int generated = 0; Vector3 patch{}; Plant tree{}; Vector3 eye{}, target{}, up{}; };
Visit V[kN];
// Whether plant q stands in the view of p from eye: near the sight line, or near the eye (crowns overhang it).
bool InTheWay(const Plant& p, Vector3 eye, const Plant& q) {
    if (q.id == p.id) return false;
    const Vector3 seg = Vector3Subtract(p.base, eye);
    const float t = std::clamp(Vector3DotProduct(Vector3Subtract(q.base, eye), seg) / std::max(Vector3LengthSqr(seg), 1e-3f), 0.0f, 1.0f);
    Vector3 off = Vector3Subtract(q.base, Vector3Add(eye, Vector3Scale(seg, t)));
    off = Vector3Subtract(off, Vector3Scale(p.up, Vector3DotProduct(off, p.up))); // across the ground only
    return Vector3Length(off) < 3.0f + 4.0f * q.scale || Vector3Distance(q.base, eye) < 12.0f + 6.0f * q.scale;
}
bool cleared2[kN] = {};
double clearedAt[kN] = {}; // step time of the second clearing
Vector3 Tangent(Vector3 up) { return Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 })); }
}

SCENARIO("visual/tree-species") {
    s.tags = { "visual", "foliage", "planet" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day", "--no-hud" };
    s.env = { { "RAYLIBGAME_WEATHER_FORCE", "0" } }; // clear: rain streaks hid the subjects
    s.timeout = 600.0;
    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    for (int k = 0; k < kN; ++k) {
        const std::string name = kWoodSpeciesNames[(int)kSubjects[k].species];
        s.step("to the " + std::string(kSubjects[k].variant), 5.0, [k](ScenarioCtx& c) {
            const int v = BiomeVariants().find(kSubjects[k].variant);
            // The nearest patch to the start and to points 4 and 8 km off each way: the one growing the most trees (savanna
            // is open country with lone trees, and a patch can be a steep hillside that grows none).
            const Vector3 home = c.field->findSpawnDir(c.body->seed), h1 = Tangent(home), h2 = Vector3CrossProduct(home, h1);
            std::vector<Vector3> homes = { home };
            for (float km : { 4.0f, 8.0f })
                for (Vector3 d : { h1, Vector3Negate(h1), h2, Vector3Negate(h2) }) homes.push_back(Vector3Add(home, Vector3Scale(d, km * 1000.0f / c.field->radius)));
            Vector3 at{};
            int trees = -1, bushes = 0, stones = 0;
            TreeScatterStats why;
            // Each search can land on the same patch, so also try points 250 and 500 m around it (a hillside patch grows
            // no trees: they need gentle slopes).
            for (size_t i = 0; i < homes.size() && v >= 0 && trees < 20; ++i) {
                const Vector3 base = c.field->nearestVariantDir(Vector3Normalize(homes[i]), v), b1 = Tangent(base), b2 = Vector3CrossProduct(base, b1);
                for (int j = 0; j < 17 && trees < 20; ++j) {
                    const float a = (j - 1) * 0.785398f, r = j == 0 ? 0.0f : (j <= 8 ? 250.0f : 500.0f) / c.field->radius;
                    const Vector3 cand = Vector3Normalize(Vector3Add(base, Vector3Add(Vector3Scale(b1, cosf(a) * r), Vector3Scale(b2, sinf(a) * r))));
                    if (c.field->variantAt(cand, c.field->elevation(cand)) != v) continue;
                    TreeScatterStats st;
                    int t = 0, b = 0, sn = 0;
                    Trees_CountScatter(*c.field, c.field->foliageSeed, cand, 96.0f, t, b, sn, &st);
                    if (t > trees) { at = cand; trees = t; bushes = b; stones = sn; why = st; }
                }
            }
            const bool any = trees >= 0;
            if (!any) {
                c.note(std::string("no patch of ") + kSubjects[k].variant + " (skipped)");
                V[k].skip = true;
                return true;
            }
            V[k].patch = Vector3Scale(at, c.field->radius + c.field->elevation(at));
            V[k].generated = trees + bushes + stones;
            c.note(TextFormat("%s: dropped of %d candidates (3x3 tiles): odds %d, elevation %d, slope %d, buried %d, no ground %d, water %d, "
                              "material %d, underfoot %d, cramped %d; mean odds tree %.4f bush %.4f stone %.4f; plants read variant %s",
                              kSubjects[k].variant, why.candidates, why.odds, why.elevation, why.slope,
                              why.buried, why.noGround, why.water, why.material, why.underfoot, why.cramped, why.pTree / std::max(1, why.candidates),
                              why.pBush / std::max(1, why.candidates), why.pStone / std::max(1, why.candidates), BiomeVariantName(why.variant)));
            // The ground's materials over the patch (trees and bushes root only in grass or regolith).
            int mats[(int)VoxelMaterial::Count] = {};
            const Vector3 e1 = Tangent(at), e2 = Vector3CrossProduct(at, e1);
            for (int i = 0; i < 64; ++i) {
                const float a = i * 2.39996f, r = 96.0f * sqrtf((i + 0.5f) / 64.0f);
                const Vector3 d = Vector3Normalize(Vector3Add(at, Vector3Scale(Vector3Add(Vector3Scale(e1, cosf(a)), Vector3Scale(e2, sinf(a))), r / c.field->radius)));
                ++mats[(int)c.field->materialLod(Vector3Scale(d, c.field->radius + c.field->elevation(d) - 0.3f), 1.0f, 0.0f, d)];
            }
            std::string ground;
            for (int m = 0; m < (int)VoxelMaterial::Count; ++m) if (mats[m]) ground += TextFormat(" %d:%d", m, mats[m]);
            c.note(TextFormat("%s: patch at %.0f m elevation; generated within 96 m: %d trees, %d bushes, %d stones; ground (material:samples)%s",
                              kSubjects[k].variant, c.field->elevation(at), trees, bushes, stones, ground.c_str()));
            return true;
        });
        s.step("find a " + name, 90.0, [k, name](ScenarioCtx& c) {
            if (V[k].skip) return true;
            const Vector3 up = Vector3Normalize(V[k].patch), e1 = Tangent(up);
            c.camera(Vector3Add(V[k].patch, Vector3Scale(up, 60.0f)), Vector3Add(V[k].patch, Vector3Scale(e1, 40.0f)), up);
            c.waitingFor("the tree tiles around the patch: its plants within 96 m as generated, an idle world");
            if (!c.body->trees) { c.fail("no trees on the planet"); return true; }
            std::vector<Plant> plants;
            c.body->trees->plantsNear(V[k].patch, 400.0f, plants);
            int near = 0;
            for (const Plant& p : plants) near += Vector3Distance(p.base, V[k].patch) <= 96.0f;
            const bool loaded = near >= V[k].generated && c.body->trees->readyTilesNear(V[k].patch, 300.0f) >= 40 && WorldIdle();
            if (c.t < 2.0 || (!loaded && c.t < 60.0)) return false;
            // Mature trees (class 2) nearest the patch, else the largest class there is: up to 6 candidates for the view.
            static const int kPrefer[4] = { 2, 3, 1, 0 };
            std::vector<Plant> cands;
            for (int cls : kPrefer) {
                for (const Plant& p : plants)
                    if (p.tree() && SpeciesOfVariant(p.variant) == kSubjects[k].species && SizeClassOfVariant(p.variant) == cls) cands.push_back(p);
                if (!cands.empty()) break;
            }
            std::sort(cands.begin(), cands.end(), [&](const Plant& a, const Plant& b) { return Vector3Distance(a.base, V[k].patch) < Vector3Distance(b.base, V[k].patch); });
            if (cands.size() > 6) cands.resize(6);
            const bool found = !cands.empty();
            int count = 0, bySpecies[kWoodSpecies] = {};
            for (const Plant& p : plants) if (p.tree()) ++bySpecies[(int)SpeciesOfVariant(p.variant)];
            count = bySpecies[(int)kSubjects[k].species];
            std::string mix;
            for (int sp = 0; sp < kWoodSpecies; ++sp) if (bySpecies[sp]) mix += TextFormat(" %s %d", kWoodSpeciesNames[sp], bySpecies[sp]);
            c.note(TextFormat("%s: %d %s trees within 400 m (%d plants; trees:%s)", kSubjects[k].variant, count, name.c_str(), (int)plants.size(), mix.c_str()));
            if (!c.check(name + " grows in " + kSubjects[k].variant, found)) { V[k].skip = true; return true; }
            // The whole tree from 15-30 m at its mid-height, from whichever candidate and side (8 headings) has open ground
            // between (terrain doesn't cut the sight line) and the fewest other plants in the way. Those are felled: the ones
            // near the sight line and any near the eye.
            int bestCost = 1 << 30;
            float bestBack = 0.0f;
            for (const Plant& p : cands) {
                // (A tree under ~3 m, the krummholz, is framed closer: from 15 m it is a speck.)
                const float h = std::max(p.height, 1.0f) * 1.3f, back = std::clamp(6.0f + 2.0f * h, std::min(15.0f, 5.0f + 3.0f * h), 30.0f);
                const Vector3 e1 = Tangent(p.up), e2 = Vector3CrossProduct(p.up, e1);
                for (int d = 0; d < 8; ++d) {
                    const float a = d * 0.785398f;
                    const Vector3 side = Vector3Add(Vector3Scale(e1, cosf(a)), Vector3Scale(e2, sinf(a)));
                    Vector3 eye = Vector3Add(p.base, Vector3Add(Vector3Scale(side, back), Vector3Scale(p.up, 0.45f * h + 1.0f)));
                    const Vector3 ed = Vector3Normalize(eye);
                    const float ground = c.field->radius + c.field->elevation(ed) + 1.7f;
                    if (Vector3Length(eye) < ground) eye = Vector3Scale(ed, ground);
                    const Vector3 target = Vector3Add(p.base, Vector3Scale(p.up, 0.5f * h));
                    // Open: no rock (terrain, tors, overhangs) on the sight line, nor within 5 m of the eye (the
                    // camera-clip flag, #1442).
                    auto rock = [&](Vector3 q) { return c.field->density(q, nullptr) >= ISO_LEVEL; };
                    bool open = true;
                    // Five sight lines: the low trunk, mid-height, the crown, and either side of the crown.
                    const Vector3 across = Vector3Normalize(Vector3CrossProduct(p.up, Vector3Subtract(target, eye)));
                    const Vector3 aims[5] = { Vector3Add(p.base, Vector3Scale(p.up, 0.15f * h)), target, Vector3Add(p.base, Vector3Scale(p.up, 0.85f * h)),
                                              Vector3Add(target, Vector3Scale(across, 0.4f * h)), Vector3Subtract(target, Vector3Scale(across, 0.4f * h)) };
                    for (const Vector3& aim : aims)
                        for (int i = 1; i < 20 && open; ++i) open = !rock(Vector3Lerp(eye, aim, i / 20.0f));
                    const Vector3 dirs[5] = { e1, Vector3Negate(e1), e2, Vector3Negate(e2), p.up }; // (not down: the eye stands 1.7 m up)
                    for (const Vector3& dv : dirs)
                        for (float r : { 2.5f, 5.0f }) open = open && !rock(Vector3Add(eye, Vector3Scale(dv, r)));
                    int cost = open ? 0 : 1000;
                    cost += (int)((Vector3Length(eye) - ground) * 2.0f); // eye pushed up a slope: a worse angle
                    for (const Plant& q : plants) cost += InTheWay(p, eye, q) ? (q.stone() ? 50 : 1) : 0; // (stones stay: remove() fells only plants)
                    if (cost < bestCost) { bestCost = cost; bestBack = back; V[k].tree = p; V[k].eye = eye; V[k].target = target; V[k].up = p.up; }
                }
            }
            const Plant& p = V[k].tree;
            int cleared = 0;
            for (const Plant& q : plants) if (InTheWay(p, V[k].eye, q)) { c.body->trees->remove(q.id); ++cleared; }
            c.note(TextFormat("%s: %.0f m tall (class %d), viewed from %.0f m; %d plants cleared from the view%s", name.c_str(), p.height,
                              SizeClassOfVariant(p.variant), bestBack, cleared, bestCost >= 1000 ? " (terrain in the sight line: no open side)" : ""));
            return true;
        });
        s.step("look at the " + name, 60.0, [k, name](ScenarioCtx& c) {
            if (V[k].skip) return true;
            c.camera(V[k].eye, V[k].target, V[k].up);
            c.waitingFor("the view settled (an idle world)");
            if (c.t < 3.0 || (!WorldIdle() && c.t < 30.0)) return false;
            // Again from the eye: tiles that loaded only once the camera got here hold plants the first pass never saw.
            if (!cleared2[k]) {
                cleared2[k] = true;
                std::vector<Plant> near;
                c.body->trees->plantsNear(V[k].eye, 80.0f, near);
                int n = 0;
                for (const Plant& q : near) if (InTheWay(V[k].tree, V[k].eye, q)) { c.body->trees->remove(q.id); ++n; }
                c.note(TextFormat("%s: %d more plants cleared from the eye", name.c_str(), n));
                clearedAt[k] = n ? c.t : 0.0;
                return false;
            }
            if (c.t < clearedAt[k] + 3.0 || (!WorldIdle() && c.t < clearedAt[k] + 30.0)) return false;
            c.shot(name, 100.0); // a picture to judge by eye: the look is meant to change (#1995)
            return true;
        });
    }
    s.step("done", 1.0, [](ScenarioCtx& c) { c.keepCapture("the new tree species (#1995 phase 2): for a person to judge"); return true; });
}
