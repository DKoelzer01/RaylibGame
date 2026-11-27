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
            // The nearest patch to the start, else to points 4 km off each way: the first that grows trees (savanna is
            // open country with lone trees, so a patch can hold none).
            const Vector3 home = c.field->findSpawnDir(c.body->seed), h1 = Tangent(home), h2 = Vector3CrossProduct(home, h1);
            const float off = 4000.0f / c.field->radius;
            const Vector3 homes[5] = { home, Vector3Add(home, Vector3Scale(h1, off)), Vector3Subtract(home, Vector3Scale(h1, off)),
                                       Vector3Add(home, Vector3Scale(h2, off)), Vector3Subtract(home, Vector3Scale(h2, off)) };
            Vector3 at{};
            int trees = 0, bushes = 0, stones = 0;
            TreeScatterStats why;
            bool any = false;
            for (int i = 0; i < 5 && v >= 0; ++i) {
                const Vector3 cand = c.field->nearestVariantDir(Vector3Normalize(homes[i]), v);
                if (c.field->variantAt(cand, c.field->elevation(cand)) != v) continue;
                TreeScatterStats st;
                int t = 0, b = 0, sn = 0;
                Trees_CountScatter(*c.field, c.field->foliageSeed, cand, 96.0f, t, b, sn, &st);
                if (!any || (trees == 0 && t > 0)) { at = cand; trees = t; bushes = b; stones = sn; why = st; any = true; }
                if (t > 0) break;
            }
            if (!any) {
                c.note(std::string("no patch of ") + kSubjects[k].variant + " (skipped)");
                V[k].skip = true;
                return true;
            }
            V[k].patch = Vector3Scale(at, c.field->radius + c.field->elevation(at));
            V[k].generated = trees + bushes + stones;
            c.note(TextFormat("%s: dropped of %d candidates (3x3 tiles): odds %d, elevation %d, slope %d, buried %d, no ground %d, water %d, "
                              "material %d, underfoot %d, cramped %d", kSubjects[k].variant, why.candidates, why.odds, why.elevation, why.slope,
                              why.buried, why.noGround, why.water, why.material, why.underfoot, why.cramped));
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
            // A mature tree (class 2) nearest the patch, else the largest class there is.
            static const int kPrefer[4] = { 2, 3, 1, 0 };
            bool found = false;
            for (int cls : kPrefer) {
                float best = 1e9f;
                for (const Plant& p : plants)
                    if (p.tree() && SpeciesOfVariant(p.variant) == kSubjects[k].species && SizeClassOfVariant(p.variant) == cls) {
                        const float d = Vector3Distance(p.base, V[k].patch);
                        if (d < best) { best = d; V[k].tree = p; found = true; }
                    }
                if (found) break;
            }
            int count = 0, bySpecies[kWoodSpecies] = {};
            for (const Plant& p : plants) if (p.tree()) ++bySpecies[(int)SpeciesOfVariant(p.variant)];
            count = bySpecies[(int)kSubjects[k].species];
            std::string mix;
            for (int sp = 0; sp < kWoodSpecies; ++sp) if (bySpecies[sp]) mix += TextFormat(" %s %d", kWoodSpeciesNames[sp], bySpecies[sp]);
            c.note(TextFormat("%s: %d %s trees within 400 m (%d plants; trees:%s)", kSubjects[k].variant, count, name.c_str(), (int)plants.size(), mix.c_str()));
            if (!c.check(name + " grows in " + kSubjects[k].variant, found)) { V[k].skip = true; return true; }
            // Beside it at its mid-height, back far enough for the crown; lifted clear of the ground at the eye.
            const Plant& p = V[k].tree;
            const float h = std::max(p.height, 2.0f) * 1.3f, back = 6.0f + 1.6f * h;
            const Vector3 side = Tangent(p.up);
            Vector3 eye = Vector3Add(p.base, Vector3Add(Vector3Scale(side, back), Vector3Scale(p.up, 0.45f * h + 1.0f)));
            const Vector3 ed = Vector3Normalize(eye);
            const float ground = c.field->radius + c.field->elevation(ed) + 2.0f;
            if (Vector3Length(eye) < ground) eye = Vector3Scale(ed, ground);
            V[k].eye = eye; V[k].target = Vector3Add(p.base, Vector3Scale(p.up, 0.5f * h)); V[k].up = p.up;
            // A clear line of sight: fell every other plant whose crown could stand between the eye and the subject.
            int cleared = 0;
            const Vector3 seg = Vector3Subtract(p.base, eye);
            const float segLen2 = std::max(Vector3LengthSqr(seg), 1e-3f);
            for (const Plant& q : plants) {
                if (q.id == p.id) continue;
                const float t = std::clamp(Vector3DotProduct(Vector3Subtract(q.base, eye), seg) / segLen2, 0.0f, 1.0f);
                const Vector3 onSeg = Vector3Add(eye, Vector3Scale(seg, t));
                Vector3 off = Vector3Subtract(q.base, onSeg);
                off = Vector3Subtract(off, Vector3Scale(p.up, Vector3DotProduct(off, p.up))); // across the ground only
                if (Vector3Length(off) < 3.0f + 4.0f * q.scale) { c.body->trees->remove(q.id); ++cleared; }
            }
            c.note(TextFormat("%s: %.0f m tall (class %d), viewed from %.0f m; %d plants cleared from the view", name.c_str(), p.height, SizeClassOfVariant(p.variant), back, cleared));
            return true;
        });
        s.step("look at the " + name, 60.0, [k, name](ScenarioCtx& c) {
            if (V[k].skip) return true;
            c.camera(V[k].eye, V[k].target, V[k].up);
            c.waitingFor("the view settled (an idle world)");
            if (c.t < 3.0 || (!WorldIdle() && c.t < 30.0)) return false;
            c.shot(name, 100.0); // a picture to judge by eye: the look is meant to change (#1995)
            return true;
        });
    }
    s.step("done", 1.0, [](ScenarioCtx& c) { c.keepCapture("the new tree species (#1995 phase 2): for a person to judge"); return true; });
}
