// --test-progression (milestone 5 gate): a scripted playthrough from an empty inventory to a placed,
// powered, working LV machine, through the real systems: chopping, digging the right ores, hand
// crafting (inventory recipes), placing machines with the build tool (each costs its crafted item),
// wiring cables, feeding machines and taking their output (factory commands). The factory clock runs
// fast so the ~20 minutes of machine time take seconds.

#include "progression_test.h"
#include "buildtool.h"
#include "player.h"
#include "survival.h"
#include "terraintool.h"
#include "inventory.h"
#include "factory.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <algorithm>
#include <cstdio>
#include <functional>
#include <map>
#include <random>
#include <string>
#include <vector>

namespace {
    struct Step { std::string name; std::function<bool()> run; };
    std::vector<Step> steps;
    size_t current = 0;
    double stepStart = 0.0, testStart = 0.0;
    Planetoid* body = nullptr;
    const PlanetField* field = nullptr;
    Plot* plot = nullptr;
    Vector3 home{};                                   // body-local, the start
    std::map<std::string, std::vector<Vector3>> orePoints;
    std::map<std::string, size_t> oreNext;
    std::map<std::string, int> parts;                 // machine -> part index
    int failed = 0;

    Vector3 Surface(Vector3 dir) { return Vector3Scale(dir, field->radius + field->elevation(dir)); }
    Vector3 Offset(float east, float north) {
        const Vector3 up = Vector3Normalize(home);
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, up);
        return Vector3Normalize(Vector3Add(up, Vector3Add(Vector3Scale(t1, east / field->radius), Vector3Scale(t2, north / field->radius))));
    }

    bool Craft(const std::string& item, int want) {
        const Recipe* r = nullptr;
        for (const Recipe& rc : itemDB.recipes()) if (rc.output == item) { r = &rc; break; }
        if (!r) { printf("bench: progression: no recipe for %s\n", item.c_str()); failed = 1; return false; }
        // Like a player: craft missing hand-craftable ingredients first (wire for cables, ...).
        while (inventory.count(item) < want) {
            for (const auto& in : r->inputs)
                if (inventory.count(in.first) < in.second)
                    for (const Recipe& sub : itemDB.recipes())
                        if (sub.output == in.first) { Craft(in.first, in.second); break; }
            if (!inventory.craft(*r)) {
                printf("bench: progression: cannot craft %s:", item.c_str());
                for (const auto& in : r->inputs) printf(" %s %d/%d", in.first.c_str(), inventory.count(in.first), in.second);
                printf("\n");
                failed = 1;
                return false;
            }
        }
        return true;
    }
    // Dig points around the start at a depth until `want` of an item are held.
    bool DigFor(const std::string& item, int want, float depth) {
        static std::mt19937 rng(5);
        std::uniform_real_distribution<float> u(-40.0f, 40.0f);
        for (int tries = 0; tries < 400 && inventory.count(item) < want; ++tries) {
            const Vector3 d = Offset(u(rng), u(rng));
            const Vector3 p = Vector3Subtract(Surface(d), Vector3Scale(d, depth));
            TerrainTool_DigAt(*body, p);
        }
        return inventory.count(item) >= want;
    }
    // Sand: beaches somewhere within a few hundred metres.
    bool DigSand(int want) {
        for (float r = 20.0f; r < 600.0f && inventory.count("sand") < want; r += 10.0f)
            for (int a = 0; a < 48 && inventory.count("sand") < want; ++a) {
                const float ang = a * 2.0f * PI / 48.0f;
                const Vector3 d = Offset(cosf(ang) * r, sinf(ang) * r);
                const Vector3 p = Vector3Subtract(Surface(d), Vector3Scale(d, 0.8f));
                const float base = field->baseDensity(p);
                if (base >= ISO_LEVEL && field->material(p, base) == VoxelMaterial::Sand)
                    for (int k = 0; k < 4 && inventory.count("sand") < want; ++k) TerrainTool_DigAt(*body, Vector3Subtract(p, Vector3Scale(d, 0.6f * k)));
            }
        return inventory.count("sand") >= want;
    }
    // Ore: find vein samples of each kind under the area once, then dig them.
    void FindOre() {
        for (float e = -160.0f; e <= 160.0f; e += 4.0f)
            for (float n = -160.0f; n <= 160.0f; n += 4.0f) {
                const Vector3 d = Offset(e, n);
                const Vector3 s = Surface(d);
                for (float depth = 10.0f; depth <= 90.0f; depth += 4.0f) {
                    const Vector3 p = Vector3Subtract(s, Vector3Scale(d, depth));
                    const float base = field->baseDensity(p);
                    if (base < ISO_LEVEL) continue;
                    if (field->material(p, base) != VoxelMaterial::Ore) continue;
                    orePoints[OreKindAt(p, body->seed)].push_back(p);
                }
            }
        for (auto& kv : orePoints) printf("bench: progression: %zu %s samples under the area\n", kv.second.size(), kv.first.c_str());
    }
    bool MineOre(const std::string& kind, int want) {
        auto& pts = orePoints[kind];
        size_t& i = oreNext[kind];
        while (inventory.count(kind) < want && i < pts.size()) TerrainTool_DigAt(*body, pts[i++]);
        return inventory.count(kind) >= want;
    }

    bool Place(const std::string& name, PartType t, Int3 cell) {
        const int before = plot->partCount();
        if (!BuildTool_Place(*body, *plot, { t, cell, 0 }, true)) { printf("bench: progression: could not place a %s\n", name.c_str()); failed = 1; return false; }
        parts[name] = (int)plot->parts.size() - 1;
        return plot->partCount() == before + 1;
    }
    void Give(const std::string& machine, const std::string& item, int count) {
        count = std::min(count, inventory.count(item));
        if (count > 0 && inventory.remove(item, count)) Factory::Insert(body->id, plot->id, parts[machine], item, count);
    }
    // Run a machine until `item` reaches `want` in the inventory (taking its output now and then).
    bool Collect(const std::string& machine, const std::string& item, int want) {
        static double lastTake = 0.0;
        if (GetTime() - lastTake > 0.25) { lastTake = GetTime(); Factory::Extract(body->id, plot->id, parts[machine]); }
        return inventory.count(item) >= want;
    }

    void Build() {
        steps.clear();
        auto add = [](const std::string& n, std::function<bool()> f) { steps.push_back({ n, std::move(f) }); };
        // Tier 0: wood, dirt, stone; clay, a crucible, a campfire to fire it.
        add("chop trees for 160 wood", [] { while (inventory.count("wood") < 160) if (!TerrainTool_HarvestNearest(*body, home, true)) return false; return true; });
        add("dig 84 dirt", [] { return DigFor("regolith", 84, 1.0f); });
        add("dig 72 stone", [] { return DigFor("stone", 72, 5.0f); });
        add("craft clay, a crucible mould and a campfire", [] { return Craft("clay", 28) && Craft("crucible_mold", 1) && Craft("campfire", 1); });
        add("place the campfire and fire the crucible", [] {
            static bool started = false;
            if (!started) { started = true; if (!Place("campfire", PartType::Campfire, { 0, 0, 0 })) return false; return false; }
            static bool fed = false;
            if (!fed) { fed = true; Give("campfire", "crucible_mold", 1); Give("campfire", "wood", 2); }
            return Collect("campfire", "fired_crucible", 1);
        });
        // Tier 1: a kiln; ores; ingots; bronze.
        add("find ore veins", [] { FindOre(); return true; });
        add("mine 38 copper, 8 tin, 32 iron ore and 12 coal", [] { return MineOre("copper_ore", 38) && MineOre("tin_ore", 8) && MineOre("iron_ore", 32) && MineOre("coal", 12); });
        add("craft and place a kiln", [] { return Craft("kiln", 1) && Place("kiln", PartType::Kiln, { 3, 0, 0 }); });
        add("smelt copper", [] { static bool fed = false; if (!fed) { fed = true; Give("kiln", "copper_ore", 38); Give("kiln", "wood", 22); } return Collect("kiln", "copper_ingot", 38); });
        add("smelt tin", [] { static bool fed = false; if (!fed) { fed = true; Give("kiln", "tin_ore", 8); Give("kiln", "wood", 5); } return Collect("kiln", "tin_ingot", 8); });
        add("smelt iron", [] { static bool fed = false; if (!fed) { fed = true; Give("kiln", "iron_ore", 30); Give("kiln", "wood", 24); } return Collect("kiln", "iron_ingot", 30); });
        add("alloy 28 bronze", [] {
            static bool fed = false;
            if (!fed) { fed = true; Give("kiln", "copper_ingot", 21); Give("kiln", "tin_ingot", 7); Give("kiln", "wood", 8); }
            return Collect("kiln", "bronze_ingot", 28);
        });
        // Tier 2: steam (boiler + cables + forge hammer), coke, steel.
        add("craft a boiler, a forge hammer and cables", [] {
            return Craft("bronze_plate", 10) && Craft("boiler", 1) && Craft("forge_hammer", 1) && Craft("copper_wire", 4) && Craft("cable", 4);
        });
        add("place and wire the steam machines", [] {
            if (!Place("boiler", PartType::Boiler, { 0, 0, 4 }) || !Place("forge_hammer", PartType::ForgeHammer, { 3, 0, 4 })) return false;
            for (int x = 0; x < 4; ++x) if (!BuildTool_Place(*body, *plot, { PartType::Cable, { x, 0, 3 }, 0 }, true)) { failed = 1; return false; }
            Give("boiler", "wood", 20);
            return true;
        });
        add("hammer 8 bronze plates", [] { static bool fed = false; if (!fed) { fed = true; Give("forge_hammer", "bronze_ingot", 8); } return Collect("forge_hammer", "bronze_plate", 8); });
        add("craft and place a coke oven and a blast furnace", [] {
            return Craft("coke_oven", 1) && Craft("blast_furnace", 1) && Place("coke_oven", PartType::CokeOven, { 9, 0, 0 }) && Place("blast_furnace", PartType::BlastFurnace, { 12, 0, 0 });
        });
        add("bake 6 coke", [] { static bool fed = false; if (!fed) { fed = true; Give("coke_oven", "coal", 12); } return Collect("coke_oven", "coke", 6); });
        add("make 26 steel", [] { static bool fed = false; if (!fed) { fed = true; Give("blast_furnace", "iron_ingot", 26); Give("blast_furnace", "coke", 6); } return Collect("blast_furnace", "steel_ingot", 26); });
        add("hammer 26 steel plates", [] { static bool fed = false; if (!fed) { fed = true; Give("forge_hammer", "steel_ingot", 26); Give("boiler", "wood", 40); } return Collect("forge_hammer", "steel_plate", 26); });
        // Tier 3: glass, vacuum tubes.
        add("dig 16 sand", [] { return DigSand(16); });
        add("make 4 glass", [] { static bool fed = false; if (!fed) { fed = true; Give("kiln", "sand", 16); Give("kiln", "wood", 3); } return Collect("kiln", "glass", 4); });
        add("craft 4 vacuum tubes", [] { return Craft("copper_wire", 8) && Craft("vacuum_tube", 4); });
        // Tier 4: circuits, hulls, an LV generator and the first LV machine, powered and working.
        add("craft basic circuits, hulls, a generator and an electric furnace", [] {
            return Craft("copper_wire", 8) && Craft("basic_circuit", 2) && Craft("lv_hull", 2) && Craft("copper_wire", 8) && Craft("generator", 1) && Craft("electric_furnace", 1) && Craft("cable", 4);
        });
        add("place and wire the LV machines", [] {
            if (!Place("generator", PartType::Generator, { 0, 0, 8 }) || !Place("electric_furnace", PartType::Furnace, { 3, 0, 8 })) return false;
            for (int x = 0; x < 4; ++x) if (!BuildTool_Place(*body, *plot, { PartType::Cable, { x, 0, 7 }, 0 }, true)) { failed = 1; return false; }
            Give("generator", "wood", 4);
            Give("electric_furnace", "iron_ore", 2);
            return true;
        });
        add("the electric furnace runs on LV power", [] {
            const MachineView v = Factory::View(body->id, plot->id, parts["electric_furnace"]);
            return v.valid && v.status == MachineStatus::Working;
        });
    }
}

int Progression_Update(Scene& world) {
    if (steps.empty()) {
        body = dynamic_cast<Planetoid*>(world.objects[0].get());
        field = body ? dynamic_cast<const PlanetField*>(body->field.get()) : nullptr;
        if (!field || !body->trees) { printf("bench: --test-progression needs the planet\n"); return -1; }
        inventory.clear();
        QuestBook_SetTier(0);
        Factory::SetSpeed(60);
        home = body->toLocal((player.position - world.renderOrigin).toFloat());
        if (!body->plots) body->plots = std::make_unique<PlotSet>();
        // The workshop: a plot on the ground a few metres ahead.
        const Vector3 up = Vector3Normalize(home);
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 }));
        const Vector3 d = Vector3Normalize(Vector3Add(up, Vector3Scale(t1, 8.0f / field->radius)));
        plot = body->plots->plotFor(Vector3Scale(d, field->radius + field->elevation(d)), 0.0f, true);
        Build();
        current = 0;
        testStart = stepStart = GetTime();
        printf("bench: progression: %zu steps, inventory empty, quest tier 0\n", steps.size());
    }
    if (failed) return -1;
    if (current >= steps.size()) return 1;
    if (steps[current].run()) {
        printf("bench: progression: [%5.1f s] %s (quest tier %d)\n", GetTime() - testStart, steps[current].name.c_str(), QuestBook_Tier());
        ++current;
        stepStart = GetTime();
        if (current >= steps.size()) {
            Factory::SetSpeed(1);
            printf("bench: progression: bare hands to a working LV machine in %.0f s, quest tier %d -> %s\n", GetTime() - testStart, QuestBook_Tier(), QuestBook_Tier() >= 5 ? "PASS" : "FAIL");
            return QuestBook_Tier() >= 5 ? 1 : -1;
        }
    } else if (GetTime() - stepStart > 120.0 || failed) {
        printf("bench: progression: FAIL at \"%s\" (wood %d, copper %d/%d ore/ingot, tin %d/%d, iron %d/%d, coal %d, bronze %d, steel %d/%d)\n", steps[current].name.c_str(),
               inventory.count("wood"), inventory.count("copper_ore"), inventory.count("copper_ingot"), inventory.count("tin_ore"), inventory.count("tin_ingot"),
               inventory.count("iron_ore"), inventory.count("iron_ingot"), inventory.count("coal"), inventory.count("bronze_ingot"), inventory.count("steel_ingot"), inventory.count("steel_plate"));
        Factory::SetSpeed(1);
        return -1;
    }
    return 0;
}
