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
#include "animals.h"
#include "farming.h"
#include "navigation.h"
#include "suit.h"

#include <raymath.h>
#include <algorithm>
#include <cstdio>
#include <cstring>
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
    Scene* scene = nullptr;
    Vector3 home{};                                   // body-local, the start
    std::map<std::string, std::vector<Vector3>> orePoints;
    std::map<std::string, size_t> oreNext;
    std::map<std::string, int> parts;                 // machine -> part index
    int failed = 0;
    bool full = false;          // --test-progression full: on through MV and HV to titanium alloy
    void BuildFull();

    Vector3 Surface(Vector3 dir) { return Vector3Scale(dir, field->radius + field->elevation(dir)); }
    Vector3 Offset(float east, float north) {
        const Vector3 up = Vector3Normalize(home);
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, up);
        return Vector3Normalize(Vector3Add(up, Vector3Add(Vector3Scale(t1, east / field->radius), Vector3Scale(t2, north / field->radius))));
    }

    // Like a player putting building stone, logs and gems in a chest: digging and chopping now give
    // those too (M21, M22), and a full 45-slot grid refuses a craft's output.
    void Stash() {
        static const char* kStash[] = { "granite", "basalt", "sandstone", "limestone", "marble", "oak_log", "birch_log", "pine_log",
                                        "quartz", "amethyst", "ruby", "sapphire", "emerald", "grain_seeds", "root_seeds", "flax_seeds", "berries" };
        for (const char* id : kStash) { const int n = inventory.count(id); if (n > 0) inventory.remove(id, n); }
    }
    bool Craft(const std::string& item, int want) {
        Stash();
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
                    orePoints[OreKindOf(*body, p)].push_back(p);
                }
            }
        // Deposits of the rarer kinds are worth a trip: prospect further and deeper for any missing.
        for (const char* rare : { "bauxite_ore", "ilmenite_ore", "coal" }) {
            if (orePoints[rare].size() >= 40) continue;
            for (float e = -800.0f; e <= 800.0f && orePoints[rare].size() < 400; e += 8.0f)
                for (float n = -800.0f; n <= 800.0f && orePoints[rare].size() < 400; n += 8.0f) {
                    const Vector3 d = Offset(e, n);
                    const Vector3 s = Surface(d);
                    for (float depth = 10.0f; depth <= 300.0f; depth += 4.0f) {
                        const Vector3 p = Vector3Subtract(s, Vector3Scale(d, depth));
                        const float base = field->baseDensity(p);
                        if (base < ISO_LEVEL || field->material(p, base) != VoxelMaterial::Ore) continue;
                        if (strcmp(OreKindOf(*body, p), rare) == 0) orePoints[rare].push_back(p);
                    }
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
        static std::map<std::string, double> lastTake; // per machine (one shared timer starved the second of two)
        if (GetTime() - lastTake[machine] > 0.25) { lastTake[machine] = GetTime(); Factory::Extract(body->id, plot->id, parts[machine]); }
        return inventory.count(item) >= want;
    }

    void Build() {
        steps.clear();
        auto add = [](const std::string& n, std::function<bool()> f) { steps.push_back({ n, std::move(f) }); };
        // Tier 0: wood, dirt, stone; clay, a crucible, a campfire to fire it.
        add("chop trees for 160 wood", [] { while (inventory.count("wood") < 160) if (!TerrainTool_HarvestNearest(*body, home, true)) return false; return true; });
        add("dig 84 dirt", [] { return DigFor("dirt", 84, 1.0f); });
        add("harvest bushes for fibre, pick up loose stones; craft and hold a stone pickaxe", [] {
            while (inventory.count("fibre") < 2) if (!TerrainTool_HarvestNearest(*body, home, false)) return false;
            while (inventory.count("stone") < 4) if (!TerrainTool_HarvestNearest(*body, home, false, true)) return false;
            if (!Craft("stone_pickaxe", 1)) return false;
            terrainTool.heldTool = "stone_pickaxe";
            return true;
        });
        add("dig 72 stone", [] { return DigFor("stone", 72, 5.0f); });
        add("craft clay, a crucible mould and a campfire", [] { return Craft("clay", 28) && Craft("crucible_mold", 1) && Craft("campfire", 1); });
        add("place the campfire and fire the crucible", [] {
            static bool started = false;
            if (!started) { started = true; if (!Place("campfire", PartType::Campfire, { 0, 0, 0 })) return false; return false; }
            static bool fed = false;
            if (!fed) { fed = true; Give("campfire", "crucible_mold", 1); Give("campfire", "wood", 2); }
            return Collect("campfire", "fired_crucible", 1);
        });
        // Tier 1: eat; a kiln; ores; ingots; bronze.
        add("harvest bushes for berries and eat them", [] {
            while (inventory.count("berries") < 1) if (!TerrainTool_HarvestNearest(*body, home, false)) return false;
            survival.food = 90.0f; // a morning's work
            return Survival_UseItem("berries") && QuestBook_CountOf("eat") >= 1.0;
        });
        add("find ore veins", [] { FindOre(); return true; });
        add("mine 38 copper, 8 tin, 52 iron ore and 16 coal", [] { return MineOre("copper_ore", 38) && MineOre("tin_ore", 8) && MineOre("iron_ore", 52) && MineOre("coal", 16); });
        add("craft and place a kiln", [] { return Craft("kiln", 1) && Place("kiln", PartType::Kiln, { 3, 0, 0 }); });
        add("smelt copper", [] { static bool fed = false; if (!fed) { fed = true; Give("kiln", "copper_ore", 38); Give("kiln", "wood", 22); } return Collect("kiln", "copper_ingot", 38); });
        add("smelt tin", [] { static bool fed = false; if (!fed) { fed = true; Give("kiln", "tin_ore", 8); Give("kiln", "wood", 5); } return Collect("kiln", "tin_ingot", 8); });
        add("smelt iron", [] { static bool fed = false; if (!fed) { fed = true; Give("kiln", "iron_ore", 38); Give("kiln", "wood", 30); } return Collect("kiln", "iron_ingot", 38); });
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
        add("bake 8 coke", [] { static bool fed = false; if (!fed) { fed = true; Give("coke_oven", "coal", 16); } return Collect("coke_oven", "coke", 8); });
        add("make 34 steel", [] { static bool fed = false; if (!fed) { fed = true; Give("blast_furnace", "iron_ingot", 34); Give("blast_furnace", "coke", 8); } return Collect("blast_furnace", "steel_ingot", 34); });
        add("hammer 34 steel plates", [] { static bool fed = false; if (!fed) { fed = true; Give("forge_hammer", "steel_ingot", 34); Give("boiler", "wood", 50); } return Collect("forge_hammer", "steel_plate", 34); });
        add("pen two hornhides in a ring of fences (then take their hide)", [] {
            static int phase = 0;
            if (phase == 0) {
                while (inventory.count("wood") < 60) if (!TerrainTool_HarvestNearest(*body, home, true)) return false;
                if (!Craft("fence", 34)) return false;
                // An 8 m square of fences 40 m east of the workshop, following the ground.
                const Vector3 c = Surface(Offset(40.0f, 0.0f));
                Plot* pen = body->plots->plotFor(c, 0.0f, true);
                const uint8_t across = (uint8_t)RotateOrientation(0, 1); // turned fences sit a cell toward -z
                for (int a = -8; a <= 8; a += 2)
                    for (int side = 0; side < 4; ++side) {
                        const bool alongX = side < 2;
                        if (!alongX && a == 8) continue; // the rows reach the corners
                        const int cx = alongX ? a : side == 2 ? -8 : 8, cz = alongX ? (side == 0 ? -8 : 9) : a;
                        const Vector3 col = Vector3Normalize(pen->toBody({ (cx + 0.5f) * kCellSize, 0.0f, (cz + 0.5f) * kCellSize }));
                        const int g = pen->cellAt(pen->toPlot(Surface(col))).y;
                        if (!BuildTool_Place(*body, *pen, { PartType::Fence, { cx, g, cz }, alongX ? across : (uint8_t)0 }, true)) { printf("bench: progression: a fence would not go down\n"); failed = 1; return false; }
                    }
                // Two hornhides found grazing nearby, fed and led in (spawned inside: the AI is --test-animals' job).
                Animals_Spawn(*scene, *body, "hornhide", c, 2, false, 1.5f);
                phase = 1;
                return false;
            }
            if (QuestBook_CountOf("penned") < 2.0) return false; // the pen check runs every 2 s
            // Slaughtered for their hide (6) for a food pouch and the field suit.
            for (int k = 0; k < 4 && inventory.count("hide") < 4; ++k) { const int i = Animals_Nearest(Surface(Offset(40.0f, 0.0f)), 20.0f); if (i < 0) break; Animals_Hit(i, 999.0f); }
            return inventory.count("hide") >= 4;
        });
        // Tier 3: glass, vacuum tubes; a field; two suit modules.
        add("dig 24 sand", [] { return DigSand(24); });
        add("make 6 glass", [] { static bool fed = false; if (!fed) { fed = true; Give("kiln", "sand", 24); Give("kiln", "wood", 4); } return Collect("kiln", "glass", 6); });
        add("craft 4 vacuum tubes", [] { return Craft("copper_wire", 8) && Craft("vacuum_tube", 4); });
        add("till a field and plant seeds from the bushes", [] {
            while (inventory.count("grain_seeds") + inventory.count("root_seeds") + inventory.count("flax_seeds") < 1)
                if (!TerrainTool_HarvestNearest(*body, home, false)) return false;
            const Vector3 g = Surface(Offset(-20.0f, 10.0f));
            if (!Farm_Till(*scene, *body, g)) { printf("bench: progression: could not till\n"); failed = 1; return false; }
            for (auto& pl : body->plots->plots)
                for (int i = (int)pl->parts.size() - 1; i >= 0; --i)
                    if (pl->parts[i].type == PartType::Farmland && Farm_Plant(*body, *pl, i, "")) return QuestBook_CountOf("plant") >= 1.0;
            return false;
        });
        add("craft and fit two suit modules (a food pouch and an oxygen tank)", [] {
            while (inventory.count("fibre") < 40) if (!TerrainTool_HarvestNearest(*body, home, false)) return false;
            if (!Craft("food_pouch", 1) || !Craft("oxygen_tank_module", 1)) return false;
            Suit_Reset();
            return Suit_Fit("food_pouch") && Suit_Fit("oxygen_tank_module") && QuestBook_CountOf("modules") >= 2.0;
        });
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
        add("mark the workshop on the map; craft and wear a field suit", [] {
            Nav_AddWaypoint(body->name, Vector3Normalize(home), "Workshop");
            if (!Craft("field_suit", 1) || !Suit_Wear("field_suit")) return false;
            return Suit_Modules().size() == 4 && QuestBook_CountOf("waypoint") >= 1.0;
        });
        add("the quest book opens the MV tier", [] { return QuestBook_Tier() >= 5; });
        if (full) BuildFull();
    }

    // A machine run: feed it once, then collect until `want` of `item` is held.
    // Machines hold 4 stacks of 64, so it is topped up in batches the way a player would; what is
    // not meant for it (the inventory beyond `count`) is kept back.
    std::function<bool()> Run(const std::string& machine, std::vector<std::pair<std::string, int>> feed, const std::string& item, int want) {
        auto keep = std::make_shared<std::map<std::string, int>>();
        auto last = std::make_shared<double>(-1.0);
        return [=] {
            if (*last < 0.0) for (const auto& f : feed) (*keep)[f.first] = std::max(0, inventory.count(f.first) - f.second);
            if (GetTime() - *last > 0.3) {
                *last = GetTime();
                for (const auto& f : feed) {
                    const int spare = inventory.count(f.first) - (*keep)[f.first];
                    if (spare > 0) Give(machine, f.first, std::min(spare, 24));
                }
            }
            return Collect(machine, item, want);
        };
    }
    bool Cables(int z, int x0, int x1) {
        for (int x = x0; x <= x1; ++x) if (!BuildTool_Place(*body, *plot, { PartType::Cable, { x, 0, z }, 0 }, true)) { printf("bench: progression: no cable at %d,%d\n", x, z); failed = 1; return false; }
        return true;
    }

    // milestone 8: the rest of the ladder, through LV ore processing, MV and HV to titanium alloy
    // (rocketry's tier), and a miner working its deposit.
    void BuildFull() {
        auto add = [](const std::string& n, std::function<bool()> f) { steps.push_back({ n, std::move(f) }); };
        add("mine for the MV and HV ladder (iron, coal, copper) and dig stone", [] {
            return MineOre("iron_ore", 185) && MineOre("coal", 530) && MineOre("copper_ore", 160) && DigFor("stone", 40, 6.0f);
        });
        add("smelt 170 iron in the kiln", Run("kiln", { { "iron_ore", 170 }, { "coal", 34 } }, "iron_ingot", 170));
        add("smelt 145 copper in the kiln", Run("kiln", { { "copper_ore", 145 }, { "coal", 29 } }, "copper_ingot", 145));
        add("bake 172 coke", Run("coke_oven", { { "coal", 344 } }, "coke", 172));
        add("make 163 steel", Run("blast_furnace", { { "iron_ingot", 163 }, { "coke", 163 } }, "steel_ingot", 163));
        // Bauxite needs an iron pickaxe, ilmenite a steel one.
        add("craft iron and steel pickaxes; mine bauxite and ilmenite", [] {
            if (!Craft("iron_pickaxe", 1)) return false;
            terrainTool.heldTool = "iron_pickaxe";
            if (!MineOre("bauxite_ore", 14)) return false;
            if (!Craft("steel_pickaxe", 1)) return false;
            terrainTool.heldTool = "steel_pickaxe";
            return MineOre("ilmenite_ore", 4);
        });
        add("fuel the boiler again", [] { Give("boiler", "coal", 20); return true; });
        add("hammer 160 steel plates", Run("forge_hammer", { { "steel_ingot", 160 } }, "steel_plate", 160));
        add("dig 184 sand", [] { return DigSand(184); });
        add("make 46 glass", Run("kiln", { { "sand", 184 }, { "coal", 12 } }, "glass", 46));
        add("craft 40 vacuum tubes and 12 basic circuits", [] { return Craft("vacuum_tube", 40) && Craft("basic_circuit", 12); });
        // LV ore processing line (its own generator and cables).
        add("craft the LV ore line: macerator, ore washer, centrifuge, water pump, a generator", [] {
            return Craft("macerator", 1) && Craft("ore_washer", 1) && Craft("centrifuge", 1) && Craft("water_pump", 1) && Craft("generator", 1) && Craft("cable", 16);
        });
        add("place and wire the LV ore line", [] {
            if (!Place("macerator", PartType::Macerator, { 0, 0, 12 }) || !Place("ore_washer", PartType::OreWasher, { 3, 0, 12 }) || !Place("centrifuge", PartType::Centrifuge, { 6, 0, 12 }) ||
                !Place("water_pump", PartType::WaterPump, { 9, 0, 12 }) || !Place("generator2", PartType::Generator, { 12, 0, 12 })) return false;
            if (!Cables(11, 0, 13)) return false;
            Give("generator2", "coal", 12);
            return true;
        });
        add("pump 60 water", Run("water_pump", {}, "water", 60));
        add("crush 14 bauxite, 4 ilmenite and 12 copper ore", [] {
            static bool fed = false;
            if (!fed) { fed = true; Give("macerator", "bauxite_ore", 14); Give("macerator", "ilmenite_ore", 4); Give("macerator", "copper_ore", 12); }
            return Collect("macerator", "bauxite_crushed", 28) && inventory.count("ilmenite_crushed") >= 8 && inventory.count("copper_crushed") >= 24;
        });
        add("wash them (purified ore)", [] {
            static bool fed = false;
            if (!fed) { fed = true; Give("ore_washer", "bauxite_crushed", 28); Give("ore_washer", "ilmenite_crushed", 8); Give("ore_washer", "copper_crushed", 24); Give("ore_washer", "water", 60); }
            return Collect("ore_washer", "bauxite_purified", 28) && inventory.count("ilmenite_purified") >= 8 && inventory.count("copper_purified") >= 24;
        });
        add("centrifuge them (dust; gold from copper)", [] {
            static bool fed = false;
            if (!fed) { fed = true; Give("centrifuge", "bauxite_purified", 28); Give("centrifuge", "ilmenite_purified", 8); Give("centrifuge", "copper_purified", 24); }
            return Collect("centrifuge", "bauxite_dust", 28) && inventory.count("ilmenite_dust") >= 8 && inventory.count("gold_dust") >= 6;
        });
        add("refuel the first LV generator", [] { Give("generator", "coal", 8); return true; });
        add("smelt 6 gold", Run("electric_furnace", { { "gold_dust", 6 } }, "gold_ingot", 6));
        add("smelt the 24 copper dust (the processing pays: 12 ore gave 24 ingots and gold)", Run("electric_furnace", { { "copper_dust", 24 } }, "copper_ingot", 24));
        // MV.
        add("craft the MV row: hulls, an MV generator, an electrolyzer, a chemical reactor", [] {
            return Craft("mv_hull", 3) && Craft("mv_generator", 1) && Craft("electrolyzer", 1) && Craft("chemical_reactor", 1) && Craft("cable", 12);
        });
        add("place and wire the MV row (128 V)", [] {
            if (!Place("mv_generator", PartType::MvGenerator, { 0, 0, 16 }) || !Place("electrolyzer", PartType::Electrolyzer, { 3, 0, 16 }) ||
                !Place("chemical_reactor", PartType::ChemicalReactor, { 6, 0, 16 })) return false;
            if (!Cables(15, 0, 11)) return false;
            Give("mv_generator", "coal", 8);
            return true;
        });
        add("pump 20 more water", Run("water_pump", {}, "water", 20));
        add("electrolyse water into 10 hydrogen", Run("electrolyzer", { { "water", 10 } }, "hydrogen", 10));
        add("make 5 polyethylene", Run("chemical_reactor", { { "coal", 5 }, { "hydrogen", 10 } }, "polyethylene", 5));
        add("electrolyse bauxite dust into 12 aluminium (MV gate)", Run("electrolyzer", { { "bauxite_dust", 24 } }, "aluminium_ingot", 12));
        add("refuel the boiler", [] { Give("boiler", "coal", 6); return true; });
        add("hammer 11 aluminium plates", Run("forge_hammer", { { "aluminium_ingot", 11 } }, "aluminium_plate", 11));
        add("a miner on the MV row works its deposit", [] {
            static bool placed = false;
            if (!placed) {
                placed = true;
                if (!Craft("basic_circuit", 5) || !Craft("good_circuit", 2) || !Craft("miner", 1) || !Place("miner", PartType::Miner, { 9, 0, 16 })) { failed = 1; return false; }
            }
            const MachineView v = Factory::View(body->id, plot->id, parts["miner"]);
            static bool reported = false;
            if (v.valid && v.output.find(" x") != std::string::npos && !reported) { reported = true; printf("bench: progression: the miner makes %s\n", v.output.c_str()); }
            return v.valid && (v.output.find("_ore") != std::string::npos || v.output.find("coal") != std::string::npos);
        });
        // HV.
        add("craft good and advanced circuits, HV hulls, an HV generator, an arc furnace", [] {
            return Craft("good_circuit", 3) && Craft("advanced_circuit", 3) && Craft("hv_hull", 2) && Craft("hv_generator", 1) && Craft("arc_furnace", 1) && Craft("cable", 8);
        });
        add("place and wire the HV row (512 V)", [] {
            if (!Place("hv_generator", PartType::HvGenerator, { 0, 0, 20 }) || !Place("arc_furnace", PartType::ArcFurnace, { 3, 0, 20 })) return false;
            if (!Cables(19, 0, 5)) return false;
            Give("hv_generator", "coke", 3);
            return true;
        });
        add("arc furnace: 3 titanium from ilmenite and coke", Run("arc_furnace", { { "ilmenite_dust", 6 }, { "coke", 3 } }, "titanium_ingot", 3));
        add("arc furnace: titanium alloy (HV gate)", Run("arc_furnace", { { "titanium_ingot", 3 }, { "aluminium_ingot", 1 } }, "titanium_alloy_ingot", 4));
        // The other HV machines (placed without charging and seeded: a recipe check, not more grind).
        add("HV recipe check: the refinery makes rocket fuel, the precision assembler alloy plates", [] {
            static bool placed = false;
            if (!placed) {
                placed = true;
                if (!BuildTool_Place(*body, *plot, { PartType::Refinery, { 6, 0, 20 }, 0 }, false) || !BuildTool_Place(*body, *plot, { PartType::PrecisionAssembler, { 9, 0, 20 }, 0 }, false)) { failed = 1; return false; }
                parts["refinery"] = (int)plot->parts.size() - 2;
                parts["precision_assembler"] = (int)plot->parts.size() - 1;
                for (int x = 6; x <= 11; ++x) BuildTool_Place(*body, *plot, { PartType::Cable, { x, 0, 19 }, 0 }, false);
                return false;
            }
            // Feed once the new layout has reached the simulation (items for machines it does not know yet are dropped).
            static bool fed = false;
            if (!fed && GetTime() - stepStart > 0.5) {
                fed = true;
                Give("hv_generator", "coke", 2);
                Factory::Seed(body->id, plot->id, parts["refinery"], "coal", 1);
                Factory::Seed(body->id, plot->id, parts["refinery"], "hydrogen", 4);
                Factory::Seed(body->id, plot->id, parts["refinery"], "oxygen", 2);
                Give("precision_assembler", "titanium_alloy_ingot", 4);
            }
            Collect("refinery", "rocket_fuel", 4);
            return Collect("precision_assembler", "titanium_alloy_plate", 4) && inventory.count("rocket_fuel") >= 4;
        });
    }
}

int Progression_Update(Scene& world, bool fullLadder) {
    full = fullLadder;
    scene = &world;
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
            const int need = full ? 7 : 5;
            printf("bench: progression: bare hands to %s in %.0f s, quest tier %d -> %s\n", full ? "titanium alloy (rocketry)" : "a working LV machine",
                   GetTime() - testStart, QuestBook_Tier(), QuestBook_Tier() >= need ? "PASS" : "FAIL");
            return QuestBook_Tier() >= need ? 1 : -1;
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
