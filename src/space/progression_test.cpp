// --test-progression (milestone 5 gate): a scripted playthrough from an empty inventory to a placed,
// powered, working LV machine, through the real systems: chopping, digging the right ores, hand
// crafting (inventory recipes), placing machines with the build tool (each costs its crafted item),
// wiring cables, feeding machines and taking their output (factory commands). The factory clock runs
// fast so the ~20 minutes of machine time take seconds.

#include "progression_test.h"
#include "fluid.h"
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
#include "tools.h"
#include "navigation.h"
#include "suit.h"
#include "rotationplot.h"
#include "rotationsources.h"
#include "researchgate.h"
#include "research.h"
#include "researchtab.h"

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
    std::vector<std::string> stepHistory; // the current step, a line a second (#372: printed if it fails)
    struct PendingGive { std::string machine, item; int count; };
    std::vector<PendingGive> pendingGives; // deliveries waiting for the factory to take a machine in (#372)
    int frameNo = 0, placedFrame = -100; // test frames; the last Place()
    void FlushGives();
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
    bool ev = false;            // --test-progression ev: tiers 9-10 (M24)
    void BuildFull();

    Vector3 Surface(Vector3 dir) { return Vector3Scale(dir, field->radius + field->elevation(dir)); }
    Vector3 Offset(float east, float north) {
        const Vector3 up = Vector3Normalize(home);
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, up);
        return Vector3Normalize(Vector3Add(up, Vector3Add(Vector3Scale(t1, east / field->radius), Vector3Scale(t2, north / field->radius))));
    }

    // Like a player putting building stone, logs and gems in a chest: digging and chopping now give
    // those too (M21, M22), and a full grid refuses a craft's output.
    void Stash() {
        static const char* kStash[] = { "granite", "basalt", "sandstone", "limestone", "marble",
                                        "quartz", "amethyst", "ruby", "sapphire", "emerald", "grain_seeds", "root_seeds", "flax_seeds", "berries" };
        for (const char* id : kStash) { const int n = inventory.count(id); if (n > 0) inventory.remove(id, n); }
    }
    // Saw one log (of whatever wood is carried) into planks; false with no logs.
    bool SawLog() {
        for (const Recipe& rc : itemDB.recipes())
            if (rc.inputs.size() == 1 && ItemMatches("any_log", rc.inputs[0].first) && inventory.count(rc.inputs[0].first) > 0) return inventory.craft(rc);
        return false;
    }
    bool Craft(const std::string& item, int want) {
        Stash();
        const Recipe* r = nullptr;
        for (const Recipe& rc : itemDB.recipes()) if (rc.output == item) { r = &rc; break; }
        if (!r) { printf("bench: progression: no recipe for %s\n", item.c_str()); failed = 1; return false; }
        // Like a player: craft missing hand-craftable ingredients first (wire for cables, ...).
        while (inventory.count(item) < want) {
            for (const auto& in : r->inputs) {
                // Wildcard planks: saw whichever logs are carried.
                while (in.first == "any_planks" && inventory.count(in.first) < in.second && SawLog()) {}
                if (inventory.count(in.first) < in.second)
                    for (const Recipe& sub : itemDB.recipes())
                        if (sub.output == in.first) { Craft(in.first, in.second); break; }
            }
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
    // Clay (#206): a terrain material now (river and pond beds, swamps, rare pockets in the dirt): search outward
    // from the start through the soil layer and dig what turns up. Prints how far the nearest clay was.
    bool DigClay(int want) {
        float first = -1.0f;
        for (float r = 0.0f; r < 600.0f && inventory.count("clay") < want; r += 4.0f)
            for (int a = 0; a < (r < 1.0f ? 1 : 64) && inventory.count("clay") < want; ++a) {
                const float ang = a * 2.0f * PI / 64.0f;
                const Vector3 d = Offset(cosf(ang) * r, sinf(ang) * r);
                for (float depth = 0.8f; depth <= 2.2f && inventory.count("clay") < want; depth += 0.7f) {
                    const Vector3 p = Vector3Subtract(Surface(d), Vector3Scale(d, depth));
                    const float base = field->baseDensity(p);
                    if (base < ISO_LEVEL || field->material(p, base) != VoxelMaterial::Clay) continue;
                    if (first < 0.0f) first = r;
                    for (int k = 0; k < 4 && inventory.count("clay") < want; ++k) TerrainTool_DigAt(*body, Vector3Subtract(p, Vector3Scale(d, 0.5f * k)));
                }
            }
        printf("bench: progression: nearest clay %.0f m from the start; %d clay\n", first, inventory.count("clay"));
        return inventory.count("clay") >= want;
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
        // (#532/#593: every kind the ladder mines - a start can sit where the near ground has no copper or iron)
        for (const char* rare : { "bauxite_ore", "ilmenite_ore", "coal", "iron_ore", "copper_ore", "tin_ore" }) {
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
    // The pickaxe the script last held: when it wears out, a player crafts another (#587's harder strata wear the stone
    // pickaxe out before the MV ladder; without one the digs failed silently and used up the ore samples).
    std::string lastPickaxe;
    bool KeepPickaxe() {
        if (!terrainTool.heldTool.empty() && inventory.count(terrainTool.heldTool) > 0) { lastPickaxe = terrainTool.heldTool; return true; }
        if (lastPickaxe.empty() || !Craft(lastPickaxe, 1)) return false;
        printf("bench: progression: crafted a new %s (the last one wore out)\n", lastPickaxe.c_str());
        terrainTool.heldTool = lastPickaxe;
        return true;
    }
    bool MineOre(const std::string& kind, int want) {
        auto& pts = orePoints[kind];
        size_t& i = oreNext[kind];
        while (inventory.count(kind) < want && i < pts.size()) {
            if (!KeepPickaxe()) break;
            TerrainTool_DigAt(*body, pts[i++]);
        }
        return inventory.count(kind) >= want;
    }

    bool Place(const std::string& name, PartType t, Int3 cell) {
        const int before = plot->partCount();
        if (!BuildTool_Place(*body, *plot, { t, cell, 0 }, true)) { printf("bench: progression: could not place a %s\n", name.c_str()); failed = 1; return false; }
        parts[name] = (int)plot->parts.size() - 1;
        placedFrame = frameNo; // gives wait a few frames: the layout goes to the factory with the next Sync (#372)
        return plot->partCount() == before + 1;
    }
    // (#1271 L5) A turned part (a shaft, a windmill, a press about x) whose footprint starts at `cell`.
    bool PlaceO(const std::string& name, PartType t, Int3 cell, Int3 axis) {
        if (!BuildTool_Place(*body, *plot, RotationPlot_PartAt(t, cell, RotationPlot_OrientFor(axis)), true)) { printf("bench: progression: could not place a %s\n", name.c_str()); failed = 1; return false; }
        parts[name] = (int)plot->parts.size() - 1;
        placedFrame = frameNo;
        return true;
    }
    void Give(const std::string& machine, const std::string& item, int count) {
        if (item.rfind("any_", 0) == 0) { // a wildcard: the real items carried, the most plentiful first
            std::vector<std::pair<int, std::string>> have;  // (fuel slots hold one kind: the biggest pile fuels best)
            for (const auto& kv : inventory.all()) if (ItemMatches(item, kv.first) && kv.second > 0) have.push_back({ kv.second, kv.first });
            std::sort(have.rbegin(), have.rend());
            // What its fuel slot already burns goes first (it takes one kind at a time; a player sees the slot).
            Factory::View(body->id, plot->id, parts[machine]);
            Factory::WaitTicks(2);
            const MachineView v = Factory::View(body->id, plot->id, parts[machine]);
            if (v.hasFuelSlot && v.fuelSlot.second > 0)
                std::stable_partition(have.begin(), have.end(), [&](const std::pair<int, std::string>& h) { return h.second == v.fuelSlot.first; });
            for (const auto& h : have)
                if (count > 0) { const int k = std::min(count, h.first); Give(machine, h.second, k); count -= k; }
            return;
        }
        count = std::min(count, inventory.count(item));
        if (count <= 0) return;
        // The factory takes a just-placed machine's layout on a later frame: an insert for a part it doesn't know is
        // refunded to the inventory, and the step then waited on a machine with nothing in it (#372; found by the
        // lighting session). The layout goes over with the next frame's Factory::Sync: deliveries within 3 frames of a Place() wait.
        if (frameNo - placedFrame < 3) { pendingGives.push_back({ machine, item, count }); return; }
        if (inventory.remove(item, count)) Factory::Insert(body->id, plot->id, parts[machine], item, count);
    }
    void FlushGives() {
        if (pendingGives.empty() || frameNo - placedFrame < 3) return;
        Factory::WaitTicks(2); // the factory has taken the layout in
        std::vector<PendingGive> again;
        for (const PendingGive& g : pendingGives) {
            if (frameNo - placedFrame < 3) { again.push_back(g); continue; }
            const int c = std::min(g.count, inventory.count(g.item));
            if (c > 0 && inventory.remove(g.item, c)) Factory::Insert(body->id, plot->id, parts[g.machine], g.item, c);
        }
        pendingGives.swap(again);
    }
    // Run a machine until `item` reaches `want` in the inventory (taking its output now and then).
    bool Collect(const std::string& machine, const std::string& item, int want) {
        static std::map<std::string, double> lastTake; // per machine (one shared timer starved the second of two)
        if (GetTime() - lastTake[machine] > 0.25) { lastTake[machine] = GetTime(); Factory::Extract(body->id, plot->id, parts[machine]); }
        // The boiler needs water (#229): a bucket whenever its tank has room for one, as a player keeping it filled
        // from the river would. Not a bucket a real second (#372): the factory runs at 60x here (~1200 ticks a second at
        // 1 L a tick of burning), so the tank ran dry now and then; a dry boiler cools as fast as it warms, and each
        // refill burnt logs warming it again, until the 7 logs were gone with the hammer unpowered.
        static double lastWater = 0.0;
        if (parts.count("boiler") && GetTime() - lastWater > 0.05) {
            const MachineView b = Factory::View(body->id, plot->id, parts["boiler"]); // watched: fills from the next tick
            if (!b.valid || b.water <= b.waterMax - 1000) {
                lastWater = GetTime();
                if (inventory.count("water_bucket") == 0) inventory.add("water_bucket", 1);
                if (inventory.remove("water_bucket", 1)) Factory::Insert(body->id, plot->id, parts["boiler"], "water_bucket", 1);
                if (inventory.count("bucket") > 0) inventory.remove("bucket", inventory.count("bucket")); // the empties (refilled at the river)
            }
        }
        // A second-by-second record of the machine (and the boiler, the generators), printed if the step fails (#372).
        static double lastNote = 0.0;
        if (GetTime() - lastNote > 1.0) {
            lastNote = GetTime();
            std::string line = TextFormat("t=%3.0f s", GetTime() - stepStart);
            for (const std::string& name : { machine, std::string("boiler"), std::string("generator"), std::string("generator2") }) {
                if (!parts.count(name) || (name != machine && line.find(" " + name + ":") != std::string::npos)) continue;
                const MachineView v = Factory::View(body->id, plot->id, parts[name]);
                if (!v.valid) continue;
                line += TextFormat(" | %s: %s %d/%d en %d", name.c_str(), MachineStatusName(v.status), v.progress, v.duration, v.energy);
                if (v.waterMax > 0) line += TextFormat(" water %d heat %d", v.water, v.heat);
                if (v.hasFuelSlot) line += TextFormat(" fuel x%d", v.fuelSlot.second);
            }
            stepHistory.push_back(line);
            if (stepHistory.size() > 40) stepHistory.erase(stepHistory.begin());
        }
        return inventory.count(item) >= want;
    }

    void Build() {
        steps.clear();
        auto add = [](const std::string& n, std::function<bool()> f) { steps.push_back({ n, std::move(f) }); };
        // Tier 0: wood, dirt, stone; clay, a crucible, a campfire to fire it.
        // (#1271 L6b-3) the fibre for every rope first: once 340 logs are felled the bushes near home ran out
        add("harvest bushes for 70 fibre", [] { while (inventory.count("fibre") < 70) if (!TerrainTool_HarvestNearest(*body, home, false)) { printf("bench: progression: no bush left near home (fibre %d)\n", inventory.count("fibre")); return false; } return true; });
        add("chop trees for 340 logs", [] { while (inventory.count("any_log") < 340) if (!TerrainTool_HarvestNearest(*body, home, true)) return false; return true; });
        add("dig 84 dirt", [] { return DigFor("dirt", 84, 1.0f); });
        add("harvest bushes for fibre, pick up loose stones; craft and hold a stone pickaxe", [] {
            while (inventory.count("fibre") < 2) if (!TerrainTool_HarvestNearest(*body, home, false)) return false;
            while (inventory.count("stone") < 4) if (!TerrainTool_HarvestNearest(*body, home, false, true)) return false;
            if (!Craft("stone_pickaxe", 1)) return false;
            terrainTool.heldTool = "stone_pickaxe";
            return true;
        });
        add("dig 90 stone (any rock, #587)", [] { return DigFor("any_stone", 90, 5.0f); });
        add("dig 44 clay (#206: found, not crafted); craft a crucible mould and a campfire", [] { return DigClay(44) && Craft("crucible_mold", 1) && Craft("campfire", 1); });
        add("place the campfire and fire the crucible", [] {
            static bool started = false;
            if (!started) { started = true; if (!Place("campfire", PartType::Campfire, { 0, 0, 0 })) return false; return false; }
            static bool fed = false;
            if (!fed) { fed = true; Give("campfire", "crucible_mold", 1); Give("campfire", "any_log", 1); }
            return Collect("campfire", "fired_crucible", 1);
        });
        // Tier 1: eat; a kiln; ores; ingots; bronze.
        add("harvest bushes for berries and eat them", [] {
            while (inventory.count("berries") < 1) if (!TerrainTool_HarvestNearest(*body, home, false)) return false;
            survival.food = 90.0f; // a morning's work
            return Survival_UseItem("berries") && QuestBook_CountOf("eat") >= 1.0;
        });
        add("find ore veins", [] { FindOre(); return true; });
        add("mine 80 copper, 12 tin, 200 iron ore and 140 coal", [] { return MineOre("copper_ore", 80) && MineOre("tin_ore", 12) && MineOre("iron_ore", 200) && MineOre("coal", 140); });
        add("craft and place a kiln", [] { return Craft("kiln", 1) && Place("kiln", PartType::Kiln, { 3, 0, 0 }); });
        add("smelt copper", [] { static bool fed = false; if (!fed) { fed = true; Give("kiln", "copper_ore", 80); Give("kiln", "any_log", 16); } return Collect("kiln", "copper_ingot", 80); });
        add("smelt tin", [] { static bool fed = false; if (!fed) { fed = true; Give("kiln", "tin_ore", 12); Give("kiln", "any_log", 3); } return Collect("kiln", "tin_ingot", 12); });
        add("alloy 44 bronze", [] {
            static bool fed = false;
            if (!fed) { fed = true; Give("kiln", "copper_ingot", 33); Give("kiln", "tin_ingot", 11); Give("kiln", "any_log", 4); }
            return Collect("kiln", "bronze_ingot", 44);
        });
        // Tier 1 (#1271): a fired mould, a cast head, the anvil, a bronze pickaxe.
        add("mould clay and fire it in the kiln", [] {
            static bool fed = false;
            if (!fed) { fed = true; if (!Craft("clay_mould", 1)) return false; Give("kiln", "clay_mould", 1); Give("kiln", "any_log", 1); }
            return Collect("kiln", "fired_mould", 1);
        });
        add("cast a rough pick head in the kiln", [] {
            static bool fed = false;
            if (!fed) { fed = true; Give("kiln", "fired_mould", 1); Give("kiln", "bronze_ingot", 2); Give("kiln", "any_log", 1); }
            return Collect("kiln", "rough_bronze_pick_head", 1);
        });
        add("craft and place an anvil", [] { return Craft("anvil", 1) && Place("anvil", PartType::Anvil, { 0, 0, 12 }); });
        add("hammer the head at the anvil; haft a bronze pickaxe", [] {
            static bool fed = false;
            if (!fed) { fed = true; Give("anvil", "rough_bronze_pick_head", 1); }
            Factory::Strike(body->id, plot->id, parts["anvil"]); // (a player holding RMB on it)
            return Collect("anvil", "bronze_pick_head", 1) && Craft("bronze_pickaxe", 1);
        });
        // Tier 2 (#1271): charcoal, a bloomery, ingots at the anvil; a windmill turns a press.
        add("burn 198 logs into charcoal in the campfire (coal under them)", [] {
            static bool fed = false;
            if (!fed) { fed = true; Give("campfire", "coal", 30); Give("campfire", "any_log", 198); }
            return Collect("campfire", "charcoal", 198);
        });
        add("craft and place a bloomery", [] { return Craft("bloomery", 1) && Place("bloomery", PartType::Bloomery, { 3, 0, 12 }); });
        add("smelt 198 iron blooms", [] { static bool fed = false; if (!fed) { fed = true; Give("bloomery", "iron_ore", 198); Give("bloomery", "charcoal", 198); } return Collect("bloomery", "iron_bloom", 198); }); // 2 ore stay for the electric furnace
        add("hammer 16 blooms into iron at the anvil", [] {
            static bool fed = false;
            if (!fed) { fed = true; Give("anvil", "iron_bloom", 16); }
            Factory::Strike(body->id, plot->id, parts["anvil"]);
            return Collect("anvil", "iron_ingot", 16);
        });
        add("craft a press, a windmill and shafts; set them up in the wind", [] {
            while (inventory.count("fibre") < 24) if (!TerrainTool_HarvestNearest(*body, home, false)) { printf("bench: progression: no bush left near home (fibre %d)\n", inventory.count("fibre")); return false; }
            if (!Craft("rope", 8) || !Craft("mechanical_press", 1) || !Craft("windmill", 1) || !Craft("gearbox", 1) || !Craft("shaft", 2)) return false; // (shafts last: the press, windmill and cogs take theirs)
            RotationSources_SetWindOverride(true, 10.0f); // (a good wind: the weather's own varies)
            return Place("press", PartType::Press, { 6, 0, 12 }) && PlaceO("shaft", PartType::Shaft, { 8, 0, 12 }, { 1, 0, 0 })
                && Place("gearbox", PartType::Gearbox, { 9, 0, 12 }) && PlaceO("windmill", PartType::Windmill, { 10, 0, 12 }, { 1, 0, 0 });
        });
        add("press 182 blooms into iron with the windmill (its 5 minutes of tier 2)", [] {
            static bool fed = false;
            if (!fed) { fed = true; Give("press", "iron_bloom", 182); }
            return Collect("press", "iron_ingot", 182 + 12) && QuestBook_CountOf("wheel_press") >= 300.0;
        });
        add("press 16 bronze plates", [] { static bool fed = false; if (!fed) { fed = true; Give("press", "bronze_ingot", 16); } return Collect("press", "bronze_plate", 16); });
        // Glass first: the steam labs take 2 each.
        add("dig 96 sand", [] { return DigSand(96); });
        add("make 24 glass", [] { static bool fed = false; if (!fed) { fed = true; Give("kiln", "sand", 96); Give("kiln", "any_log", 7); } return Collect("kiln", "glass", 24); });
        // (#1271 L6b-3) Tier 3 with the research gate on: the steam lab and its packs come with tier 3; a steam lab over the
        // windmill's gearbox researches steam power, a second on the steam engine's shaft wire drawing (tier 3's working build).
        add("tier 3: a steam lab over the windmill's gearbox", [] {
            return Craft("steam_lab", 1) && Craft("shaft", 1) && PlaceO("labshaft", PartType::Shaft, { 9, 1, 12 }, { 0, 1, 0 }) && Place("steam_lab", PartType::SteamLab, { 9, 2, 12 });
        });
        add("press 40 iron plates for steam packs", [] { static bool fed = false; if (!fed) { fed = true; Give("press", "iron_ingot", 40); } return Collect("press", "iron_plate", 40); });
        add("craft 20 steam packs; the steam lab researches steam power", [] {
            static bool fed = false;
            ResearchState& rs = Research_Team(ResearchTab_Team());
            if (!fed) { fed = true; if (!Craft("steam_pack", 20)) return false; rs.queue = { "steam_power" }; Give("steam_lab", "steam_pack", 20); }
            return rs.Done("steam_power");
        });
        // Tier 3: a boiler and a steam engine, coke, steel.
        add("craft a boiler, a steam engine, a second steam lab and a gearbox", [] {
            return Craft("bronze_plate", 10) && Craft("boiler", 1) && Craft("iron_plate", 6) && Craft("cogwheel", 2) && Craft("steam_engine", 1)
                && Craft("steam_lab", 1) && Craft("gearbox", 1) && Craft("shaft", 2);
        });
        add("place the boiler, the steam engine against it and the lab on its shaft", [] {
            if (!Place("boiler", PartType::Boiler, { 0, 0, 16 }) || !Place("steam_engine", PartType::SteamEngine, { 2, 0, 16 })) return false;
            if (!PlaceO("engineshaft", PartType::Shaft, { 4, 0, 16 }, { 1, 0, 0 }) || !Place("gearbox2", PartType::Gearbox, { 5, 0, 16 })
                || !PlaceO("labshaft2", PartType::Shaft, { 5, 1, 16 }, { 0, 1, 0 }) || !Place("steam_lab2", PartType::SteamLab, { 5, 2, 16 })) return false;
            Give("boiler", "any_log", 12);
            return true;
        });
        add("the steam lab on the steam engine researches wire drawing (tier 3's working build)", [] {
            static bool fed = false;
            ResearchState& rs = Research_Team(ResearchTab_Team());
            if (!fed) { fed = true; if (!Craft("steam_pack", 20)) return false; rs.queue = { "wire_drawing" }; Give("steam_lab2", "steam_pack", 20); }
            Collect("boiler", "-", 0); // (keeps the boiler's water up; its empty buckets back)
            return rs.Done("wire_drawing");
        });
        add("the other tier 3 core nodes (scripted: 100 more steam packs of the same)", [] {
            ResearchState& rs = Research_Team(ResearchTab_Team());
            for (const ResearchNode& n : Research_Nodes()) if (n.tier == 3 && n.core) rs.done.insert(n.id);
            printf("bench: progression: tier 3's other core nodes granted (the lab path ran for two of them)\n");
            return true;
        });
        add("craft and place a coke oven and a blast furnace", [] {
            return Craft("coke_oven", 1) && Craft("blast_furnace", 1) && Place("coke_oven", PartType::CokeOven, { 9, 0, 0 }) && Place("blast_furnace", PartType::BlastFurnace, { 12, 0, 0 });
        });
        add("bake 28 coke", [] { static bool fed = false; if (!fed) { fed = true; Give("coke_oven", "coal", 56); } return Collect("coke_oven", "coke", 28); });
        add("make 110 steel", [] { static bool fed = false; if (!fed) { fed = true; Give("blast_furnace", "iron_ingot", 110); Give("blast_furnace", "coke", 28); } return Collect("blast_furnace", "steel_ingot", 110); });
        add("press 110 steel plates", [] { static bool fed = false; if (!fed) { fed = true; Give("press", "steel_ingot", 110); } return Collect("press", "steel_plate", 110); });
        add("pen two hornhides in a ring of fences (then take their hide)", [] {
            static int phase = 0;
            if (phase == 0) {
                while (inventory.count("any_log") < 20) if (!TerrainTool_HarvestNearest(*body, home, true)) return false;
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
        add("craft 18 vacuum tubes", [] { return Craft("copper_wire", 36) && Craft("vacuum_tube", 18); });
        add("tier 3's milestone: hand in 50 steel plates, 20 cogwheels, 10 vacuum tubes and a find at a lab", [] {
            if (!Craft("cogwheel", 20)) return false;
            if (inventory.count("ruin_schematic") == 0) inventory.add("ruin_schematic", 1); // (scripted: the finds step is to come)
            ResearchGate_HandIn();
            const bool met = ResearchGate_MilestoneMet(Research_Team(ResearchTab_Team()), 3);
            if (!met) { printf("bench: progression: tier 3 milestone not met:"); if (const TierMilestone* m = ResearchGate_Milestone(3)) for (int l = 0; l < (int)m->lines.size(); ++l) printf(" %d=%d", l, (int)Research_Team(ResearchTab_Team()).milestoneDone.count({ 3, l })); printf("\n"); }
            return met;
        });
        add("the quest book reaches tier 4 (the research milestone met)", [] { return QuestBook_Tier() >= 4; });
        add("the tier 4 core nodes (scripted: LV packs come with the LV assembler, the next step)", [] {
            ResearchState& rs = Research_Team(ResearchTab_Team());
            for (const ResearchNode& n : Research_Nodes()) if (n.tier == 4 && n.core) rs.done.insert(n.id);
            return true;
        });
        add("till a field and plant seeds from the bushes", [] {
            while (inventory.count("grain_seeds") + inventory.count("root_seeds") + inventory.count("flax_seeds") < 1)
                if (!TerrainTool_HarvestNearest(*body, home, false)) return false;
            // (#1085) a till paints the ground's surface band: the real ground in that column (the density surface, from 3 m over
            // the height map down; the mesh strays from the height map), not the height map's point
            const Vector3 col = Offset(-20.0f, 10.0f);
            auto ground = [&] {
                Vector3 g = Surface(col);
                for (float r = Vector3Length(g) + 3.0f, r0 = Vector3Length(g) - 3.0f; r > r0; r -= 0.05f)
                    if (body->densityAt(Vector3Scale(col, r)) >= ISO_LEVEL) { g = Vector3Scale(col, r); break; }
                return g;
            };
            Vector3 g = ground();
            if (!IsTillable(body->materialAt(g))) { // frozen ground or rock: a cold start places dirt and tills that
                TerrainEdit fill{}; fill.fill = true; fill.radius = 1.0f; fill.paint = (int)VoxelMaterial::Regolith;
                fill.center = Vector3Add(g, Vector3Scale(col, 0.3f));
                ApplyTerrainEdit(*body, fill);
                g = ground();
            }
            if (!Farm_Till(*scene, *body, g)) { printf("bench: progression: could not till\n"); failed = 1; return false; }
            return Farm_Plant(*body, g, "") && QuestBook_CountOf("plant") >= 1.0; // (#1085) on the tilled ground itself
        });
        add("water the field: a bucket of water poured beside it", [] {
            // The bucket's pour (TerrainTool_Bucket counts "pour" through the game's hook), done directly:
            // the quest task needs the event, the field gets the water.
            Fluid::Pour(*body, Surface(Offset(-22.0f, 10.0f)), 1.0f);
            QuestBook_Count("pour");
            return QuestBook_CountOf("pour") >= 1.0;
        });
        add("craft and fit two suit modules (a food pouch and an oxygen tank)", [] {
            while (inventory.count("fibre") < 40) if (!TerrainTool_HarvestNearest(*body, home, false)) { printf("bench: progression: no bush left near home (fibre %d)\n", inventory.count("fibre")); return false; }
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
            Give("generator", "any_log", 2);
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
        add("craft and place a macerator (#996: a tier 4 core node now, the user's decision)", [] {
            return Craft("lv_hull", 1) && Craft("copper_wire", 8) && Craft("basic_circuit", 2) && Craft("macerator", 1) && Place("macerator", PartType::Macerator, { 6, 0, 8 });
        });
        add("craft and place an LV assembler, cabled to the generator", [] {
            static bool fed = false;
            if (!fed) {
                // (scripted: its parts, more than the base run makes: 12 steel plates, 4 vacuum tubes, 16 copper wire)
                inventory.add("steel_plate", 12); inventory.add("vacuum_tube", 4); inventory.add("copper_wire", 16);
                if (!Craft("lv_assembler", 1) || !Place("lv_assembler", PartType::LvAssembler, { 9, 0, 8 })) return false;
                for (int x = 4; x < 11; ++x) BuildTool_Place(*body, *plot, { PartType::Cable, { x, 0, 7 }, 0 }, false); // (the cable is scripted too)
                fed = true;
            }
            return true;
        });
        add("the LV assembler assembles a basic circuit on the generator's power", [] {
            static bool fed = false;
            if (!fed) {
                Give("generator", "any_log", 4);
                Factory::SetPick(body->id, plot->id, parts["lv_assembler"], "basic_circuit"); // (its recipe picked, as on its screen: an LV pack takes the same parts)
                inventory.add("vacuum_tube", 2); inventory.add("copper_wire", 4); inventory.add("steel_plate", 1); // (a circuit's parts, scripted too)
                Give("lv_assembler", "vacuum_tube", 2); Give("lv_assembler", "copper_wire", 4); Give("lv_assembler", "steel_plate", 1);
                fed = true;
            }
            return Collect("lv_assembler", "basic_circuit", 1);
        });
        add("tier 4's tasks all done (its research milestone waits for the LV assembler)", [] {
            static double told = 0.0;
            bool all = true;
            for (int k = 0; k < QuestBook_TaskCount(4); ++k)
                if (!QuestBook_TaskDone(4, k)) {
                    all = false;
                    if (GetTime() - told > 10.0) printf("bench: progression: tier 4 task not done yet: %s\n", QuestBook_TaskText(4, k));
                }
            if (!all && GetTime() - told > 10.0) told = GetTime();
            return all;
        });
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

    // M24: tiers 9 and 10 on Terra. The player arrives with what a moon and a ring trip bring home
    // (anorthite, regolith, nickel-iron, platinum and iridium ore) and the HV-era basics; the ladder
    // works the real machines: helium-3 in the centrifuge, ring metals and aluminium in the arc
    // furnace, superconductors in the precision assembler, then a fusion reactor powering a mass
    // driver that launches a cargo capsule.
    void BuildEv() {
        auto add = [](const std::string& n, std::function<bool()> f) { steps.push_back({ n, std::move(f) }); };
        add("arrive home with moon and ring materials (and HV-era machines and parts)", [] {
            const std::pair<const char*, int> carried[] = {
                { "anorthite", 16 }, { "regolith", 64 }, { "nickel_iron_ore", 8 }, { "platinum_ore", 32 }, { "iridium_ore", 4 }, { "ice", 28 },
                { "generator", 1 }, { "centrifuge", 1 }, { "hv_generator", 1 }, { "arc_furnace", 1 }, { "precision_assembler", 1 }, { "cable", 40 },
                { "oak_log", 40 }, { "coke", 24 }, { "copper_wire", 56 }, { "titanium_alloy_plate", 8 }, { "advanced_circuit", 3 },
                { "steel_plate", 8 }, { "aluminium_plate", 4 }, { "iron_plate", 2 } };
            for (const auto& c : carried) inventory.add(c.first, c.second);
            return true;
        });
        add("place the LV row (centrifuge) and the HV row (arc furnace, precision assembler)", [] {
            if (!Place("generator", PartType::Generator, { 0, 0, 0 }) || !Place("centrifuge", PartType::Centrifuge, { 3, 0, 0 }) || !Cables(-1, 0, 4)) return false;
            if (!Place("hv_generator", PartType::HvGenerator, { 0, 0, 6 }) || !Place("arc_furnace", PartType::ArcFurnace, { 3, 0, 6 }) ||
                !Place("precision_assembler", PartType::PrecisionAssembler, { 6, 0, 6 }) || !Cables(5, 0, 7)) return false;
            return true;
        });
        add("fuel both generators", [] {
            if (GetTime() - stepStart < 0.5) return false; // the layout reaches the simulation first
            Give("generator", "oak_log", 40); Give("hv_generator", "coke", 24);
            return true;
        });
        add("centrifuge: 2 helium-3 from moon regolith", Run("centrifuge", { { "regolith", 64 } }, "helium3", 2));
        add("arc furnace: aluminium from anorthite", Run("arc_furnace", { { "anorthite", 16 } }, "aluminium_ingot", 4));
        add("arc furnace: nickel and iron from nickel-iron", Run("arc_furnace", { { "nickel_iron_ore", 8 } }, "nickel_ingot", 4));
        add("arc furnace: 16 platinum", Run("arc_furnace", { { "platinum_ore", 32 } }, "platinum_ingot", 16));
        add("arc furnace: 2 iridium", Run("arc_furnace", { { "iridium_ore", 4 } }, "iridium_ingot", 2));
        add("precision assembler: 28 superconductors", [] {
            static auto run = Run("precision_assembler", { { "platinum_ingot", 14 }, { "copper_wire", 56 }, { "ice", 28 } }, "superconductor", 28);
            return run();
        });
        add("craft two EV hulls, a fusion reactor, a mass driver and a cargo capsule", [] {
            return Craft("ev_hull", 2) && Craft("fusion_reactor", 1) && Craft("mass_driver", 1) && Craft("cargo_capsule", 1);
        });
        add("place the EV row (2048 V): fusion reactor and mass driver", [] {
            return Place("fusion_reactor", PartType::FusionReactor, { 0, 0, 12 }) && Place("mass_driver", PartType::MassDriver, { 3, 0, 12 }) && Cables(11, 0, 4);
        });
        add("fuel the fusion reactor with helium-3", [] {
            if (GetTime() - stepStart < 0.5) return false;
            Give("fusion_reactor", "helium3", 2);
            return true;
        });
        add("mass driver: launch a cargo capsule", Run("mass_driver", { { "cargo_capsule", 1 } }, "launched_capsule", 1));
    }

    // milestone 8: the rest of the ladder, through LV ore processing, MV and HV to titanium alloy
    // (rocketry's tier), and a miner working its deposit.
    void BuildFull() {
        ResearchGate_SetOn(false); // (#1271 L6b-3) the full run beyond tier 4 isn't scripted for research yet (LV packs need the LV assembler)
        auto add = [](const std::string& n, std::function<bool()> f) { steps.push_back({ n, std::move(f) }); };
        add("mine for the MV and HV ladder (iron, coal, copper) and dig stone", [] {
            const bool ok = MineOre("iron_ore", 185) && MineOre("coal", 530) && MineOre("copper_ore", 160) && DigFor("any_stone", 40, 6.0f);
            static bool told = false; // (the step retries every frame)
            if (!ok && !told && (told = true)) printf("bench: progression: the MV ladder dig stopped at iron %d/185, coal %d/530, copper %d/160, any_stone %d/40 "
                            "(samples left: iron %zu, coal %zu, copper %zu; holding %s x%d)\n", inventory.count("iron_ore"), inventory.count("coal"),
                            inventory.count("copper_ore"), inventory.count("any_stone"), orePoints["iron_ore"].size() - oreNext["iron_ore"],
                            orePoints["coal"].size() - oreNext["coal"], orePoints["copper_ore"].size() - oreNext["copper_ore"],
                            terrainTool.heldTool.empty() ? "nothing" : terrainTool.heldTool.c_str(), terrainTool.heldTool.empty() ? 0 : inventory.count(terrainTool.heldTool));
            return ok;
        });
        // (#1271 L2, L5) iron from blooms: charcoal from logs, the bloomery, the windmill-driven press
        add("chop 180 more logs for charcoal", [] { while (inventory.count("any_log") < 180) if (!TerrainTool_HarvestNearest(*body, home, true)) return false; return true; });
        add("burn 170 charcoal", Run("campfire", { { "coal", 22 }, { "any_log", 170 } }, "charcoal", 170));
        add("smelt 170 iron blooms", Run("bloomery", { { "iron_ore", 170 }, { "charcoal", 170 } }, "iron_bloom", 170));
        add("press 170 blooms into iron", Run("press", { { "iron_bloom", 170 } }, "iron_ingot", 170));
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
        add("press 160 steel plates", Run("press", { { "steel_ingot", 160 } }, "steel_plate", 160));
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
        // (#1427 L7b) tier 5's MV assembler (its own recipe; scripted parts and placement: the full run checks the ladder, not more grind)
        add("place an MV assembler (tier 5's task)", [] { return BuildTool_Place(*body, *plot, { PartType::LvAssembler, { 12, 0, 16 }, 0, (uint8_t)1 }, false); });
        add("pump 20 more water", Run("water_pump", {}, "water", 20));
        add("electrolyse water into 10 hydrogen", Run("electrolyzer", { { "water", 10 } }, "hydrogen", 10));
        add("make 5 polyethylene", Run("chemical_reactor", { { "coal", 5 }, { "hydrogen", 10 } }, "polyethylene", 5));
        add("electrolyse bauxite dust into 12 aluminium (MV gate)", Run("electrolyzer", { { "bauxite_dust", 24 } }, "aluminium_ingot", 12));
        add("refuel the boiler", [] { Give("boiler", "coal", 6); return true; });
        add("press 11 aluminium plates", Run("press", { { "aluminium_ingot", 11 } }, "aluminium_plate", 11));
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
        add("place an HV assembler (tier 6's task)", [] { return BuildTool_Place(*body, *plot, { PartType::LvAssembler, { 12, 0, 20 }, 0, (uint8_t)2 }, false); });
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

int Progression_Update(Scene& world, bool fullLadder, bool evLadder) {
    full = fullLadder;
    scene = &world;
    if (steps.empty()) {
        body = dynamic_cast<Planetoid*>(world.objects[0].get());
        field = body ? dynamic_cast<const PlanetField*>(body->field.get()) : nullptr;
        if (!field || !body->trees) { printf("bench: --test-progression needs the planet\n"); return -1; }
        inventory.clear();
        QuestBook_SetTier(0);
        Factory::SetSpeed(60);
        home = body->toLocalOf(player.position, world.renderOrigin);
        if (!body->plots) body->plots = std::make_unique<PlotSet>();
        // The workshop: a plot on the ground a few metres ahead.
        const Vector3 up = Vector3Normalize(home);
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 }));
        const Vector3 d = Vector3Normalize(Vector3Add(up, Vector3Scale(t1, 8.0f / field->radius)));
        plot = body->plots->plotFor(Vector3Scale(d, field->radius + field->elevation(d)), 0.0f, true);
        Build();
        if (evLadder) { steps.clear(); BuildEv(); QuestBook_SetTier(9); }
        ev = evLadder;
        current = 0;
        testStart = stepStart = GetTime();
        printf("bench: progression: %zu steps, inventory empty, quest tier 0\n", steps.size());
    }
    if (failed) return -1;
    if (current >= steps.size()) return 1;
    ++frameNo;
    FlushGives(); // deliveries to machines the factory had not taken in yet (#372)
    if (steps[current].run()) {
        printf("bench: progression: [%5.1f s] %s (quest tier %d)\n", GetTime() - testStart, steps[current].name.c_str(), QuestBook_Tier());
        ++current;
        stepStart = GetTime();
        stepHistory.clear();
        if (current >= steps.size()) {
            Factory::SetSpeed(1);
            if (ev) {
                bool done = QuestBook_Tier() > 10;
                if (!done && QuestBook_Tier() == 10) { done = true; for (int k = 0; k < QuestBook_TaskCount(10); ++k) done = done && QuestBook_TaskDone(10, k); }
                printf("bench: progression: moon and ring haul to a launched capsule in %.0f s, quest tier %d -> %s\n", GetTime() - testStart, QuestBook_Tier(), done ? "PASS" : "FAIL");
                return done ? 1 : -1;
            }
            const int need = full ? 7 : 5;
            // (#1271 L6b-3) the base run, gate on: tier 4's tasks all done; its research milestone (40 circuits, the LV packs'
            // throughput) needs the LV assembler (L6c), so the base run passes with only the milestone left
            bool ok = QuestBook_Tier() >= need;
            if (!full && !ok && QuestBook_Tier() == 4 && ResearchGate_On()) {
                ok = true;
                for (int k = 0; k < QuestBook_TaskCount(4); ++k) ok = ok && QuestBook_TaskDone(4, k);
            }
            printf("bench: progression: bare hands to %s in %.0f s, quest tier %d -> %s\n", full ? "titanium alloy (rocketry)" : "a working LV machine",
                   GetTime() - testStart, QuestBook_Tier(), ok ? "PASS" : "FAIL");
            return ok ? 1 : -1;
        }
    } else if (GetTime() - stepStart > 120.0 || failed) {
        printf("bench: progression: FAIL at \"%s\" (logs %d, copper %d/%d ore/ingot, tin %d/%d, iron %d/%d, coal %d, bronze %d, steel %d/%d)\n", steps[current].name.c_str(),
               inventory.count("any_log"), inventory.count("copper_ore"), inventory.count("copper_ingot"), inventory.count("tin_ore"), inventory.count("tin_ingot"),
               inventory.count("iron_ore"), inventory.count("iron_ingot"), inventory.count("coal"), inventory.count("bronze_ingot"), inventory.count("steel_ingot"), inventory.count("steel_plate"));
        // Every placed machine's state (#372: the electrolyse step failed now and then with no clue why).
        for (const auto& [name, part] : parts) {
            Factory::View(body->id, plot->id, part); // watch it; the view fills from the next tick
            Factory::WaitTicks(2);
            const MachineView v = Factory::View(body->id, plot->id, part);
            if (!v.valid) continue;
            printf("bench: progression:   %s: %s, progress %d/%d, energy %d/%d, in [%s], out [%s]%s\n", name.c_str(), MachineStatusName(v.status), v.progress, v.duration,
                   v.energy, v.energyMax, v.input.c_str(), v.output.c_str(), v.hasFuelSlot ? TextFormat(", fuel %s x%d", v.fuelSlot.first.c_str(), v.fuelSlot.second) : "");
        }
        printf("bench: progression:   inventory: water %d, hydrogen %d, oxygen %d\n", inventory.count("water"), inventory.count("hydrogen"), inventory.count("oxygen"));
        for (const std::string& h : stepHistory) printf("bench: progression:   history %s\n", h.c_str()); // the step's last 40 s (#372)
        Factory::SetSpeed(1);
        return -1;
    }
    return 0;
}
