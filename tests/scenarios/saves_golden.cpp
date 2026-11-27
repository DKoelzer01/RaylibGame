// saves/golden-write and saves/golden-load (#1693): golden saves, one per save version, kept in git under tests/saves/v<N>/.
//
// golden-write builds a small known world (seed 42, the sandbox plane): a plot with a chest holding known items and a campfire
// holding coal, a carved pit in the ground, the player at a known spot; saves it and copies the save (the world file, the
// player records, the region files; not the backups or the erosion cache) into <run dir>/golden-v<kSaveVersion>/ together with
// facts.txt (what it built, read back by golden-load). The command that produces a new one:
//   $GAME --seed 42 --scenario saves/golden-write --scenario-out <dir>   ->  <dir>/saves/golden-write/golden-v13/
// (the run folder is <dir>/<scenario name>); commit that folder as tests/saves/v13/.
//
// golden-load loads every tests/saves/v*/ through the migration chain (a copy beside the run's save: the load never touches
// the committed files) and checks facts.txt against the loaded world. No directory yet: it notes so and passes.
#include "scenario.h"
#include "factory.h"
#include "gamesave.h"
#include "player.h"
#include "plots.h"
#include "saveversion.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>

namespace fs = std::filesystem;

namespace {
// facts.txt: "key value..." lines.
using Facts = std::map<std::string, std::string>;
Facts ReadFacts(const fs::path& p) {
    Facts f;
    std::ifstream in(p);
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const size_t sp = line.find(' ');
        if (sp != std::string::npos && f.find(line.substr(0, sp)) == f.end()) f[line.substr(0, sp)] = line.substr(sp + 1);
    }
    return f;
}
std::string Num(double v) { char b[64]; snprintf(b, sizeof b, "%.6f", v); return b; }

Planetoid* Body(ScenarioCtx& c) { return c.world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(c.world.objects[0].get()); }
// The fuel in a machine's slot: View() is asked for, the sim publishes it, then it is read (as the machine scenarios do).
int FuelOf(ScenarioCtx& c, uint64_t bodyId, uint32_t plotId, int part) {
    Factory::View(bodyId, plotId, part);
    Factory::WaitTicks(2);
    Factory::Sync(c.world);
    return Factory::View(bodyId, plotId, part).fuelSlot.second;
}

// The world file and what sits beside it, into `to` (not backups, temporaries or the erosion cache: rebuilt from the seed).
void CopySave(const fs::path& from, const fs::path& to) {
    std::error_code ec;
    fs::create_directories(to, ec);
    for (const auto& e : fs::directory_iterator(from, ec)) {
        const std::string n = e.path().filename().string();
        if (n.rfind("world.sav", 0) != 0 || n.find(".backups") != std::string::npos || n.find(".tmp") != std::string::npos ||
            n.find(".erosion") != std::string::npos || n.find(".migrating") != std::string::npos) continue;
        fs::copy(e.path(), to / n, fs::copy_options::recursive | fs::copy_options::overwrite_existing, ec);
    }
}

struct Built { uint32_t plotId = 0; int chest = -1, fire = -1; Vector3 dig{}; } B;
}

SCENARIO("saves/golden-write") {
    s.tags = { "saves", "local" }; // (a producer, not a gate: run on the sandbox when a save version is bumped)
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
    s.step("build the known world", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const Vector3 dir = c.spawnDirLocal();
        const float ground = c.field->radius + c.field->elevation(dir);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        Plot* plot = c.body->plots->plotFor(Vector3Scale(dir, ground + 1.0f), 0.0f, true);
        B.plotId = plot->id;
        B.chest = plot->place({ PartType::Chest, { 0, 0, 0 }, 0 });
        B.fire = plot->place({ PartType::Campfire, { 4, 0, 0 }, 0 });
        if (!c.check("chest and campfire placed", B.chest >= 0 && B.fire >= 0)) return true;
        Factory::Sync(c.world);
        Factory::WaitTicks(3);
        Factory::Insert(c.body->id, B.plotId, B.chest, "iron_ingot", 40);
        Factory::Insert(c.body->id, B.plotId, B.chest, "copper_ingot", 7);
        Factory::InsertFuel(c.body->id, B.plotId, B.fire, "coal", 3);
        Factory::WaitTicks(5);
        // A pit 12 m along the ground from the start.
        const Vector3 side = Vector3Normalize(Vector3CrossProduct(dir, { 0.0f, 1.0f, 0.0f }));
        const Vector3 along = Vector3Normalize(Vector3Add(dir, Vector3Scale(side, 12.0f / c.field->radius)));
        B.dig = Vector3Scale(along, c.field->radius + c.field->elevation(along) - 0.5f);
        const TerrainEdit pit{ B.dig, 2.5f, -1.0f, true };
        ApplyTerrainEdits(*c.body, &pit, 1);
        c.camera(Vector3Add(B.dig, Vector3Scale(dir, 6.0f)), B.dig, dir);
        return true;
    });
    s.step("save and copy the save directory", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("the terrain edit applied");
        if (c.T < 1.0) return false;
        Factory::Sync(c.world);
        Factory::WaitTicks(3);
        const bool saved = SaveWorld(c.world);
        if (!c.check("saved", saved)) return true;
        const fs::path out = fs::path(Scenario_RunDir()) / ("golden-v" + std::to_string(kSaveVersion));
        std::error_code ec;
        fs::remove_all(out, ec);
        CopySave(fs::path(WorldSavePath()).parent_path(), out);
        const DVec3 local = c.body->toLocalD(player.position);
        std::ostringstream f;
        f << "version " << kSaveVersion << "\nseed 42\n";
        f << "player_local " << Num(local.x) << " " << Num(local.y) << " " << Num(local.z) << "\n";
        f << "strokes " << c.body->edits->strokes << "\n";
        f << "dig " << Num(B.dig.x) << " " << Num(B.dig.y) << " " << Num(B.dig.z) << " " << Num(c.body->densityAt(B.dig)) << "\n";
        f << "plot " << B.plotId << " chest " << B.chest << " fire " << B.fire << "\n";
        for (const auto& [item, n] : Factory::ChestContents(c.body->id, B.plotId, B.chest)) f << "chest_item " << item << " " << n << "\n";
        f << "fire_coal " << FuelOf(c, c.body->id, B.plotId, B.fire) << "\n";
        std::ofstream(out / "facts.txt") << f.str();
        c.note("golden save written to " + out.string());
        c.check("the copy holds the world file", fs::exists(out / "world.sav"));
        c.check("the copy holds the regions", fs::exists(out / "world.sav.regions"));
        c.check("the copy holds the player record", fs::exists(out / "world.sav.players"));
        return true;
    });
}

namespace {
struct LoadState { std::vector<fs::path> dirs; size_t i = 0; int phase = 0; long f0 = 0; bool listed = false; fs::path copy; Facts facts; std::string oldPath; } L;
}

SCENARIO("saves/golden-load") {
    s.tags = { "saves" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 180.0;
    s.step("every tests/saves/v*/ loads through the chain with its facts", 120.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        if (!L.listed) {
            L.listed = true;
            std::error_code ec;
            // The runner's run folder links only resources/ and assets/ from the checkout: find tests/ beside the real
            // resources/ when it isn't in the working directory (the first runner run found none and checked nothing).
            fs::path saves = "tests/saves";
            if (!fs::is_directory(saves, ec)) saves = fs::canonical("resources", ec).parent_path() / "tests" / "saves";
            for (const auto& e : fs::directory_iterator(saves, ec))
                if (e.is_directory() && e.path().filename().string().rfind("v", 0) == 0) L.dirs.push_back(e.path());
            std::sort(L.dirs.begin(), L.dirs.end());
            L.oldPath = WorldSavePath();
            // The goldens are in git (tests/saves/v13 on): none found is a broken run, not a pass.
            if (!c.check("tests/saves/v*/ found", !L.dirs.empty(), saves.string())) return true;
        }
        if (L.i >= L.dirs.size()) { SetWorldSavePath(L.oldPath.c_str()); return true; }
        const fs::path& dir = L.dirs[L.i];
        const std::string name = dir.filename().string();
        if (L.phase == 0) { // a copy beside the run's save (loading an older save backs it up and writes beside it)
            L.copy = fs::path(L.oldPath).parent_path() / ("golden_" + name);
            std::error_code ec;
            fs::remove_all(L.copy, ec);
            fs::create_directories(L.copy, ec);
            // Entry by entry: a directory copied onto an existing one fails on MinGW's libstdc++ ("File exists") and the
            // load found no file there (#1881: the scenario had only run on Linux).
            for (const auto& e : fs::directory_iterator(dir, ec))
                fs::copy(e.path(), L.copy / e.path().filename(), fs::copy_options::recursive | fs::copy_options::overwrite_existing, ec);
            L.facts = ReadFacts(L.copy / "facts.txt");
            SetWorldSavePath((L.copy / "world.sav").string().c_str());
            const bool loaded = LoadWorldSave(c.world);
            c.check(name + ": loads", loaded, LastLoadError());
            if (!loaded) { ++L.i; return false; }
            L.phase = 1; L.f0 = c.frame;
            return false;
        }
        if (c.frame - L.f0 < 20) return false; // the layouts reach the factory
        Factory::Sync(c.world);
        Factory::WaitTicks(5);
        Planetoid* body = Body(c);
        if (!c.check(name + ": the body is there", body != nullptr)) { ++L.i; L.phase = 0; return false; }
        auto vals = [&](const std::string& key) { std::istringstream in(L.facts[key]); std::vector<double> v; double x; while (in >> x) v.push_back(x); return v; };
        const std::vector<double> pl = vals("player_local");
        if (pl.size() == 3) {
            const DVec3 now = body->toLocalD(player.position);
            c.near(name + ": the player's spot", (now - DVec3(pl[0], pl[1], pl[2])).length(), 0.0, 0.05);
        }
        c.check(name + ": the player record loaded", GameSave_PlayerLoaded());
        c.eq(name + ": terrain strokes", (double)body->edits->strokes, std::atof(L.facts["strokes"].c_str()));
        const std::vector<double> dig = vals("dig");
        if (dig.size() == 4) c.near(name + ": the carved pit's density", body->densityAt(Vector3{ (float)dig[0], (float)dig[1], (float)dig[2] }), dig[3], 0.01);
        std::istringstream pf(L.facts["plot"]);
        std::string w1, w2;
        unsigned plotId = 0; int chest = -1, fire = -1;
        pf >> plotId >> w1 >> chest >> w2 >> fire;
        c.check(name + ": the plot is there", body->plots && body->plots->byId(plotId) != nullptr);
        // The chest holds what was put in: every chest_item fact (several lines: read the file again).
        std::ifstream in(L.copy / "facts.txt");
        std::string line;
        std::map<std::string, int> want;
        while (std::getline(in, line)) {
            std::istringstream ls(line);
            std::string key, item;
            int n = 0;
            if ((ls >> key >> item >> n) && key == "chest_item") want[item] = n;
        }
        std::map<std::string, int> got;
        for (const auto& [item, n] : Factory::ChestContents(body->id, plotId, chest)) got[item] += n;
        c.check(name + ": the chest's items", got == want, std::to_string(got.size()) + " kinds, wanted " + std::to_string(want.size()));
        c.eq(name + ": the campfire's coal", FuelOf(c, body->id, plotId, fire), std::atof(L.facts["fire_coal"].c_str()));
        ++L.i; L.phase = 0;
        return false;
    });
}
