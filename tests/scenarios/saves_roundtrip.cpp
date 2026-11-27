// saves/roundtrip (#1881): a save loads back to the same save. A small known world (the sandbox plane: a plot with a chest of
// known items and a campfire with coal, a carved pit, the player at a spot) is saved, loaded, and saved again at once; the
// two save directories (the world file, the player records, the region files) must be byte-equal. A difference is printed
// as the file and its first differing line (text) or byte offset (binary): something the save writes that the load doesn't
// read back, or reads back differently. Lines that legitimately move between two saves (the wall-clock save stamp) are
// listed in kMoves with the reason.
#include "scenario.h"
#include "factory.h"
#include "player.h"
#include "plots.h"
#include "world.h"
#include <raymath.h>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {
// Line prefixes that may differ between a save and the save of its load, and why. (Empty until a run shows one that must.)
const std::vector<std::pair<std::string, std::string>> kMoves = {};

struct State { uint32_t plotId = 0; int chest = -1, fire = -1; fs::path a, b; long loadFrame = -1; } R;

// The world file and what sits beside it (not backups, temporaries or the erosion cache: rebuilt from the seed).
void CopySave(const fs::path& from, const fs::path& to) {
    std::error_code ec;
    fs::remove_all(to, ec);
    fs::create_directories(to, ec);
    for (const auto& e : fs::directory_iterator(from, ec)) {
        const std::string n = e.path().filename().string();
        if (n.rfind("world.sav", 0) != 0 || n.find(".backups") != std::string::npos || n.find(".tmp") != std::string::npos ||
            n.find(".erosion") != std::string::npos || n.find(".migrating") != std::string::npos) continue;
        fs::copy(e.path(), to / n, fs::copy_options::recursive | fs::copy_options::overwrite_existing, ec);
    }
}
std::string ReadAll(const fs::path& p) { std::ifstream in(p, std::ios::binary); return std::string(std::istreambuf_iterator<char>(in), {}); }
bool Moves(const std::string& line) {
    for (const auto& [prefix, why] : kMoves) if (line.rfind(prefix, 0) == 0) return true;
    return false;
}
// "" when equal (allowing kMoves lines); else where they first differ.
std::string FirstDifference(const std::string& x, const std::string& y) {
    const bool text = x.find('\0') == std::string::npos && y.find('\0') == std::string::npos;
    if (!text) {
        if (x == y) return "";
        size_t i = 0;
        while (i < x.size() && i < y.size() && x[i] == y[i]) ++i;
        return TextFormat("binary: %zu vs %zu bytes, first difference at byte %zu", x.size(), y.size(), i);
    }
    std::istringstream sx(x), sy(y);
    std::string lx, ly;
    for (int n = 1;; ++n) {
        const bool gx = (bool)std::getline(sx, lx), gy = (bool)std::getline(sy, ly);
        if (!gx && !gy) return "";
        if (gx && gy && (lx == ly || (Moves(lx) && Moves(ly)))) continue;
        return TextFormat("line %d: \"%s\" vs \"%s\"", n, gx ? lx.substr(0, 120).c_str() : "(end)", gy ? ly.substr(0, 120).c_str() : "(end)");
    }
}
}

SCENARIO("saves/roundtrip") {
    s.tags = { "saves", "world" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 120.0;
    s.step("build a known world", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const Vector3 dir = c.spawnDirLocal();
        const float ground = c.field->radius + c.field->elevation(dir);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        Plot* plot = c.body->plots->plotFor(Vector3Scale(dir, ground + 1.0f), 0.0f, true);
        R.plotId = plot->id;
        R.chest = plot->place({ PartType::Chest, { 0, 0, 0 }, 0 });
        R.fire = plot->place({ PartType::Campfire, { 4, 0, 0 }, 0 });
        if (!c.check("chest and campfire placed", R.chest >= 0 && R.fire >= 0)) return true;
        Factory::Sync(c.world);
        Factory::WaitTicks(3);
        Factory::Insert(c.body->id, R.plotId, R.chest, "iron_ingot", 40);
        Factory::Insert(c.body->id, R.plotId, R.chest, "copper_ingot", 7);
        Factory::InsertFuel(c.body->id, R.plotId, R.fire, "coal", 3);
        Factory::WaitTicks(5);
        const Vector3 side = Vector3Normalize(Vector3CrossProduct(dir, { 0.0f, 1.0f, 0.0f }));
        const Vector3 along = Vector3Normalize(Vector3Add(dir, Vector3Scale(side, 12.0f / c.field->radius)));
        const Vector3 dig = Vector3Scale(along, c.field->radius + c.field->elevation(along) - 0.5f);
        const TerrainEdit pit{ dig, 2.5f, -1.0f, true };
        ApplyTerrainEdits(*c.body, &pit, 1);
        c.camera(Vector3Add(dig, Vector3Scale(dir, 6.0f)), dig, dir);
        return true;
    });
    s.step("save, load, save again", 60.0, [](ScenarioCtx& c) {
        if (R.loadFrame < 0) {
            c.waitingFor("the terrain edit applied");
            if (c.t < 1.0) return false;
            Factory::Sync(c.world);
            Factory::WaitTicks(3);
            if (!c.check("the first save", SaveWorld(c.world))) return true;
            const fs::path dir = fs::path(WorldSavePath()).parent_path();
            R.a = fs::path(Scenario_RunDir()) / "save-a"; R.b = fs::path(Scenario_RunDir()) / "save-b";
            CopySave(dir, R.a);
            if (!c.check("the load", LoadWorldSave(c.world), LastLoadError())) return true;
            // Saved again at once, before any sim step: nothing in the world has moved since the load.
            if (!c.check("the second save (right after the load)", SaveWorld(c.world))) return true;
            CopySave(dir, R.b);
            R.loadFrame = c.frame;
        }
        int files = 0, differ = 0;
        std::error_code ec;
        for (const auto& e : fs::recursive_directory_iterator(R.a, ec)) {
            if (!e.is_regular_file()) continue;
            const fs::path rel = fs::relative(e.path(), R.a), other = R.b / rel;
            ++files;
            if (!fs::exists(other)) { ++differ; c.check("in both saves: " + rel.generic_string(), false, "missing after the load"); continue; }
            const std::string d = FirstDifference(ReadAll(e.path()), ReadAll(other));
            if (!d.empty()) { ++differ; c.check("the same after a load: " + rel.generic_string(), false, d); }
        }
        for (const auto& e : fs::recursive_directory_iterator(R.b, ec))
            if (e.is_regular_file() && !fs::exists(R.a / fs::relative(e.path(), R.b))) {
                ++differ; c.check("in both saves: " + fs::relative(e.path(), R.b).generic_string(), false, "only after the load");
            }
        c.note(TextFormat("%d files compared, %d differ (save-a and save-b in the run folder)", files, differ));
        c.check("the save has its world file and regions", fs::exists(R.a / "world.sav") && fs::exists(R.a / "world.sav.regions"));
        c.eq("files that differ between a save and the save of its load", differ, 0);
        return true;
    });
}
