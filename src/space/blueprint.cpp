#include "blueprint.h"
#include "buildtool.h"
#include "gameui.h"
#include "inventory.h"
#include "audio.h"
#include "worldtypes.h"

#include <raylib.h>
#include <raymath.h>
#include <algorithm>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <map>

namespace fs = std::filesystem;

namespace {
    std::vector<Blueprint> list;
    bool loaded = false;
    int active = -1;

    void LoadDir(const std::string& dir) {
        std::error_code ec;
        if (!fs::exists(dir, ec)) return;
        std::vector<std::string> files;
        for (const auto& e : fs::directory_iterator(dir, ec))
            if (e.is_regular_file() && e.path().extension() == ".bp") files.push_back(e.path().string());
        std::sort(files.begin(), files.end());
        for (const std::string& f : files) {
            Blueprint bp;
            if (Blueprint_LoadFile(f, bp)) list.push_back(std::move(bp));
        }
    }
    void Ensure() {
        if (loaded) return;
        loaded = true;
        LoadDir("resources/blueprints");
        LoadDir("blueprints");
    }
    // Items a blueprint costs.
    std::map<std::string, int> Cost(const Blueprint& bp) {
        std::map<std::string, int> c;
        for (const Part& p : bp.parts) c[PartItem(p)] += PartCostCount(p.type);
        return c;
    }
}

const std::vector<Blueprint>& Blueprint_List() { Ensure(); return list; }
void Blueprint_Reload() { list.clear(); loaded = false; active = -1; Ensure(); }

bool Blueprint_SaveFile(const Plot& plot, const std::string& path, const std::string& name) {
    // Relative to the lowest corner of the parts' origins.
    Int3 lo{ 1 << 30, 1 << 30, 1 << 30 };
    int n = 0;
    for (const Part& p : plot.parts) {
        if (p.type == PartType::Count) continue;
        lo.x = std::min(lo.x, p.origin.x); lo.y = std::min(lo.y, p.origin.y); lo.z = std::min(lo.z, p.origin.z);
        ++n;
    }
    if (n == 0) return false;
    std::error_code ec;
    fs::create_directories(fs::path(path).parent_path(), ec);
    FILE* f = fopen(path.c_str(), "w");
    if (!f) return false;
    fprintf(f, "# Blueprint: part type, cell x y z (from its lowest corner), orientation, material\nname %s\nparts %d\n", name.c_str(), n);
    for (const Part& p : plot.parts)
        if (p.type != PartType::Count) fprintf(f, "%d %d %d %d %d %d\n", (int)p.type, p.origin.x - lo.x, p.origin.y - lo.y, p.origin.z - lo.z, p.orient, p.material);
    fclose(f);
    return true;
}

bool Blueprint_LoadFile(const std::string& path, Blueprint& out) {
    FILE* f = fopen(path.c_str(), "r");
    if (!f) return false;
    out = Blueprint{};
    out.path = path;
    out.name = fs::path(path).stem().string();
    char line[256];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char name[128];
        int n;
        if (sscanf(line, "name %127[^\n]", name) == 1) { out.name = name; continue; }
        if (sscanf(line, "parts %d", &n) == 1) continue;
        int t, x, y, z, o, m;
        if (sscanf(line, "%d %d %d %d %d %d", &t, &x, &y, &z, &o, &m) == 6 && t >= 0 && t < (int)PartType::Count)
            out.parts.push_back({ (PartType)t, { x, y, z }, (uint8_t)o, (uint8_t)std::clamp(m, 0, BuildMaterialCount() - 1) });
    }
    fclose(f);
    // Bottom up, so supports go in before what rests on them.
    std::stable_sort(out.parts.begin(), out.parts.end(), [](const Part& a, const Part& b) { return a.origin.y < b.origin.y; });
    return !out.parts.empty();
}

const Blueprint* Blueprint_Active() { Ensure(); return active >= 0 && active < (int)list.size() ? &list[active] : nullptr; }
void Blueprint_SetActive(int i) { Ensure(); active = i >= 0 && i < (int)list.size() ? i : -1; }
void Blueprint_Cycle() {
    Ensure();
    if (list.empty()) { GameUI_Toast("No blueprints (B saves the plot you look at)"); return; }
    active = active + 1 >= (int)list.size() ? -1 : active + 1;
    GameUI_Toast(active >= 0 ? TextFormat("Blueprint: %s (RMB places what you can afford; N next)", list[active].name.c_str()) : "Blueprint off");
}

int Blueprint_Place(Planetoid& body, Plot& plot, const Blueprint& bp, Int3 o, bool charge, std::string& missing) {
    int placed = 0;
    std::map<std::string, int> short_;
    for (const Part& p : bp.parts) {
        Part q = p;
        q.origin = { o.x + p.origin.x, o.y + p.origin.y, o.z + p.origin.z };
        if (!plot.canPlace(q)) continue; // already there (filling in) or blocked
        const char* item = PartItem(q);
        if (charge && inventory.count(item) < PartCostCount(q.type)) { short_[item] += PartCostCount(q.type); continue; }
        if (BuildTool_Place(body, plot, q, charge)) ++placed;
    }
    missing.clear();
    for (const auto& kv : short_) missing += TextFormat("%s%d %s", missing.empty() ? "" : ", ", kv.second, kv.first.c_str());
    return placed;
}

bool Blueprint_Keys(Planetoid* body, Plot* lookingAt) {
    (void)body;
    if (IsKeyPressed(KEY_N)) { Blueprint_Cycle(); return true; }
    if (IsKeyPressed(KEY_B) && lookingAt) {
        const std::time_t now = std::time(nullptr);
        char stamp[32];
        std::strftime(stamp, sizeof stamp, "%Y%m%d_%H%M%S", std::localtime(&now));
        const std::string path = std::string("blueprints/plot_") + stamp + ".bp";
        if (Blueprint_SaveFile(*lookingAt, path, std::string("plot ") + stamp)) {
            Blueprint_Reload();
            GameUI_Toast(TextFormat("Saved %s (N to place it)", path.c_str()));
            Audio_Play(Sfx::Build, 0.4f, 1.6f);
        }
        return true;
    }
    return false;
}

void Blueprint_DrawGhost(Planetoid& body, const Plot& frame, Int3 origin, bool ok) {
    const Blueprint* bp = Blueprint_Active();
    if (!bp) return;
    std::vector<std::pair<Int3, CellShape>> fp;
    const Color col = ok ? Fade(SKYBLUE, 0.8f) : Fade(RED, 0.8f);
    for (const Part& p : bp->parts) {
        Part q = p;
        q.origin = { origin.x + p.origin.x, origin.y + p.origin.y, origin.z + p.origin.z };
        Plot::Footprint(q, fp);
        for (const auto& c : fp) {
            // Edges of the cell (12 lines).
            Vector3 v[8];
            for (int i = 0; i < 8; ++i)
                v[i] = body.toWorld(frame.toBody({ (c.first.x + (i & 1)) * kCellSize, (c.first.y + ((i >> 1) & 1)) * kCellSize, (c.first.z + ((i >> 2) & 1)) * kCellSize }));
            const int e[12][2] = { { 0, 1 }, { 2, 3 }, { 4, 5 }, { 6, 7 }, { 0, 2 }, { 1, 3 }, { 4, 6 }, { 5, 7 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
            for (auto& ed : e) DrawLine3D(v[ed[0]], v[ed[1]], col);
        }
    }
}

void Blueprint_DrawHUD() {
    const Blueprint* bp = Blueprint_Active();
    if (!bp) return;
    const auto cost = Cost(*bp);
    int x = GetScreenWidth() - 360, y = 160;
    DrawRectangle(x - 8, y - 8, 350, 40 + 20 * (int)cost.size(), Fade(BLACK, 0.55f));
    DrawText(TextFormat("BLUEPRINT  %s  (%zu parts)", bp->name.c_str(), bp->parts.size()), x, y, 18, SKYBLUE);
    y += 24;
    for (const auto& kv : cost) {
        const int have = inventory.count(kv.first);
        DrawText(TextFormat("%s  %d / %d", kv.first.c_str(), std::min(have, kv.second), kv.second), x, y, 16, have >= kv.second ? RAYWHITE : ORANGE);
        y += 20;
    }
}
