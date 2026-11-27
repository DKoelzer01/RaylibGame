#include "profiler.h"
#include "uipalette.h"
#include "blueprint.h"
#include "playerinput.h"
#include "buildtool.h"
#include "gameui.h"
#include "inventory.h"
#include "audio.h"
#include "terraintool.h"
#include "worldtypes.h"
#include "playerrecord.h" // (#1886) the local player's team
#include "researchtab.h"

#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>

namespace fs = std::filesystem;

namespace {
    std::vector<Blueprint> list;
    bool loaded = false;
    int active = -1;
    std::map<int, std::vector<std::string>> teamLibs; // (#1886) team -> its blueprints' .bp texts, in the order added
    int clientTeam = 0;                               // (#1886) a client's team, from the host's store (0: not a client)

    std::string Normalised(std::string text) { // LF line ends and a final newline: one blueprint compares equal however it was written
        text.erase(std::remove(text.begin(), text.end(), '\r'), text.end());
        if (!text.empty() && text.back() != '\n') text += '\n';
        return text;
    }

    // skipSame: leave out a file whose text is already in the list (a blueprint saved to the team's and the player's)
    void LoadDir(const std::string& dir, bool skipSame) {
        std::error_code ec;
        if (!fs::exists(dir, ec)) return;
        std::vector<std::string> files;
        for (const auto& e : fs::directory_iterator(dir, ec))
            if (e.is_regular_file() && e.path().extension() == ".bp") files.push_back(e.path().string());
        std::sort(files.begin(), files.end());
        for (const std::string& f : files) {
            Blueprint bp;
            if (!Blueprint_LoadFile(f, bp)) continue;
            if (skipSame && std::any_of(list.begin(), list.end(), [&](const Blueprint& b) { return b.text == bp.text; })) continue;
            list.push_back(std::move(bp));
        }
    }
    void Ensure() {
        if (loaded) return;
        loaded = true;
        if (const int team = Blueprint_LocalTeam()) // (#1886) the team's library first, then the player's own, then the presets
            if (auto it = teamLibs.find(team); it != teamLibs.end())
                for (const std::string& text : it->second) {
                    Blueprint bp;
                    bp.name = "team blueprint";
                    if (Blueprint_ParseText(text, bp)) { bp.team = team; list.push_back(std::move(bp)); }
                }
        LoadDir("blueprints", true);
        LoadDir("resources/blueprints", false);
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

std::string Blueprint_Text(const Plot& plot, const std::string& name) {
    // Relative to the lowest corner of the parts' origins.
    Int3 lo{ 1 << 30, 1 << 30, 1 << 30 };
    int n = 0;
    for (const Part& p : plot.parts) {
        if (p.type == PartType::Count) continue;
        lo.x = std::min(lo.x, p.origin.x); lo.y = std::min(lo.y, p.origin.y); lo.z = std::min(lo.z, p.origin.z);
        ++n;
    }
    if (n == 0) return std::string();
    std::string s = TextFormat("# Blueprint: part type, cell x y z (from its lowest corner), orientation, material\nname %s\nparts %d\n", name.c_str(), n);
    for (const Part& p : plot.parts)
        if (p.type != PartType::Count) s += TextFormat("%d %d %d %d %d %d\n", (int)p.type, p.origin.x - lo.x, p.origin.y - lo.y, p.origin.z - lo.z, p.orient, p.material);
    return s;
}

bool Blueprint_SaveFile(const Plot& plot, const std::string& path, const std::string& name) {
    const std::string text = Blueprint_Text(plot, name);
    if (text.empty()) return false;
    std::error_code ec;
    fs::create_directories(fs::path(path).parent_path(), ec);
    FILE* f = fopen(path.c_str(), "w");
    if (!f) return false;
    fputs(text.c_str(), f);
    fclose(f);
    return true;
}

bool Blueprint_ParseText(const std::string& text, Blueprint& out) {
    out.parts.clear();
    out.text = Normalised(text);
    std::istringstream in(text);
    for (std::string l; std::getline(in, l);) {
        if (!l.empty() && l.back() == '\r') l.pop_back();
        if (l.empty() || l[0] == '#') continue;
        const char* line = l.c_str();
        char name[128];
        int n;
        if (sscanf(line, "name %127[^\n]", name) == 1) { out.name = name; continue; }
        if (sscanf(line, "parts %d", &n) == 1) continue;
        int t, x, y, z, o, m;
        if (sscanf(line, "%d %d %d %d %d %d", &t, &x, &y, &z, &o, &m) == 6 && t >= 0 && t < (int)PartType::Count)
            out.parts.push_back({ (PartType)t, { x, y, z }, (uint8_t)o, (uint8_t)std::clamp(m, 0, BuildMaterialCount() - 1) });
    }
    // Bottom up, so supports go in before what rests on them.
    std::stable_sort(out.parts.begin(), out.parts.end(), [](const Part& a, const Part& b) { return a.origin.y < b.origin.y; });
    return !out.parts.empty();
}

bool Blueprint_LoadFile(const std::string& path, Blueprint& out) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    std::ostringstream text;
    text << f.rdbuf();
    out = Blueprint{};
    out.path = path;
    out.name = fs::path(path).stem().string();
    return Blueprint_ParseText(text.str(), out);
}

// ---- Team libraries (#1886) ----
void Blueprint_AddToTeam(int team, const std::string& text) {
    if (!team) return;
    const std::string t = Normalised(text);
    Blueprint probe;
    if (!Blueprint_ParseText(t, probe)) return; // no parts: nothing to keep
    std::vector<std::string>& lib = teamLibs[team];
    if (std::find(lib.begin(), lib.end(), t) == lib.end()) { lib.push_back(t); loaded = false; }
}

int Blueprint_TeamCount(int team) { auto it = teamLibs.find(team); return it == teamLibs.end() ? 0 : (int)it->second.size(); }

std::string Blueprint_TeamText(int team) {
    std::string s;
    if (auto it = teamLibs.find(team); it != teamLibs.end())
        for (const std::string& text : it->second) s += "blueprint\n" + text + "blueprint_end\n";
    return s;
}

void Blueprint_SetTeamText(int team, const std::string& text) {
    teamLibs.erase(team);
    loaded = false;
    std::istringstream in(text);
    std::string bp;
    bool inside = false;
    for (std::string l; std::getline(in, l);) {
        if (!l.empty() && l.back() == '\r') l.pop_back();
        if (!inside) { inside = l == "blueprint"; bp.clear(); continue; }
        if (l == "blueprint_end") { Blueprint_AddToTeam(team, bp); inside = false; continue; }
        bp += l; bp += '\n';
    }
}

void Blueprint_ApplyTeamText(int team, const std::string& text) {
    clientTeam = team; // a client applies its own team's store
    Blueprint_SetTeamText(team, text);
    active = -1;
}

void Blueprint_ResetTeams() { teamLibs.clear(); clientTeam = 0; loaded = false; list.clear(); active = -1; }

int Blueprint_LocalTeam() { return clientTeam ? clientTeam : Teams_Of(PlayerId_Local()); }

const Blueprint* Blueprint_Active() { Ensure(); return active >= 0 && active < (int)list.size() ? &list[active] : nullptr; }
void Blueprint_SetActive(int i) { Ensure(); active = i >= 0 && i < (int)list.size() ? i : -1; }
void Blueprint_Cycle() {
    Ensure();
    if (list.empty()) { GameUI_Toast("No blueprints (B saves the plot you look at)"); return; }
    active = active + 1 >= (int)list.size() ? -1 : active + 1;
    GameUI_Toast(active >= 0 ? TextFormat("Blueprint: %s (RMB places what you can afford; N next)", list[active].name.c_str()) : "Blueprint off");
}

// ---- Ghosts (#1018) ----
namespace {
struct Ghost { std::string body; uint32_t plot = 0; Part part{}; };
std::vector<Ghost> ghosts;
double nextGhostBuild = 0.0;
constexpr double kGhostEvery = 0.25; // s between two ghosts built (4 a second)
bool SamePart(const Part& a, const Part& b) { return a.type == b.type && a.origin.x == b.origin.x && a.origin.y == b.origin.y && a.origin.z == b.origin.z && a.orient == b.orient && a.material == b.material; }
Part At(const Part& p, Int3 o) { Part q = p; q.origin = { o.x + p.origin.x, o.y + p.origin.y, o.z + p.origin.z }; return q; }
void AddGhost(const Planetoid& body, const Plot& plot, const Part& q) {
    for (const Ghost& g : ghosts) if (g.plot == plot.id && g.body == body.name && SamePart(g.part, q)) return;
    ghosts.push_back({ body.name, plot.id, q });
}
// The 12 edges of each cell of a part's footprint, in the plot's cell space (under its model matrix).
void DrawCellEdges(const Part& q, Color col) {
    static std::vector<std::pair<Int3, CellShape>> fp;
    Plot::Footprint(q, fp);
    for (const auto& c : fp) {
        Vector3 v[8];
        for (int i = 0; i < 8; ++i)
            v[i] = { (c.first.x + (i & 1)) * kCellSize, (c.first.y + ((i >> 1) & 1)) * kCellSize, (c.first.z + ((i >> 2) & 1)) * kCellSize };
        const int e[12][2] = { { 0, 1 }, { 2, 3 }, { 4, 5 }, { 6, 7 }, { 0, 2 }, { 1, 3 }, { 4, 6 }, { 5, 7 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
        for (auto& ed : e) DrawLine3D(v[ed[0]], v[ed[1]], col);
    }
}
}

int Blueprint_Place(Planetoid& body, Plot& plot, const Blueprint& bp, Int3 o, bool charge, std::string& missing, bool leaveGhosts) {
    PROFILE_SCOPE("blueprint place");
    int placed = 0;
    std::map<std::string, int> short_;
    BuildTool_BeginFoundationBatch(); // (#1167) the parts' foundations as one terrain edit batch, at the end
    for (const Part& p : bp.parts) {
        const Part q = At(p, o);
        if (!plot.canPlace(q)) continue; // already there (filling in) or blocked
        const char* item = PartItem(q);
        if (charge && !terrainTool.creative && inventory.count(item) < PartCostCount(q.type)) {
            short_[item] += PartCostCount(q.type);
            if (leaveGhosts) AddGhost(body, plot, q); // (#1018) built later, from the inventory, when you come by
            continue;
        }
        if (BuildTool_Place(body, plot, q, charge)) ++placed;
    }
    BuildTool_EndFoundationBatch(body);
    missing.clear();
    for (const auto& kv : short_) missing += TextFormat("%s%d %s", missing.empty() ? "" : ", ", kv.second, kv.first.c_str());
    return placed;
}

int Blueprint_GhostCount() { return (int)ghosts.size(); }
int Blueprint_GhostsOn(const Plot& plot, const Blueprint& bp, Int3 o) {
    int n = 0;
    for (const Part& p : bp.parts) { const Part q = At(p, o); for (const Ghost& g : ghosts) if (g.plot == plot.id && SamePart(g.part, q)) { ++n; break; } }
    return n;
}
int Blueprint_RemoveGhosts(const Plot& plot, const Blueprint& bp, Int3 o) {
    const size_t before = ghosts.size();
    for (const Part& p : bp.parts) {
        const Part q = At(p, o);
        ghosts.erase(std::remove_if(ghosts.begin(), ghosts.end(), [&](const Ghost& g) { return g.plot == plot.id && SamePart(g.part, q); }), ghosts.end());
    }
    return (int)(before - ghosts.size());
}

int Blueprint_UpdateGhosts(Planetoid* body, Vector3 feetLocal, double now) {
    if (!body || !body->plots || ghosts.empty() || now < nextGhostBuild) return 0;
    for (size_t i = 0; i < ghosts.size(); ++i) {
        const Ghost& g = ghosts[i];
        if (g.body != body->name) continue;
        Plot* plot = body->plots->byId(g.plot);
        if (!plot) { ghosts.erase(ghosts.begin() + (long)i); return 0; } // its plot is gone (it emptied)
        if (!plot->canPlace(g.part)) { ghosts.erase(ghosts.begin() + (long)i); return 0; } // built, or the cell is taken: dropped
        const Vector3 at = plot->toBody({ (g.part.origin.x + 0.5f) * kCellSize, (g.part.origin.y + 0.5f) * kCellSize, (g.part.origin.z + 0.5f) * kCellSize });
        if (Vector3Distance(at, feetLocal) > kGhostReach * body->scale) continue;
        if (!terrainTool.creative && inventory.count(PartItem(g.part)) < PartCostCount(g.part.type)) continue; // not yet
        const Part part = g.part;
        if (!BuildTool_Place(*body, *plot, part, true)) continue;
        ghosts.erase(ghosts.begin() + (long)i);
        nextGhostBuild = now + kGhostEvery;
        Audio_Play(Sfx::Build, 0.3f, 1.3f);
        return 1;
    }
    return 0;
}

void Blueprint_DrawGhosts(Planetoid& body, const DVec3& renderOrigin, Vector3 cameraLocal) {
    if (ghosts.empty() || !body.plots) return;
    const Color col = Fade(SKYBLUE, 0.45f);
    for (auto& plot : body.plots->plots) {
        bool any = false;
        for (const Ghost& g : ghosts) any = any || (g.plot == plot->id && g.body == body.name);
        if (!any || Vector3Distance(plot->toBody({ 0, 0, 0 }), cameraLocal) > 120.0f) continue;
        rlDrawRenderBatchActive();
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(plot->cellModel(body, renderOrigin))); // (#271: under the plot's double-precision model)
        for (const Ghost& g : ghosts) if (g.plot == plot->id && g.body == body.name) DrawCellEdges(g.part, col);
        rlDrawRenderBatchActive();
        rlPopMatrix();
    }
}

void Blueprint_ResetGhosts() { ghosts.clear(); nextGhostBuild = 0.0; }
void Blueprint_SaveGhosts(FILE* f) {
    for (const Ghost& g : ghosts)
        fprintf(f, "ghost %s %u %d %d %d %d %d %d\n", g.body.c_str(), g.plot, (int)g.part.type, g.part.origin.x, g.part.origin.y, g.part.origin.z, (int)g.part.orient, (int)g.part.material);
}
bool Blueprint_LoadGhostLine(const char* tok, FILE* f) {
    if (strcmp(tok, "ghost") != 0) return false;
    char body[64] = {}; unsigned plot = 0; int type = 0, x = 0, y = 0, z = 0, orient = 0, mat = 0;
    if (fscanf(f, " %63s %u %d %d %d %d %d %d", body, &plot, &type, &x, &y, &z, &orient, &mat) == 8 && type >= 0 && type < (int)PartType::Count) {
        Ghost g; g.body = body; g.plot = plot;
        g.part.type = (PartType)type; g.part.origin = { x, y, z }; g.part.orient = (uint8_t)orient; g.part.material = (uint8_t)mat;
        ghosts.push_back(g);
    }
    return true;
}

bool Blueprint_Keys(Planetoid* body, Plot* lookingAt) {
    (void)body;
    if (PlayerInput_BindingPressed(Act::BlueprintCycle)) { Blueprint_Cycle(); return true; }
    if (PlayerInput_BindingPressed(Act::BlueprintSave) && lookingAt) {
        const std::time_t now = std::time(nullptr);
        char stamp[32];
        std::strftime(stamp, sizeof stamp, "%Y%m%d_%H%M%S", std::localtime(&now));
        const std::string path = std::string("blueprints/plot_") + stamp + ".bp";
        const std::string name = std::string("plot ") + stamp;
        if (Blueprint_SaveFile(*lookingAt, path, name)) {
            // (#1886) and into the team's library (the user: both). A client's addition stays local until the net
            // carries it to the host (#1886 step 4).
            const int team = clientTeam ? clientTeam : ResearchTab_Team();
            Blueprint_AddToTeam(team, Blueprint_Text(*lookingAt, name));
            Blueprint_Reload();
            GameUI_Toast(TextFormat("Saved %s, and to your team's blueprints (N to place it)", path.c_str()));
            Audio_Play(Sfx::Build, 0.4f, 1.6f);
        }
        return true;
    }
    return false;
}

void Blueprint_DrawGhost(Planetoid& body, const Plot& frame, Int3 origin, bool ok, const DVec3& renderOrigin) {
    const Blueprint* bp = Blueprint_Active();
    if (!bp) return;
    std::vector<std::pair<Int3, CellShape>> fp;
    const Color col = ok ? Fade(UiPalette_Mode() == ColorVision::Normal ? SKYBLUE : Ui().good, 0.8f) : Fade(Ui().bad, 0.8f); // (#1432)
    // In plot coordinates under the plot's double-precision model, as its blocks are drawn (#271).
    rlDrawRenderBatchActive();
    rlPushMatrix();
    rlMultMatrixf(MatrixToFloat(frame.cellModel(body, renderOrigin)));
    struct PopAtEnd { ~PopAtEnd() { rlDrawRenderBatchActive(); rlPopMatrix(); } } popAtEnd;
    for (const Part& p : bp->parts) {
        Part q = p;
        q.origin = { origin.x + p.origin.x, origin.y + p.origin.y, origin.z + p.origin.z };
        Plot::Footprint(q, fp);
        for (const auto& c : fp) {
            // Edges of the cell (12 lines).
            Vector3 v[8];
            for (int i = 0; i < 8; ++i)
                v[i] = { (c.first.x + (i & 1)) * kCellSize, (c.first.y + ((i >> 1) & 1)) * kCellSize, (c.first.z + ((i >> 2) & 1)) * kCellSize };
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
