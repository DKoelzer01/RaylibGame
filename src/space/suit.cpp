#include "suit.h"
#include "playersession.h"
#include "skills.h"
#include "audio.h"
#include "gameui.h"
#include "inventory.h"
#include "player.h"
#include "survival.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"
#include "imgui.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>

namespace {

struct SuitDef { std::string id, name; int slots = 2; float oxygen = 1.0f, power = 1.0f; bool vacuum = false; };
struct ModuleDef { std::string id, effect, text; float value = 0.0f; };
std::vector<SuitDef> suits;
std::vector<ModuleDef> modules;
bool loaded = false;

SuitState& S() { return LocalSession().suit; } // the local player's worn suit, modules and scanner (#313)

// Player and survival numbers before the suit touches them.
const float kWalk = 6.0f, kJetFuel = 0.1f, kJetAccel = 9.0f, kWalkableCos = 0.55f;


std::vector<std::string> Columns(const std::string& line) {
    std::vector<std::string> cols;
    size_t start = 0, bar;
    while ((bar = line.find('|', start)) != std::string::npos) { cols.push_back(line.substr(start, bar - start)); start = bar + 1; }
    cols.push_back(line.substr(start));
    for (auto& c : cols) { const size_t a = c.find_first_not_of(" \t"), b = c.find_last_not_of(" \t\r\n"); c = a == std::string::npos ? "" : c.substr(a, b - a + 1); }
    return cols;
}

void Load() {
    if (loaded) return;
    loaded = true;
    char buf[1024];
    if (FILE* f = fopen("resources/data/suits.txt", "r")) {
        while (fgets(buf, sizeof buf, f)) {
            if (buf[0] == '#' || buf[0] == '\n' || buf[0] == '\r') continue;
            const auto c = Columns(buf);
            if (c.size() < 5) continue;
            suits.push_back({ c[0], c[1], atoi(c[2].c_str()), (float)atof(c[3].c_str()), (float)atof(c[4].c_str()), c.size() >= 6 && c[5] == "1" });
        }
        fclose(f);
    }
    if (FILE* f = fopen("resources/data/suit_modules.txt", "r")) {
        while (fgets(buf, sizeof buf, f)) {
            if (buf[0] == '#' || buf[0] == '\n' || buf[0] == '\r') continue;
            const auto c = Columns(buf);
            if (c.size() < 4) continue;
            modules.push_back({ c[0], c[1], c[3], (float)atof(c[2].c_str()) });
        }
        fclose(f);
    }
}
const SuitDef* FindSuit(const std::string& id) { Load(); for (const auto& s : suits) if (s.id == id) return &s; return nullptr; }
const ModuleDef* FindModule(const std::string& id) { Load(); for (const auto& m : modules) if (m.id == id) return &m; return nullptr; }

SuitStats StatsOf(const SuitState& st);
void Recompute() { S().stats = StatsOf(S()); }
SuitStats StatsOf(const SuitState& st) { // (#337: any player's suit, not only the local one)
    SuitStats s;
    const SuitDef* suit = FindSuit(st.worn);
    if (suit) { s.slots = suit->slots; s.oxygen = suit->oxygen; s.power = suit->power; s.vacuum = suit->vacuum; }
    for (const std::string& id : st.fitted) {
        const ModuleDef* m = id.empty() ? nullptr : FindModule(id);
        if (!m) continue;
        const std::string& e = m->effect;
        if (e == "oxygen") s.oxygen += m->value;
        else if (e == "power") { s.power += m->value; s.powered = true; } // #217
        else if (e == "jetpack") { s.jetFuel += m->value; s.jetThrust += m->value; s.jetpack = true; } // #217
        else if (e == "walk") { s.walk += m->value; s.servo = true; }
        else if (e == "lamp") s.lamp += m->value;
        else if (e == "food") s.food = std::max(0.0f, s.food - m->value);
        else if (e == "underwater") s.underwater = std::min(s.underwater, m->value);
        else if (e == "insulation") s.cold = std::min(s.cold, m->value);
        else if (e == "thermal") { s.cold = 0.0f; s.heatProof = true; }
        else if (e == "grip") s.grip = true;
        else if (e == "scanner") s.scanner = true;
        else if (e == "hover") s.hover = true;
    }
    return s;
}

// Fit the slot list to the suit: extra modules go back to the inventory.
void Resize(int slots) {
    for (size_t i = slots; i < S().fitted.size(); ++i) if (!S().fitted[i].empty()) inventory.add(S().fitted[i], 1);
    std::vector<std::string> keep;
    for (const std::string& id : S().fitted) if (!id.empty() && (int)keep.size() < slots) keep.push_back(id);
    // Modules past the new count that were packed forward stay fitted; the rest went back above.
    S().fitted.assign(slots, "");
    for (size_t i = 0; i < keep.size(); ++i) S().fitted[i] = keep[i];
}

void Scan(Scene& world) {
    S().scanDist = -1.0f; S().scanKind.clear();
    Planetoid* body = player.nearest;
    const PlanetField* pf = body ? dynamic_cast<const PlanetField*>(body->field.get()) : nullptr;
    if (!pf) return;
    const Vector3 me = body->toLocal((player.position - world.renderOrigin).toFloat());
    const Vector3 up = Vector3Normalize(me);
    Vector3 east = Vector3CrossProduct({ 0, 1, 0 }, up);
    east = Vector3Length(east) < 1e-3f ? Vector3{ 1, 0, 0 } : Vector3Normalize(east);
    const Vector3 north = Vector3CrossProduct(up, east);
    float best = 64.0f;
    for (float a = -64.0f; a <= 64.0f; a += 8.0f)
        for (float b = -64.0f; b <= 64.0f; b += 8.0f)
            for (float d = 2.0f; d <= 64.0f; d += 6.0f) {
                const Vector3 p = Vector3Add(me, Vector3Add(Vector3Add(Vector3Scale(east, a), Vector3Scale(north, b)), Vector3Scale(up, -d)));
                const float dist = Vector3Distance(p, me);
                if (dist >= best) continue;
                const Vector3 dir = Vector3Normalize(p);
                const float depth = pf->radius + pf->elevation(dir) - Vector3Length(p);
                if (depth < 1.0f) continue;
                const int k = pf->oreKind(p, depth);
                if (k < 0) continue;
                best = dist; S().scanDist = dist; S().scanKind = PlanetField::OreKindName(k); S().scanLocal = p;
            }
}

} // namespace

const SuitStats& Suit_Stats() { return S().stats; }
const std::string& Suit_Worn() { return S().worn; }
const std::vector<std::string>& Suit_Modules() { return S().fitted; }
bool Suit_IsSuit(const std::string& id) { return FindSuit(id) != nullptr; }
bool Suit_IsModule(const std::string& id) { return FindModule(id) != nullptr; }
const char* Suit_ModuleText(const std::string& id) { const ModuleDef* m = FindModule(id); return m ? m->text.c_str() : ""; }

bool Suit_Wear(const std::string& id) {
    const SuitDef* s = FindSuit(id);
    if (!s || id == S().worn || inventory.count(id) <= 0) return false;
    inventory.remove(id, 1);
    if (!S().worn.empty()) inventory.add(S().worn, 1);
    S().worn = id;
    Resize(s->slots);
    Recompute();
    Audio_Play(Sfx::Craft, 0.5f, 0.8f);
    GameUI_Toast(TextFormat("Wearing the %s (%d module slots)", s->name.c_str(), s->slots));
    return true;
}

bool Suit_Fit(const std::string& id) {
    if (!FindModule(id) || inventory.count(id) <= 0) return false;
    for (std::string& slot : S().fitted)
        if (slot.empty()) {
            inventory.remove(id, 1);
            slot = id;
            { int nfit = 0; for (const std::string& m : S().fitted) nfit += m.empty() ? 0 : 1; QuestBook_AtLeast("modules", nfit); }
            Recompute();
            Audio_Play(Sfx::Build, 0.4f, 1.4f);
            return true;
        }
    GameUI_Toast("No free module slot (a better suit has more)");
    return false;
}

bool Suit_Unfit(int i) {
    if (i < 0 || i >= (int)S().fitted.size() || S().fitted[i].empty()) return false;
    inventory.add(S().fitted[i], 1);
    S().fitted[i].clear();
    Recompute();
    return true;
}

void Suit_Reset() {
    Load();
    S().worn = "basic_suit";
    const SuitDef* s = FindSuit(S().worn);
    S().fitted.assign(s ? s->slots : 2, "");
    Recompute();
}

void Suit_ApplyStatsTo(Player& p, const SuitStats& st, float walkMul) {
    p.walkSpeed = kWalk * st.walk * walkMul;
    p.jetpackFuelPerSec = kJetFuel / st.jetFuel;
    p.canHover = st.hover; // #163
    p.jetpackAccel = st.jetpack && st.powered ? kJetAccel * st.jetThrust : 0.0f; // #217: no thruster pack (or no power), no jetpack
    if (!st.jetpack || !st.powered) p.jetpackFuel = 0.0f;
    p.walkableCos = st.grip ? 0.3f : kWalkableCos; // ~73 deg with grip boots
}

void Suit_Apply(Scene& world, float dt) {
    (void)dt;
    Load();
    if (S().fitted.empty() && !FindSuit(S().worn)) Suit_Reset();
    Suit_ApplyStatsTo(player, S().stats, (float)Skills_Stat(LocalSession().skills, "walk_speed")); // (#919: Athletics stacks with the suit)
    if (S().stats.scanner && SimTime() - S().lastScan > 1.0) { S().lastScan = SimTime(); Scan(world); }
    if (!S().stats.scanner) S().scanDist = -1.0f;
}

void Suit_DrawHUD() {
    if (!S().stats.scanner || survival.dead) return;
    const int x = GetScreenWidth() - 300, y = 380;
    if (S().scanDist < 0.0f) { DrawText("ORE SCANNER: nothing within 64 m", x, y, 16, Fade(RAYWHITE, 0.6f)); return; }
    DrawText(TextFormat("ORE SCANNER: %s %.0f m", S().scanKind.c_str(), S().scanDist), x, y, 16, Color{ 255, 200, 120, 255 });
}

// (#822) A carried module into a given slot (a fitted one there goes back to the inventory first): true if fitted.
bool Suit_FitInto(int slot, const std::string& id) {
    Load();
    if (slot < 0 || slot >= (int)S().fitted.size() || !Suit_IsModule(id) || inventory.count(id) <= 0) return false;
    if (!S().fitted[slot].empty()) Suit_Unfit(slot);
    inventory.remove(id, 1);
    S().fitted[slot] = id;
    Recompute();
    int nfit = 0; for (const std::string& m : S().fitted) nfit += m.empty() ? 0 : 1;
    QuestBook_AtLeast("modules", nfit);
    Audio_Play(Sfx::Build, 0.4f, 1.4f);
    return true;
}

// (#822) The suits and modules you carry, each with Wear / Fit (the inventory screen, under the player).
void Suit_DrawCarried() {
    Load();
    bool any = false;
    for (const auto& m : modules) {
        const int n = inventory.count(m.id);
        if (n <= 0) continue;
        any = true;
        ImGui::PushID(m.id.c_str());
        const ItemDef* d = itemDB.find(m.id);
        if (ImGui::SmallButton("Fit")) Suit_Fit(m.id);
        ImGui::SameLine();
        ImGui::TextWrapped("%s x%d - %s", d ? d->name.c_str() : m.id.c_str(), n, m.text.c_str());
        ImGui::PopID();
    }
    for (const auto& su : suits) {
        if (su.id == S().worn || inventory.count(su.id) <= 0) continue;
        any = true;
        ImGui::PushID(su.id.c_str());
        if (ImGui::SmallButton("Wear")) Suit_Wear(su.id);
        ImGui::SameLine();
        ImGui::TextWrapped("%s - %d slots, oxygen and power x%.1f", su.name.c_str(), su.slots, su.oxygen);
        ImGui::PopID();
    }
    if (!any) ImGui::TextDisabled("No other suits or modules carried (craft them on the Crafting screen)");
}

void Suit_Save(FILE* f) {
    fprintf(f, "suit %s %d", S().worn.c_str(), (int)S().fitted.size());
    for (const std::string& m : S().fitted) fprintf(f, " %s", m.empty() ? "-" : m.c_str());
    fprintf(f, "\n");
}

bool Suit_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "suit") != 0) return false;
    char id[64] = {};
    int n = 0;
    if (fscanf(f, " %63s %d", id, &n) != 2) return true;
    S().worn = id;
    S().fitted.assign(std::clamp(n, 0, 16), "");
    for (int i = 0; i < n; ++i) {
        char m[64] = {};
        if (fscanf(f, " %63s", m) != 1) break;
        if (i < (int)S().fitted.size() && strcmp(m, "-") != 0) S().fitted[i] = m;
    }
    const SuitDef* s = FindSuit(S().worn);
    if (s && (int)S().fitted.size() != s->slots) Resize(s->slots);
    Recompute();
    return true;
}
SuitStats Suit_StatsOf(const SuitState& st) { return StatsOf(st); } // (#337)

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"
TEST(SuitStats, NoThrusterPackMeansNoJetpack) { // (#1087: a server's remote players kept the Player defaults: a full, working jetpack)
    Player p; // (a fresh player: jetpack thrust and a full tank by default)
    SuitStats basic; // the starting suit (#217): no thruster pack, no power
    Suit_ApplyStatsTo(p, basic);
    EXPECT_FLOAT_EQ(p.jetpackAccel, 0.0f) << "no thrust";
    EXPECT_FLOAT_EQ(p.jetpackFuel, 0.0f) << "and an empty tank";
    SuitStats packed; packed.jetpack = true; packed.powered = true; packed.jetThrust = 1.5f;
    Player q; q.jetpackFuel = 0.4f;
    Suit_ApplyStatsTo(q, packed);
    EXPECT_GT(q.jetpackAccel, 0.0f) << "a thruster pack with power flies";
    EXPECT_FLOAT_EQ(q.jetpackFuel, 0.4f) << "and keeps its fuel";
    SuitStats unpowered = packed; unpowered.powered = false;
    Suit_ApplyStatsTo(q, unpowered);
    EXPECT_FLOAT_EQ(q.jetpackAccel, 0.0f) << "no power, no jetpack";
}
#endif
