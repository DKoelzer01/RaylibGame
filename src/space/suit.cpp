#include "suit.h"
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

std::string worn = "basic_suit";
std::vector<std::string> fitted(2);
SuitStats stats;

// Player and survival numbers before the suit touches them.
const float kWalk = 6.0f, kJetFuel = 0.1f, kJetAccel = 9.0f, kWalkableCos = 0.55f;

// The ore scanner: the nearest deposit within 64 m, looked for once a second.
std::string scanKind;
float scanDist = -1.0f;
Vector3 scanLocal{};
double lastScan = -10.0;

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

void Recompute() {
    SuitStats s;
    const SuitDef* suit = FindSuit(worn);
    if (suit) { s.slots = suit->slots; s.oxygen = suit->oxygen; s.power = suit->power; s.vacuum = suit->vacuum; }
    for (const std::string& id : fitted) {
        const ModuleDef* m = id.empty() ? nullptr : FindModule(id);
        if (!m) continue;
        const std::string& e = m->effect;
        if (e == "oxygen") s.oxygen += m->value;
        else if (e == "power") s.power += m->value;
        else if (e == "jetpack") { s.jetFuel += m->value; s.jetThrust += m->value; }
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
    stats = s;
}

// Fit the slot list to the suit: extra modules go back to the inventory.
void Resize(int slots) {
    for (size_t i = slots; i < fitted.size(); ++i) if (!fitted[i].empty()) inventory.add(fitted[i], 1);
    std::vector<std::string> keep;
    for (const std::string& id : fitted) if (!id.empty() && (int)keep.size() < slots) keep.push_back(id);
    // Modules past the new count that were packed forward stay fitted; the rest went back above.
    fitted.assign(slots, "");
    for (size_t i = 0; i < keep.size(); ++i) fitted[i] = keep[i];
}

void Scan(Scene& world) {
    scanDist = -1.0f; scanKind.clear();
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
                best = dist; scanDist = dist; scanKind = PlanetField::OreKindName(k); scanLocal = p;
            }
}

} // namespace

const SuitStats& Suit_Stats() { return stats; }
const std::string& Suit_Worn() { return worn; }
const std::vector<std::string>& Suit_Modules() { return fitted; }
bool Suit_IsSuit(const std::string& id) { return FindSuit(id) != nullptr; }
bool Suit_IsModule(const std::string& id) { return FindModule(id) != nullptr; }
const char* Suit_ModuleText(const std::string& id) { const ModuleDef* m = FindModule(id); return m ? m->text.c_str() : ""; }

bool Suit_Wear(const std::string& id) {
    const SuitDef* s = FindSuit(id);
    if (!s || id == worn || inventory.count(id) <= 0) return false;
    inventory.remove(id, 1);
    if (!worn.empty()) inventory.add(worn, 1);
    worn = id;
    Resize(s->slots);
    Recompute();
    Audio_Play(Sfx::Craft, 0.5f, 0.8f);
    GameUI_Toast(TextFormat("Wearing the %s (%d module slots)", s->name.c_str(), s->slots));
    return true;
}

bool Suit_Fit(const std::string& id) {
    if (!FindModule(id) || inventory.count(id) <= 0) return false;
    for (std::string& slot : fitted)
        if (slot.empty()) {
            inventory.remove(id, 1);
            slot = id;
            { int nfit = 0; for (const std::string& m : fitted) nfit += m.empty() ? 0 : 1; QuestBook_AtLeast("modules", nfit); }
            Recompute();
            Audio_Play(Sfx::Build, 0.4f, 1.4f);
            return true;
        }
    GameUI_Toast("No free module slot (a better suit has more)");
    return false;
}

bool Suit_Unfit(int i) {
    if (i < 0 || i >= (int)fitted.size() || fitted[i].empty()) return false;
    inventory.add(fitted[i], 1);
    fitted[i].clear();
    Recompute();
    return true;
}

void Suit_Reset() {
    Load();
    worn = "basic_suit";
    const SuitDef* s = FindSuit(worn);
    fitted.assign(s ? s->slots : 2, "");
    Recompute();
}

void Suit_Apply(Scene& world, float dt) {
    (void)dt;
    Load();
    if (fitted.empty() && !FindSuit(worn)) Suit_Reset();
    player.walkSpeed = kWalk * stats.walk;
    player.jetpackFuelPerSec = kJetFuel / stats.jetFuel;
    player.canHover = stats.hover; // #163
    player.jetpackAccel = kJetAccel * stats.jetThrust;
    player.walkableCos = stats.grip ? 0.3f : kWalkableCos; // ~73 deg with grip boots
    if (stats.scanner && GetTime() - lastScan > 1.0) { lastScan = GetTime(); Scan(world); }
    if (!stats.scanner) scanDist = -1.0f;
}

void Suit_DrawHUD() {
    if (!stats.scanner || survival.dead) return;
    const int x = GetScreenWidth() - 300, y = 380;
    if (scanDist < 0.0f) { DrawText("ORE SCANNER: nothing within 64 m", x, y, 16, Fade(RAYWHITE, 0.6f)); return; }
    DrawText(TextFormat("ORE SCANNER: %s %.0f m", scanKind.c_str(), scanDist), x, y, 16, Color{ 255, 200, 120, 255 });
}

void Suit_DrawTab() {
    Load();
    const SuitDef* s = FindSuit(worn);
    ImGui::Text("Worn: %s  (%d module slots)", s ? s->name.c_str() : worn.c_str(), (int)fitted.size());
    ImGui::Separator();
    for (int i = 0; i < (int)fitted.size(); ++i) {
        ImGui::PushID(i);
        // Each slot is a drop target: drag a carried module here (a fitted one swaps back out).
        const ItemDef* d = fitted[i].empty() ? nullptr : itemDB.find(fitted[i]);
        const std::string label = fitted[i].empty() ? TextFormat("Slot %d: empty (drag a module here)", i + 1)
                                                    : TextFormat("Slot %d: %s", i + 1, d ? d->name.c_str() : fitted[i].c_str());
        ImGui::Selectable(label.c_str(), false, 0, ImVec2(ImGui::GetContentRegionAvail().x - 70.0f, 0.0f));
        if (!fitted[i].empty() && ImGui::BeginDragDropSource()) {
            ImGui::SetDragDropPayload("SUIT_SLOT", &i, sizeof(int));
            ImGui::Text("Remove %s", d ? d->name.c_str() : fitted[i].c_str());
            ImGui::EndDragDropSource();
        }
        if (ImGui::BeginDragDropTarget()) {
            if (const ImGuiPayload* pl = ImGui::AcceptDragDropPayload("SUIT_MODULE")) {
                const std::string id((const char*)pl->Data);
                if (inventory.count(id) > 0) {
                    if (!fitted[i].empty()) Suit_Unfit(i);
                    inventory.remove(id, 1);
                    fitted[i] = id;
                    Recompute();
                    int nfit = 0; for (const std::string& m : fitted) nfit += m.empty() ? 0 : 1;
                    QuestBook_AtLeast("modules", nfit);
                    Audio_Play(Sfx::Build, 0.4f, 1.4f);
                }
            }
            ImGui::EndDragDropTarget();
        }
        if (!fitted[i].empty()) {
            ImGui::SameLine();
            if (ImGui::SmallButton("Remove")) Suit_Unfit(i);
            ImGui::TextDisabled("    %s", Suit_ModuleText(fitted[i]));
        }
        ImGui::PopID();
    }
    ImGui::Separator();
    ImGui::TextUnformatted("Carried (drag onto a slot, or Fit):");
    if (ImGui::BeginDragDropTarget()) { // a fitted module dropped back here comes out
        if (const ImGuiPayload* pl = ImGui::AcceptDragDropPayload("SUIT_SLOT")) Suit_Unfit(*(const int*)pl->Data);
        ImGui::EndDragDropTarget();
    }
    bool any = false;
    for (const auto& m : modules) {
        const int n = inventory.count(m.id);
        if (n <= 0) continue;
        any = true;
        ImGui::PushID(m.id.c_str());
        const ItemDef* d = itemDB.find(m.id);
        if (ImGui::SmallButton("Fit")) Suit_Fit(m.id);
        ImGui::SameLine();
        ImGui::Selectable(TextFormat("%s x%d - %s", d ? d->name.c_str() : m.id.c_str(), n, m.text.c_str()));
        if (ImGui::BeginDragDropSource()) {
            ImGui::SetDragDropPayload("SUIT_MODULE", m.id.c_str(), m.id.size() + 1);
            ImGui::Text("%s", d ? d->name.c_str() : m.id.c_str());
            ImGui::EndDragDropSource();
        }
        ImGui::PopID();
    }
    for (const auto& su : suits) {
        if (su.id == worn || inventory.count(su.id) <= 0) continue;
        any = true;
        ImGui::PushID(su.id.c_str());
        if (ImGui::SmallButton("Wear")) Suit_Wear(su.id);
        ImGui::SameLine();
        ImGui::Text("%s - %d slots, oxygen and power x%.1f", su.name.c_str(), su.slots, su.oxygen);
        ImGui::PopID();
    }
    if (!any) ImGui::TextDisabled("No modules or suits carried (craft them: Crafting tab)");
    ImGui::Separator();
    const SuitStats& st = stats;
    ImGui::Text("Oxygen x%.2f   power x%.2f   jetpack fuel x%.2f, thrust x%.2f", st.oxygen, st.power, st.jetFuel, st.jetThrust);
    ImGui::Text("Walk x%.2f   lamp x%.2f   hunger x%.2f   underwater oxygen x%.1f", st.walk, st.lamp, st.food, st.underwater);
    ImGui::Text("Cold felt %.0f%%   heat %s   grip %s   ore scanner %s", st.cold * 100.0f, st.heatProof ? "regulated" : "hurts above 35 C",
                st.grip ? "on" : "off", st.scanner ? "on" : "off");
}

void Suit_Save(FILE* f) {
    fprintf(f, "suit %s %d", worn.c_str(), (int)fitted.size());
    for (const std::string& m : fitted) fprintf(f, " %s", m.empty() ? "-" : m.c_str());
    fprintf(f, "\n");
}

bool Suit_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "suit") != 0) return false;
    char id[64] = {};
    int n = 0;
    if (fscanf(f, " %63s %d", id, &n) != 2) return true;
    worn = id;
    fitted.assign(std::clamp(n, 0, 16), "");
    for (int i = 0; i < n; ++i) {
        char m[64] = {};
        if (fscanf(f, " %63s", m) != 1) break;
        if (i < (int)fitted.size() && strcmp(m, "-") != 0) fitted[i] = m;
    }
    const SuitDef* s = FindSuit(worn);
    if (s && (int)fitted.size() != s->slots) Resize(s->slots);
    Recompute();
    return true;
}
