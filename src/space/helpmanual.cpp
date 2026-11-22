#include "helpmanual.h"
#include "structures.h"
#include "inventory.h"
#include "machinetiers.h"
#include "factory.h"
#include "plots.h"
#include "itemicons.h"
#include "recipebook.h"
#include "imgui.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

// The in-game database (#354, user 2026-11-14: like Factorio's Factoriopedia, with an option to show only what has been
// discovered). Left: search, "Discovered only", the hand-written guides (resources/data/help.txt) and every item by kind.
// Right: a page. An item's page is built from the data files, so it stays current: what makes it (by hand, in which
// machine), what it is used in, what a machine makes and burns, what burns it. Every item named on a page is a link;
// Back walks the pages visited.

namespace {
    struct Topic { std::string title, text; };
    std::vector<Topic> topics;
    bool loaded = false;
    char search[64] = "";
    bool discoveredOnly = false; // saved with the world (HelpManual_Save)

    struct Page { int topic = 0; std::string item; }; // a guide (topic index) or an item page (item id)
    Page page;
    std::vector<Page> history;

    void Load() {
        loaded = true;
        topics.clear();
        if (FILE* f = fopen("resources/data/help.txt", "r")) {
            char line[2048];
            while (fgets(line, sizeof(line), f)) {
                std::string s(line);
                while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
                if (s.rfind("##", 0) == 0) continue;                       // comments
                if (s.rfind("# ", 0) == 0) { topics.push_back({ s.substr(2), "" }); continue; }
                if (topics.empty()) continue;
                std::string& t = topics.back().text;
                if (s.empty()) { if (!t.empty() && t.back() != '\n') t += "\n"; continue; } // paragraph break
                if (!t.empty()) t += "\n";
                t += s;
            }
            fclose(f);
        }
    }

    std::string Lower(std::string s) { for (char& c : s) c = (char)std::tolower((unsigned char)c); return s; }
    bool Contains(const std::string& hay, const std::string& needleLower) { return needleLower.empty() || Lower(hay).find(needleLower) != std::string::npos; }

    bool Visible(const std::string& id) { return !discoveredOnly || RecipeBook_Discovered(id); }
    std::string NameOf(const std::string& id) { const ItemDef* d = itemDB.find(id); return d ? d->name : id; }

    // The machine part an item places (machines, generators), or Count.
    PartType MachineOfItem(const std::string& id) {
        for (int t = 0; t < (int)PartType::Count; ++t)
            if (PartMachineName((PartType)t) && id == PartCostItem((PartType)t)) return (PartType)t;
        return PartType::Count;
    }
    std::string MachineItem(PartType t) { const char* c = PartCostItem(t); return c ? c : ""; }

    void Open(const Page& p) { history.push_back(page); page = p; }

    void Icon(const std::string& id, float size) {
        if (!ItemIcons_Has(id)) { ImGui::Dummy(ImVec2(size, size)); return; }
        float uv0[2], uv1[2];
        ItemIcons_UV(id, uv0, uv1);
        ImGui::Image((ImTextureID)(intptr_t)ItemIcons_Texture().id, ImVec2(size, size), ImVec2(uv0[0], uv0[1]), ImVec2(uv1[0], uv1[1]));
    }

    // "3 x Iron ingot" as a link to its page (or ??? when hidden by "Discovered only").
    void Link(const std::string& id, int count, int uid) {
        ImGui::PushID(uid);
        if (IsWildcardItem(id)) { ImGui::TextUnformatted((std::to_string(count) + " x " + (id == "any_planks" ? "any planks" : id == "any_log" ? "any log" : id == "any_stone" ? "any stone" : id)).c_str()); ImGui::PopID(); return; }
        if (!Visible(id)) { ImGui::TextDisabled("%d x ???", count); ImGui::PopID(); return; }
        Icon(id, 20.0f);
        ImGui::SameLine();
        const std::string label = (count > 0 ? std::to_string(count) + " x " : std::string()) + NameOf(id);
        if (ImGui::SmallButton(label.c_str())) Open({ -1, id });
        ImGui::PopID();
    }

    // One recipe line: inputs -> outputs, each a link.
    void RecipeLine(const std::vector<std::pair<std::string, int>>& in, const std::vector<std::pair<std::string, int>>& out, int& uid) {
        for (size_t i = 0; i < in.size(); ++i) { if (i) { ImGui::SameLine(); ImGui::TextDisabled("+"); ImGui::SameLine(); } Link(in[i].first, in[i].second, uid++); }
        ImGui::SameLine(); ImGui::TextDisabled("->"); ImGui::SameLine();
        for (size_t i = 0; i < out.size(); ++i) { if (i) { ImGui::SameLine(); ImGui::TextDisabled(","); ImGui::SameLine(); } Link(out[i].first, out[i].second, uid++); }
    }

    void Heading(const char* s) { ImGui::Spacing(); ImGui::TextColored(ImVec4(0.55f, 0.8f, 1.0f, 1.0f), "%s", s); }

    void DrawItemPage(const std::string& id) {
        const ItemDef* d = itemDB.find(id);
        if (!d) { ImGui::TextDisabled("Unknown item"); return; }
        if (!Visible(id)) { ImGui::TextDisabled("Not discovered yet."); return; }
        Icon(id, 48.0f);
        ImGui::SameLine();
        ImGui::BeginGroup();
        ImGui::TextUnformatted(d->name.c_str());
        ImGui::TextDisabled("%s, %d per slot", d->kind == ItemKind::Tool ? "tool" : d->kind == ItemKind::Block ? "block" : "resource", StackLimit(id));
        ImGui::EndGroup();
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextUnformatted(d->description.c_str());
        ImGui::PopTextWrapPos();
        int uid = 0;
        const std::vector<Factory::RecipeInfo> mr = Factory::RecipeList();
        bool any = false;
        for (const Recipe& r : itemDB.recipes()) {
            if (r.output != id) continue;
            if (!any) { Heading("Made by hand"); any = true; }
            RecipeLine(r.inputs, { { r.output, r.outputCount } }, uid);
        }
        any = false;
        for (const auto& r : mr) {
            bool makes = false;
            for (const auto& o : r.out) makes |= o.first == id;
            if (!makes) continue;
            if (!any) { Heading("Made in machines"); any = true; }
            Link(MachineItem(r.machine), 0, uid++);
            ImGui::SameLine(); ImGui::TextDisabled("(%.1f s%s)", r.ticks / 20.0f, r.eu > 0 ? TextFormat(", %d EU/t", r.eu) : "");
            ImGui::Indent(); RecipeLine(r.in, r.out, uid); ImGui::Unindent(); // on its own line: machine recipes ran off the page
        }
        any = false;
        for (const Recipe& r : itemDB.recipes()) {
            bool uses = false;
            for (const auto& in : r.inputs) uses |= ItemMatches(in.first, id);
            if (!uses) continue;
            if (!any) { Heading("Used by hand in"); any = true; }
            RecipeLine(r.inputs, { { r.output, r.outputCount } }, uid);
        }
        any = false;
        for (const auto& r : mr) {
            bool uses = false;
            for (const auto& in : r.in) uses |= in.first == id;
            if (!uses) continue;
            if (!any) { Heading("Used in machines"); any = true; }
            Link(MachineItem(r.machine), 0, uid++);
            ImGui::Indent(); RecipeLine(r.in, r.out, uid); ImGui::Unindent();
        }
        any = false; // fuel for
        for (int t = 0; t < (int)PartType::Count; ++t) {
            if (!PartMachineName((PartType)t)) continue;
            const std::vector<std::string> fuels = Factory::FuelsOf((PartType)t);
            if (std::find(fuels.begin(), fuels.end(), id) == fuels.end()) continue;
            if (!any) { Heading("Fuel for"); any = true; }
            Link(MachineItem((PartType)t), 0, uid++);
        }
        PartType m = MachineOfItem(id);
        int step = 0;
        { PartType vt; int vs = 0; if (MachineTiers_ItemPart(id, vt, vs)) { m = vt; step = vs; } } // (#1427 L7d) a tiered variant
        if (m != PartType::Count) { // a machine: what it makes and burns
            const int myTier = MachineTiers_Tiered(m) ? MachineTiers_PartTier(m, step, 0) : -1;
            if (myTier >= 0) { // (#1427 L7d) its tiers, each its own item
                Heading("Tiers");
                ImGui::TextWrapped("%s (%d V). Runs its tier's recipes and every one below; those under its tier are overclocked "
                                   "(x%d speed for x%d EU a tick per tier).", MachineTiers_Name(myTier), MachineTiers_Volts(myTier),
                                   MachineTiers_Overclock().speed, MachineTiers_Overclock().eu);
                for (int s = 0; s <= MachineTiers_MaxStep(m); ++s) { if (s) ImGui::SameLine(); Link(MachineTiers_VariantItem(m, s), 0, uid++); }
            }
            if (const PartType next = MachineTiers_Successor(m); next != PartType::Count) { // (#1427 L7e) its job from LV on
                Heading("From LV on");
                ImGui::TextWrapped("An electric machine does this job, with tiers:");
                Link(MachineItem(next), 0, uid++);
            }
            any = false;
            for (const auto& r : mr) {
                if (r.machine != m) continue;
                if (!any) { Heading("Recipes in this machine"); any = true; }
                const int steps = myTier >= 0 && r.tier >= 0 ? MachineTiers_Steps(myTier, r.tier) : 0;
                if (steps < 0) ImGui::TextDisabled("(%s: a %s machine or better)", MachineTiers_Name(r.tier), MachineTiers_Name(r.tier));
                else ImGui::TextDisabled("(%.1f s%s%s)", MachineTiers_Ticks(r.ticks, steps) / 20.0f, r.eu > 0 ? TextFormat(", %d EU/t", MachineTiers_Eu(r.eu, steps)) : "",
                                         steps > 0 ? ", overclocked" : "");
                ImGui::SameLine();
                RecipeLine(r.in, r.out, uid);
            }
            const std::vector<std::string> fuels = Factory::FuelsOf(m);
            if (!fuels.empty()) {
                Heading("Burns");
                for (size_t i = 0; i < fuels.size(); ++i) { if (i) ImGui::SameLine(); Link(fuels[i], 0, uid++); }
            }
        }
    }
}

int HelpManual_TopicCount() { if (!loaded) Load(); return (int)topics.size() + 1; } // + the item database

void HelpManual_OpenItem(const std::string& id) { if (!loaded) Load(); Open({ -1, id }); }

void HelpManual_DrawTab() {
    if (!loaded) Load();
    ImGui::BeginDisabled(history.empty());
    if (ImGui::Button("< Back") && !history.empty()) { page = history.back(); history.pop_back(); }
    ImGui::EndDisabled();
    ImGui::SameLine();
    ImGui::SetNextItemWidth(240.0f);
    ImGui::InputTextWithHint("##helpsearch", "Search guides and items", search, sizeof(search));
    ImGui::SameLine();
    ImGui::Checkbox("Discovered only", &discoveredOnly);
    const std::string q = Lower(search);
    ImGui::BeginChild("helptopics", ImVec2(230.0f, 0.0f), true);
    if (ImGui::CollapsingHeader("Guides", ImGuiTreeNodeFlags_DefaultOpen))
        for (int i = 0; i < (int)topics.size(); ++i) {
            const Topic& t = topics[i];
            if (!Contains(t.title, q) && !Contains(t.text, q)) continue;
            ImGui::PushID(i);
            if (ImGui::Selectable(t.title.c_str(), page.item.empty() && page.topic == i)) Open({ i, "" });
            ImGui::PopID();
        }
    const struct { ItemKind kind; const char* title; } kinds[3] = { { ItemKind::Tool, "Tools" }, { ItemKind::Block, "Machines and blocks" }, { ItemKind::Resource, "Resources" } };
    for (const auto& k : kinds) {
        if (!ImGui::CollapsingHeader(k.title, q.empty() ? 0 : ImGuiTreeNodeFlags_DefaultOpen)) continue;
        for (const ItemDef& d : itemDB.items()) {
            if (d.kind != k.kind || IsWildcardItem(d.id) || !Visible(d.id)) continue;
            if (!Contains(d.name, q) && !Contains(d.description, q) && !Contains(d.id, q)) continue;
            ImGui::PushID(d.id.c_str());
            Icon(d.id, 18.0f);
            ImGui::SameLine();
            if (ImGui::Selectable(d.name.c_str(), page.item == d.id)) Open({ -1, d.id });
            ImGui::PopID();
        }
    }
    ImGui::EndChild();
    ImGui::SameLine();
    ImGui::BeginChild("helptext", ImVec2(0.0f, 0.0f), true);
    if (!page.item.empty()) DrawItemPage(page.item);
    else if (page.topic >= 0 && page.topic < (int)topics.size()) {
        const Topic& t = topics[page.topic];
        ImGui::TextUnformatted(t.title.c_str());
        ImGui::Separator();
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextUnformatted(t.text.c_str());
        if (t.title == "Places") { ImGui::Spacing(); ImGui::TextUnformatted(Structures_PlacesText(discoveredOnly).c_str()); } // each kind of site (#164)
        ImGui::PopTextWrapPos();
    }
    ImGui::EndChild();
}

void HelpManual_Save(FILE* f) { fprintf(f, "wiki_discovered %d\n", discoveredOnly ? 1 : 0); }
bool HelpManual_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "wiki_discovered") != 0) return false;
    int v = 0;
    if (fscanf(f, " %d", &v) == 1) discoveredOnly = v != 0;
    return true;
}
const std::string& HelpManual_PageItem() { return page.item; }
