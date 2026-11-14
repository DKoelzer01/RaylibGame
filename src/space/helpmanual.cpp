#include "helpmanual.h"
#include "inventory.h"
#include "imgui.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <string>
#include <vector>

namespace {
    struct Topic { std::string title, text; };
    std::vector<Topic> topics;
    bool loaded = false;
    int selected = 0;
    char search[64] = "";
    const char* kItemsTitle = "Items";

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
        topics.push_back({ kItemsTitle, "" }); // generated when drawn (the item list can change between loads)
    }

    std::string Lower(std::string s) { for (char& c : s) c = (char)std::tolower((unsigned char)c); return s; }
    bool Contains(const std::string& hay, const std::string& needleLower) { return needleLower.empty() || Lower(hay).find(needleLower) != std::string::npos; }

    // An item's line: its description and, if it is crafted by hand, what it takes.
    std::string ItemLine(const ItemDef& d) {
        std::string s = d.description;
        for (const Recipe& r : itemDB.recipes()) {
            if (r.output != d.id) continue;
            s += "\n    Crafted: ";
            for (size_t i = 0; i < r.inputs.size(); ++i) {
                const ItemDef* in = itemDB.find(r.inputs[i].first);
                s += (i ? ", " : "") + std::to_string(r.inputs[i].second) + " " + (in ? in->name : r.inputs[i].first);
            }
            if (r.outputCount > 1) s += " (makes " + std::to_string(r.outputCount) + ")";
            break;
        }
        return s;
    }
}

int HelpManual_TopicCount() { if (!loaded) Load(); return (int)topics.size(); }

void HelpManual_DrawTab() {
    if (!loaded) Load();
    ImGui::SetNextItemWidth(260.0f);
    ImGui::InputTextWithHint("##helpsearch", "Search the manual and items", search, sizeof(search));
    const std::string q = Lower(search);
    ImGui::BeginChild("helptopics", ImVec2(200.0f, 0.0f), true);
    for (int i = 0; i < (int)topics.size(); ++i) {
        const Topic& t = topics[i];
        bool match = Contains(t.title, q) || Contains(t.text, q);
        if (!match && t.title == kItemsTitle)
            for (const ItemDef& d : itemDB.items()) if (Contains(d.name, q) || Contains(d.description, q)) { match = true; break; }
        if (!match) continue;
        if (ImGui::Selectable(t.title.c_str(), selected == i)) selected = i;
    }
    ImGui::EndChild();
    ImGui::SameLine();
    ImGui::BeginChild("helptext", ImVec2(0.0f, 0.0f), true);
    if (selected >= 0 && selected < (int)topics.size()) {
        const Topic& t = topics[selected];
        ImGui::TextUnformatted(t.title.c_str());
        ImGui::Separator();
        if (t.title == kItemsTitle) {
            for (const ItemDef& d : itemDB.items()) {
                if (!Contains(d.name, q) && !Contains(d.description, q) && !Contains(d.id, q)) continue;
                ImGui::TextColored(ImVec4(0.55f, 0.8f, 1.0f, 1.0f), "%s", d.name.c_str());
                ImGui::PushTextWrapPos(0.0f);
                ImGui::TextUnformatted(ItemLine(d).c_str());
                ImGui::PopTextWrapPos();
                ImGui::Spacing();
            }
        } else {
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextUnformatted(t.text.c_str());
            ImGui::PopTextWrapPos();
        }
    }
    ImGui::EndChild();
}
