#include "questtree.h"
#include "audio.h"
#include "gameui.h"
#include "inventory.h"
#include "survival.h"
#include "imgui.h"

#include <raylib.h>
#include <algorithm>
#include <cstring>
#include <map>
#include <set>
#include <vector>

namespace {
    struct Line { std::string id, name; Color color = WHITE; };
    struct Node {
        std::string id, line, kind, target, title;
        double n = 1.0;
        std::vector<std::string> reqNodes;
        int reqTier = -1;
        std::vector<std::pair<std::string, int>> reward;
        float x = 0.0f, y = 0.0f; // layout (columns, rows)
    };
    std::vector<Line> lines;
    std::vector<Node> nodes;
    std::set<std::string> done;
    bool loaded = false;

    std::vector<std::string> Split(const std::string& s, char sep) {
        std::vector<std::string> out;
        size_t start = 0, at;
        while ((at = s.find(sep, start)) != std::string::npos) { out.push_back(s.substr(start, at - start)); start = at + 1; }
        out.push_back(s.substr(start));
        for (auto& c : out) { const size_t a = c.find_first_not_of(" \t"), b = c.find_last_not_of(" \t\r\n"); c = a == std::string::npos ? "" : c.substr(a, b - a + 1); }
        return out;
    }

    const Node* Find(const std::string& id) { for (const Node& n : nodes) if (n.id == id) return &n; return nullptr; }
    const Line* LineOf(const std::string& id) { for (const Line& l : lines) if (l.id == id) return &l; return nullptr; }

    void Load() {
        if (loaded) return;
        loaded = true;
        FILE* f = fopen("resources/data/quests.txt", "r");
        if (!f) return;
        char buf[1024];
        while (fgets(buf, sizeof buf, f)) {
            if (buf[0] == '#' || buf[0] == '\n' || buf[0] == '\r') continue;
            const auto c = Split(buf, '|');
            if (c[0] == "line" && c.size() >= 4) {
                Line l; l.id = c[1]; l.name = c[2];
                int r, g, b;
                if (sscanf(c[3].c_str(), "%d,%d,%d", &r, &g, &b) == 3) l.color = { (unsigned char)r, (unsigned char)g, (unsigned char)b, 255 };
                lines.push_back(l);
            } else if (c[0] == "node" && c.size() >= 7) {
                Node n; n.id = c[1]; n.line = c[2]; n.title = c[5];
                for (const std::string& r : Split(c[3], ',')) {
                    if (r.empty() || r == "-") continue;
                    int t;
                    if (sscanf(r.c_str(), "tier %d", &t) == 1) n.reqTier = std::max(n.reqTier, t);
                    else n.reqNodes.push_back(r);
                }
                char kind[32] = {}, target[64] = {};
                double cnt = 1.0;
                const int got = sscanf(c[4].c_str(), "%31s %63s %lf", kind, target, &cnt);
                if (got < 1) continue;
                n.kind = kind; n.target = got >= 2 ? target : ""; n.n = got >= 3 ? cnt : 1.0;
                for (const std::string& r : Split(c[6], ',')) {
                    char id[64] = {}; int k = 1;
                    if (r != "-" && sscanf(r.c_str(), "%63s %d", id, &k) >= 1) n.reward.push_back({ id, k });
                }
                nodes.push_back(n);
            }
        }
        fclose(f);
        // Layout: the main line in column 0 (row = tier); each side line its own column, a node one row
        // below its deepest requirement (a tier requirement: that tier's row); siblings side by side.
        std::map<std::string, float> row;
        std::map<std::pair<std::string, int>, int> used; // (line, row) -> nodes there
        std::map<std::string, int> width;                // line -> columns it needs (most siblings in a row)
        std::vector<int> sibOf(nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i) {
            Node& n = nodes[i];
            float y = n.reqTier >= 0 ? (float)n.reqTier : 0.0f;
            for (const std::string& r : n.reqNodes) if (row.count(r)) y = std::max(y, row[r] + 1.0f);
            row[n.id] = y;
            n.y = y;
            sibOf[i] = used[{ n.line, (int)y }]++;
            width[n.line] = std::max(width[n.line], sibOf[i] + 1);
        }
        std::map<std::string, float> start;
        float x = 1.0f;
        for (const Line& l : lines) { start[l.id] = x; x += (float)std::max(1, width[l.id]) + 0.15f; }
        for (size_t i = 0; i < nodes.size(); ++i) nodes[i].x = start[nodes[i].line] + (float)sibOf[i];
    }

    bool ReqsDone(const Node& n) {
        if (n.reqTier >= 0 && QuestBook_Tier() < n.reqTier) return false;
        for (const std::string& r : n.reqNodes) if (!done.count(r)) return false;
        return true;
    }
    std::string TaskText(const Node& n) {
        if (n.kind == "get") return TextFormat("Have %g %s (ever)", n.n, n.target.c_str());
        if (n.kind == "place") return TextFormat("Place a %s", n.target.c_str());
        if (n.kind == "do") return TextFormat("%s: %g", n.target.c_str(), n.n);
        return n.kind + " " + n.target;
    }
    std::string RewardText(const Node& n) {
        std::string s;
        for (const auto& r : n.reward) s += TextFormat("%s%d %s", s.empty() ? "" : ", ", r.second, r.first.c_str());
        return s.empty() ? "-" : s;
    }
}

void QuestTree_Update() {
    Load();
    for (const Node& n : nodes) {
        if (done.count(n.id) || !ReqsDone(n) || !QuestBook_Met(n.kind, n.target, n.n)) continue;
        done.insert(n.id);
        for (const auto& r : n.reward) if (itemDB.find(r.first)) inventory.add(r.first, r.second);
        const Line* l = LineOf(n.line);
        GameUI_Toast(TextFormat("%s: %s done  (+%s)", l ? l->name.c_str() : "Quest", n.title.c_str(), RewardText(n).c_str()));
        Audio_Play(Sfx::Craft, 0.6f, 1.2f);
    }
}

void QuestTree_Reset() { Load(); done.clear(); }
int QuestTree_NodeCount() { Load(); return (int)nodes.size(); }
bool QuestTree_IsDone(const std::string& id) { return done.count(id) > 0; }
bool QuestTree_IsOpen(const std::string& id) { const Node* n = Find(id); return n && !done.count(id) && ReqsDone(*n); }

void QuestTree_Save(FILE* f) { for (const std::string& id : done) fprintf(f, "side_done %s\n", id.c_str()); }
bool QuestTree_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "side_done") != 0) return false;
    char id[64] = {};
    if (fscanf(f, " %63s", id) == 1) done.insert(id);
    return true;
}

void QuestTree_DrawTab() {
    Load();
    const float cw = 232.0f, ch = 74.0f, bw = 210.0f, bh = 46.0f, pad = 12.0f;
    float maxX = 0.0f, maxY = (float)QuestBook_TierCount();
    for (const Node& n : nodes) { maxX = std::max(maxX, n.x); maxY = std::max(maxY, n.y + 1.0f); }
    // Legend: the lines' colours.
    for (const Line& l : lines) {
        ImGui::TextColored(ImVec4(l.color.r / 255.0f, l.color.g / 255.0f, l.color.b / 255.0f, 1.0f), "%s", l.name.c_str());
        ImGui::SameLine();
    }
    ImGui::TextDisabled("  green: done   white: open   grey: locked (hover for details)");
    ImGui::BeginChild("questtree", ImVec2(0, 0), true, ImGuiWindowFlags_HorizontalScrollbar);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 o = ImGui::GetCursorScreenPos();
    ImGui::Dummy(ImVec2((maxX + 1.0f) * cw + pad * 2.0f, maxY * ch + pad * 2.0f));
    auto boxAt = [&](float x, float y) { return ImVec2(o.x + pad + x * cw, o.y + pad + y * ch); };
    const int cur = QuestBook_Tier();
    // Edges first.
    for (const Node& n : nodes) {
        const ImVec2 b = boxAt(n.x, n.y);
        const ImVec2 to(b.x + bw * 0.5f, b.y);
        if (n.reqTier >= 0) {
            const ImVec2 t = boxAt(0.0f, (float)n.reqTier);
            dl->AddLine(ImVec2(t.x + bw, t.y + bh * 0.5f), to, IM_COL32(150, 150, 150, 120), 1.5f);
        }
        for (const std::string& r : n.reqNodes)
            if (const Node* p = Find(r)) {
                const ImVec2 pb = boxAt(p->x, p->y);
                dl->AddLine(ImVec2(pb.x + bw * 0.5f, pb.y + bh), to, done.count(r) ? IM_COL32(120, 210, 120, 200) : IM_COL32(150, 150, 150, 120), 2.0f);
            }
    }
    // The main line: one box per tier.
    for (int t = 0; t < QuestBook_TierCount(); ++t) {
        const ImVec2 b = boxAt(0.0f, (float)t);
        const ImU32 fill = t < cur ? IM_COL32(40, 90, 40, 255) : t == cur ? IM_COL32(70, 70, 30, 255) : IM_COL32(35, 35, 40, 255);
        dl->AddRectFilled(b, ImVec2(b.x + bw, b.y + bh), fill, 6.0f);
        dl->AddRect(b, ImVec2(b.x + bw, b.y + bh), IM_COL32(230, 210, 140, 255), 6.0f, 0, t == cur ? 2.5f : 1.0f);
        dl->AddText(ImVec2(b.x + 8, b.y + 5), IM_COL32(240, 230, 200, 255), TextFormat("Tier %d", t));
        dl->AddText(ImVec2(b.x + 8, b.y + 24), IM_COL32(220, 220, 220, 255), QuestBook_TierName(t));
    }
    // Side nodes.
    for (const Node& n : nodes) {
        const ImVec2 b = boxAt(n.x, n.y);
        const bool isDone = done.count(n.id) > 0, open = !isDone && ReqsDone(n);
        const Line* l = LineOf(n.line);
        const ImU32 lc = l ? IM_COL32(l->color.r, l->color.g, l->color.b, 255) : IM_COL32_WHITE;
        const ImU32 fill = isDone ? IM_COL32(40, 90, 40, 255) : open ? IM_COL32(45, 45, 55, 255) : IM_COL32(28, 28, 32, 255);
        dl->AddRectFilled(b, ImVec2(b.x + bw, b.y + bh), fill, 6.0f);
        dl->AddRect(b, ImVec2(b.x + bw, b.y + bh), lc, 6.0f, 0, open ? 2.5f : 1.0f);
        const ImU32 tc = isDone ? IM_COL32(170, 240, 170, 255) : open ? IM_COL32(240, 240, 240, 255) : IM_COL32(130, 130, 130, 255);
        dl->AddText(ImVec2(b.x + 8, b.y + 5), tc, n.title.c_str());
        dl->AddText(ImVec2(b.x + 8, b.y + 24), IM_COL32(170, 170, 170, 255), isDone ? "done" : open ? TaskText(n).c_str() : "locked");
        if (ImGui::IsMouseHoveringRect(b, ImVec2(b.x + bw, b.y + bh))) {
            ImGui::BeginTooltip();
            ImGui::Text("%s  (%s)", n.title.c_str(), l ? l->name.c_str() : n.line.c_str());
            ImGui::Text("Task: %s", TaskText(n).c_str());
            ImGui::Text("Reward: %s", RewardText(n).c_str());
            std::string req;
            if (n.reqTier >= 0) req += TextFormat("tier %d", n.reqTier);
            for (const std::string& r : n.reqNodes) { const Node* p = Find(r); req += (req.empty() ? "" : ", ") + (p ? p->title : r); }
            ImGui::Text("Needs: %s", req.empty() ? "-" : req.c_str());
            ImGui::EndTooltip();
        }
    }
    ImGui::EndChild();
}
