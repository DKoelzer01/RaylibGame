#include "questtree.h"
#include "audio.h"
#include "gameui.h"
#include "inventory.h"
#include "survival.h"
#include "imgui.h"
#include "imgui_internal.h" // (#1077) BringWindowToDisplayFront

#include <raylib.h>
#include <algorithm>
#include <cstring>
#include <map>
#include <set>
#include <cfloat>
#include <cmath>
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
    // (#338) the selected quest book's side-line progress
    std::set<std::string>& Done() { return QuestBook_Selected().treeDone; }
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
        for (const std::string& r : n.reqNodes) if (!Done().count(r)) return false;
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
        if (Done().count(n.id) || !ReqsDone(n) || !QuestBook_Met(n.kind, n.target, n.n)) continue;
        Done().insert(n.id);
        // (#816) the reward goes to the completing player: here the local one, or on a server (gQuestReward) the team member
        // whose event finished it, through the server's inventory and its deltas
        for (const auto& r : n.reward) if (itemDB.find(r.first)) { if (gQuestReward) gQuestReward(r.first, r.second); else inventory.add(r.first, r.second); }
        const Line* l = LineOf(n.line);
        const char* news = TextFormat("%s: %s done  (+%s)", l ? l->name.c_str() : "Quest", n.title.c_str(), RewardText(n).c_str());
        if (gQuestNews) gQuestNews(news);
        else { GameUI_Toast(news); Audio_Play(Sfx::Craft, 0.6f, 1.2f); }
    }
}

void QuestTree_Reset() { Load(); Done().clear(); }
int QuestTree_NodeCount() { Load(); return (int)nodes.size(); }
bool QuestTree_IsDone(const std::string& id) { return Done().count(id) > 0; }
bool QuestTree_IsOpen(const std::string& id) { const Node* n = Find(id); return n && !Done().count(id) && ReqsDone(*n); }

void QuestTree_Save(FILE* f) { for (const std::string& id : Done()) fprintf(f, "side_done %s\n", id.c_str()); }
bool QuestTree_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "side_done") != 0) return false;
    char id[64] = {};
    if (fscanf(f, " %63s", id) == 1) Done().insert(id);
    return true;
}

static int openTier = -1; // (#890) the tier whose objectives are open in the pop-out (-1: none)
static int focusedTier = -1; // (#1077) the tier the pop-out was last focused for (-1 after a close: a reopen focuses it again)
void QuestTree_TestOpenTier(int tier) { openTier = tier; } // (shots and tests)
bool QuestTree_ClosePopout() { if (openTier < 0) return false; openTier = -1; return true; }

void QuestTree_DrawTab() {
    Load();
    // (#891) the boxes in font heights (they were sized for a 13 px font: at 2x the text ran out of them)
    const float fh = ImGui::GetFontSize();
    const float cw = fh * 17.8f, ch = fh * 5.7f, bw = fh * 16.2f, bh = fh * 3.6f, pad = fh;
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
    // Edges first. (#914 review) Routed along the gaps between rows and columns, never across a box: out of the parent into
    // a gap, along the column gap left of the child to the row gap above it, across to the child's top centre.
    const float gx = (cw - bw) * 0.5f, gy = (ch - bh) * 0.5f;
    auto route = [&](ImVec2 from, float alongX, ImVec2 b, ImU32 col, float thick) {
        const ImVec2 to(b.x + bw * 0.5f, b.y);
        const float yc = b.y - gy, xc = b.x - gx;
        std::vector<ImVec2> pts{ from, ImVec2(alongX, from.y) };
        if (std::fabs(from.y - yc) > 0.5f) { pts.push_back(ImVec2(xc, from.y)); pts.push_back(ImVec2(xc, yc)); }
        pts.push_back(ImVec2(to.x, yc));
        pts.push_back(to);
        dl->AddPolyline(pts.data(), (int)pts.size(), col, ImDrawFlags_None, thick);
    };
    for (const Node& n : nodes) {
        const ImVec2 b = boxAt(n.x, n.y);
        if (n.reqTier >= 0) { // from the tier box's right side into the gap beside it
            const ImVec2 t = boxAt(0.0f, (float)n.reqTier);
            const ImVec2 s(t.x + bw, t.y + bh * 0.5f);
            route(s, s.x + gx, b, IM_COL32(150, 150, 150, 120), 1.5f);
        }
        for (const std::string& r : n.reqNodes)
            if (const Node* p = Find(r)) { // from the parent's bottom centre down into the row gap below it
                const ImVec2 pb = boxAt(p->x, p->y);
                const ImVec2 s(pb.x + bw * 0.5f, pb.y + bh + gy);
                dl->AddLine(ImVec2(s.x, pb.y + bh), s, Done().count(r) ? IM_COL32(120, 210, 120, 200) : IM_COL32(150, 150, 150, 120), 2.0f);
                route(s, s.x, b, Done().count(r) ? IM_COL32(120, 210, 120, 200) : IM_COL32(150, 150, 150, 120), 2.0f);
            }
    }
    // (#914 review) A node's text fits its box: whole; else shrunk to 70%; else cut at a word with "..." (never mid-word;
    // the tooltip has it all). The UI font is monospaced and its glyph advances snap to whole pixels, so a size is chosen
    // by whole-pixel advances (a 4% shrink drew exactly as wide as the full size and still ran out of the box).
    auto fitText = [&](ImVec2 at, float room, ImU32 col, const std::string& s) {
        ImFont* font = ImGui::GetFont();
        const float fd = ImGui::GetFontSize();
        const int len = std::max(1, (int)s.size());
        const float adv = font->CalcTextSizeA(fd, FLT_MAX, 0.0f, s.c_str()).x / (float)len; // one glyph at fd
        auto drawn = [&](float f, int n) { return std::ceil(adv * f / fd) * (float)n; };     // width at size f, snapped
        if (drawn(fd, len) <= room) { dl->AddText(font, fd, at, col, s.c_str()); return; }
        const float px = std::floor(room / (float)len);                                       // whole pixels a glyph may take
        if (px >= std::ceil(adv * 0.7f)) {
            const float f = fd * px / adv;
            dl->AddText(font, f, ImVec2(at.x, at.y + (fd - f) * 0.5f), col, s.c_str());
            return;
        }
        const float f = fd * 0.7f;
        std::string cut = s;
        while (drawn(f, (int)cut.size() + 3) > room) {
            const size_t sp = cut.find_last_of(' ');
            if (sp == std::string::npos) break;
            cut.resize(sp);
            while (!cut.empty() && (cut.back() == ' ' || cut.back() == ',' || cut.back() == ':' || cut.back() == '(')) cut.pop_back();
        }
        cut += "...";
        dl->AddText(font, f, ImVec2(at.x, at.y + (fd - f) * 0.5f), col, cut.c_str());
    };
    // The main line: one box per tier (done / now / locked). A click opens its objectives in a pop-out beside it (#890).
    bool tierClicked = false;
    for (int t = 0; t < QuestBook_TierCount(); ++t) {
        const ImVec2 b = boxAt(0.0f, (float)t);
        const ImVec2 e(b.x + bw, b.y + bh);
        const bool hover = ImGui::IsMouseHoveringRect(b, e) && ImGui::IsWindowHovered();
        const ImU32 fill = t < cur ? IM_COL32(40, 90, 40, 255) : t == cur ? IM_COL32(70, 70, 30, 255) : IM_COL32(35, 35, 40, 255);
        dl->AddRectFilled(b, e, hover ? IM_COL32(80, 80, 60, 255) : fill, 6.0f);
        dl->AddRect(b, e, IM_COL32(230, 210, 140, 255), 6.0f, 0, t == cur || t == openTier ? 2.5f : 1.0f);
        dl->AddText(ImVec2(b.x + 8, b.y + fh * 0.35f), IM_COL32(240, 230, 200, 255), TextFormat("Tier %d", t));
        const char* state = t < cur ? "done" : t == cur ? "now" : "locked";
        dl->AddText(ImVec2(e.x - ImGui::CalcTextSize(state).x - 8, b.y + fh * 0.35f), t < cur ? IM_COL32(150, 230, 150, 255) : t == cur ? IM_COL32(240, 220, 120, 255) : IM_COL32(140, 140, 140, 255), state);
        dl->AddText(ImVec2(b.x + 8, b.y + fh * 1.85f), IM_COL32(220, 220, 220, 255), QuestBook_TierName(t));
        if (hover && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            openTier = openTier == t ? -1 : t; // the same tier again closes it
            tierClicked = true;
        }
    }
    // Side nodes.
    for (const Node& n : nodes) {
        const ImVec2 b = boxAt(n.x, n.y);
        const bool isDone = Done().count(n.id) > 0, open = !isDone && ReqsDone(n);
        const Line* l = LineOf(n.line);
        const ImU32 lc = l ? IM_COL32(l->color.r, l->color.g, l->color.b, 255) : IM_COL32_WHITE;
        const ImU32 fill = isDone ? IM_COL32(40, 90, 40, 255) : open ? IM_COL32(45, 45, 55, 255) : IM_COL32(28, 28, 32, 255);
        dl->AddRectFilled(b, ImVec2(b.x + bw, b.y + bh), fill, 6.0f);
        dl->AddRect(b, ImVec2(b.x + bw, b.y + bh), lc, 6.0f, 0, open ? 2.5f : 1.0f);
        const ImU32 tc = isDone ? IM_COL32(170, 240, 170, 255) : open ? IM_COL32(240, 240, 240, 255) : IM_COL32(130, 130, 130, 255);
        fitText(ImVec2(b.x + 8, b.y + fh * 0.35f), bw - 16.0f, tc, n.title);
        fitText(ImVec2(b.x + 8, b.y + fh * 1.85f), bw - 16.0f, IM_COL32(170, 170, 170, 255), isDone ? std::string("done") : open ? TaskText(n) : std::string("locked"));
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
    // (#890) The open tier's objectives: a semi-transparent pop-out beside its box, each with its progress and a pin box
    // (up to kMaxPins pinned to the HUD, the player's own). Esc, a click outside it or the tier again closes it.
    if (openTier >= 0 && openTier < QuestBook_TierCount()) {
        { const ImVec2 b = boxAt(0.0f, (float)openTier); ImGui::SetNextWindowPos(ImVec2(b.x + bw + 8.0f, b.y), ImGuiCond_Always); } // (beside its box, wherever the tree is scrolled)
        ImGui::SetNextWindowBgAlpha(0.93f); // (semi-transparent, but the tree behind must not show through the text)
        if (focusedTier != openTier) { focusedTier = openTier; ImGui::SetNextWindowFocus(); } // (#1077) focus it when a tier opens (it opened behind the Tab window)
        bool keep = true;
        const ImGuiWindowFlags f = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize |
                                   ImGuiWindowFlags_NoSavedSettings;
        if (ImGui::Begin("##tierobjectives", nullptr, f)) {
            // (#1077) and keep it in front of the tree it belongs to: a click on the tree (or anywhere in the Tab window) would
            // otherwise bring that window over it again
            ImGui::BringWindowToDisplayFront(ImGui::GetCurrentWindow());
            const int t = openTier;
            ImGui::TextColored(ImVec4(0.95f, 0.85f, 0.55f, 1.0f), "Tier %d: %s", t, QuestBook_TierName(t));
            ImGui::PushTextWrapPos(ImGui::GetFontSize() * 26.0f);
            ImGui::TextDisabled("%s", QuestBook_TierText(t));
            ImGui::PopTextWrapPos();
            ImGui::Separator();
            ImGui::PushTextWrapPos(ImGui::GetFontSize() * 30.0f); // (long objectives wrap: one ran off the screen)
            for (int k = 0; k < QuestBook_TaskCount(t); ++k) {
                const bool kd = QuestBook_TaskDone(t, k);
                ImGui::PushID(k);
                if (!kd && t <= cur) {
                    bool pinned = QuestBook_Pinned(t, k);
                    const bool full = !pinned && (int)QuestBook_Pins().size() >= kMaxPins;
                    ImGui::BeginDisabled(full);
                    if (ImGui::Checkbox("##pin", &pinned)) QuestBook_Pin(t, k, pinned);
                    ImGui::EndDisabled();
                    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) ImGui::SetTooltip(full ? "Three are pinned: unpin one first" : pinned ? "Pinned to the HUD" : "Pin to the HUD");
                    ImGui::SameLine();
                } else { ImGui::Dummy(ImVec2(ImGui::GetFrameHeight(), 0.0f)); ImGui::SameLine(); }
                ImGui::TextColored(kd ? ImVec4(0.45f, 0.85f, 0.45f, 1.0f) : t <= cur ? ImVec4(0.95f, 0.95f, 0.95f, 1.0f) : ImVec4(0.6f, 0.6f, 0.62f, 1.0f),
                                   "%s%s", QuestBook_TaskText(t, k), kd ? "   (done)" : "");
                ImGui::PopID();
            }
            ImGui::PopTextWrapPos();
            if (t > cur) ImGui::TextDisabled("Locked: finish tier %d first", t - 1);
            const bool outside = !ImGui::IsWindowHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem) && ImGui::IsMouseClicked(ImGuiMouseButton_Left) && !tierClicked;
            if (outside || ImGui::IsKeyPressed(ImGuiKey_Escape)) keep = false;
        }
        ImGui::End();
        if (!keep) openTier = -1;
    } else focusedTier = -1;
}
