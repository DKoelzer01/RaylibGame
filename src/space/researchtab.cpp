#include "researchtab.h"
#include "gameui.h"
#include "menus.h"
#include "net.h"
#include "playerrecord.h"
#include "research.h"
#include "researchgate.h"
#include "survival.h"
#include "tome.h"          // (#1251) tome pages per team
#include "playersession.h" // the local player's magic state (found)
#include "imgui.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <map>

namespace {
std::map<std::string, ImVec2> drawnCenters; // (tests) this frame's node boxes
float zoom = 1.0f;

const char* TierName(int t) {
    static const char* kNames[] = { "Steam and steel", "LV", "MV", "HV", "Rocketry" };
    return t >= 3 && t <= 7 ? kNames[t - 3] : "";
}
// A name that fits `room` px at `size`: whole, else cut at a word with "..." (the tooltip has it all).
std::string Fit(const std::string& s, float size, float room) {
    ImFont* font = ImGui::GetFont();
    if (font->CalcTextSizeA(size, FLT_MAX, 0.0f, s.c_str()).x <= room) return s;
    std::string cut = s;
    while (!cut.empty() && font->CalcTextSizeA(size, FLT_MAX, 0.0f, (cut + "...").c_str()).x > room) {
        const size_t sp = cut.find_last_of(' ');
        cut.resize(sp == std::string::npos || sp == 0 ? cut.size() - 1 : sp);
    }
    return cut + "...";
}
ImU32 LineColour(const std::string& line) {
    if (line == "combat") return IM_COL32(220, 110, 100, 255);
    if (line == "logistics") return IM_COL32(110, 180, 230, 255);
    if (line == "efficiency") return IM_COL32(120, 210, 140, 255);
    if (line == "magic") return IM_COL32(190, 130, 230, 255);
    if (line == "mining") return IM_COL32(210, 170, 100, 255);
    if (line == "exploration") return IM_COL32(230, 220, 120, 255);
    return IM_COL32(200, 200, 200, 255);
}
}

// The team the tab shows: the player's (single player makes a team of one the first time, as claims do: buildtool's LocalTeam).
int ResearchTab_Team() {
    int team = Teams_Of(PlayerId_Local());
    if (!team) team = Teams_CreateSolo(PlayerId_Local(), settings.mpName);
    return team;
}

bool ResearchTab_TestNodeCenter(const std::string& id, Vector2& center) {
    auto it = drawnCenters.find(id);
    if (it == drawnCenters.end()) return false;
    center = { it->second.x, it->second.y };
    return true;
}

void ResearchTab_Draw() {
    drawnCenters.clear();
    if (Research_Nodes().empty()) Research_Load();
    const bool joined = gNetRole == NetRole::Client; // (the team's research is the server's: synced in a later step)
    const int team = ResearchTab_Team();
    ResearchState& st = Research_Team(team);
    const int tier = QuestBook_Tier();
    ImGui::TextWrapped("Research from tier 3: labs turn research packs into the nodes below; a tier's core nodes are part of its milestone.");
    ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
    ImGui::TextWrapped("%s%s", joined ? "Joined: your team's research is the server's (shown here in a later step). "
                                      : ResearchGate_On() ? "Machines, recipes and machine recipes unlock with their node; labs research the queue from the top (I at a lab hands in what the milestone wants). Click a node to queue it. "
                                      : "Research gating is off (RAYLIBGAME_LADDER=0): nothing is locked. Click a node to queue it. ", // (#1271: the old "not active yet" note)
                       joined ? "" : "Drag to move, the wheel scrolls, Ctrl + wheel zooms.");
    ImGui::PopStyleColor();

    const float avail = ImGui::GetContentRegionAvail().x;
    const float sideW = std::clamp(avail * 0.28f, 220.0f, 420.0f);
    // ---- The graph ----
    ImGui::BeginChild("researchgraph", ImVec2(avail - sideW - 8.0f, 0), true, ImGuiWindowFlags_HorizontalScrollbar);
    {
        const float fh = ImGui::GetFontSize() * zoom;
        const float bw = fh * 13.0f, bh = fh * 3.2f, gapX = fh * 1.2f, gapY = fh * 0.8f, colW = 2.0f * bw + 3.0f * gapX;
        ImDrawList* dl = ImGui::GetWindowDrawList();
        const ImVec2 o = ImGui::GetCursorScreenPos();
        // Place: per tier column, core nodes down the left (the spine), branch nodes down the right in line order.
        std::map<std::string, ImVec2> at; // top-left
        int rows = 0;
        for (int t = 3; t <= 7; ++t) {
            int coreRow = 0, branchRow = 0;
            std::vector<const ResearchNode*> branch;
            for (const ResearchNode& n : Research_Nodes()) if (n.tier == t) { if (n.core) at[n.id] = ImVec2(o.x + (t - 3) * colW + gapX, o.y + fh * 2.0f + coreRow++ * (bh + gapY)); else branch.push_back(&n); }
            std::stable_sort(branch.begin(), branch.end(), [](const ResearchNode* a, const ResearchNode* b) { return a->line < b->line; });
            for (const ResearchNode* n : branch) at[n->id] = ImVec2(o.x + (t - 3) * colW + 2.0f * gapX + bw, o.y + fh * 2.0f + branchRow++ * (bh + gapY));
            rows = std::max(rows, std::max(coreRow, branchRow));
        }
        // The canvas: a press in it is the graph's (#1078: a Dummy let ImGui drag the Tab window from it).
        ImGui::InvisibleButton("canvas", ImVec2(5.0f * colW, fh * 2.0f + rows * (bh + gapY) + gapY), ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonMiddle);
        const bool hovered = ImGui::IsWindowHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);
        // Tier headings: current tier and the one after it by name; later ones hidden.
        for (int t = 3; t <= 7; ++t) {
            const ImVec2 h(o.x + (t - 3) * colW + gapX, o.y + fh * 0.4f);
            const bool shown = t <= tier + 1;
            dl->AddText(ImGui::GetFont(), fh, h, t <= tier ? IM_COL32(240, 220, 150, 255) : IM_COL32(150, 150, 160, 255),
                        shown ? TextFormat("Tier %d: %s", t, TierName(t)) : TextFormat("Tier %d", t));
            if (!shown) dl->AddText(ImGui::GetFont(), fh * 0.85f, ImVec2(h.x, h.y + fh * 1.8f), IM_COL32(110, 110, 120, 255), "(reach the tier before it)");
        }
        // Lines to the prerequisites (only between drawn nodes).
        for (const ResearchNode& n : Research_Nodes()) {
            if (Research_Sight(st, n, tier) == ResearchSight::Hidden) continue;
            for (const std::string& d : n.needs) {
                const ResearchNode* p = Research_Find(d);
                if (!p || Research_Sight(st, *p, tier) == ResearchSight::Hidden || !at.count(d)) continue;
                const ImU32 col = st.Done(d) ? IM_COL32(120, 210, 120, 200) : IM_COL32(150, 150, 150, 110);
                if (at[d].x == at[n.id].x) { // the same column (the spine): along its left side, never across a box
                    const float x = at[d].x - gapX * 0.45f, ya = at[d].y + bh * 0.5f, yb = at[n.id].y + bh * 0.5f;
                    const ImVec2 pts[4] = { ImVec2(at[d].x, ya), ImVec2(x, ya), ImVec2(x, yb), ImVec2(at[n.id].x, yb) };
                    dl->AddPolyline(pts, 4, col, ImDrawFlags_None, 1.5f);
                } else dl->AddLine(ImVec2(at[d].x + bw, at[d].y + bh * 0.5f), ImVec2(at[n.id].x, at[n.id].y + bh * 0.5f), col, 1.5f);
            }
        }
        // The nodes.
        for (const ResearchNode& n : Research_Nodes()) {
            const ResearchSight sight = Research_Sight(st, n, tier);
            if (sight == ResearchSight::Hidden) continue;
            const ImVec2 p = at[n.id], q(p.x + bw, p.y + bh);
            drawnCenters[n.id] = ImVec2((p.x + q.x) * 0.5f, (p.y + q.y) * 0.5f);
            const bool done = st.Done(n.id), queued = std::find(st.queue.begin(), st.queue.end(), n.id) != st.queue.end();
            const bool hover = hovered && ImGui::IsMouseHoveringRect(p, q);
            const ImU32 fill = done ? IM_COL32(40, 90, 40, 255) : sight == ResearchSight::NameOnly ? IM_COL32(26, 26, 30, 255) : queued ? IM_COL32(70, 64, 30, 255) : IM_COL32(40, 42, 52, 255);
            dl->AddRectFilled(p, q, hover ? IM_COL32(80, 80, 70, 255) : fill, 6.0f);
            dl->AddRect(p, q, n.core ? IM_COL32(230, 210, 140, 255) : LineColour(n.line), 6.0f, 0, n.core ? 2.0f : 1.2f);
            // The text stays inside its box (a long name is cut there; the tooltip has it all); queued shows as the box's colour.
            dl->PushClipRect(ImVec2(p.x + 4, p.y), ImVec2(q.x - 6, q.y), true);
            dl->AddText(ImGui::GetFont(), fh, ImVec2(p.x + 8, p.y + fh * 0.3f), sight == ResearchSight::NameOnly ? IM_COL32(140, 140, 150, 255) : IM_COL32(235, 235, 235, 255), Fit(n.name, fh, bw - 16.0f).c_str()); // (cut at a word, #1195 review)
            if (sight == ResearchSight::Full) {
                const int units = n.units, have = st.progress.count(n.id) ? st.progress.at(n.id) : 0;
                const char* sub = done ? "done" : have ? TextFormat("%d / %d units, %.0f s each", have, units, n.unitSeconds) : TextFormat("%d units, %.0f s each", units, n.unitSeconds);
                dl->AddText(ImGui::GetFont(), fh * 0.8f, ImVec2(p.x + 8, p.y + fh * 1.7f), IM_COL32(170, 170, 170, 255), sub);
            } else dl->AddText(ImGui::GetFont(), fh * 0.8f, ImVec2(p.x + 8, p.y + fh * 1.7f), IM_COL32(110, 110, 120, 255), "next tier");
            dl->PopClipRect();
            if (hover) {
                ImGui::BeginTooltip();
                ImGui::Text("%s  (tier %d, %s)", n.name.c_str(), n.tier, n.core ? "core: part of the milestone" : n.line.c_str());
                if (sight == ResearchSight::Full) {
                    std::string packs; for (const std::string& k : n.packs) packs += (packs.empty() ? "" : ", ") + k;
                    ImGui::Text("Cost: %d units, each %s and %.0f s in a lab", n.units, packs.c_str(), n.unitSeconds);
                    std::string needs; for (const std::string& d : n.needs) { const ResearchNode* r = Research_Find(d); needs += (needs.empty() ? "" : ", ") + (r ? r->name : d); }
                    ImGui::Text("Needs: %s", needs.empty() ? "-" : needs.c_str());
                    if (!n.tome.empty()) { // (#1251) the magic cross link: the team's tome page, or a hint before magic is found
                        const TomePage* page = Tome_Find(n.tome);
                        if (Tome_Team(ResearchTab_Team()).Has(n.tome)) ImGui::Text("Tome page %s: your team has it", page ? page->name.c_str() : n.tome.c_str());
                        else if (page && LocalSession().magic.found) ImGui::Text("Also needs the tome page: %s", page->name.c_str());
                        else ImGui::TextDisabled("Locked: something in the world resonates. Find it first.");
                    }
                    if (!done && !joined) ImGui::TextDisabled(queued ? "Click: take it out of the queue" : "Click: queue it");
                } else ImGui::TextDisabled("Its cost shows once your team reaches tier %d", n.tier);
                ImGui::EndTooltip();
            }
            if (hover && ImGui::IsMouseClicked(ImGuiMouseButton_Left) && sight == ResearchSight::Full && !done && !joined) {
                if (queued) st.queue.erase(std::remove(st.queue.begin(), st.queue.end(), n.id), st.queue.end());
                else st.queue.push_back(n.id);
            }
        }
        // Pan (middle button, or the left from an empty spot) and Ctrl + wheel zoom around the cursor, as the quest tree (#1078).
        static bool panning = false;
        const ImGuiIO& io = ImGui::GetIO();
        bool overNode = false;
        for (const auto& [id, c] : drawnCenters) overNode = overNode || (fabsf(io.MousePos.x - c.x) < bw * 0.5f && fabsf(io.MousePos.y - c.y) < bh * 0.5f);
        if (hovered && ((ImGui::IsMouseClicked(ImGuiMouseButton_Left) && !overNode) || ImGui::IsMouseClicked(ImGuiMouseButton_Middle))) panning = true;
        if (!ImGui::IsMouseDown(ImGuiMouseButton_Left) && !ImGui::IsMouseDown(ImGuiMouseButton_Middle)) panning = false;
        if (panning) { ImGui::SetScrollX(ImGui::GetScrollX() - io.MouseDelta.x); ImGui::SetScrollY(ImGui::GetScrollY() - io.MouseDelta.y); }
        if (hovered && io.KeyCtrl && io.MouseWheel != 0.0f) {
            const float z1 = std::clamp(zoom * (io.MouseWheel > 0.0f ? 1.15f : 1.0f / 1.15f), 0.5f, 2.0f);
            const ImVec2 base(o.x + ImGui::GetScrollX(), o.y + ImGui::GetScrollY());
            const ImVec2 u((io.MousePos.x - o.x) / zoom, (io.MousePos.y - o.y) / zoom);
            ImGui::SetScrollX(std::max(0.0f, base.x - io.MousePos.x + u.x * z1));
            ImGui::SetScrollY(std::max(0.0f, base.y - io.MousePos.y + u.y * z1));
            zoom = z1;
        }
    }
    ImGui::EndChild();
    ImGui::SameLine();
    // ---- The queue and the milestone checklist ----
    ImGui::BeginChild("researchside", ImVec2(0, 0), true);
    ImGui::PushTextWrapPos(0.0f); // (the panel is narrow: every line wraps in it)
    ImGui::TextUnformatted("Queue");
    ImGui::Separator();
    if (st.queue.empty()) ImGui::TextDisabled("Empty: click nodes to queue them (labs take the first they have packs for)");
    for (size_t i = 0; i < st.queue.size(); ++i) {
        ImGui::PushID((int)i);
        const ResearchNode* n = Research_Find(st.queue[i]);
        ImGui::Text("%zu. %s", i + 1, n ? n->name.c_str() : st.queue[i].c_str());
        if (!joined) { // (on their own line: beside a long name they were pushed out of the panel)
            if (i > 0) { if (ImGui::SmallButton("Up")) std::swap(st.queue[i], st.queue[i - 1]); ImGui::SameLine(); }
            if (ImGui::SmallButton("Remove")) { st.queue.erase(st.queue.begin() + (long)i); ImGui::PopID(); break; }
        }
        ImGui::PopID();
    }
    ImGui::Spacing();
    const int mt = std::clamp(tier < 3 ? 3 : tier, 3, 7); // the milestone the team works on (research from tier 3)
    ImGui::Text("Tier %d milestone", mt);
    ImGui::Separator();
    int core = 0, coreDone = 0;
    for (const ResearchNode& n : Research_Nodes()) if (n.tier == mt && n.core) { ++core; coreDone += st.Done(n.id); }
    ImGui::Text("Core research: %d / %d", coreDone, core);
    // (#1271 L6b-2) the milestone's lines, in any order: deliveries (handed in at a lab: I), a working build, finds
    if (const TierMilestone* m = ResearchGate_Milestone(mt))
        for (int l = 0; l < (int)m->lines.size(); ++l) {
            const bool done = st.milestoneDone.count({ mt, l }) > 0;
            const std::string text = m->lines[l].text.empty() ? ResearchGate_LineLabel(m->lines[l]) : m->lines[l].text;
            if (done) ImGui::TextColored(ImVec4(0.55f, 0.85f, 0.55f, 1.0f), "[x] %s", text.c_str());
            else ImGui::Text("[ ] %s", text.c_str());
        }
    ImGui::TextDisabled("Deliveries and finds: look at a lab and press I with them in your inventory (a ruin's schematic stands in for any find for now).");
    ImGui::PopTextWrapPos();
    ImGui::EndChild();
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1192): the tab's clicks (a headless ImGui context, the mouse injected) ----
#include "unit_helpers.h"

TEST(ResearchTab, AClickQueuesAndUnqueuesOnlyWhatTheTeamMaySee) {
    ImGuiContext* prev = ImGui::GetCurrentContext();
    ImGuiContext* ctx = ImGui::CreateContext();
    ImGui::SetCurrentContext(ctx);
    ImGuiIO& io = ImGui::GetIO();
    io.DisplaySize = ImVec2(1600.0f, 900.0f); io.DeltaTime = 1.0f / 60.0f; io.IniFilename = nullptr;
    io.BackendFlags |= ImGuiBackendFlags_RendererHasTextures;
    io.Fonts->AddFontDefault();
    ASSERT_TRUE(Research_Load());
    Teams_Reset();
    const int tier0 = QuestBook_Tier();
    QuestBook_SetTier(3);
    auto frame = [&](Vector2 m, bool down) {
        io.AddMousePosEvent(m.x, m.y); io.AddMouseButtonEvent(0, down);
        ImGui::NewFrame();
        ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f)); ImGui::SetNextWindowSize(ImVec2(1600.0f, 900.0f));
        ImGui::Begin("game window", nullptr, ImGuiWindowFlags_NoDecoration);
        ResearchTab_Draw();
        ImGui::End();
        ImGui::Render();
    };
    auto click = [&](Vector2 m) { frame(m, false); frame(m, true); frame(m, false); };
    frame({ 5, 5 }, false); frame({ 5, 5 }, false);
    const ResearchState& st = Research_Team(ResearchTab_Team());
    Vector2 steam{}, circuits{}, aluminium{};
    ASSERT_TRUE(ResearchTab_TestNodeCenter("steam_power", steam)) << "tier 3, the team's tier: drawn";
    ASSERT_TRUE(ResearchTab_TestNodeCenter("basic_circuits", circuits)) << "tier 4, the next: drawn (names only)";
    EXPECT_FALSE(ResearchTab_TestNodeCenter("aluminium", aluminium)) << "tier 5: hidden";
    click(steam);
    EXPECT_EQ(st.queue, std::vector<std::string>{ "steam_power" }) << "a click queues a node of the team's tier";
    click(circuits);
    EXPECT_EQ(st.queue.size(), 1u) << "a next-tier node (names only) can't be queued";
    click(steam);
    EXPECT_TRUE(st.queue.empty()) << "a second click takes it out";
    QuestBook_SetTier(tier0);
    Teams_Reset();
    ImGui::DestroyContext(ctx);
    ImGui::SetCurrentContext(prev);
}
#endif
