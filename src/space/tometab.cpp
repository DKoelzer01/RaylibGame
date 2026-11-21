#include "tometab.h"

#include "gameui.h"
#include "inventory.h"
#include "magic.h"
#include "playersession.h"
#include "research.h"
#include "researchtab.h"
#include "tome.h"
#include "imgui.h"

#include <algorithm>
#include <map>
#include <vector>

namespace {
std::map<std::string, ImVec2> drawnCenters; // (tests) this frame's page boxes

ImVec4 Colour(unsigned rgb, float a = 1.0f) { return ImVec4(((rgb >> 16) & 255) / 255.0f, ((rgb >> 8) & 255) / 255.0f, (rgb & 255) / 255.0f, a); }
const char* BranchName(const std::string& b) {
    if (b == "alchemy") return "Alchemy";
    if (b == "artifice") return "Artifice";
    if (b == "wards") return "Wards and rituals";
    if (b == "attunement") return "Attunement";
    return b.c_str();
}
// A name that is wider than `room` breaks at the last space that fits (the user's rule from #1192: whole names at a
// readable size, never a smaller font or "..."). A single word wider than the box stays whole.
std::string WrapName(const std::string& name, float room) {
    if (ImGui::CalcTextSize(name.c_str()).x <= room) return name;
    size_t cut = std::string::npos;
    for (size_t sp = name.find(' '); sp != std::string::npos; sp = name.find(' ', sp + 1)) {
        if (ImGui::CalcTextSize(name.substr(0, sp).c_str()).x > room) break;
        cut = sp;
    }
    return cut == std::string::npos ? name : name.substr(0, cut) + "\n" + name.substr(cut + 1);
}
std::string CostText(const TomePage& p) {
    std::string s;
    for (const auto& [id, n] : p.cost) { const AspectDef* a = Aspects_Find(id); s += (s.empty() ? "" : ", ") + std::to_string(n) + " " + (a ? a->name : id); }
    return s.empty() ? "nothing" : s;
}
}

bool TomeTab_Shown() { return Magic_Enabled() && LocalSession().magic.found; }

bool TomeTab_TestPageCenter(const std::string& id, Vector2& center) {
    auto it = drawnCenters.find(id);
    if (it == drawnCenters.end()) return false;
    center = { it->second.x, it->second.y };
    return true;
}

void TomeTab_Draw() {
    drawnCenters.clear();
    if (Tome_Pages().empty()) Tome_Load();
    MagicPlayerState& me = LocalSession().magic;
    const int team = ResearchTab_Team();
    TomeState& tome = Tome_Team(team);
    const ResearchState& research = Research_Team(team);
    // The reader's aspects.
    ImGui::TextUnformatted("Your aspects:");
    bool any = false;
    for (const AspectDef& a : Aspects_All()) {
        const auto it = me.aspects.find(a.id);
        if (it == me.aspects.end() || it->second < 0.05f) continue;
        ImGui::SameLine();
        ImGui::TextColored(Colour(a.colour), "%s %.1f", a.name.c_str(), it->second);
        if (ImGui::IsItemHovered() && !a.latin.empty()) ImGui::SetTooltip("%s (%s)", a.name.c_str(), a.latin.c_str());
        any = true;
    }
    if (!any) { ImGui::SameLine(); ImGui::TextDisabled("none yet: study things with magic sight (hold Y on them)"); }
    ImGui::Separator();
    // The pages in four branches.
    static const char* kBranches[] = { "attunement", "alchemy", "artifice", "wards" };
    const float colW = ImGui::GetContentRegionAvail().x / 4.0f;
    const float boxH = ImGui::GetTextLineHeightWithSpacing() * 3.4f; // room for a name on two lines and the state
    ImGui::Columns(4, "tome", false);
    for (const char* branch : kBranches) {
        ImGui::TextUnformatted(BranchName(branch));
        for (const TomePage& p : Tome_Pages()) {
            if (p.branch != branch) continue;
            const bool owned = tome.Has(p.id);
            std::string why;
            const bool ready = !owned && Tome_CanUnlock(tome, me, research, p, &why);
            ImGui::PushStyleColor(ImGuiCol_Button, owned ? ImVec4(0.16f, 0.36f, 0.16f, 1.0f) : ready ? ImVec4(0.30f, 0.22f, 0.45f, 1.0f) : ImVec4(0.12f, 0.12f, 0.14f, 1.0f));
            ImGui::PushStyleColor(ImGuiCol_Text, owned ? ImVec4(0.7f, 0.95f, 0.7f, 1.0f) : ready ? ImVec4(1, 1, 1, 1) : ImVec4(0.55f, 0.55f, 0.58f, 1.0f));
            // (capture 2026-11-21: "The Transmutation Circ...") a long name wraps to a second line at a word, never shrinks
            const std::string label = WrapName(p.name, colW - 12.0f - 2.0f * ImGui::GetStyle().FramePadding.x) + "\n" + (owned ? "written" : ready ? "ready: click to write" : "locked") + "##" + p.id;
            const bool clicked = ImGui::Button(label.c_str(), ImVec2(colW - 12.0f, boxH));
            const ImVec2 lo = ImGui::GetItemRectMin(), hi = ImGui::GetItemRectMax();
            drawnCenters[p.id] = ImVec2((lo.x + hi.x) * 0.5f, (lo.y + hi.y) * 0.5f);
            ImGui::PopStyleColor(2);
            if (ImGui::IsItemHovered()) {
                ImGui::BeginTooltip();
                ImGui::Text("%s%s%s", p.name.c_str(), p.latin.empty() ? "" : "  -  ", p.latin.c_str());
                if (!p.text.empty()) ImGui::TextWrapped("%s", p.text.c_str());
                ImGui::Text("Costs: %s", CostText(p).c_str());
                if (!owned && !ready) ImGui::TextDisabled("Not yet: %s", why.c_str());
                ImGui::EndTooltip();
            }
            if (clicked && ready && Tome_Unlock(tome, me, research, p.id)) {
                Tome_ApplyEffects(tome, me);
                GameUI_Toast(TextFormat("Written in the tome: %s", p.name.c_str()));
            }
        }
        ImGui::NextColumn();
    }
    ImGui::Columns(1);
    ImGui::Separator();
    // Spending items (the research table's job).
    if (!tome.Has("research_table")) { ImGui::TextDisabled("Write The Research Table to spend things for their aspects."); return; }
    ImGui::TextUnformatted("Spend a carried item for its aspects:");
    int shown = 0;
    for (const ItemDef& d : itemDB.items()) {
        const auto* asp = Aspects_Of("item:" + d.id);
        if (!asp || inventory.count(d.id) <= 0) continue;
        std::string what;
        for (const auto& [id, n] : *asp) { const AspectDef* a = Aspects_Find(id); what += (what.empty() ? "" : ", ") + std::to_string(n) + " " + (a ? a->name : id); }
        if (shown++ % 3) ImGui::SameLine();
        if (ImGui::Button(TextFormat("%s x%d (%s)##spend_%s", d.name.c_str(), inventory.count(d.id), what.c_str(), d.id.c_str())))
            Tome_SpendItem(tome, me, inventory, d.id);
    }
    if (!shown) ImGui::TextDisabled("Nothing you carry is made of anything you can name yet.");
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1251): writing a page with a click spends the aspects and the team owns it ----
#include "unit_helpers.h"
#include "playerrecord.h"

TEST(TomeTab, AClickWritesAReadyPage) {
    UnitTest::EnsureItemDB();
    ASSERT_TRUE(Research_Load());
    ASSERT_TRUE(Tome_Load());
    ImGuiContext* prev = ImGui::GetCurrentContext();
    ImGuiContext* ctx = ImGui::CreateContext();
    ImGui::SetCurrentContext(ctx);
    ImGuiIO& io = ImGui::GetIO();
    io.DisplaySize = ImVec2(1600.0f, 900.0f); io.DeltaTime = 1.0f / 60.0f; io.IniFilename = nullptr;
    io.BackendFlags |= ImGuiBackendFlags_RendererHasTextures;
    io.Fonts->AddFontDefault();
    Teams_Reset();
    MagicPlayerState saved = LocalSession().magic;
    MagicPlayerState& me = LocalSession().magic;
    me = MagicPlayerState{};
    me.found = true;
    TomeState& tome = Tome_Team(ResearchTab_Team());
    tome.pages.insert("first_sight");
    me.aspects = { { "air", 6.0f }, { "light", 4.0f }, { "order", 3.0f } };
    auto frame = [&](Vector2 m, bool down) {
        io.AddMousePosEvent(m.x, m.y); io.AddMouseButtonEvent(0, down);
        ImGui::NewFrame();
        ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f)); ImGui::SetNextWindowSize(ImVec2(1600.0f, 900.0f));
        ImGui::Begin("game window", nullptr, ImGuiWindowFlags_NoDecoration);
        TomeTab_Draw();
        ImGui::End();
        ImGui::Render();
    };
    frame({ 5, 5 }, false); frame({ 5, 5 }, false);
    Vector2 wider{}, far{};
    ASSERT_TRUE(TomeTab_TestPageCenter("wider_sight", wider));
    ASSERT_TRUE(TomeTab_TestPageCenter("far_sight", far));
    frame(far, false); frame(far, true); frame(far, false);
    EXPECT_FALSE(Tome_Team(ResearchTab_Team()).Has("far_sight")) << "a locked page doesn't write";
    frame(wider, false); frame(wider, true); frame(wider, false);
    EXPECT_TRUE(Tome_Team(ResearchTab_Team()).Has("wider_sight")) << "a ready page writes on a click";
    EXPECT_FLOAT_EQ(me.aspects.count("air") ? me.aspects["air"] : 0.0f, 0.0f) << "its cost spent";
    EXPECT_EQ(me.sightTier, 1) << "its effect applies at once";
    LocalSession().magic = saved;
    Teams_Reset();
    ImGui::DestroyContext(ctx);
    ImGui::SetCurrentContext(prev);
}
#endif
