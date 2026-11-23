// ui/quest-tier-popup (#1558): the quest tree's tier pop-out sat beside its tier's box at the box's top, so a tier near the
// bottom of the screen ran it off the bottom (playtest 2026-11-22). Each tier's pop-out is opened in turn on the Quests
// screen; its window must lie inside the display (it moves up, or to the box's left, to stay there).
#include "scenario.h"
#include "gameui.h"
#include "survival.h"
#include "questtree.h"
#include "world.h"
#include "imgui.h"
#include "imgui_internal.h"

SCENARIO("ui/quest-tier-popup") {
    s.tags = { "world", "ui" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.step("settle", 30.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T >= 1.0 && WorldIdle(); });
    s.step("every tier's pop-out on screen", 60.0, [](ScenarioCtx& c) {
        static int tier = -1;
        static double at = 0.0;
        const int count = QuestBook_TierCount();
        if (tier < 0) { GameUI_ShowTab("Quests"); tier = 0; at = c.t; QuestTree_TestOpenTier(tier); }
        if (c.t - at < 0.4) return false; // (it sizes itself, then places itself by that size: a few frames)
        const ImVec2 screen = ImGui::GetIO().DisplaySize;
        if (ImGuiWindow* w = ImGui::FindWindowByName("##tierobjectives")) {
            const bool inside = w->Pos.x >= 0.0f && w->Pos.y >= 0.0f && w->Pos.x + w->Size.x <= screen.x + 0.5f && w->Pos.y + w->Size.y <= screen.y + 0.5f;
            c.check(TextFormat("tier %d's pop-out inside the %.0fx%.0f display", tier, screen.x, screen.y), inside,
                    TextFormat("at (%.0f, %.0f), %.0fx%.0f", w->Pos.x, w->Pos.y, w->Size.x, w->Size.y));
            if (tier == count - 1) c.grab("last-tier"); // (RAYLIBGAME_SLIDE_SHOTS: a picture to look at)
        } else c.check(TextFormat("tier %d's pop-out opened", tier), false);
        if (++tier >= count) { QuestTree_TestOpenTier(-1); return true; }
        QuestTree_TestOpenTier(tier);
        at = c.t;
        return false;
    });
}
