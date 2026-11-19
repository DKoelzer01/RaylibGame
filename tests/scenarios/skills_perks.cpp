// skills/perks-act (#919): every family of slice 2a's perks changes the number it says, measured through the game's own
// code (not the stat table): with no perks, then with ranks set on the local player's skills.
//   Athletics: walk speed (Suit_Apply) and the safe landing speed (a 9 m/s landing hurts without, not with);
//   Woodcutting / Foraging / Farming: yields through SkillScaled (the mean of many rolls);
//   Digging / Woodcutting: tool wear (SkillPlay_WearStrokes on a shovel and an axe);
//   Building: the hold to take a part down (BuildTool_BreakHold);
//   Survival: hunger over a minute (Survival_StepFor on a copy of the player's survival state).
#include "scenario.h"
#include "buildtool.h"
#include "playersession.h"
#include "playerinput.h"
#include "skillplay.h"
#include "skills.h"
#include "suit.h"
#include "survival.h"
#include "terraintool.h"
#include <string>

namespace {
double MeanScaled(int n, const char* stat, int rolls = 2000) { double s = 0.0; for (int k = 0; k < rolls; ++k) s += SkillScaled(n, stat); return s / rolls; }
float HungerOverAMinute(ScenarioCtx& c) {
    SurvivalState s = LocalSession().survival; // (a copy: the player's own state is untouched)
    s.food = 80.0f; s.dead = false;
    for (int k = 0; k < 60; ++k) Survival_StepFor(c.world, 1.0f, PlayerInput::None(), s, LocalSession().player, LocalSession().inventory, Suit_Stats(), true);
    return 80.0f - s.food;
}
}

SCENARIO("skills/perks-act") {
    s.tags = { "skills" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 60.0;
    s.record = false; // (numbers only: nothing to watch)
    s.step("each perk family changes its number", 20.0, [](ScenarioCtx& c) {
        c.waitingFor("the player standing in the world");
        if (c.T < 1.0) return false;
        terrainTool.creative = false;
        SkillState& sk = LocalSession().skills;
        const SkillState saved = sk;
        sk.ranks.clear();
        // Baselines: no perks.
        Suit_Apply(c.world, 0.0f);
        const float walk0 = LocalSession().player.walkSpeed;
        SurvivalState fall = LocalSession().survival; fall.health = 100.0f; fall.dead = false;
        Survival_OnLandingFor(fall, Suit_Stats(), 9.0f, true);
        const float fallHurt0 = 100.0f - fall.health;
        const double logs0 = MeanScaled(10, "log_yield");
        const int shovel0 = SkillPlay_WearStrokes("stone_shovel", 100), axe0 = SkillPlay_WearStrokes("stone_axe", 100);
        const float hold0 = BuildTool_BreakHold();
        const float hunger0 = HungerOverAMinute(c);
        // With perks.
        sk.ranks = { { "long_stride", 3 }, { "marathon", 1 }, { "light_landing", 2 }, { "cats_feet", 1 }, { "lumberjack", 3 }, { "timber", 1 },
                     { "nimble_fingers", 3 }, { "green_thumb", 3 }, { "steady_haft", 2 }, { "oiled_haft", 2 }, { "quick_hands", 3 },
                     { "iron_stomach", 3 }, { "hardy", 1 } };
        Suit_Apply(c.world, 0.0f);
        c.near("walk speed x (Long Stride 3, Marathon)", LocalSession().player.walkSpeed / walk0, 1.25, 0.01);
        fall = LocalSession().survival; fall.health = 100.0f; fall.dead = false;
        Survival_OnLandingFor(fall, Suit_Stats(), 9.0f, true);
        c.gt("a 9 m/s landing hurts without perks (health lost)", fallHurt0, 15.0);
        c.eq("... and not with Light Landing 2 and Cat's Feet (health lost)", 100.0f - fall.health, 0.0);
        c.near("logs per 10 (Lumberjack 3, Timber: x1.55)", MeanScaled(10, "log_yield"), 15.5, 0.3);
        c.near("logs per 10 without perks", logs0, 10.0, 0.01);
        c.near("fibre per 10 (Nimble Fingers 3: x1.45)", MeanScaled(10, "fibre_yield"), 14.5, 0.3);
        c.near("crop per 10 (Green Thumb 3: x1.3)", MeanScaled(10, "crop_yield"), 13.0, 0.3);
        c.eq("shovel wear for 100 strokes without perks", shovel0, 100);
        c.near("... with Steady Haft 2 (x0.8)", SkillPlay_WearStrokes("stone_shovel", 100), 80, 1);
        c.eq("axe wear for 100 strokes without perks", axe0, 100);
        c.near("... with Oiled Haft 2 (x0.8)", SkillPlay_WearStrokes("stone_axe", 100), 80, 1);
        c.near("break hold x (Quick Hands 3: 1 / 1.45)", BuildTool_BreakHold() / hold0, 1.0 / 1.45, 0.01);
        const float hunger1 = HungerOverAMinute(c);
        c.gt("hunger over a minute without perks", hunger0, 0.0);
        c.near("hunger x (Iron Stomach 3, Hardy: 0.55)", hunger0 > 0.0f ? hunger1 / hunger0 : 0.0f, 0.55, 0.02);
        sk = saved;
        Suit_Apply(c.world, 0.0f);
        return true;
    });
}
