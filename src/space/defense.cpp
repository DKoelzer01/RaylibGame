// Combat step 4 (#1218): see defense.h.
#include "defense.h"

#include <cstdlib>

bool Combat4_On() {
    static const bool on = !(getenv("RAYLIBGAME_COMBAT4") && getenv("RAYLIBGAME_COMBAT4")[0] == '0');
    return on;
}

// The combat doc's shields (section 3, guesses).
// (#1218 playtest 2026-11-22: "could not parry") Parry windows doubled, GUESS: raising on the wind-up's telegraph (0.4-0.6 s, plus
// ~0.25 s to react) landed ~0.25 s before the strike, outside 0.15 s: an ordinary block every time. Was 0.22 / 0.15 / 0.10.
static const ShieldSpec kBuckler = { "buckler", 0.55f, 45.0f, 0.45f, 1.0f };
static const ShieldSpec kRound = { "round shield", 0.75f, 60.0f, 0.33f, 1.0f };
static const ShieldSpec kTower = { "tower shield", 0.95f, 75.0f, 0.22f, 0.7f };

const ShieldSpec* Shield_Of(const std::string& id) {
    if (id.find("buckler") != std::string::npos) return &kBuckler;
    if (id.find("tower_shield") != std::string::npos) return &kTower;
    if (id.find("shield") != std::string::npos && id.find("heat_shield") == std::string::npos) return &kRound;
    return nullptr;
}

bool Parry_InWindow(double blockStartedAt, double now, const ShieldSpec& s) {
    const double since = now - blockStartedAt;
    return since >= 0.0 && since <= s.parryWindow + 1e-6;
}

int DoubleTap::press(int dir, double now) {
    if (dir < 0 || dir > 3) return -1;
    if (now - last[dir] <= DefenseTune::kDoubleTap) { last[dir] = -1e9; return dir; }
    last[dir] = now;
    return -1;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>
TEST(Defense, Shields) {
    EXPECT_EQ(Shield_Of("stone"), nullptr);
    EXPECT_EQ(Shield_Of("heat_shield_tile"), nullptr);
    ASSERT_NE(Shield_Of("wooden_shield"), nullptr);
    EXPECT_FLOAT_EQ(Shield_Of("wooden_shield")->stops, 0.75f);  // today's
    EXPECT_FLOAT_EQ(Shield_Of("iron_buckler")->parryWindow, 0.45f);
    EXPECT_FLOAT_EQ(Shield_Of("steel_tower_shield")->raisedWalk, 0.7f);
}

TEST(Defense, ParryWindow) {
    const ShieldSpec& round = *Shield_Of("wooden_shield");
    EXPECT_TRUE(Parry_InWindow(10.0, 10.1, round));
    EXPECT_TRUE(Parry_InWindow(10.0, 10.3, round));              // (#1218: a reaction to the telegraph)
    EXPECT_FALSE(Parry_InWindow(10.0, 10.4, round));             // too early a raise: an ordinary block
    EXPECT_TRUE(Parry_InWindow(10.0, 10.4, *Shield_Of("buckler"))); // a buckler's window is wider
    EXPECT_FALSE(Parry_InWindow(10.0, 9.9, round));              // raised after the hit
}

TEST(Defense, DoubleTap) {
    DoubleTap d;
    EXPECT_EQ(d.press(2, 1.0), -1);
    EXPECT_EQ(d.press(2, 1.2), 2);   // within 0.25 s
    EXPECT_EQ(d.press(2, 1.3), -1);  // the pair was used up
    EXPECT_EQ(d.press(3, 2.0), -1);
    EXPECT_EQ(d.press(3, 2.4), -1);  // too slow
    EXPECT_EQ(d.press(0, 3.0), -1);
    EXPECT_EQ(d.press(1, 3.1), -1);  // a different direction
}
#endif
