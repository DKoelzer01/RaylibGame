#include "burden.h"

#include "magicsight.h"

#include <algorithm>

// All numbers are the doc's guesses (section 6), to tune after playtests.

BurdenBand Burden_Band(float b) {
    if (b >= 80.0f) return BurdenBand::Overburdened;
    if (b >= 60.0f) return BurdenBand::Frayed;
    if (b >= 40.0f) return BurdenBand::Marked;
    if (b >= 20.0f) return BurdenBand::Unsettled;
    return BurdenBand::Clear;
}

const char* Burden_BandName(BurdenBand b) {
    switch (b) {
        case BurdenBand::Clear: return "Clear";
        case BurdenBand::Unsettled: return "Unsettled";
        case BurdenBand::Marked: return "Marked";
        case BurdenBand::Frayed: return "Frayed";
        case BurdenBand::Overburdened: return "Overburdened";
    }
    return "";
}

const char* Burden_BandSymptoms(BurdenBand b) {
    switch (b) {
        case BurdenBand::Clear: return "No symptoms.";
        case BurdenBand::Unsettled: return "A shimmer at the edges, whispers in quiet places; mana regenerates 10% slower.";
        case BurdenBand::Marked: return "A harmless wisp may follow you (a cast shoos it); mana costs 10% more; magic sight flickers now and then.";
        case BurdenBand::Frayed: return "Big spells (5+ parts) may fizzle, refunding half; strain drains 25% slower; animals shy away.";
        case BurdenBand::Overburdened: return "Once a day a weak shade comes for you (telegraphed, killable); backlash comes sooner.";
    }
    return "";
}

float Burden_RegenFactor(float b) { return Burden_Band(b) >= BurdenBand::Unsettled ? 0.9f : 1.0f; }
float Burden_CostFactor(float b) { return Burden_Band(b) >= BurdenBand::Marked ? 1.1f : 1.0f; }
float Burden_MisfireChance(float b, int parts) { return parts >= 5 && Burden_Band(b) >= BurdenBand::Frayed ? 0.15f : 0.0f; }
float Burden_BacklashAt(float b) { return Burden_Band(b) >= BurdenBand::Overburdened ? 85.0f : 100.0f; }

namespace { float CharmCut(const MagicPlayerState& s) { return s.wardingCharm ? 0.75f : 1.0f; } }

void Burden_Add(MagicPlayerState& s, float amount) {
    if (amount <= 0.0f) return;
    s.burden = std::min(100.0f, s.burden + amount * CharmCut(s));
}

bool Strain_Add(MagicPlayerState& s, float amount) {
    if (amount <= 0.0f) return false;
    // A quarter of the strain gained while over half strain becomes burden (only the part above half counts).
    const float before = s.strain, after = std::min(100.0f, s.strain + amount);
    const float overHalf = std::max(0.0f, after - std::max(before, 50.0f));
    s.strain = after;
    if (overHalf > 0.0f) Burden_Add(s, 0.25f * overHalf);
    if (s.strain >= Burden_BacklashAt(s.burden)) { // backlash: the pool empties, the strain resets to half, 3 burden
        s.pool = 0.0f;
        s.strain = 50.0f;
        Burden_Add(s, 3.0f);
        return true;
    }
    return false;
}

void Strain_Drain(MagicPlayerState& s, float dt, bool asleep) {
    float rate = 100.0f / 300.0f;                     // a full bar in 5 minutes
    if (asleep) rate *= 3.0f;                         // faster asleep in a bed
    if (Burden_Band(s.burden) >= BurdenBand::Frayed) rate *= 0.75f;
    s.strain = std::max(0.0f, s.strain - rate * std::max(0.0f, dt));
}

float Burden_BacklashHealth(float health, float maxHealth) {
    // Damage that never kills: it stops at a quarter of the player's health (the doc).
    const float floorHp = 0.25f * maxHealth;
    if (health <= floorHp) return health;
    return std::max(floorHp, health - 0.4f * maxHealth);
}

void Burden_DayPassed(MagicPlayerState& s, bool castToday) { if (!castToday) s.burden = std::max(0.0f, s.burden - 1.0f); }
void Burden_Cure(MagicPlayerState& s, float amount) { s.burden = std::max(0.0f, s.burden - std::max(0.0f, amount)); }
bool Burden_LatePagesShown(const MagicPlayerState& s) { return s.burden >= 20.0f; }

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1255): bands, strain into burden, backlash, drains, cures, the charm, the save lines ----
#include <gtest/gtest.h>
#include <sstream>

TEST(Burden, BandsAndTheirEffects) {
    EXPECT_EQ(Burden_Band(0), BurdenBand::Clear);
    EXPECT_EQ(Burden_Band(19.9f), BurdenBand::Clear);
    EXPECT_EQ(Burden_Band(20), BurdenBand::Unsettled);
    EXPECT_EQ(Burden_Band(45), BurdenBand::Marked);
    EXPECT_EQ(Burden_Band(60), BurdenBand::Frayed);
    EXPECT_EQ(Burden_Band(100), BurdenBand::Overburdened);
    EXPECT_FLOAT_EQ(Burden_RegenFactor(10), 1.0f);
    EXPECT_FLOAT_EQ(Burden_RegenFactor(25), 0.9f);
    EXPECT_FLOAT_EQ(Burden_CostFactor(39), 1.0f);
    EXPECT_FLOAT_EQ(Burden_CostFactor(40), 1.1f);
    EXPECT_EQ(Burden_MisfireChance(70, 4), 0.0f) << "only big spells";
    EXPECT_GT(Burden_MisfireChance(70, 5), 0.0f);
    EXPECT_EQ(Burden_MisfireChance(50, 7), 0.0f) << "only Frayed and up";
    EXPECT_FLOAT_EQ(Burden_BacklashAt(79), 100.0f);
    EXPECT_FLOAT_EQ(Burden_BacklashAt(85), 85.0f);
    for (int b = 0; b < 5; ++b) EXPECT_STRNE(Burden_BandSymptoms((BurdenBand)b), "") << "the HUD always says what a band does";
}

TEST(Burden, AQuarterOfTheStrainOverHalfStays) {
    MagicPlayerState s;
    EXPECT_FALSE(Strain_Add(s, 40.0f));
    EXPECT_FLOAT_EQ(s.burden, 0.0f) << "under half strain: no burden";
    EXPECT_FALSE(Strain_Add(s, 20.0f));                 // 40 -> 60: 10 of it over half
    EXPECT_FLOAT_EQ(s.burden, 2.5f);
    s.wardingCharm = true;
    EXPECT_FALSE(Strain_Add(s, 20.0f));                 // 60 -> 80: all over half, cut by the charm
    EXPECT_FLOAT_EQ(s.burden, 2.5f + 5.0f * 0.75f);
}

TEST(Burden, BacklashEmptiesThePoolAndNeverKills) {
    MagicPlayerState s;
    s.strain = 90.0f; s.pool = 70.0f;
    EXPECT_TRUE(Strain_Add(s, 15.0f));
    EXPECT_FLOAT_EQ(s.pool, 0.0f);
    EXPECT_FLOAT_EQ(s.strain, 50.0f);
    EXPECT_GT(s.burden, 3.0f) << "the backlash itself leaves burden";
    EXPECT_FLOAT_EQ(Burden_BacklashHealth(100.0f, 100.0f), 60.0f);
    EXPECT_FLOAT_EQ(Burden_BacklashHealth(30.0f, 100.0f), 25.0f) << "it stops at a quarter";
    EXPECT_FLOAT_EQ(Burden_BacklashHealth(10.0f, 100.0f), 10.0f) << "never lower than it was";
    MagicPlayerState heavy; heavy.burden = 90.0f; heavy.strain = 80.0f;
    EXPECT_TRUE(Strain_Add(heavy, 6.0f)) << "Overburdened: backlash at 85";
}

TEST(Burden, StrainDrainsBurdenLingersCuresWork) {
    MagicPlayerState s;
    s.strain = 100.0f;
    Strain_Drain(s, 150.0f, false);
    EXPECT_NEAR(s.strain, 50.0f, 0.01f) << "five minutes for a full bar";
    Strain_Drain(s, 50.0f, true);
    EXPECT_NEAR(s.strain, 0.0f, 0.01f) << "three times faster asleep";
    MagicPlayerState frayed; frayed.burden = 65.0f; frayed.strain = 100.0f;
    Strain_Drain(frayed, 150.0f, false);
    EXPECT_NEAR(frayed.strain, 62.5f, 0.01f) << "Frayed: 25% slower";
    s.burden = 30.0f;
    Burden_DayPassed(s, true);
    EXPECT_FLOAT_EQ(s.burden, 30.0f) << "a day of casting: none of it goes";
    Burden_DayPassed(s, false);
    EXPECT_FLOAT_EQ(s.burden, 29.0f);
    EXPECT_TRUE(Burden_LatePagesShown(s));
    Burden_Cure(s, 15.0f);
    EXPECT_FLOAT_EQ(s.burden, 14.0f);
    EXPECT_FALSE(Burden_LatePagesShown(s)) << "the late pages show at 20 and up";
    Burden_Cure(s, 50.0f);
    EXPECT_FLOAT_EQ(s.burden, 0.0f);
}

TEST(Burden, SavedWithThePlayer) {
    MagicPlayerState a; a.strain = 12.5f; a.burden = 33.0f; a.wardingCharm = true;
    MagicPlayerState b;
    std::istringstream in(MagicSight_SaveText(a));
    for (std::string line; std::getline(in, line);) MagicSight_LoadLine(b, line.c_str());
    EXPECT_NEAR(b.strain, 12.5f, 1e-3f);
    EXPECT_NEAR(b.burden, 33.0f, 1e-3f);
    EXPECT_TRUE(b.wardingCharm);
}

#include "research.h"
#include "tome.h"
TEST(Burden, TheLatePagesNeedBurdenTwenty) {
    ASSERT_TRUE(Research_Load());
    ASSERT_TRUE(Tome_Load());
    EXPECT_EQ(Tome_Validate(), "");
    const TomePage* p = Tome_Find("void_mastery");
    ASSERT_TRUE(p);
    TomeState team; team.pages = { "first_sight", "wider_sight" };
    MagicPlayerState me; me.aspects["void"] = 99; me.aspects["order"] = 99;
    ResearchState research;
    std::string why;
    me.burden = 10.0f;
    EXPECT_FALSE(Tome_CanUnlock(team, me, research, *p, &why)) << "hidden under burden 20";
    me.burden = 25.0f;
    EXPECT_TRUE(Tome_CanUnlock(team, me, research, *p, &why)) << why;
}
#endif
