// Combat step 1 (#1215): see damage.h.
#include "damage.h"

#include <algorithm>
#include <cstdlib>

using namespace StatusTune;

bool Combat1_On() {
    static const bool on = !(getenv("RAYLIBGAME_COMBAT1") && getenv("RAYLIBGAME_COMBAT1")[0] == '0');
    return on;
}

const char* DamageType_Name(DamageType t) {
    switch (t) {
        case DamageType::Blunt: return "blunt";
        case DamageType::Slash: return "slash";
        case DamageType::Pierce: return "pierce";
        case DamageType::Heat: return "heat";
        case DamageType::Cold: return "cold";
        case DamageType::Energy: return "energy";
        default: return "?";
    }
}

const char* StatusKind_Name(StatusKind s) {
    switch (s) {
        case StatusKind::Burn: return "burn";
        case StatusKind::Chill: return "chill";
        case StatusKind::Bleed: return "bleed";
        case StatusKind::Poison: return "poison";
        case StatusKind::Shock: return "shock";
        default: return "none";
    }
}

float Damage_Apply(float amount, DamageType type, const Resistances& res, bool ignoresArmour) {
    if (amount <= 0.0f) return 0.0f;
    const int ti = std::clamp((int)type, 0, (int)DamageType::Count - 1);
    const float resisted = amount * (1.0f - std::clamp(res.r[ti], 0.0f, 1.0f));
    const float afterArmour = ignoresArmour ? resisted : resisted - std::max(0.0f, res.armour);
    return std::max(amount * kMinDamageShare, afterArmour);
}

StatusKind Status_FromHit(DamageType type, StatusKind explicitStatus, bool heavy, bool fromBehind, bool bite) {
    if (explicitStatus != StatusKind::None) return explicitStatus;
    switch (type) {
        case DamageType::Heat: return StatusKind::Burn;
        case DamageType::Cold: return StatusKind::Chill;
        case DamageType::Energy: return StatusKind::Shock;
        case DamageType::Slash: if (heavy || fromBehind) return StatusKind::Bleed; break;
        default: break;
    }
    return bite ? StatusKind::Poison : StatusKind::None;
}

void Status_Apply(StatusState& s, StatusKind k) {
    switch (k) {
        case StatusKind::Burn: s.burn = kBurnTime; break;
        case StatusKind::Chill:
            s.chill = kChillTime;
            if (++s.chillStacks >= kChillFreezeStacks) { s.frozen = kFreezeTime; s.chillStacks = 0; } // three stacks: frozen
            break;
        case StatusKind::Bleed: s.bleed = kBleedTime; s.bleedStacks = std::min(kBleedMaxStacks, s.bleedStacks + 1); break;
        case StatusKind::Poison: s.poison = kPoisonTime; break;
        case StatusKind::Shock: s.interrupted = true; break;
        default: break;
    }
}

float Status_Step(StatusState& s, float dt, bool inWater, const Resistances& res) {
    if (dt <= 0.0f) return 0.0f;
    float dmg = 0.0f;
    if (s.burn > 0.0f) {
        if (inWater) s.burn = 0.0f;                                              // water puts it out
        else { dmg += Damage_Apply(kBurnDps, DamageType::Heat, res) * dt; s.burn = std::max(0.0f, s.burn - dt); }
    }
    if (s.bleed > 0.0f) {
        dmg += Damage_Apply(kBleedDps * (float)s.bleedStacks, DamageType::Slash, res) * dt;
        s.bleed = std::max(0.0f, s.bleed - dt);
        if (s.bleed <= 0.0f) s.bleedStacks = 0;
    }
    if (s.poison > 0.0f) { dmg += kPoisonDps * dt; s.poison = std::max(0.0f, s.poison - dt); } // ignores armour (and resistance: no type)
    if (s.chill > 0.0f) { s.chill = std::max(0.0f, s.chill - dt); if (s.chill <= 0.0f) s.chillStacks = 0; }
    if (s.frozen > 0.0f) s.frozen = std::max(0.0f, s.frozen - dt);
    s.interrupted = false;
    return dmg;
}

float Status_SpeedFactor(const StatusState& s) {
    if (s.frozen > 0.0f) return 0.0f;
    return s.chill > 0.0f ? 1.0f - kChillSlow : 1.0f;
}

void Status_CurePoison(StatusState& s) { s.poison = 0.0f; }

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>
// #1215: the damage formula and each status effect, with the decided numbers.
TEST(Damage, FormulaResistanceArmourAndFloor) {
    Resistances res;
    EXPECT_FLOAT_EQ(Damage_Apply(20.0f, DamageType::Blunt, res), 20.0f);            // nothing in the way: today's damage
    res.r[(int)DamageType::Slash] = 0.5f; res.armour = 3.0f;
    EXPECT_FLOAT_EQ(Damage_Apply(20.0f, DamageType::Slash, res), 7.0f);             // 20 x 0.5 - 3
    EXPECT_FLOAT_EQ(Damage_Apply(20.0f, DamageType::Blunt, res), 17.0f);            // other types: armour only
    EXPECT_FLOAT_EQ(Damage_Apply(20.0f, DamageType::Slash, res, true), 10.0f);      // ignoring armour
    res.armour = 100.0f;
    EXPECT_FLOAT_EQ(Damage_Apply(20.0f, DamageType::Blunt, res), 2.0f);             // never below 10%
    res.r[(int)DamageType::Heat] = 1.0f;
    EXPECT_FLOAT_EQ(Damage_Apply(20.0f, DamageType::Heat, res), 2.0f);              // immune still takes the floor
}

TEST(Damage, StatusFromHit) {
    EXPECT_EQ(Status_FromHit(DamageType::Heat, StatusKind::None, false, false, false), StatusKind::Burn);
    EXPECT_EQ(Status_FromHit(DamageType::Cold, StatusKind::None, false, false, false), StatusKind::Chill);
    EXPECT_EQ(Status_FromHit(DamageType::Energy, StatusKind::None, false, false, false), StatusKind::Shock);
    EXPECT_EQ(Status_FromHit(DamageType::Slash, StatusKind::None, false, false, false), StatusKind::None);
    EXPECT_EQ(Status_FromHit(DamageType::Slash, StatusKind::None, true, false, false), StatusKind::Bleed);   // a heavy
    EXPECT_EQ(Status_FromHit(DamageType::Slash, StatusKind::None, false, true, false), StatusKind::Bleed);   // from behind
    EXPECT_EQ(Status_FromHit(DamageType::Pierce, StatusKind::None, false, false, true), StatusKind::Poison); // a bite
    EXPECT_EQ(Status_FromHit(DamageType::Blunt, StatusKind::Burn, false, false, false), StatusKind::Burn);  // explicit wins
}

static float RunFor(StatusState& s, float seconds, bool inWater, const Resistances& res = {}) {
    float total = 0.0f;
    for (float t = 0.0f; t < seconds - 1e-4f; t += 1.0f / 60.0f) total += Status_Step(s, 1.0f / 60.0f, inWater, res);
    return total;
}

TEST(Damage, BurnTicksAndWaterPutsItOut) {
    StatusState s; Status_Apply(s, StatusKind::Burn);
    EXPECT_NEAR(RunFor(s, 5.0f, false), StatusTune::kBurnDps * StatusTune::kBurnTime, 0.1f); // 3/s for 4 s = 12
    EXPECT_FALSE(s.any());
    StatusState w; Status_Apply(w, StatusKind::Burn);
    EXPECT_NEAR(RunFor(w, 1.0f, true), 0.0f, 1e-4f);                                       // in water: out at once
    EXPECT_FLOAT_EQ(w.burn, 0.0f);
}

TEST(Damage, ChillSlowsAndThreeStacksFreeze) {
    StatusState s; Status_Apply(s, StatusKind::Chill);
    EXPECT_FLOAT_EQ(Status_SpeedFactor(s), 0.7f);
    Status_Apply(s, StatusKind::Chill); Status_Apply(s, StatusKind::Chill);
    EXPECT_FLOAT_EQ(Status_SpeedFactor(s), 0.0f);                                           // frozen
    RunFor(s, 1.1f, false);
    EXPECT_FLOAT_EQ(Status_SpeedFactor(s), 0.7f);                                           // thawed, still chilled
    RunFor(s, 3.0f, false);
    EXPECT_FLOAT_EQ(Status_SpeedFactor(s), 1.0f);
}

TEST(Damage, BleedStacksToThree) {
    StatusState s;
    for (int i = 0; i < 5; ++i) Status_Apply(s, StatusKind::Bleed);
    EXPECT_EQ(s.bleedStacks, 3);
    EXPECT_NEAR(RunFor(s, 7.0f, false), StatusTune::kBleedDps * 3 * StatusTune::kBleedTime, 0.2f); // 2/s x 3 for 6 s = 36
    EXPECT_EQ(s.bleedStacks, 0);
}

TEST(Damage, PoisonIgnoresArmourAndHerbCures) {
    Resistances heavy; heavy.armour = 50.0f;
    StatusState s; Status_Apply(s, StatusKind::Poison);
    EXPECT_NEAR(RunFor(s, 4.0f, false, heavy), StatusTune::kPoisonDps * 4.0f, 0.1f);        // 1.5/s through any armour
    Status_CurePoison(s);
    EXPECT_NEAR(RunFor(s, 4.0f, false, heavy), 0.0f, 1e-4f);
}

TEST(Damage, ShockInterruptsForOneStep) {
    StatusState s; Status_Apply(s, StatusKind::Shock);
    EXPECT_TRUE(s.interrupted);
    Status_Step(s, 1.0f / 60.0f, false, {});
    EXPECT_FALSE(s.interrupted);
}
#endif
