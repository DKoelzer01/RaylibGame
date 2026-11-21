#include "potions.h"

#include "alchemy.h"
#include "burden.h"
#include "magicsight.h"
#include "suit.h"
#include "survival.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

bool Potion_Drink(const PotionDef& p, PotionState& ps, MagicPlayerState& magic, SurvivalState& survival, long long day, std::string* why) {
    auto refuse = [&](const char* w) { if (why) *why = w; return false; };
    switch (p.effect) {
    case PotionEffect::Heal:
        if (survival.health >= 99.9f) return refuse("You are not hurt");
        survival.health = std::min(100.0f, survival.health + p.amount);
        return true;
    case PotionEffect::CleanseBurden:
        if (p.onceADay && ps.cleansedDay == day) return refuse("The draught works only once a day");
        if (magic.burden <= 0.0f) return refuse("You carry no burden");
        Burden_Cure(magic, p.amount);
        ps.cleansedDay = day;
        return true;
    default:
        ps.left[(int)p.effect] = p.seconds; // drunk again: its time starts over
        return true;
    }
}

void Potion_Tick(PotionState& ps, float dt) {
    for (auto it = ps.left.begin(); it != ps.left.end();) {
        it->second -= std::max(0.0f, dt);
        if (it->second <= 0.0f) it = ps.left.erase(it); else ++it;
    }
}
float Potion_Left(const PotionState& ps, int effect) { auto it = ps.left.find(effect); return it == ps.left.end() ? 0.0f : it->second; }

SuitStats Potion_Suit(const PotionState& ps, const SuitStats& suit) {
    SuitStats s = suit;
    if (Potion_Left(ps, (int)PotionEffect::WaterBreathing) > 0.0f) s.underwater = 0.0f;
    if (Potion_Left(ps, (int)PotionEffect::ColdResist) > 0.0f) s.cold = 0.0f;
    return s;
}

float Potion_StrokeScale(const PotionState& ps) {
    if (Potion_Left(ps, (int)PotionEffect::MiningHaste) <= 0.0f) return 1.0f;
    for (const PotionDef& p : Alchemy_Potions()) if (p.effect == PotionEffect::MiningHaste && p.amount > 0.0f) return p.amount;
    return 1.5f;
}

std::string Potion_SaveText(const PotionState& ps) {
    std::string out;
    char buf[96];
    for (const auto& [e, t] : ps.left) if (t > 0.0f) { snprintf(buf, sizeof buf, "magic_potion %d %.2f\n", e, t); out += buf; }
    if (ps.cleansedDay >= 0) { snprintf(buf, sizeof buf, "magic_potion_cleansed %lld\n", ps.cleansedDay); out += buf; }
    return out;
}
bool Potion_LoadLine(PotionState& ps, const char* line) {
    int e = 0; float t = 0.0f; long long d = 0;
    if (sscanf(line, "magic_potion_cleansed %lld", &d) == 1) { ps.cleansedDay = d; return true; }
    if (sscanf(line, "magic_potion %d %f", &e, &t) == 2) { if (e >= 0 && e <= (int)PotionEffect::CleanseBurden && t > 0.0f) ps.left[e] = t; return true; }
    return false;
}

// ---- The game ----
#include "audio.h"
#include "gameui.h"
#include "inventory.h"
#include "magic.h"
#include "player.h"
#include "playersession.h"
#include "terraintool.h"
#include "world.h"
#include "worldtypes.h"

namespace { PotionState local; }

PotionState& Potion_Local() { return local; }
bool Potion_IsPotion(const std::string& item) { return Magic_Enabled() && Alchemy_PotionOfItem(item) != nullptr; }

long long Potion_Day() {
    const Planetoid* b = player.nearest;
    const float spin = b ? std::max({ fabsf(b->spin.x), fabsf(b->spin.y), fabsf(b->spin.z) }) : 0.0f;
    const double dayLength = spin > 1e-6f ? 360.0 / spin : 1200.0;
    return (long long)std::floor(WorldTime() / dayLength);
}

bool Potion_UseLocal(const std::string& item) {
    const PotionDef* p = Alchemy_PotionOfItem(item);
    if (!p || inventory.count(item) <= 0) return false;
    std::string why;
    if (!Potion_Drink(*p, local, LocalSession().magic, survival, Potion_Day(), &why)) { GameUI_Toast(why.c_str()); Audio_Play(Sfx::Error, 0.4f); return false; }
    inventory.remove(item, 1);
    inventory.add("phial", 1); // the empty phial comes back
    GameUI_Toast(TextFormat("%s: %s", p->name.c_str(), p->text.c_str()));
    Audio_Play(Sfx::Pickup, 0.5f, 0.7f);
    return true;
}

void Potion_Update(float dt) {
    if (!Magic_Enabled()) return;
    static bool hooked = false;
    if (!hooked) { hooked = true; gTerrainStrokeScale = [] { return Potion_StrokeScale(local); }; }
    Potion_Tick(local, dt);
}

bool Potion_NightVision() { return Potion_Left(local, (int)PotionEffect::NightVision) > 0.0f; }

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1254): what each potion does ----
#include <gtest/gtest.h>
#include "tome.h"

namespace {
const PotionDef& Shipped(const char* id) {
    static bool loaded = false;
    if (!loaded) { loaded = Tome_Load() && Alchemy_Load(); }
    const PotionDef* p = Alchemy_PotionOfItem(std::string("potion_") + id);
    if (!p) { static PotionDef none; ADD_FAILURE() << "no potion " << id; return none; }
    return *p;
}
}

TEST(Potions, TheCleansingDraughtWorksOnceADay) {
    const PotionDef& c = Shipped("cleansing");
    PotionState ps; MagicPlayerState m; SurvivalState s;
    m.burden = 40.0f;
    std::string why;
    EXPECT_TRUE(Potion_Drink(c, ps, m, s, 3, &why));
    EXPECT_FLOAT_EQ(m.burden, 25.0f) << "15 off";
    EXPECT_FALSE(Potion_Drink(c, ps, m, s, 3, &why)) << "the same day";
    EXPECT_FLOAT_EQ(m.burden, 25.0f);
    EXPECT_FALSE(why.empty());
    EXPECT_TRUE(Potion_Drink(c, ps, m, s, 4, &why)) << "the next day";
    EXPECT_FLOAT_EQ(m.burden, 10.0f);
    m.burden = 0.0f;
    EXPECT_FALSE(Potion_Drink(c, ps, m, s, 9, &why)) << "no burden: not wasted";
}

TEST(Potions, HealingAndTimedEffects) {
    PotionState ps; MagicPlayerState m; SurvivalState s;
    s.health = 100.0f;
    EXPECT_FALSE(Potion_Drink(Shipped("healing"), ps, m, s, 0)) << "not hurt: not wasted";
    s.health = 30.0f;
    EXPECT_TRUE(Potion_Drink(Shipped("healing"), ps, m, s, 0));
    EXPECT_FLOAT_EQ(s.health, 70.0f);
    s.health = 90.0f;
    EXPECT_TRUE(Potion_Drink(Shipped("healing"), ps, m, s, 0));
    EXPECT_FLOAT_EQ(s.health, 100.0f) << "never above full";

    SuitStats suit; suit.underwater = 2.0f; suit.cold = 1.0f;
    EXPECT_FLOAT_EQ(Potion_Suit(ps, suit).underwater, 2.0f);
    EXPECT_TRUE(Potion_Drink(Shipped("water_breathing"), ps, m, s, 0));
    EXPECT_TRUE(Potion_Drink(Shipped("cold_resist"), ps, m, s, 0));
    EXPECT_FLOAT_EQ(Potion_Suit(ps, suit).underwater, 0.0f) << "no oxygen lost under water";
    EXPECT_FLOAT_EQ(Potion_Suit(ps, suit).cold, 0.0f);
    EXPECT_FLOAT_EQ(Potion_StrokeScale(ps), 1.0f);
    EXPECT_TRUE(Potion_Drink(Shipped("mining_haste"), ps, m, s, 0));
    EXPECT_FLOAT_EQ(Potion_StrokeScale(ps), 1.5f);
    // Time runs out: water breathing lasts 120 s.
    Potion_Tick(ps, 119.0f);
    EXPECT_FLOAT_EQ(Potion_Suit(ps, suit).underwater, 0.0f);
    Potion_Tick(ps, 2.0f);
    EXPECT_FLOAT_EQ(Potion_Suit(ps, suit).underwater, 2.0f);
    EXPECT_FLOAT_EQ(Potion_StrokeScale(ps), 1.0f) << "mining haste (120 s) gone too";
    EXPECT_GT(Potion_Left(ps, (int)PotionEffect::ColdResist), 0.0f) << "cold resistance lasts 300 s";
}

TEST(Potions, StateSavesAndLoads) {
    PotionState ps; ps.left[(int)PotionEffect::NightVision] = 42.5f; ps.cleansedDay = 17;
    PotionState back;
    const std::string text = Potion_SaveText(ps);
    size_t a = 0;
    while (a < text.size()) { const size_t b = text.find('\n', a); EXPECT_TRUE(Potion_LoadLine(back, text.substr(a, b - a).c_str())); a = b + 1; }
    EXPECT_NEAR(Potion_Left(back, (int)PotionEffect::NightVision), 42.5f, 0.01f);
    EXPECT_EQ(back.cleansedDay, 17);
    EXPECT_FALSE(Potion_LoadLine(back, "magic_pool 1 2"));
}
#endif
