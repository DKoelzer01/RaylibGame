// Combat step 2 (#1216): see melee.h.
#include "melee.h"
#include "composed.h"
#include "tools.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>

using namespace MeleeTune;

bool Combat2_On() {
    static const bool on = !(getenv("RAYLIBGAME_COMBAT2") && getenv("RAYLIBGAME_COMBAT2")[0] == '0');
    return on;
}

// The combat doc's table (2.2, light timings in ms: wind-up / hit window / recovery). Poise: low 5, medium 10, medium-high
// 14, high 20, very high 30 (guess; #1217 tunes them). Arcs are a guess: generous, the doc's "first-person aim".
static const WeaponSpec kWeapons[(int)WeaponClass::Count] = {
    { "tools",      0.25f, 0.10f, 0.30f, 1.6f, DamageType::Blunt,  6.0f,  5.0f,  70.0f, false },
    { "knife",      0.18f, 0.09f, 0.22f, 1.5f, DamageType::Slash,  6.0f,  5.0f,  70.0f, false },
    { "sword",      0.30f, 0.12f, 0.35f, 2.0f, DamageType::Slash, 10.0f, 10.0f, 100.0f, false },
    { "axe",        0.38f, 0.12f, 0.45f, 2.0f, DamageType::Slash, 12.0f, 14.0f,  90.0f, false },
    { "mace",       0.45f, 0.13f, 0.55f, 1.9f, DamageType::Blunt, 14.0f, 20.0f,  80.0f, false },
    { "spear",      0.35f, 0.11f, 0.40f, 2.8f, DamageType::Pierce, 10.0f, 10.0f, 20.0f, false },
    { "two-handed", 0.60f, 0.15f, 0.70f, 2.4f, DamageType::Slash, 20.0f, 30.0f, 120.0f, true },
};

const WeaponSpec& Weapon_Spec(WeaponClass c) { return kWeapons[std::clamp((int)c, 0, (int)WeaponClass::Count - 1)]; }

WeaponClass Weapon_ClassOf(const std::string& id) {
    auto has = [&](const char* w) { return id.find(w) != std::string::npos; };
    if (has("greatsword") || has("maul")) return WeaponClass::TwoHanded;
    if (has("knife") || has("dagger")) return WeaponClass::Knife;
    if (has("sword")) return WeaponClass::Sword;
    if (has("axe") && !has("pickaxe")) return WeaponClass::Axe;
    if (has("mace") || has("hammer") || has("club")) return WeaponClass::Mace;
    if (has("spear")) return WeaponClass::Spear;
    return WeaponClass::Tools;
}

float Weapon_TimeScale(const std::string& id, const std::string& data) {
    if (!Composed_Is(id)) return 1.0f;
    const ComposedStats s = Composed_Stats(id, data), ref = Composed_Stats(id, Composed_MakeData("iron"));
    if (!s.valid || !ref.valid || ref.weight <= 0.0f) return 1.0f;
    return std::clamp(powf(s.weight / ref.weight, 0.35f), 0.85f, 1.2f); // GUESS: the exponent and the clamp (#1219)
}

WeaponSpec Weapon_SpecOf(const std::string& id, const std::string& data) {
    WeaponSpec w = Weapon_Spec(Weapon_ClassOf(id));
    const float k = Weapon_TimeScale(id, data);
    w.windUp *= k; w.hit *= k; w.recovery *= k;
    return w;
}

float Weapon_BaseDamage(const std::string& id, const std::string& data) {
    if (Composed_Is(id)) { // (#1219) a composed copy: from its head (3 + 3 x hardness) and quality
        const ComposedStats s = Composed_Stats(id, data);
        if (s.valid) return s.damage;
    }
    const ToolSpec t = ToolOf(id, data);
    return t.type == ToolType::Hands ? 2.0f : 3.0f + 3.0f * t.tier;
}

DamageType Weapon_TypeOf(const std::string& id) {
    if (id.find("pickaxe") != std::string::npos) return DamageType::Pierce;
    return Weapon_Spec(Weapon_ClassOf(id)).type;
}

float Melee_WindUpTime(const MeleeState& s, const WeaponSpec& w) { return s.heavy ? w.windUp * kHeavyWindUp : w.windUp; }
bool Melee_HyperArmour(const MeleeState& s, const WeaponSpec& w) { return s.phase == MeleePhase::WindUp && (s.heavy || w.twoHanded); }
float Melee_Damage(float base, const MeleeState& s) { return s.heavy ? base * kHeavyDamage : base; }
float Melee_Poise(const MeleeState& s, const WeaponSpec& w) { return s.heavy ? w.poise * kHeavyPoise : w.poise; }
bool Melee_Struck(const MeleeState& s, uint64_t id) { return std::find(s.struck.begin(), s.struck.end(), id) != s.struck.end(); }

float Melee_SweepAngle(const MeleeState& s, const WeaponSpec& w, float f) {
    const float half = 0.5f * w.arcDeg * DEG2RAD;
    const float a = half - 2.0f * half * std::clamp(f, 0.0f, 1.0f); // right to left
    return s.chain == 1 ? -a : a;                                     // the chained light comes back
}

MeleeStep Melee_Step(MeleeState& s, const WeaponSpec& w, const MeleeInput& in, float dt) {
    MeleeStep o;
    auto start = [&](int chain) {
        if (in.stamina - o.cost < w.stamina) { o.tired = true; s.phase = MeleePhase::Idle; s.chain = 0; s.queued = false; return false; }
        s.phase = MeleePhase::WindUp; s.t = 0.0f; s.heavy = false; s.chain = chain; s.queued = false; s.sweep = 0.0f;
        s.struck.clear();
        o.started = true; o.cost += w.stamina; ++s.swings;
        return true;
    };
    auto cancel = [&]() { s.phase = MeleePhase::Idle; s.t = 0.0f; s.chain = 0; s.queued = false; s.heavy = false; o.cancelled = true; };

    if (s.phase == MeleePhase::Idle) {
        if (!in.pressed || !start(0)) return o;
    } else {
        // Cancels and interrupts, judged as the step begins.
        if (s.phase == MeleePhase::WindUp) {
            if ((in.block || in.dodge) && s.t < 0.5f * Melee_WindUpTime(s, w)) { cancel(); return o; }
            if (in.interrupted && !Melee_HyperArmour(s, w)) { cancel(); return o; }
        }
        if (s.phase == MeleePhase::Recovery && in.dodge && s.t >= kDodgeCancelFrom * w.recovery) { cancel(); o.finished = true; return o; }
        // A light pressed during the first light's hit window or recovery chains (one chain only; heavies don't chain).
        if (in.pressed && (s.phase == MeleePhase::Hit || s.phase == MeleePhase::Recovery) && !s.heavy && s.chain == 0) s.queued = true;
    }

    float left = dt;
    for (int guard = 0; left > 1e-7f && guard < 8; ++guard) {
        if (s.phase == MeleePhase::WindUp) {
            const float end = Melee_WindUpTime(s, w);
            const float d = std::min(left, end - s.t);
            s.t += d; left -= d;
            if (s.t < end - 1e-6f) break;
            // The light wind-up ends: still held, it grows into a heavy (when the stamina for it is there).
            const float extra = w.stamina * (kHeavyStamina - 1.0f);
            if (!s.heavy && in.held && in.stamina - o.cost >= extra) {
                s.heavy = true; o.becameHeavy = true; o.cost += extra; ++s.heavies;
                continue;
            }
            s.phase = MeleePhase::Hit; s.t = 0.0f; s.sweep = 0.0f;
        } else if (s.phase == MeleePhase::Hit) {
            const float from = s.t / w.hit;
            const float d = std::min(left, w.hit - s.t);
            s.t += d; left -= d;
            const float to = std::min(1.0f, s.t / w.hit);
            if (!o.sweeping) o.sweepFrom = from;
            o.sweeping = true; o.sweepTo = to; s.sweep = to;
            if (s.t < w.hit - 1e-6f) break;
            s.phase = MeleePhase::Recovery; s.t = 0.0f;
        } else if (s.phase == MeleePhase::Recovery) {
            const float end = s.queued ? kChainAt * w.recovery : w.recovery;
            const float d = std::min(left, std::max(0.0f, end - s.t));
            s.t += d; left -= d;
            if (s.t < end - 1e-6f) break;
            if (s.queued) { if (!start(1)) { o.finished = true; break; } continue; }
            s.phase = MeleePhase::Idle; s.t = 0.0f; s.chain = 0; o.finished = true;
            break;
        } else break;
    }
    return o;
}

// Closest points of two segments (Ericson, Real-Time Collision Detection 5.1.9).
float Melee_SegmentDistance(Vector3 p0, Vector3 p1, Vector3 q0, Vector3 q1, Vector3* onP) {
    const Vector3 d1 = Vector3Subtract(p1, p0), d2 = Vector3Subtract(q1, q0), r = Vector3Subtract(p0, q0);
    const float a = Vector3DotProduct(d1, d1), e = Vector3DotProduct(d2, d2), f = Vector3DotProduct(d2, r);
    float s = 0.0f, t = 0.0f;
    if (a <= 1e-9f && e <= 1e-9f) { s = t = 0.0f; }
    else if (a <= 1e-9f) { s = 0.0f; t = std::clamp(f / e, 0.0f, 1.0f); }
    else {
        const float c = Vector3DotProduct(d1, r);
        if (e <= 1e-9f) { t = 0.0f; s = std::clamp(-c / a, 0.0f, 1.0f); }
        else {
            const float b = Vector3DotProduct(d1, d2), den = a * e - b * b;
            s = den > 1e-9f ? std::clamp((b * f - c * e) / den, 0.0f, 1.0f) : 0.0f;
            t = (b * s + f) / e;
            if (t < 0.0f) { t = 0.0f; s = std::clamp(-c / a, 0.0f, 1.0f); }
            else if (t > 1.0f) { t = 1.0f; s = std::clamp((b - c) / a, 0.0f, 1.0f); }
        }
    }
    const Vector3 cp = Vector3Add(p0, Vector3Scale(d1, s)), cq = Vector3Add(q0, Vector3Scale(d2, t));
    if (onP) *onP = cp;
    return Vector3Distance(cp, cq);
}

bool Melee_SweepCapsule(Vector3 eye, Vector3 forward, Vector3 up, float reach, float a0, float a1,
                        Vector3 q0, Vector3 q1, float r, Vector3* contact, float* dist) {
    Vector3 right = Vector3CrossProduct(forward, up);
    if (Vector3Length(right) < 1e-5f) return false;
    right = Vector3Normalize(right);
    const int n = std::max(2, (int)ceilf(fabsf(a1 - a0) / (kSweepStepDeg * DEG2RAD)) + 1);
    bool hit = false;
    float best = 1e9f;
    for (int i = 0; i < n; ++i) {
        const float a = a0 + (a1 - a0) * (float)i / (float)(n - 1);
        const Vector3 dir = Vector3Normalize(Vector3Add(Vector3Scale(forward, cosf(a)), Vector3Scale(right, sinf(a))));
        const Vector3 p0 = Vector3Add(eye, Vector3Scale(dir, kBladeStart)), p1 = Vector3Add(eye, Vector3Scale(dir, reach));
        Vector3 on;
        if (Melee_SegmentDistance(p0, p1, q0, q1, &on) > r) continue;
        const float d = Vector3Distance(eye, on);
        if (d < best) { best = d; hit = true; if (contact) *contact = on; }
    }
    if (hit && dist) *dist = best;
    return hit;
}

Vector3 Melee_SoftTurn(Vector3 forward, Vector3 up, Vector3 toTarget, float dt) {
    if (Vector3Length(toTarget) < 1e-4f || dt <= 0.0f) return forward;
    const float k = 1.0f - expf(-dt / kLockTurnTime);
    const Vector3 t = Vector3Normalize(toTarget);
    const Vector3 f = Vector3Normalize(Vector3Lerp(forward, t, k));
    return fabsf(Vector3DotProduct(f, up)) < 0.985f ? f : forward; // the look's pitch clamp
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>
namespace {
constexpr float kDt = 1.0f / 60.0f;
// Steps until `pred` holds (at most `max`); the seconds it took.
template <class P> float RunUntil(MeleeState& s, const WeaponSpec& w, MeleeInput in, P pred, int max = 600) {
    for (int i = 1; i <= max; ++i) { Melee_Step(s, w, in, kDt); in.pressed = false; if (pred(s)) return i * kDt; }
    return -1.0f;
}
}

TEST(Melee, LightTimingsFollowTheTable) {
    const WeaponSpec& w = Weapon_Spec(WeaponClass::Sword);
    MeleeState s; MeleeInput in; in.pressed = true; in.held = true;
    const MeleeStep first = Melee_Step(s, w, in, kDt);
    EXPECT_TRUE(first.started); EXPECT_FLOAT_EQ(first.cost, w.stamina);
    in.pressed = false; in.held = false; // a tap
    const float toHit = kDt + RunUntil(s, w, in, [](const MeleeState& m) { return m.phase == MeleePhase::Hit; });
    EXPECT_NEAR(toHit, w.windUp, kDt + 1e-4f);
    EXPECT_FALSE(s.heavy);
    const float rest = RunUntil(s, w, in, [](const MeleeState& m) { return m.phase == MeleePhase::Idle; });
    EXPECT_NEAR(toHit + rest, w.windUp + w.hit + w.recovery, kDt + 1e-4f);
    EXPECT_EQ(s.swings, 1);
}

TEST(Melee, HoldingThroughTheWindUpMakesAHeavy) {
    const WeaponSpec& w = Weapon_Spec(WeaponClass::Tools);
    MeleeState s; MeleeInput in; in.pressed = true; in.held = true;
    float cost = 0.0f, t = 0.0f; bool heavy = false;
    for (int i = 0; i < 120 && s.phase != MeleePhase::Hit; ++i) {
        const MeleeStep o = Melee_Step(s, w, in, kDt); in.pressed = false; cost += o.cost; t += kDt; heavy |= o.becameHeavy;
    }
    EXPECT_TRUE(heavy); EXPECT_TRUE(s.heavy);
    EXPECT_NEAR(t, w.windUp * MeleeTune::kHeavyWindUp, kDt + 1e-4f);
    EXPECT_NEAR(cost, w.stamina * MeleeTune::kHeavyStamina, 1e-4f);
    EXPECT_FLOAT_EQ(Melee_Damage(10.0f, s), 18.0f);
    EXPECT_FLOAT_EQ(Melee_Poise(s, w), w.poise * 2.0f);
    // Not enough stamina for the heavy part: it stays a light.
    MeleeState l; MeleeInput lo; lo.pressed = true; lo.held = true; lo.stamina = w.stamina + 1.0f;
    Melee_Step(l, w, lo, kDt); lo.pressed = false; lo.stamina = 1.0f; // the swing paid: 1 left
    RunUntil(l, w, lo, [](const MeleeState& m) { return m.phase == MeleePhase::Hit; });
    EXPECT_FALSE(l.heavy);
}

TEST(Melee, TwoLightsChainButNotThree) {
    const WeaponSpec& w = Weapon_Spec(WeaponClass::Knife);
    MeleeState s; MeleeInput in; in.pressed = true;
    Melee_Step(s, w, in, kDt);
    in.pressed = false;
    RunUntil(s, w, in, [](const MeleeState& m) { return m.phase == MeleePhase::Hit; });
    in.pressed = true; Melee_Step(s, w, in, kDt); in.pressed = false; // pressed in the hit window: queued
    EXPECT_TRUE(s.queued);
    RunUntil(s, w, in, [](const MeleeState& m) { return m.phase == MeleePhase::Recovery; });
    const float toSecond = RunUntil(s, w, in, [](const MeleeState& m) { return m.phase == MeleePhase::WindUp; });
    EXPECT_NEAR(toSecond, MeleeTune::kChainAt * w.recovery, kDt + 1e-4f);
    EXPECT_EQ(s.chain, 1);
    EXPECT_FLOAT_EQ(Melee_SweepAngle(s, w, 0.0f), -Melee_SweepAngle(MeleeState{}, w, 0.0f)); // swings back
    RunUntil(s, w, in, [](const MeleeState& m) { return m.phase == MeleePhase::Hit; });
    in.pressed = true; Melee_Step(s, w, in, kDt); in.pressed = false; // a third press: no chain
    EXPECT_FALSE(s.queued);
    RunUntil(s, w, in, [](const MeleeState& m) { return m.phase == MeleePhase::Idle; });
    EXPECT_EQ(s.swings, 2);
}

TEST(Melee, CancelsAndHyperArmour) {
    const WeaponSpec& w = Weapon_Spec(WeaponClass::Axe);
    MeleeInput press; press.pressed = true;
    { MeleeState s; Melee_Step(s, w, press, kDt); MeleeInput b; b.block = true; EXPECT_TRUE(Melee_Step(s, w, b, kDt).cancelled); }
    { MeleeState s; Melee_Step(s, w, press, kDt); MeleeInput n;
      RunUntil(s, w, n, [&](const MeleeState& m) { return m.t > 0.6f * w.windUp; });
      MeleeInput b; b.block = true; EXPECT_FALSE(Melee_Step(s, w, b, kDt).cancelled); }           // too late: committed
    { MeleeState s; Melee_Step(s, w, press, kDt); MeleeInput sh; sh.interrupted = true;
      EXPECT_TRUE(Melee_Step(s, w, sh, kDt).cancelled); }                                           // a shock stops a light
    { MeleeState s; MeleeInput h = press; h.held = true; Melee_Step(s, w, h, kDt); h.pressed = false;
      RunUntil(s, w, h, [](const MeleeState& m) { return m.heavy; });
      MeleeInput sh; sh.held = true; sh.interrupted = true;
      EXPECT_TRUE(Melee_HyperArmour(s, w)); EXPECT_FALSE(Melee_Step(s, w, sh, kDt).cancelled); }    // not a heavy
    { const WeaponSpec& two = Weapon_Spec(WeaponClass::TwoHanded); MeleeState s; Melee_Step(s, two, press, kDt);
      MeleeInput sh; sh.interrupted = true; EXPECT_FALSE(Melee_Step(s, two, sh, kDt).cancelled); }   // nor a two-hander
    // A dodge cancels only the last third of the recovery.
    MeleeState s; Melee_Step(s, w, press, kDt); MeleeInput n;
    RunUntil(s, w, n, [](const MeleeState& m) { return m.phase == MeleePhase::Recovery; });
    MeleeInput d; d.dodge = true;
    EXPECT_FALSE(Melee_Step(s, w, d, kDt).cancelled);
    RunUntil(s, w, n, [&](const MeleeState& m) { return m.t >= 0.7f * w.recovery; });
    EXPECT_TRUE(Melee_Step(s, w, d, kDt).cancelled);
    EXPECT_EQ(s.phase, MeleePhase::Idle);
}

TEST(Melee, TooTiredToSwing) {
    const WeaponSpec& w = Weapon_Spec(WeaponClass::Mace);
    MeleeState s; MeleeInput in; in.pressed = true; in.stamina = w.stamina - 1.0f;
    const MeleeStep o = Melee_Step(s, w, in, kDt);
    EXPECT_TRUE(o.tired); EXPECT_FALSE(o.started); EXPECT_EQ(s.phase, MeleePhase::Idle);
}

TEST(Melee, SweepFindsCapsulesInTheArcAndReach) {
    const Vector3 eye{ 0, 1.6f, 0 }, fwd{ 0, 0, 1 }, up{ 0, 1, 0 };
    const float half = 35.0f * DEG2RAD;
    auto capsuleAt = [&](float deg, float d, Vector3& q0, Vector3& q1) { // an upright capsule d m away at deg (+ right)
        const Vector3 right = Vector3Normalize(Vector3CrossProduct(fwd, up));
        const Vector3 c = Vector3Add(eye, Vector3Add(Vector3Scale(fwd, d * cosf(deg * DEG2RAD)), Vector3Scale(right, d * sinf(deg * DEG2RAD))));
        q0 = { c.x, 0.4f, c.z }; q1 = { c.x, 1.6f, c.z };
    };
    Vector3 q0, q1; float dist = 0.0f;
    capsuleAt(30.0f, 1.2f, q0, q1);
    EXPECT_TRUE(Melee_SweepCapsule(eye, fwd, up, 1.6f, half, -half, q0, q1, 0.3f, nullptr, &dist));
    EXPECT_NEAR(dist, 1.2f, 0.1f); // the blade crosses its axis
    EXPECT_FALSE(Melee_SweepCapsule(eye, fwd, up, 1.6f, 0.0f, -half, q0, q1, 0.05f));             // the left half misses it
    capsuleAt(-60.0f, 1.2f, q0, q1);
    EXPECT_FALSE(Melee_SweepCapsule(eye, fwd, up, 1.6f, half, -half, q0, q1, 0.3f));             // outside the arc
    capsuleAt(0.0f, 2.5f, q0, q1);
    EXPECT_FALSE(Melee_SweepCapsule(eye, fwd, up, 1.6f, half, -half, q0, q1, 0.3f));             // out of reach
    EXPECT_TRUE(Melee_SweepCapsule(eye, fwd, up, 2.8f, half, -half, q0, q1, 0.3f));              // a spear's reach
}

TEST(Melee, WeaponClasses) {
    EXPECT_EQ(Weapon_ClassOf(""), WeaponClass::Tools);
    EXPECT_EQ(Weapon_ClassOf("iron_pickaxe"), WeaponClass::Tools);
    EXPECT_EQ(Weapon_TypeOf("iron_pickaxe"), DamageType::Pierce);
    EXPECT_EQ(Weapon_ClassOf("stone_axe"), WeaponClass::Axe);
    EXPECT_EQ(Weapon_ClassOf("steel_greatsword"), WeaponClass::TwoHanded);
    EXPECT_EQ(Weapon_ClassOf("bronze_sword"), WeaponClass::Sword);
    EXPECT_EQ(Weapon_ClassOf("flint_knife"), WeaponClass::Knife);
    EXPECT_EQ(Weapon_ClassOf("war_hammer"), WeaponClass::Mace);
    EXPECT_EQ(Weapon_ClassOf("wooden_spear"), WeaponClass::Spear);
    // (#1219) the composed classes
    EXPECT_EQ(Weapon_ClassOf("weapon_knife"), WeaponClass::Knife);
    EXPECT_EQ(Weapon_ClassOf("weapon_sword"), WeaponClass::Sword);
    EXPECT_EQ(Weapon_ClassOf("weapon_greatsword"), WeaponClass::TwoHanded);
    EXPECT_EQ(Weapon_ClassOf("weapon_mace"), WeaponClass::Mace);
    EXPECT_EQ(Weapon_ClassOf("weapon_spear"), WeaponClass::Spear);
    EXPECT_EQ(Weapon_ClassOf("tool_axe"), WeaponClass::Axe);
    EXPECT_EQ(Weapon_ClassOf("tool_pickaxe"), WeaponClass::Tools);
}

TEST(Melee, AHeavierHeadSwingsSlower) { // (#1219)
    ASSERT_TRUE(Composed_Load());
    EXPECT_FLOAT_EQ(Weapon_TimeScale("weapon_sword", Composed_MakeData("iron")), 1.0f) << "iron is the reference";
    EXPECT_FLOAT_EQ(Weapon_TimeScale("stone_pickaxe", ""), 1.0f) << "not composed";
    const float bronze = Weapon_TimeScale("weapon_sword", Composed_MakeData("bronze")), titanium = Weapon_TimeScale("weapon_sword", Composed_MakeData("titanium"));
    const float stone = Weapon_TimeScale("weapon_sword", Composed_MakeData("stone"));
    EXPECT_GT(bronze, 1.0f) << "bronze is denser than iron";
    EXPECT_LT(titanium, 1.0f) << "titanium is light";
    EXPECT_GE(stone, 0.85f) << "clamped";
    EXPECT_LE(bronze, 1.2f);
    const WeaponSpec base = Weapon_Spec(WeaponClass::Sword), b = Weapon_SpecOf("weapon_sword", Composed_MakeData("bronze"));
    EXPECT_FLOAT_EQ(b.windUp, base.windUp * bronze);
    EXPECT_FLOAT_EQ(b.recovery, base.recovery * bronze);
    EXPECT_FLOAT_EQ(b.reach, base.reach) << "reach stays the class's";
    EXPECT_FLOAT_EQ(Weapon_SpecOf("", "").windUp, Weapon_Spec(WeaponClass::Tools).windUp) << "hands as before";
}

TEST(Melee, ComposedDamageComesFromTheHead) { // (#1219)
    ASSERT_TRUE(Composed_Load());
    EXPECT_FLOAT_EQ(Weapon_BaseDamage("", ""), 2.0f) << "hands";
    EXPECT_FLOAT_EQ(Weapon_BaseDamage("weapon_sword", Composed_MakeData("bronze")), 9.0f);  // 3 + 3 x 2.0
    EXPECT_FLOAT_EQ(Weapon_BaseDamage("weapon_sword", Composed_MakeData("steel")), 12.0f);  // 3 + 3 x 3.0
    EXPECT_FLOAT_EQ(Weapon_BaseDamage("weapon_knife", Composed_MakeData("flint")), 6.6f);   // 3 + 3 x 1.2
    EXPECT_FLOAT_EQ(Weapon_BaseDamage("tool_pickaxe", Composed_MakeData("stone")), 6.0f) << "the same as an old stone pickaxe (3 + 3 x tier 1)";
    EXPECT_NEAR(Weapon_BaseDamage("weapon_sword", Composed_MakeData("steel", {}, {}, 110)), 13.2f, 1e-4f) << "quality";
}

TEST(Melee, SoftTurnClosesOnTheTarget) {
    const Vector3 up{ 0, 1, 0 }, to{ 1, 0, 0 };
    Vector3 f{ 0, 0, 1 };
    for (int i = 0; i < 60; ++i) f = Melee_SoftTurn(f, up, to, kDt); // a second
    EXPECT_GT(Vector3DotProduct(f, to), 0.99f);
}
#endif
