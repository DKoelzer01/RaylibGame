#include "combat.h"
#include "playersession.h"
#include "simclock.h"
#include "worldtypes.h"
#include "voxel.h"
#include "scene.h"
#include "tools.h"
#include "hotbar.h"
#include "survival.h"
#include "navigation.h"
#include "gameui.h"
#include "net.h"
#include "itemdata.h"

#include "raymath.h"
#include "rlgl.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <vector>

// Tunables (#165 first slice; the numbers agreed on the issue).
static constexpr float kStaminaRegenPerSec = 25.0f;   // after kRegenDelay of not spending
static constexpr double kRegenDelay = 1.0;            // s
static constexpr float kBlockCostPerDamage = 1.5f;    // stamina per health point a block stops
static constexpr float kBlockArcCos = 0.5f;           // a block covers +-60 deg around the view (level)
static constexpr float kShieldBlockFraction = 0.75f;  // share of the damage a shield stops
static constexpr float kBlockedKnockback = 0.4f;      // share of the push left through a block
static constexpr float kKnockLift = 2.0f;             // m/s up with a push on the ground: off it, a short hop
static constexpr double kParryRecovery = 0.4;         // (#1776) s after a parry window before the next tap can parry (GUESS)

static bool testHoldBlock = false; // --test-combat: RMB held
static void PickUpGraves(PlayerSession& s);

// (#1218) Bucklers, round and tower shields (defense.h); RAYLIBGAME_COMBAT4=0: the wooden shield only, as before.
static bool IsShield(const std::string& id) { return Combat4_On() ? Shield_Of(id) != nullptr : id == "wooden_shield"; }
static const ShieldSpec& ShieldFor(const std::string& id) { static const ShieldSpec* round = Shield_Of("wooden_shield"); const ShieldSpec* s = Shield_Of(id); return s ? *s : *round; }

// Does RMB already do something with what is in the main hand? Then it keeps doing that, and a shield in the off
// hand does not block: food and feed (eat, #201), fillers (place), buckets, shovels and hoes (till), torches and
// snares (set down).
static bool RmbUsedByHand(const PlayerSession& s) {
    const std::string& held = s.inventory.grid()[Inventory::kGrid + s.hotbar.selected].id;
    if (held.empty()) return false;
    if (s.terrainTool.holdingItem || TerrainTool_IsFiller(held) || held == "bucket" || held == "water_bucket") return true;
    if (Hotbar_IsLight(held) || Hotbar_IsPlaceable(held)) return true;
    const ToolType t = ToolOf(held).type;
    return t == ToolType::Hoe; // (its RMB tills; #1085: a shovel no longer does)
}

void Combat_Step(PlayerSession& s, float dt, const PlayerInput& in) {
    SurvivalState& sv = s.survival;
    // The block: RMB with a shield in the off hand (and RMB not taken by the main hand or build mode), while stamina lasts.
    s.terrainTool.rmbBlocks = IsShield(s.inventory.offHand().id) && !s.buildTool.active && !RmbUsedByHand(s);
    if (Combat4_On()) { const ShieldSpec& sh = ShieldFor(s.inventory.offHand().id); sv.shieldStops = sh.stops; sv.shieldArcDeg = sh.arcDeg; sv.shieldParryWindow = sh.parryWindow; }
    // (#1776, user) An RMB tap parries, a hold blocks. The press (or the shield coming into the hand with RMB down) opens
    // the shield's parry window, whether or not RMB stays down; held past the window, or past a parry that landed, it is
    // the block. After each tap a short recovery before the next can parry (no spamming the window).
    const double now = SimTime();
    const bool held = s.terrainTool.rmbBlocks && (in.down(Act::Secondary) || testHoldBlock) && sv.stamina > 0.0f && !sv.dead;
    if (held && !sv.rmbHeld) {
        sv.rmbDownAt = now;
        if (Combat4_On() && now >= sv.parryReadyAt) {
            sv.parryStartedAt = now; sv.parryLanded = false;
            sv.parryReadyAt = now + sv.shieldParryWindow + kParryRecovery;
        }
    }
    sv.rmbHeld = held;
    sv.parrying = Combat4_On() && !sv.parryLanded && !sv.dead && now - sv.parryStartedAt <= sv.shieldParryWindow + 1e-6;
    const bool wasBlocking = sv.blocking;
    sv.blocking = held && (!Combat4_On() || sv.parryLanded || now - sv.rmbDownAt > sv.shieldParryWindow + 1e-6);
    if (sv.blocking && !wasBlocking) sv.blockStartedAt = now;
    // Stamina comes back a second after the last spend.
    if (SimTime() - sv.staminaUsedAt >= kRegenDelay) sv.stamina = std::min(100.0f, sv.stamina + kStaminaRegenPerSec * dt);
    // (#1215) Status effects: their damage over time, and chill's slow (frozen: no moving) on the walk speed.
    if (Combat1_On() && !sv.dead) {
        if (sv.status.interrupted) ++sv.shocks; // (#1216) seen by the melee wind-up before Status_Step clears it
        const float dot = Status_Step(sv.status, dt, s.player.inWater, sv.resist);
        if (dot > 0.0f && !sv.invulnerable) { // (#1561)
            sv.health -= dot;
            sv.lastHarm = sv.status.poison > 0.0f ? "poison" : sv.status.burn > 0.0f ? "burn" : "bleed"; // (#844 logs)
        }
        s.player.statusSpeed = Status_SpeedFactor(sv.status);
    } else s.player.statusSpeed = 1.0f;
    // (#1217) Poise comes back after a rest; a stagger slows the walk.
    if (Combat3_On()) {
        Poise_Step(sv.poise, dt);
        if (Poise_Staggered(sv.poise)) s.player.statusSpeed *= PoiseTune::kPlayerStaggerSlow;
    }
    // (#1218) Sprinting spends stamina; run dry, the player is winded; a raised tower shield slows; the dodge.
    if (Combat4_On() && !sv.dead) {
        Player& p = s.player;
        const bool moving = in.down(Act::Forward) || in.down(Act::Back) || in.down(Act::Left) || in.down(Act::Right);
        if (in.down(Act::Sprint) && moving && p.grounded && !p.flyMode) { sv.stamina = std::max(0.0f, sv.stamina - DefenseTune::kSprintStamina * dt); sv.staminaUsedAt = SimTime(); }
        if (sv.stamina <= 0.0f && !sv.staminaEmpty) sv.windedLeft = DefenseTune::kWindedTime;
        sv.staminaEmpty = sv.stamina <= 0.0f;
        if (sv.windedLeft > 0.0f) { sv.windedLeft = std::max(0.0f, sv.windedLeft - dt); p.statusSpeed *= DefenseTune::kWindedWalk; }
        if (sv.blocking) p.statusSpeed *= ShieldFor(s.inventory.offHand().id).raisedWalk;
        sv.dodgeInvuln = std::max(0.0f, sv.dodgeInvuln - dt);
        // A dodge: double-tap a direction (Settings, on by default), or the dodge key toward the held direction (back if none).
        static const Act kDirs[4] = { Act::Forward, Act::Back, Act::Left, Act::Right };
        int dir = -1;
        if (gDoubleTapDodge) for (int k = 0; k < 4 && dir < 0; ++k) if (in.hit(kDirs[k])) dir = sv.taps.press(k, SimTime());
        if (dir < 0 && in.hit(Act::Dodge)) { dir = 1; for (int k = 0; k < 4; ++k) if (in.down(kDirs[k])) { dir = k; break; } }
        if (dir >= 0 && p.grounded && !p.flyMode && p.dodgeLeft <= 0.0f && sv.windedLeft <= 0.0f && sv.stamina >= DefenseTune::kDodgeStamina) {
            const Vector3 flat = Vector3Normalize(Vector3Subtract(p.forward, Vector3Scale(p.up, Vector3DotProduct(p.forward, p.up))));
            const Vector3 side = Vector3CrossProduct(flat, p.up); // the walk's right
            const Vector3 dirs[4] = { flat, Vector3Negate(flat), Vector3Negate(side), side };
            p.dodgeDir = dirs[dir]; p.dodgeLeft = DefenseTune::kDodgeTime; p.dodgeSpeed = DefenseTune::kDodgeDistance / DefenseTune::kDodgeTime;
            sv.dodgeInvuln = DefenseTune::kDodgeInvuln;
            sv.stamina -= DefenseTune::kDodgeStamina; sv.staminaUsedAt = SimTime();
            ++sv.dodges;
        }
    }
    PickUpGraves(s);
}

void Combat_CurePoison(PlayerSession& s) { Status_CurePoison(s.survival.status); } // (#1215) a herb, when it exists

DamageResult Combat_Damage(PlayerSession& s, const DamageEvent& e) { return Combat_DamageFor(s.survival, s.player, e); }
// (#337) any player's state: the local session's, or a remote player's on a server (its block from its held input).
DamageResult Combat_DamageFor(SurvivalState& sv, Player& p, const DamageEvent& e) {
    DamageResult r;
    if (sv.dead || sv.invulnerable || e.amount <= 0.0f) return r; // (#1561: creative takes no hits, no knockback)
    // (#1218) A dodge's invulnerability: the hit passes through.
    if (Combat4_On() && sv.dodgeInvuln > 0.0f && e.kind != DamageKind::Other) { r.dodged = true; return r; }
    float amount = e.amount, push = e.knockback;
    // Away from the source, level with the ground where the player stands.
    Vector3 away = { 0, 0, 0 };
    if (e.hasSource) {
        const Vector3 d = (p.position - e.source).toFloat(); // a short difference: fine in float
        away = Vector3Subtract(d, Vector3Scale(p.up, Vector3DotProduct(d, p.up)));
        const float l = Vector3Length(away);
        away = l > 1e-4f ? Vector3Scale(away, 1.0f / l) : Vector3{ 0, 0, 0 };
    }
    // Facing the source: within the shield's arc around the level view.
    auto facing = [&] {
        if (!e.hasSource || Vector3LengthSqr(away) <= 0.0f) return false;
        const Vector3 fwd = Vector3Normalize(Vector3Subtract(p.forward, Vector3Scale(p.up, Vector3DotProduct(p.forward, p.up))));
        return Vector3DotProduct(fwd, Vector3Negate(away)) >= (Combat4_On() ? cosf(sv.shieldArcDeg * DEG2RAD) : kBlockArcCos);
    };
    // (#1218, #1776) A parry: the tap's window overlaps the attack's live span [t - hitWindow, t] (the window may have just
    // closed as the blow lands), facing it. A melee hit, a bite or a projectile is stopped whole, for no stamina (a melee
    // attacker staggers: animals.cpp). One parry per tap.
    if (Combat4_On() && e.parryable && e.kind != DamageKind::Other && !sv.parryLanded && facing()) {
        const double t = SimTime(), ps = sv.parryStartedAt;
        if (ps <= t + 1e-6 && ps + sv.shieldParryWindow >= t - e.hitWindow - 1e-6) {
            r.parried = true; ++sv.parries;
            sv.parryLanded = true; sv.parrying = false;
            return r;
        }
    }
    // A block stops part of a hit it faces while stamina lasts.
    if (sv.blocking && facing()) {
        {
            const float stopped = amount * (Combat4_On() ? sv.shieldStops : kShieldBlockFraction);
            const float cost = stopped * kBlockCostPerDamage;
            if (sv.stamina >= cost) {
                amount -= stopped;
                push *= kBlockedKnockback;
                r.blocked = true;
            } else {
                // Not enough stamina left: the guard breaks, the hit lands whole and the block drops.
                r.guardBroken = true;
                sv.blocking = false;
            }
            sv.stamina = std::max(0.0f, sv.stamina - cost);
            sv.staminaUsedAt = SimTime();
        }
    }
    if (Combat1_On()) { // (#1215) the hit's type through the target's resistances and armour (the 10% floor), then its effect
        if (e.type == DamageType::Energy && p.inWater) amount *= StatusTune::kShockWaterFactor; // shock: x2 in water
        amount = Damage_Apply(amount, e.type, sv.resist);
        Status_Apply(sv.status, Status_FromHit(e.type, e.status, e.heavy, e.fromBehind, e.kind == DamageKind::Bite));
    }
    // (#1217) Poise from armour: the hit wears it (a block takes its share off), and broken the player staggers (a short slow,
    // a melee wind-up broken) unless a heavy or two-handed wind-up holds it (hyper-armour).
    if (Combat3_On()) {
        const float pmax = Poise_PlayerMax(sv.resist.armour);
        if (sv.poise.max != pmax) { const float f = sv.poise.max > 0.0f ? sv.poise.value / sv.poise.max : 1.0f; sv.poise.max = pmax; sv.poise.value = pmax * f; }
        float pd = e.poise > 0.0f ? e.poise : e.amount * PoiseTune::kAttackPoisePerDamage;
        if (r.blocked) pd *= 1.0f - (Combat4_On() ? sv.shieldStops : kShieldBlockFraction);
        if (Poise_Hit(sv.poise, pd, e.heavy, false, sv.hyperArmour, PoiseTune::kPlayerStaggerTime).staggered) { ++sv.staggers; r.staggered = true; }
    }
    sv.health -= amount;
    sv.hitDamageTaken += amount;
    sv.lastHarm = "hit"; // (#844)
    sv.lastHitAt = SimTime();
    sv.lastHitStrength = std::min(1.0f, amount / 20.0f);
    r.dealt = amount;
    r.killed = sv.health <= 0.0f; // Survival_Update turns this into a death on its next step
    if (push > 0.0f && Vector3LengthSqr(away) > 0.0f && !p.flyMode)
        p.pendingPush = Vector3Add(p.pendingPush, Vector3Add(Vector3Scale(away, push), Vector3Scale(p.up, kKnockLift)));
    return r;
}

// Death losses by difficulty (#165). Peaceful and normal keep the inventory with a small loss; hard leaves all of it in
// a grave at the death spot. Graves are world state (any player who died), kept until emptied, and saved.
struct Grave {
    std::string body;
    DVec3 local{};                   // feet at death, body-local (double: planet scale)
    std::string name;                // its waypoint ("Grave 1", ...)
    std::vector<InvSlot> items;      // what it holds; off-hand items go back to the off hand
    int offHand = -1;                // index in items of the off hand's, -1 none
    // A pile of items the player dropped (#564: a stack let go outside the inventory window) rather than a grave: no
    // waypoint, drawn as a small bundle, and picked up only once the player has walked away from it and comes back
    // (`armed`), so dropping something to make room doesn't take it straight back.
    bool dropped = false, armed = true;
    // Multiplayer (#804): whose it is ("" single player: anyone's), when they died and when anyone of their team was last
    // online (wall clock, s), a stable id (the replication's) and its stone (0 cross, 1 headstone, 2 monolith).
    std::string owner;
    int64_t diedAt = 0, ownerSeen = 0;
    uint32_t id = 0;
    uint8_t style = 0;
};
static std::vector<Grave> graves;
static int graveCounter = 0;
static uint32_t graveVersion = 1;          // bumps on every change of the grave set (the server resends its list)
static constexpr double kGraveReach = 2.5; // m: walking up to it takes everything back
static constexpr double kDropMerge = 1.0;  // m: a drop this close to an unarmed pile joins it
static int64_t WallNow() { return (int64_t)time(nullptr); }

uint8_t Combat_GraveStyle(const std::string& body, const DVec3& local, uint32_t id) {
    // Seeded by where and which it is: every machine (and every reload) picks the same stone.
    uint64_t h = 1469598103934665603ull;
    auto mix = [&](uint64_t v) { h ^= v; h *= 1099511628211ull; };
    for (char c : body) mix((uint8_t)c);
    mix((uint64_t)llround(local.x * 10.0)); mix((uint64_t)llround(local.y * 10.0)); mix((uint64_t)llround(local.z * 10.0)); mix(id);
    h ^= h >> 29;
    return (uint8_t)(h % 3);
}

bool Combat_SameTeam(const std::string& a, const std::string& b) {
    if (a == b) return true;
    return gSameTeam ? gSameTeam(a, b) : false; // teams: #343
}

bool Combat_CanLoot(const std::string& owner, int64_t diedAt, const std::string& who, int64_t now, bool inOwnerClaim) {
    if (owner.empty() || Combat_SameTeam(owner, who)) return true;
    return !inOwnerClaim && now - diedAt >= kGraveOpenToAllSeconds;
}

void Combat_OnDeath(PlayerSession& s) { Combat_OnDeathFor(s.inventory, &s.toolWear, s.player, std::string(), DeathDrops::Difficulty, true); }

bool Combat_OnDeathFor(Inventory& inv, ToolWearState* wear, const Player& p, const std::string& owner, DeathDrops mode, bool local) {
    auto& cells = inv.grid();
    if (mode == DeathDrops::Keep) { if (local) GameUI_Toast("You kept everything you carried"); return false; }
    if ((mode == DeathDrops::Grave || worldDifficulty >= 2) && p.nearest) {
        Grave g;
        g.body = p.nearest->name;
        g.local = p.nearest->toLocalD(p.position);
        for (int i = 0; i < Inventory::kCells; ++i) {
            if (cells[i].empty()) continue;
            if (i == Inventory::kOffHand) g.offHand = (int)g.items.size();
            g.items.push_back(cells[i]);
            if (cells[i].hasData()) inv.takeInstance(i); // (#1446) logged for a remote player's client, record kept in the grave
            else cells[i] = InvSlot{};
        }
        if (g.items.empty()) return false;
        g.id = (uint32_t)++graveCounter;
        g.name = owner.empty() ? std::string(TextFormat("Grave %d", graveCounter)) : owner + "'s grave " + std::to_string(graveCounter);
        g.owner = owner;
        g.diedAt = g.ownerSeen = WallNow();
        g.style = Combat_GraveStyle(g.body, g.local, g.id);
        if (local) {
            const DVec3 l = g.local;
            Nav_AddWaypoint(g.body, Vector3{ (float)l.x, (float)l.y, (float)l.z }, g.name, false); // normalised there
            GameUI_Toast(TextFormat("Everything you carried lies in your grave (%s, marked on the compass)", g.name.c_str()));
        }
        graves.push_back(std::move(g));
        ++graveVersion;
        return true;
    }
    // Peaceful and normal: each tool wears by a tenth of its life (never past breaking); on normal each resource stack
    // also loses a tenth, rounded down to whole units.
    int lost = 0;
    for (int i = 0; i < Inventory::kCells; ++i) {
        InvSlot& c = cells[i];
        if (c.empty()) continue;
        const ItemKind k = KindOf(c.id);
        if (k == ItemKind::Tool) {
            const int max = ToolDurability(ToolOf(c.id));
            if (max > 0 && wear) { int& u = wear->used[c.id]; u = std::min(u + max / 10, max - 1); }
        } else if (k == ItemKind::Resource && worldDifficulty >= 1) {
            const int64_t units = (int64_t)(c.amount() * 0.1);
            if (units > 0) { c.centi -= units * 100; lost += (int)units; }
        }
    }
    if (local) GameUI_Toast(lost > 0 ? TextFormat("You lost %d items and your tools took wear", lost) : "Your tools took wear");
    return false;
}

static void PickUpGraves(PlayerSession& s) {
    if (gNetRole == NetRole::Client) return; // (#804) the server picks up for us; our list is only its mirror
    if (s.player.nearest) Combat_PickUpGravesFor(s.player.nearest->name, s.player.nearest->toLocalD(s.player.position), s.inventory, std::string(), s.survival.dead, true);
}

bool Combat_PickUpGravesFor(const std::string& body, const DVec3& me, Inventory& inv, const std::string& who, bool dead, bool local) {
    if (body.empty() || dead) return false;
    const int64_t now = WallNow();
    bool moved = false;
    for (size_t gi = 0; gi < graves.size();) {
        Grave& g = graves[gi];
        const bool near = g.body == body && (g.local - me).length() <= kGraveReach;
        // A dropped pile waits until its dropper has walked away from it (the dropper's own walk arms it: on a server another
        // player's pass must not, #865); anyone else may take it 2 s after the drop, as in Minecraft.
        if (g.dropped && !g.armed) {
            if (who == g.owner) { if (!near) g.armed = true; ++gi; continue; }
            if (now - g.diedAt < 2) { ++gi; continue; }
        }
        if (!near) { ++gi; continue; }
        // (#804) whose: the owner and their team at once, others after 15 minutes and never inside the owner's team claim
        if (!g.dropped && !Combat_CanLoot(g.owner, g.diedAt, who, now, gInOwnerTeamClaim && gInOwnerTeamClaim(g.owner, g.body, g.local))) { ++gi; continue; }
        auto& cells = inv.grid();
        std::vector<InvSlot> left;
        for (size_t i = 0; i < g.items.size(); ++i) {
            const InvSlot& it = g.items[i];
            if ((int)i == g.offHand && cells[Inventory::kOffHand].empty()) {
                cells[Inventory::kOffHand] = it; moved = true;
                if (it.hasData() && inv.journal) inv.instanceLog.push_back({ it.id, +1, it.data, {} }); // (#1446)
                continue;
            }
            if (it.hasData()) { // (#1446) that one item back, record and all; no room: it stays in the grave
                if (!g.dropped) ++inv.muted;
                const int cell = inv.addInstance(it.id, it.data, it.heat);
                if (!g.dropped) --inv.muted;
                inv.overflowed = false;
                if (cell >= 0) moved = true; else left.push_back(it);
                continue;
            }
            if (!g.dropped) ++inv.muted; // a grave: no pickup feed line per stack (a dropped pile shows them)
            const double rest = inv.addAmount(it.id, it.amount());
            if (!g.dropped) --inv.muted;
            inv.overflowed = false;
            if (rest < it.amount()) moved = true;
            if (rest > 0.0) left.push_back(InvSlot{ it.id, (int64_t)llround(rest * 100.0) });
        }
        if (left.empty()) {
            if (!g.dropped && local) {
                Nav_RemoveWaypoint(g.body, g.name);
                GameUI_Toast(TextFormat("You took everything back from %s", g.name.c_str()));
            }
            ++graveVersion; // (a dropped pile too: clients see it go, #865)
            graves.erase(graves.begin() + gi);
            continue;
        }
        if (g.items.size() != left.size()) ++graveVersion;
        g.items = std::move(left); g.offHand = -1;
        if (local) {
            if (g.dropped) GameUI_Toast("No room: the rest stays on the ground");
            else GameUI_Toast(TextFormat("No room for the rest: %d stacks stay in %s", (int)g.items.size(), g.name.c_str()));
        }
        if (g.dropped) g.armed = false; // again after walking away
        ++gi;
    }
    return moved;
}

int Combat_UpdateGraveOwners(const std::vector<std::string>& online, int64_t now) {
    int removed = 0;
    for (size_t gi = 0; gi < graves.size();) {
        Grave& g = graves[gi];
        if (g.dropped || g.owner.empty()) { ++gi; continue; }
        for (const std::string& n : online) if (Combat_SameTeam(g.owner, n)) { g.ownerSeen = now; break; }
        if (now - g.ownerSeen > kGraveDecaySeconds) { graves.erase(graves.begin() + gi); ++removed; ++graveVersion; continue; }
        ++gi;
    }
    return removed;
}

std::vector<GraveInfo> Combat_GraveList() {
    std::vector<GraveInfo> out;
    for (const Grave& g : graves) {
        GraveInfo gi{ g.id, g.body, g.local, g.owner, g.name, g.style, g.diedAt };
        if (g.dropped) { // (#865) dropped piles too, with what they hold for their bundle's colours
            gi.dropped = true;
            for (size_t k = 0; k < g.items.size() && k < 3; ++k) gi.items.push_back(g.items[k].id);
        }
        out.push_back(std::move(gi));
    }
    return out;
}

uint32_t Combat_GraveVersion() { return graveVersion; }

void Combat_TestAddGrave(const std::string& body, const DVec3& local, uint8_t style, const std::string& owner) {
    Grave g;
    g.body = body; g.local = local; g.owner = owner; g.style = style % 3;
    g.id = (uint32_t)++graveCounter; g.name = TextFormat("Grave %d", graveCounter);
    g.diedAt = g.ownerSeen = WallNow();
    g.items.push_back(InvSlot{ "stone", 100 });
    graves.push_back(std::move(g));
    ++graveVersion;
}

void Combat_SetGraveMirror(const std::vector<GraveInfo>& list, const std::string& me) {
    // A client: the server's graves and dropped piles (#865: everyone's, drawn here; the server picks them up) replace ours.
    // New graves of our own or our team's get a waypoint (and our own a toast); ones gone lose theirs.
    std::vector<Grave> before;
    for (Grave& g : graves) if (!g.dropped) before.push_back(std::move(g));
    graves.clear();
    auto mine = [&](const std::string& owner) { return !owner.empty() && Combat_SameTeam(owner, me); };
    for (const GraveInfo& gi : list) {
        Grave g;
        g.id = gi.id; g.body = gi.body; g.local = gi.local; g.owner = gi.owner; g.name = gi.name; g.style = gi.style % 3; g.diedAt = gi.diedAt;
        if (gi.dropped) {
            g.dropped = true;
            for (const std::string& id : gi.items) g.items.push_back(InvSlot{ id, 100 }); // (for its colours: the server has the amounts)
            graves.push_back(std::move(g));
            continue;
        }
        const bool known = std::any_of(before.begin(), before.end(), [&](const Grave& b) { return b.id == gi.id; });
        if (!known && mine(g.owner)) {
            Nav_AddWaypoint(g.body, Vector3{ (float)g.local.x, (float)g.local.y, (float)g.local.z }, g.name, false);
            if (g.owner == me) GameUI_Toast(TextFormat("Everything you carried lies in your grave (%s, marked on the compass)", g.name.c_str()));
        }
        graves.push_back(std::move(g));
    }
    for (const Grave& b : before) {
        const bool still = std::any_of(list.begin(), list.end(), [&](const GraveInfo& gi) { return gi.id == b.id; });
        if (!still && mine(b.owner)) Nav_RemoveWaypoint(b.body, b.name);
    }
    ++graveVersion;
}

int Combat_GraveCount() { return (int)std::count_if(graves.begin(), graves.end(), [](const Grave& g) { return !g.dropped; }); }

bool Combat_GraveMark(const std::string& body, const std::string& name, std::string& label) {
    for (const Grave& g : graves) {
        if (g.dropped || g.body != body || g.name != name) continue;
        label = g.name;
        if (!g.owner.empty()) { // (#828) multiplayer: the owner's team alone until kGraveOpenToAllSeconds after the death
            const int64_t left = g.diedAt + kGraveOpenToAllSeconds - WallNow();
            label += left > 0 ? TextFormat(", open to all in %lld min", (long long)((left + 59) / 60)) : ", open to all";
        }
        return true;
    }
    return false;
}


// (#1225) Thieves: the nearest dropped pile (not a grave) within `radius` of `at` on `body`; and taking one (its items into `out`).
bool Combat_PileNear(const std::string& body, Vector3 at, float radius, Vector3& pileAt) {
    double best = radius;
    bool found = false;
    for (const Grave& g : graves) {
        if (!g.dropped || g.body != body) continue;
        const double d = (g.local - DVec3(at)).length();
        if (d < best) { best = d; pileAt = g.local.toFloat(); found = true; }
    }
    return found;
}

bool Combat_TakePileNear(const std::string& body, Vector3 at, float radius, std::vector<InvSlot>& out) {
    int bestI = -1; double best = radius;
    for (size_t i = 0; i < graves.size(); ++i) {
        const Grave& g = graves[i];
        if (!g.dropped || g.body != body) continue;
        const double d = (g.local - DVec3(at)).length();
        if (d < best) { best = d; bestI = (int)i; }
    }
    if (bestI < 0) return false;
    out.insert(out.end(), graves[bestI].items.begin(), graves[bestI].items.end());
    graves.erase(graves.begin() + bestI);
    ++graveVersion;
    return true;
}

bool Combat_DropFor(const std::string& body, const DVec3& at, const InvSlot& stack, const std::string& who) {
    if (stack.empty() || body.empty()) return false;
    ++graveVersion; // (#865: a server's clients see it)
    for (Grave& g : graves) // onto a pile this player just dropped here
        if (g.dropped && !g.armed && g.owner == who && g.body == body && (g.local - at).length() < kDropMerge) {
            auto it = std::find_if(g.items.begin(), g.items.end(), [&](const InvSlot& c) { return c.id == stack.id && !c.hasData() && !stack.hasData(); }); // (#1730) a record never merges
            if (it != g.items.end()) it->centi += stack.centi; else g.items.push_back(stack);
            g.diedAt = WallNow();
            return true;
        }
    static uint32_t dropCounter = 0;
    Grave g;
    g.body = body;
    g.local = at;
    g.items.push_back(stack);
    g.dropped = true; g.armed = false;
    g.owner = who; g.diedAt = WallNow();          // who dropped it, and when (others may take it 2 s later)
    g.id = 0x80000000u | ++dropCounter;           // (apart from graves' ids, which also number their names)
    graves.push_back(std::move(g));
    return true;
}

bool Combat_Drop(PlayerSession& s, const InvSlot& stack) {
    if (stack.empty() || !s.player.nearest) return false;
    if (gNetRole == NetRole::Client) { // (#865) the server drops it from the inventory it keeps; its delta takes it from ours
        if (!gDropRequest) return false;
        ++s.inventory.muted; // (back first, as MachineAct_Insert; #1730: an item with a record goes back as itself)
        if (stack.hasData()) s.inventory.addInstance(stack.id, stack.data, stack.heat); else s.inventory.addAmount(stack.id, stack.amount(), stack.heat);
        --s.inventory.muted;
        gDropRequest(stack.id, stack.centi, stack.data);
        return true;
    }
    return Combat_DropFor(s.player.nearest->name, s.player.nearest->toLocalD(s.player.position), stack, std::string()); // (at the feet)
}

int Combat_DroppedPiles() { return (int)graves.size() - Combat_GraveCount(); }

// Hit feel (#165), the local view only: the camera shakes for 0.3 s after a hit, up to 6 cm at 20 damage, and the
// held tool freezes for 70 ms (a hit-stop the simulation never sees). Reduce motion turns both off.
static constexpr double kShakeTime = 0.3, kHitStopTime = 0.07;
static constexpr float kShakeAmplitude = 0.06f; // m

void Combat_ShakeCamera(Scene& world, const PlayerSession& s, bool reduceMotion) {
    const double t = SimTime() - s.survival.lastHitAt;
    if (reduceMotion || t < 0.0 || t > kShakeTime || s.survival.lastHitStrength <= 0.0f) return;
    const float a = kShakeAmplitude * s.survival.lastHitStrength * (float)(1.0 - t / kShakeTime);
    const Vector3 fwd = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
    const Vector3 right = Vector3Normalize(Vector3CrossProduct(fwd, world.camera.up));
    const float tt = (float)t;
    const Vector3 off = Vector3Add(Vector3Scale(right, a * sinf(tt * 71.0f)), Vector3Scale(world.camera.up, a * 0.7f * sinf(tt * 53.0f + 1.3f)));
    world.camera.position = Vector3Add(world.camera.position, off);
    world.camera.target = Vector3Add(world.camera.target, off);
}

bool Combat_HitStop(const PlayerSession& s, bool reduceMotion) {
    const double t = SimTime() - s.survival.lastHitAt, dealt = SimTime() - s.survival.lastDealtAt; // (#1216) and hits landed
    return !reduceMotion && ((t >= 0.0 && t < kHitStopTime && s.survival.lastHitStrength > 0.0f) || (dealt >= 0.0 && dealt < kHitStopTime));
}

void Combat_DrawGraves(Scene& world) {
    // A small cairn of stones with a red marker on top, drawn in the body's frame (double anchor: no jitter at 60 km).
    for (const Grave& g : graves) {
        Planetoid* b = nullptr;
        for (auto& o : world.objects) if (auto* pb = dynamic_cast<Planetoid*>(o.get()); pb && pb->name == g.body) { b = pb; break; }
        if (!b) continue;
        const Matrix m = b->localModel(g.local, world.renderOrigin);
        const Vector3 at = { m.m12, m.m13, m.m14 }; // the anchor in render space
        if (Vector3Distance(at, world.camera.position) > 300.0f) continue;
        const DVec3 l = g.local;
        const Vector3 up = Vector3Normalize(Vector3{ (float)l.x, (float)l.y, (float)l.z });
        rlDrawRenderBatchActive();
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(m));
        if (g.dropped) { // a small bundle in the colours of what it holds (up to three stacks shown)
            for (size_t i = 0; i < g.items.size() && i < 3; ++i) {
                const ItemDef* d = itemDB.find(g.items[i].id);
                const float a = (float)i * 2.1f;
                const Vector3 side = Vector3Normalize(Vector3Perpendicular(up));
                const Vector3 off = Vector3Add(Vector3Scale(up, 0.12f + 0.1f * (float)i),
                                               Vector3Scale(Vector3RotateByAxisAngle(side, up, a), i ? 0.12f : 0.0f));
                DrawCube(off, 0.22f, 0.22f, 0.22f, d ? d->color : Color{ 150, 120, 80, 255 });
                DrawCubeWires(off, 0.22f, 0.22f, 0.22f, Color{ 40, 30, 20, 255 });
            }
            rlDrawRenderBatchActive();
            rlPopMatrix();
            continue;
        }
        // (#804, user 2026-11-18: "a simple stone cross / headstone / monolith, randomly picked") - plain stone, its kind and
        // its turn seeded by the grave (every machine draws the same one). A frame with y = up, turned about it.
        const Vector3 side0 = Vector3Normalize(Vector3Perpendicular(up));
        const float turn = (float)((g.id * 2654435761u) % 360u) * DEG2RAD;
        const Vector3 side = Vector3RotateByAxisAngle(side0, up, turn), fwd = Vector3CrossProduct(side, up);
        const Matrix basis = { side.x, up.x, fwd.x, 0.0f, side.y, up.y, fwd.y, 0.0f, side.z, up.z, fwd.z, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f };
        rlMultMatrixf(MatrixToFloat(basis));
        const Color stone = { 128, 126, 122, 255 }, stoneLight = { 150, 148, 143, 255 }, edge = { 70, 68, 66, 255 };
        auto block = [&](Vector3 at, Vector3 size, Color c) { DrawCubeV(at, size, c); DrawCubeWiresV(at, size, edge); };
        switch (g.style % 3) {
            case 0: // a stone cross on a low plinth
                block({ 0.0f, 0.02f, 0.0f }, { 0.42f, 0.16f, 0.3f }, stone);
                block({ 0.0f, 0.55f, 0.0f }, { 0.14f, 1.0f, 0.12f }, stoneLight);
                block({ 0.0f, 0.78f, 0.0f }, { 0.56f, 0.13f, 0.12f }, stoneLight);
                break;
            case 1: // a headstone: a slab, stepped narrower at the top, on a base
                block({ 0.0f, 0.02f, 0.0f }, { 0.78f, 0.14f, 0.34f }, stone);
                block({ 0.0f, 0.42f, 0.0f }, { 0.62f, 0.68f, 0.14f }, stoneLight);
                block({ 0.0f, 0.81f, 0.0f }, { 0.44f, 0.12f, 0.14f }, stoneLight);
                break;
            default: // a monolith: a tall rough block with a smaller cap
                block({ 0.0f, 0.62f, 0.0f }, { 0.36f, 1.36f, 0.32f }, stone);
                block({ 0.0f, 1.35f, 0.0f }, { 0.26f, 0.12f, 0.24f }, stoneLight);
                break;
        }
        rlDrawRenderBatchActive();
        rlPopMatrix();
    }
}

void Combat_Reset() { graves.clear(); graveCounter = 0; }

// (#1730) A pile's or grave's item as one save token: "id", or "id|record" for an item with its own data (ItemData
// records hold no spaces or '|'; older saves have only the id, so they still load).
static std::string SavedItemToken(const InvSlot& it) { return it.hasData() ? it.id + "|" + it.data : it.id; }

void Combat_Save(FILE* f) {
    fprintf(f, "grave_counter %d\n", graveCounter);
    for (const Grave& g : graves) {
        if (g.dropped) { // #564
            fprintf(f, "dropped %s %.6f %.6f %.6f %d %zu", g.body.c_str(), g.local.x, g.local.y, g.local.z, g.armed ? 1 : 0, g.items.size());
            for (const InvSlot& it : g.items) fprintf(f, " %s %lld", SavedItemToken(it).c_str(), (long long)it.centi);
            fprintf(f, "\n");
            continue;
        }
        fprintf(f, "grave %s %.6f %.6f %.6f %d %zu", g.body.c_str(), g.local.x, g.local.y, g.local.z, g.offHand, g.items.size());
        for (const InvSlot& it : g.items) fprintf(f, " %s %lld", SavedItemToken(it).c_str(), (long long)it.centi);
        fprintf(f, " %s\n", g.name.c_str()); // last: the name has a space in it
        // (#804) whose, when, the stone: a line of its own after the grave's (older saves lack it: anyone's, now)
        fprintf(f, "grave_mp %u %lld %lld %d %s\n", g.id, (long long)g.diedAt, (long long)g.ownerSeen, (int)g.style, g.owner.empty() ? "-" : g.owner.c_str());
    }
}

bool Combat_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "grave_counter") == 0) { if (fscanf(f, " %d", &graveCounter) != 1) graveCounter = 0; return true; }
    if (strcmp(tok, "grave_mp") == 0) { // (#804) for the grave just read
        unsigned id = 0; long long died = 0, seen = 0; int style = 0; char owner[64] = {};
        if (fscanf(f, " %u %lld %lld %d %63s", &id, &died, &seen, &style, owner) == 5 && !graves.empty() && !graves.back().dropped) {
            Grave& g = graves.back();
            g.id = id; g.diedAt = died; g.ownerSeen = seen; g.style = (uint8_t)(style % 3);
            g.owner = strcmp(owner, "-") == 0 ? std::string() : std::string(owner);
        }
        ++graveVersion;
        return true;
    }
    const bool dropped = strcmp(tok, "dropped") == 0;
    if (strcmp(tok, "grave") != 0 && !dropped) return false;
    Grave g;
    char body[64] = {};
    size_t n = 0;
    int flag = -1; // a grave's off-hand index, a dropped pile's armed
    if (fscanf(f, " %63s %lf %lf %lf %d %zu", body, &g.local.x, &g.local.y, &g.local.z, &flag, &n) != 6) return true;
    if (dropped) { g.dropped = true; g.armed = flag != 0; } else g.offHand = flag;
    g.body = body;
    for (size_t i = 0; i < n; ++i) {
        char tok[512] = {};
        long long c = 0;
        if (fscanf(f, " %511s %lld", tok, &c) != 2) return true;
        char id[64] = {};
        std::string data; // (#1730) "id|record"
        if (const char* bar = strchr(tok, '|')) { snprintf(id, sizeof id, "%.*s", (int)(bar - tok), tok); data = bar + 1; if (!ItemData::Valid(data)) data.clear(); }
        else snprintf(id, sizeof id, "%s", tok);
        if (const char* now = RetiredItem(id); *now) snprintf(id, sizeof id, "%s", now); // #1090
        if (itemDB.find(id) && c > 0) { InvSlot s{ id, (int64_t)c }; s.data = data; g.items.push_back(s); }
        else if ((int)i == g.offHand) g.offHand = -1;
    }
    if (!dropped) { // a pile has no name
        char name[64] = {};
        if (fscanf(f, " %63[^\n]", name) != 1) return true;
        g.name = name;
    }
    if (g.offHand >= (int)g.items.size()) g.offHand = -1;
    if (!dropped) { // (#804) an older save's grave: an id of its own and a stone (a grave_mp line, when there is one, sets both)
        g.id = 0x80000000u + (uint32_t)graves.size();
        g.style = Combat_GraveStyle(g.body, g.local, g.id);
        g.diedAt = g.ownerSeen = WallNow();
    }
    if (!g.items.empty()) { graves.push_back(std::move(g)); ++graveVersion; }
    return true;
}

// --test-combat: knockback, the shield block (front, behind, guard break, stamina and its regen delay), no block
// without a shield, and a killing hit, on the local player at the spawn, in sim time.
bool Combat_TestFrame(Scene& world) {
    PlayerSession& s = LocalSession();
    SurvivalState& sv = s.survival;
    Player& p = s.player;
    InvSlot& off = s.inventory.grid()[Inventory::kOffHand];
    static int step = 0;
    static double t0 = 0.0;
    static bool ok = true;
    static DVec3 local0{};
    static float pushPlain = 0.0f;
    static int savedDifficulty = 1;
    const double now = SimTime();
    auto check = [&](bool c, const char* what) { ok &= c; printf("bench: combat: %s: %s\n", what, c ? "yes" : "NO"); };
    auto local = [&]() { return p.nearest ? p.nearest->toLocalD(p.position) : p.position; };
    // A source 2 m away along the level view (sign +1: in front, -1: behind), world space.
    auto sourceAt = [&](float sign) {
        const Vector3 f = Vector3Normalize(Vector3Subtract(p.forward, Vector3Scale(p.up, Vector3DotProduct(p.forward, p.up))));
        return p.position + Vector3Scale(f, 2.0f * sign);
    };
    auto hit = [&](float amount, float sign, float knock) {
        DamageEvent e;
        e.amount = amount; e.hasSource = true; e.source = sourceAt(sign); e.knockback = knock;
        return Combat_Damage(s, e);
    };
    // A fixed sideways direction at a body-local point, and the direction of a point `dist` m along it (spawnOn takes one).
    auto sideOf = [](Vector3 up) {
        Vector3 side = Vector3CrossProduct(up, Vector3{ 0, 1, 0 });
        if (Vector3LengthSqr(side) < 1e-6f) side = Vector3CrossProduct(up, Vector3{ 1, 0, 0 });
        return Vector3Normalize(side);
    };
    auto standOff = [&](const DVec3& l, float dist) {
        const Vector3 up = Vector3Normalize(Vector3{ (float)l.x, (float)l.y, (float)l.z });
        return Vector3Normalize(Vector3Add(Vector3Scale(up, (float)l.length()), Vector3Scale(sideOf(up), dist)));
    };
    auto finish = [&]() {
        off = InvSlot{};
        testHoldBlock = false;
        printf("bench: combat: %s\n", ok ? "PASS" : "FAIL");
        step = 99;
        return true;
    };
    switch (step) {
    case 0: // benches start in the free-fly camera: out of it, standing at the spawn, bare hands, full bars
        if (t0 == 0.0) {
            t0 = now;
            Planetoid* start = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get());
            const auto* pf = start ? dynamic_cast<const PlanetField*>(start->field.get()) : nullptr;
            if (!start || !pf) { printf("bench: combat: no planet\n"); ok = false; return finish(); }
            p.flyMode = false;
            p.spawnOn(world, *start, start->dirToWorld(pf->findSpawnDir(start->seed)));
            s.inventory.grid()[Inventory::kGrid + s.hotbar.selected] = InvSlot{};
            return false;
        }
        if (!p.grounded || !p.nearest) {
            if (now - t0 < 20.0) return false;
            printf("bench: combat: the player never stood on the ground (grounded %d, body %s, fly %d)\n", p.grounded, p.nearest ? "yes" : "none", p.flyMode);
            ok = false;
            return finish();
        }
        sv.health = sv.stamina = 100.0f;
        local0 = local();
        {
            const DamageResult r = hit(10.0f, 1.0f, 6.0f);
            check(std::fabs(r.dealt - 10.0f) < 1e-3f && !r.blocked && std::fabs(sv.health - 90.0f) < 1e-3f, "an unblocked hit takes its damage");
            check(Combat_HitStop(s, false) && !Combat_HitStop(s, true), "a hit-stop follows the hit (none with Reduce motion)");
        }
        t0 = now; step = 1;
        return false;
    case 1: // the push: away from the source
        if (now - t0 < 0.6) return false;
        pushPlain = (float)(local() - local0).length();
        printf("bench: combat: knocked back %.2f m\n", pushPlain);
        check(pushPlain > 1.0f, "knockback pushes the player away");
        // RMB held with bare hands and no shield: no block.
        testHoldBlock = true;
        t0 = now; step = 2;
        return false;
    case 2:
        if (now - t0 < 0.1) return false;
        check(!sv.blocking && !s.terrainTool.rmbBlocks, "no shield in the off hand: RMB does not block");
        off.id = "wooden_shield"; off.centi = 100;
        t0 = now; step = 3;
        return false;
    case 3: // a shield in the off hand, RMB held, landed again: a hit from in front
        if (now - t0 < 0.6 || (!p.grounded && now - t0 < 3.0)) return false; // (#1218, #1776) held past the parry window: an ordinary block
        check(sv.blocking && s.terrainTool.rmbBlocks, "a shield in the off hand: RMB blocks (and builds nothing)");
        sv.health = sv.stamina = 100.0f;
        local0 = local();
        {
            const DamageResult r = hit(20.0f, 1.0f, 6.0f);
            printf("bench: combat: blocked in front: dealt %.1f, stamina %.1f\n", r.dealt, sv.stamina);
            check(r.blocked && std::fabs(r.dealt - 5.0f) < 1e-3f && std::fabs(sv.stamina - 77.5f) < 1e-3f, "the shield stops three quarters of a hit in front, for 1.5 stamina a point");
        }
        t0 = now; step = 4;
        return false;
    case 4: {
        if (now - t0 < 0.6) return false;
        const float pushBlocked = (float)(local() - local0).length();
        printf("bench: combat: knocked back through a block %.2f m (unblocked %.2f m)\n", pushBlocked, pushPlain);
        check(pushBlocked < 0.7f * pushPlain, "a block takes most of the push");
        const double since = now - sv.staminaUsedAt;
        printf("bench: combat: stamina %.2f s after the block: %.1f\n", since, sv.stamina);
        if (since < kRegenDelay) check(std::fabs(sv.stamina - 77.5f) < 0.01f, "no stamina back in the first second");
        else printf("bench: combat: (a slow frame: the first-second check is past, skipped)\n");
        step = 5;
        return false;
    }
    case 5: { // regen 25 a second after the first second (checked against the time that really passed)
        const double since = now - sv.staminaUsedAt;
        if (since < 1.6) return false;
        const float expect = std::min(100.0f, 77.5f + 25.0f * (float)(since - kRegenDelay));
        printf("bench: combat: stamina %.2f s after the block: %.1f (expected %.1f)\n", since, sv.stamina, expect);
        check(std::fabs(sv.stamina - expect) < 1.0f, "then 25 stamina a second");
        step = 6;
        return false;
    }
    case 6: { // landed: a hit from behind, a guard break, then a killing hit
        if (!p.grounded && now - t0 < 6.0) return false;
        sv.health = sv.stamina = 100.0f;
        const DamageResult behind = hit(20.0f, -1.0f, 0.0f);
        check(!behind.blocked && std::fabs(behind.dealt - 20.0f) < 1e-3f, "a block does not cover hits from behind");
        sv.health = 100.0f; sv.stamina = 2.0f;
        const DamageResult broken = hit(20.0f, 1.0f, 0.0f);
        check(broken.guardBroken && std::fabs(broken.dealt - 20.0f) < 1e-3f && !sv.blocking, "out of stamina the guard breaks");
        off = InvSlot{};
        testHoldBlock = false;
        savedDifficulty = worldDifficulty;
        worldDifficulty = 1; // the first death: normal, carrying nothing
        sv.blocking = false; sv.health = 5.0f; sv.stamina = 100.0f;
        const DamageResult kill = hit(10.0f, 1.0f, 0.0f);
        check(kill.killed, "a hit past the last health kills");
        t0 = now; step = 7;
        return false;
    }
    case 7: // the survival step turns it into a death (the death above was on normal: nothing carried, nothing lost)
        if (!sv.dead && now - t0 < 1.0) return false;
        check(sv.dead, "the player dies");
        step = 8;
        return false;
    case 8: // respawned: on normal, die carrying an axe, 100 stone and 5 dirt
        if (sv.dead || !p.grounded) return false;
        s.inventory.clear();
        s.toolWear.used.clear();
        s.inventory.add("stone_axe", 1); s.inventory.add("stone", 100); s.inventory.add("dirt", 5);
        sv.health = 1.0f;
        hit(10.0f, 1.0f, 0.0f);
        t0 = now; step = 9;
        return false;
    case 9: {
        if (!sv.dead && now - t0 < 1.0) return false;
        const int wear = s.toolWear.used.count("stone_axe") ? s.toolWear.used["stone_axe"] : 0;
        printf("bench: combat: normal death: stone %d, dirt %d, axe wear %d of %d\n", s.inventory.count("stone"), s.inventory.count("dirt"), wear, ToolDurability(ToolOf("stone_axe")));
        check(sv.dead && s.inventory.count("stone") == 90 && s.inventory.count("dirt") == 5 && s.inventory.count("stone_axe") == 1 && wear == ToolDurability(ToolOf("stone_axe")) / 10,
              "normal: a tenth of each resource stack (rounded down) and a tenth of tool life");
        step = 10;
        return false;
    }
    case 10: // respawned: walk 20 m off (the grave must lie out of reach of the respawn)
        if (sv.dead || !p.grounded || !p.nearest) return false;
        {
            const DVec3 l = p.nearest->toLocalD(p.position);
            p.spawnOn(world, *p.nearest, p.nearest->dirToWorld(standOff(l, 20.0f)));
        }
        t0 = now; step = 100;
        return false;
    case 100: // on hard, die carrying the axe, the stone and a shield in the off hand
        if ((!p.grounded && now - t0 < 5.0) || now - t0 < 0.5) return false;
        worldDifficulty = 2;
        off.id = "wooden_shield"; off.centi = 100;
        sv.health = 1.0f;
        hit(10.0f, 1.0f, 0.0f);
        t0 = now; step = 11;
        return false;
    case 11: {
        if (!sv.dead && now - t0 < 1.0) return false;
        int held = 0;
        for (const InvSlot& c : s.inventory.grid()) if (!c.empty()) ++held;
        check(sv.dead && held == 0 && graves.size() == 1 && graves[0].items.size() == 4, "hard: everything goes into one grave");
        // The grave survives a save and a load.
        FILE* f = fopen("bench_grave.tmp", "w+"); // (tmpfile() fails on Windows without rights to the drive root)
        check(f != nullptr, "a scratch file for the grave save");
        if (f) {
            Combat_Save(f);
            const Grave before = graves[0];
            graves.clear();
            rewind(f);
            char tok[64];
            while (fscanf(f, " %63s", tok) == 1) if (!Combat_LoadLine(tok, f)) break;
            fclose(f);
            remove("bench_grave.tmp");
            const bool same = graves.size() == 1 && graves[0].name == before.name && graves[0].items.size() == before.items.size() &&
                              graves[0].offHand == before.offHand && (graves[0].local - before.local).length() < 1e-3;
            check(same, "the grave saves and loads");
            // (#804) the loot rule: the owner (and team) at once; others after 15 minutes, never inside the owner's team claim
            const int64_t t = 1795000000;
            check(Combat_CanLoot("ann", t, "ann", t, false), "the owner loots at once");
            check(!Combat_CanLoot("ann", t, "bob", t + 60, false), "another player can't within 15 minutes");
            check(Combat_CanLoot("ann", t, "bob", t + kGraveOpenToAllSeconds, false), "and can after 15 minutes");
            check(!Combat_CanLoot("ann", t, "bob", t + 3600, true), "but never inside the owner's team claim");
            check(Combat_CanLoot("", t, "bob", t, false), "single player's graves: anyone");
            gSameTeam = [](const std::string& a, const std::string& b) { return (a == "ann" && b == "cat") || (a == "cat" && b == "ann"); };
            check(Combat_CanLoot("ann", t, "cat", t, false), "a teammate loots at once");
            gSameTeam = nullptr;
            int styles[3] = { 0, 0, 0 };
            for (uint32_t i = 1; i <= 30; ++i) ++styles[Combat_GraveStyle("terra", DVec3{ i * 3.0, 60000.0, -i * 2.0 }, i) % 3];
            check(styles[0] > 0 && styles[1] > 0 && styles[2] > 0 && Combat_GraveStyle("terra", DVec3{ 3.0, 60000.0, -2.0 }, 1) == Combat_GraveStyle("terra", DVec3{ 3.0, 60000.0, -2.0 }, 1),
                  "graves pick all three stones, the same one for the same grave");
        }
        step = 12;
        return false;
    }
    case 12: // respawned: walk up to the grave (stand on it) and take it all back
        if (sv.dead || !p.grounded || graves.empty()) {
            if (graves.empty() && !sv.dead) step = 13; // the respawn was next to it already
            return false;
        }
        {
            // First stand 4 m off and look at the cairn (bench_grave.png), then step onto it.
            Planetoid* b = p.nearest;
            const DVec3 l = graves[0].local;
            if (b && b->name == graves[0].body) {
                p.spawnOn(world, *b, b->dirToWorld(standOff(l, 4.0f)));
                const Vector3 up = Vector3Normalize(Vector3{ (float)l.x, (float)l.y, (float)l.z });
                p.forward = b->dirToWorld(Vector3Normalize(Vector3Subtract(Vector3Scale(sideOf(up), -1.0f), Vector3Scale(up, 0.35f)))); // toward it, a little down
            }
        }
        t0 = now; step = 120;
        return false;
    case 120:
        if (now - t0 < 2.0) return false;
        TakeScreenshot("bench_grave.png");
        {
            Planetoid* b = p.nearest;
            const DVec3 l = graves.empty() ? DVec3{} : graves[0].local;
            if (b && !graves.empty() && b->name == graves[0].body) p.spawnOn(world, *b, b->dirToWorld(Vector3Normalize(Vector3{ (float)l.x, (float)l.y, (float)l.z })));
        }
        t0 = now; step = 13;
        return false;
    case 13:
        if (!graves.empty() && now - t0 < 5.0) return false;
        printf("bench: combat: after the grave: axe %d, stone %d, dirt %d, off hand %s\n", s.inventory.count("stone_axe"), s.inventory.count("stone"),
               s.inventory.count("dirt"), off.empty() ? "empty" : off.id.c_str());
        check(graves.empty() && s.inventory.count("stone_axe") == 1 && s.inventory.count("stone") == 90 && s.inventory.count("dirt") == 5 && off.id == "wooden_shield",
              "walking up to the grave takes everything back, the shield to the off hand");
        worldDifficulty = savedDifficulty;
        return finish();
    default:
        return true;
    }
}


#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#828): a grave's map label ----
#include "unit_helpers.h"

// A multiplayer grave reads its time until it opens to others (15 min after the death), then "open to all"; a single
// player's grave (no owner: anyone may take it) its name alone; a waypoint that is no grave is not one.
TEST(Combat, AGraveOnTheMapSaysWhenItOpensToOthers) {
    Combat_Reset();
    const int64_t now = (int64_t)time(nullptr);
    std::vector<GraveInfo> list(3);
    list[0].id = 1; list[0].body = "terra"; list[0].local = { 0, 60000, 0 }; list[0].owner = "bob-uuid"; list[0].name = "bob's grave 1"; list[0].diedAt = now - 60;
    list[1].id = 2; list[1].body = "terra"; list[1].local = { 0, 60000, 10 }; list[1].owner = "bob-uuid"; list[1].name = "bob's grave 2"; list[1].diedAt = now - 20 * 60;
    list[2].id = 3; list[2].body = "terra"; list[2].local = { 0, 60000, 20 }; list[2].name = "Grave 3"; list[2].diedAt = now;
    Combat_SetGraveMirror(list, "someone-else");
    std::string label;
    ASSERT_TRUE(Combat_GraveMark("terra", "bob's grave 1", label));
    EXPECT_EQ(label, "bob's grave 1, open to all in 14 min");
    ASSERT_TRUE(Combat_GraveMark("terra", "bob's grave 2", label));
    EXPECT_EQ(label, "bob's grave 2, open to all");
    ASSERT_TRUE(Combat_GraveMark("terra", "Grave 3", label));
    EXPECT_EQ(label, "Grave 3");
    EXPECT_FALSE(Combat_GraveMark("terra", "Copper hill", label)) << "a plain waypoint";
    EXPECT_FALSE(Combat_GraveMark("luna", "Grave 3", label)) << "on another body";
    Combat_Reset();
}
#endif
