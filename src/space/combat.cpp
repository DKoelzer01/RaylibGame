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

#include "raymath.h"
#include "rlgl.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

// Tunables (#165 first slice; the numbers agreed on the issue).
static constexpr float kStaminaRegenPerSec = 25.0f;   // after kRegenDelay of not spending
static constexpr double kRegenDelay = 1.0;            // s
static constexpr float kBlockCostPerDamage = 1.5f;    // stamina per health point a block stops
static constexpr float kBlockArcCos = 0.5f;           // a block covers +-60 deg around the view (level)
static constexpr float kShieldBlockFraction = 0.75f;  // share of the damage a shield stops
static constexpr float kBlockedKnockback = 0.4f;      // share of the push left through a block
static constexpr float kKnockLift = 2.0f;             // m/s up with a push on the ground: off it, a short hop

static bool testHoldBlock = false; // --test-combat: RMB held
static void PickUpGraves(PlayerSession& s);

static bool IsShield(const std::string& id) { return id == "wooden_shield"; }

// Does RMB already do something with what is in the main hand? Then it keeps doing that, and a shield in the off
// hand does not block: food and feed (eat, #201), fillers (place), buckets, shovels and hoes (till), torches and
// snares (set down).
static bool RmbUsedByHand(const PlayerSession& s) {
    const std::string& held = s.inventory.grid()[Inventory::kGrid + s.hotbar.selected].id;
    if (held.empty()) return false;
    if (s.terrainTool.holdingItem || TerrainTool_IsFiller(held) || held == "bucket" || held == "water_bucket") return true;
    if (Hotbar_IsLight(held) || Hotbar_IsPlaceable(held)) return true;
    const ToolType t = ToolOf(held).type;
    return t == ToolType::Shovel || t == ToolType::Hoe;
}

void Combat_Step(PlayerSession& s, float dt, const PlayerInput& in) {
    SurvivalState& sv = s.survival;
    // The block: RMB with a shield in the off hand (and RMB not taken by the main hand or build mode), while stamina lasts.
    s.terrainTool.rmbBlocks = IsShield(s.inventory.offHand().id) && !s.buildTool.active && !RmbUsedByHand(s);
    sv.blocking = s.terrainTool.rmbBlocks && (in.down(Act::Secondary) || testHoldBlock) && sv.stamina > 0.0f && !sv.dead;
    // Stamina comes back a second after the last spend.
    if (SimTime() - sv.staminaUsedAt >= kRegenDelay) sv.stamina = std::min(100.0f, sv.stamina + kStaminaRegenPerSec * dt);
    PickUpGraves(s);
}

DamageResult Combat_Damage(PlayerSession& s, const DamageEvent& e) {
    DamageResult r;
    SurvivalState& sv = s.survival;
    Player& p = s.player;
    if (sv.dead || e.amount <= 0.0f) return r;
    float amount = e.amount, push = e.knockback;
    // Away from the source, level with the ground where the player stands.
    Vector3 away = { 0, 0, 0 };
    if (e.hasSource) {
        const Vector3 d = (p.position - e.source).toFloat(); // a short difference: fine in float
        away = Vector3Subtract(d, Vector3Scale(p.up, Vector3DotProduct(d, p.up)));
        const float l = Vector3Length(away);
        away = l > 1e-4f ? Vector3Scale(away, 1.0f / l) : Vector3{ 0, 0, 0 };
    }
    // A block stops part of a hit it faces (the source within the arc around the level view) while stamina lasts.
    if (sv.blocking && e.hasSource && Vector3LengthSqr(away) > 0.0f) {
        const Vector3 fwd = Vector3Normalize(Vector3Subtract(p.forward, Vector3Scale(p.up, Vector3DotProduct(p.forward, p.up))));
        if (Vector3DotProduct(fwd, Vector3Negate(away)) >= kBlockArcCos) {
            const float stopped = amount * kShieldBlockFraction;
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
    sv.health -= amount;
    sv.hitDamageTaken += amount;
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
};
static std::vector<Grave> graves;
static int graveCounter = 0;
static constexpr double kGraveReach = 2.5; // m: walking up to it takes everything back

void Combat_OnDeath(PlayerSession& s) {
    Inventory& inv = s.inventory;
    auto& cells = inv.grid();
    Player& p = s.player;
    if (worldDifficulty >= 2 && p.nearest) {
        Grave g;
        g.body = p.nearest->name;
        g.local = p.nearest->toLocalD(p.position);
        for (int i = 0; i < Inventory::kCells; ++i) {
            if (cells[i].empty()) continue;
            if (i == Inventory::kOffHand) g.offHand = (int)g.items.size();
            g.items.push_back(cells[i]);
            cells[i] = InvSlot{};
        }
        if (g.items.empty()) return;
        g.name = TextFormat("Grave %d", ++graveCounter);
        const DVec3 l = g.local;
        Nav_AddWaypoint(g.body, Vector3{ (float)l.x, (float)l.y, (float)l.z }, g.name, false); // normalised there
        GameUI_Toast(TextFormat("Everything you carried lies in your grave (%s, marked on the compass)", g.name.c_str()));
        graves.push_back(std::move(g));
        return;
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
            if (max > 0) { int& u = s.toolWear.used[c.id]; u = std::min(u + max / 10, max - 1); }
        } else if (k == ItemKind::Resource && worldDifficulty >= 1) {
            const int64_t units = (int64_t)(c.amount() * 0.1);
            if (units > 0) { c.centi -= units * 100; lost += (int)units; }
        }
    }
    GameUI_Toast(lost > 0 ? TextFormat("You lost %d items and your tools took wear", lost) : "Your tools took wear");
}

static void PickUpGraves(PlayerSession& s) {
    Player& p = s.player;
    if (!p.nearest || s.survival.dead) return;
    const DVec3 me = p.nearest->toLocalD(p.position);
    for (size_t gi = 0; gi < graves.size();) {
        Grave& g = graves[gi];
        if (g.body != p.nearest->name || (g.local - me).length() > kGraveReach) { ++gi; continue; }
        auto& cells = s.inventory.grid();
        std::vector<InvSlot> left;
        for (size_t i = 0; i < g.items.size(); ++i) {
            const InvSlot& it = g.items[i];
            if ((int)i == g.offHand && cells[Inventory::kOffHand].empty()) { cells[Inventory::kOffHand] = it; continue; }
            ++s.inventory.muted; // no pickup feed line per stack
            const double rest = s.inventory.addAmount(it.id, it.amount());
            --s.inventory.muted;
            s.inventory.overflowed = false;
            if (rest > 0.0) left.push_back(InvSlot{ it.id, (int64_t)llround(rest * 100.0) });
        }
        if (left.empty()) {
            Nav_RemoveWaypoint(g.body, g.name);
            GameUI_Toast(TextFormat("You took everything back from %s", g.name.c_str()));
            graves.erase(graves.begin() + gi);
            continue;
        }
        g.items = std::move(left); g.offHand = -1;
        GameUI_Toast(TextFormat("No room for the rest: %d stacks stay in %s", (int)g.items.size(), g.name.c_str()));
        ++gi;
    }
}

int Combat_GraveCount() { return (int)graves.size(); }

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
    const double t = SimTime() - s.survival.lastHitAt;
    return !reduceMotion && t >= 0.0 && t < kHitStopTime && s.survival.lastHitStrength > 0.0f;
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
        DrawSphere(Vector3Scale(up, 0.15f), 0.4f, Color{ 110, 105, 100, 255 });
        DrawSphere(Vector3Scale(up, 0.55f), 0.28f, Color{ 125, 120, 115, 255 });
        DrawSphere(Vector3Scale(up, 0.85f), 0.18f, Color{ 140, 135, 128, 255 });
        DrawSphere(Vector3Scale(up, 1.1f), 0.09f, Color{ 220, 50, 40, 255 });
        rlDrawRenderBatchActive();
        rlPopMatrix();
    }
}

void Combat_Reset() { graves.clear(); graveCounter = 0; }

void Combat_Save(FILE* f) {
    fprintf(f, "grave_counter %d\n", graveCounter);
    for (const Grave& g : graves) {
        fprintf(f, "grave %s %.6f %.6f %.6f %d %zu", g.body.c_str(), g.local.x, g.local.y, g.local.z, g.offHand, g.items.size());
        for (const InvSlot& it : g.items) fprintf(f, " %s %lld", it.id.c_str(), (long long)it.centi);
        fprintf(f, " %s\n", g.name.c_str()); // last: the name has a space in it
    }
}

bool Combat_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "grave_counter") == 0) { if (fscanf(f, " %d", &graveCounter) != 1) graveCounter = 0; return true; }
    if (strcmp(tok, "grave") != 0) return false;
    Grave g;
    char body[64] = {};
    size_t n = 0;
    if (fscanf(f, " %63s %lf %lf %lf %d %zu", body, &g.local.x, &g.local.y, &g.local.z, &g.offHand, &n) != 6) return true;
    g.body = body;
    for (size_t i = 0; i < n; ++i) {
        char id[64] = {};
        long long c = 0;
        if (fscanf(f, " %63s %lld", id, &c) != 2) return true;
        if (itemDB.find(id) && c > 0) g.items.push_back(InvSlot{ id, (int64_t)c });
        else if ((int)i == g.offHand) g.offHand = -1;
    }
    char name[64] = {};
    if (fscanf(f, " %63[^\n]", name) != 1) return true;
    g.name = name;
    if (g.offHand >= (int)g.items.size()) g.offHand = -1;
    if (!g.items.empty()) graves.push_back(std::move(g));
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
        if (now - t0 < 0.1 || (!p.grounded && now - t0 < 3.0)) return false;
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
