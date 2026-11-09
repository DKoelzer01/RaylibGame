#include "survival.h"
#include "inventory.h"
#include "player.h"
#include "ship.h"
#include "audio.h"
#include "gameui.h"
#include "world.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <vector>

SurvivalState survival;

void Survival_OnLanding(float impactSpeed) {
    const float safe = 9.0f;
    if (impactSpeed <= safe || survival.dead) return;
    survival.health -= (impactSpeed - safe) * 8.0f;
    Audio_Play(Sfx::Error, 0.6f, 0.6f);
}

void Survival_Respawn(Scene& world) {
    survival.oxygen = survival.power = survival.health = 100.0f;
    survival.dead = false;
    player.flyMode = false; // otherwise the free-fly camera would put the player back where they died
    if (ship.spawned) {
        // Wake up next to the ship (or aboard it, if it is out in space).
        const Vector3 sp = ship.renderPos(world);
        const Vector3 beside = Vector3Add(sp, Vector3Scale(ship.right(), 4.0f));
        float dist; Vector3 n;
        Planetoid* body = NearestPlanetoid(world, sp);
        if (body && TerrainSurfaceDistance(*body, sp, dist, n) && dist < 6.0f) {
            player.position = world.renderOrigin + Vector3Add(beside, Vector3Scale(ship.up(), 1.0f));
            player.velocity = { 0, 0, 0 };
            ship.piloted = false;
        } else {
            ship.piloted = true;
        }
    }
    GameUI_Toast("Respawned at your ship");
}

void Survival_Update(Scene& world, float dt, bool inputEnabled) {
    if (survival.dead) {
        if (GetTime() - survival.deathTime > 3.0) Survival_Respawn(world);
        return;
    }
    const Vector3 pr = ship.piloted ? ship.renderPos(world) : (player.position - world.renderOrigin).toFloat();
    // Life support runs on ship energy: a drained ship is no refuge (charge it with a power cell).
    survival.inLifeSupport = (ship.piloted || ship.playerInRange(world, pr, ship.lifeSupportRadius)) && ship.energy > 0.0f;

    if (survival.inLifeSupport) {
        // Topping up the suit costs the ship energy; refill only as much as the ship can pay for.
        const float energyPerPoint = 0.02f;
        const float wantO2 = std::min(100.0f - survival.oxygen, survival.lifeSupportRefillPerSec * dt);
        const float wantPwr = std::min(100.0f - survival.power, survival.lifeSupportRefillPerSec * dt);
        const float want = wantO2 + wantPwr;
        const float k = want > 0.0f ? std::min(1.0f, ship.energy / (energyPerPoint * want)) : 0.0f;
        survival.oxygen += wantO2 * k;
        survival.power += wantPwr * k;
        ship.energy = std::max(0.0f, ship.energy - energyPerPoint * want * k);
        survival.health = std::min(100.0f, survival.health + 2.0f * dt);
        player.jetpackFuel = std::min(1.0f, player.jetpackFuel + 0.5f * dt); // refill the suit's thruster tank
    } else {
        survival.oxygen = std::max(0.0f, survival.oxygen - survival.oxygenDrainPerSec * dt);
        survival.power = std::max(0.0f, survival.power - (survival.powerDrainPerSec
                                                          + (survival.headlamp && !ship.piloted ? survival.headlampPowerPerSec : 0.0f)) * dt);
        // Firing the jetpack only burns fuel; recharging the tank (on the ground) is what costs suit
        // power, and only as far as the power left can pay for.
        const float refilled = player.jetpackFuel - survival.lastJetpackFuel;
        if (refilled > 0.0f) {
            const float affordable = survival.power / survival.jetpackPowerPerTank;
            const float paid = std::min(refilled, affordable);
            player.jetpackFuel = survival.lastJetpackFuel + paid;
            survival.power = std::max(0.0f, survival.power - paid * survival.jetpackPowerPerTank);
        }
    }
    survival.lastJetpackFuel = player.jetpackFuel;
    if (survival.power <= 0.0f && survival.headlamp && !ship.piloted) { survival.headlamp = false; GameUI_Toast("Headlamp out: no suit power"); }
    if (inputEnabled && IsKeyPressed(KEY_L)) {
        if (survival.power > 0.0f || ship.piloted) {
            survival.headlamp = !survival.headlamp;
            Audio_Play(Sfx::Build, 0.25f, survival.headlamp ? 1.8f : 1.5f);
        } else GameUI_Toast("No suit power for the headlamp");
    }
    if (survival.oxygen <= 0.0f) survival.health -= survival.suffocationDamagePerSec * dt;

    if (inputEnabled) {
        if (IsKeyPressed(KEY_ONE)) {
            if (survival.oxygen < 99.0f && inventory.remove("oxygen_canister", 1)) { survival.oxygen = std::min(100.0f, survival.oxygen + 50.0f); Audio_Play(Sfx::Pickup); }
            else Audio_Play(Sfx::Error, 0.4f);
        }
        if (IsKeyPressed(KEY_TWO)) {
            // Don't waste a cell on a full store (ship while piloting, otherwise the suit).
            const bool full = ship.piloted ? ship.energy >= 99.0f : survival.power >= 99.0f;
            if (!full && inventory.remove("power_cell", 1)) {
                if (ship.piloted) ship.energy = std::min(100.0f, ship.energy + 40.0f);
                else survival.power = std::min(100.0f, survival.power + 50.0f);
                Audio_Play(Sfx::Pickup);
            } else Audio_Play(Sfx::Error, 0.4f);
        }
    }

    if (survival.health <= 0.0f) {
        survival.health = 0.0f;
        survival.dead = true;
        survival.deathTime = GetTime();
        survival.deaths++;
        player.autoWalk = false;
    }
}

static void Bar(int x, int y, int w, const char* label, float value, Color c, bool warn) {
    DrawRectangle(x, y, w, 16, Fade(BLACK, 0.55f));
    const bool flash = warn && fmodf((float)GetTime(), 0.8f) < 0.4f;
    DrawRectangle(x, y, (int)(w * value / 100.0f), 16, flash ? RED : c);
    DrawRectangleLines(x, y, w, 16, Fade(RAYWHITE, 0.4f));
    DrawText(TextFormat("%s %3.0f", label, value), x + 6, y + 1, 14, RAYWHITE);
}

void Survival_DrawHUD() {
    const int x = 20, y = GetScreenHeight() - 170;
    Bar(x, y, 260, "O2", survival.oxygen, { 90, 200, 255, 255 }, survival.oxygen < 20.0f);
    Bar(x, y + 22, 260, "PWR", survival.power, { 250, 210, 70, 255 }, survival.power < 15.0f);
    Bar(x, y + 44, 260, "HP", survival.health, { 230, 70, 70, 255 }, survival.health < 30.0f);
    if (ship.spawned) Bar(x + 280, y + 44, 180, "SHIP", ship.energy, { 160, 120, 255, 255 }, ship.energy < 15.0f);
    if (survival.inLifeSupport) DrawText("LIFE SUPPORT", x + 280, y + 2, 16, SKYBLUE);
    DrawText(TextFormat("[1] O2 canister x%d   [2] power cell x%d", inventory.count("oxygen_canister"), inventory.count("power_cell")), x, y - 22, 16, Fade(RAYWHITE, 0.8f));
    if (survival.dead) {
        DrawRectangle(0, 0, GetScreenWidth(), GetScreenHeight(), Fade(BLACK, 0.6f));
        const char* msg = survival.oxygen <= 0.0f ? "You suffocated" : "You died";
        DrawText(msg, GetScreenWidth() / 2 - MeasureText(msg, 40) / 2, GetScreenHeight() / 2 - 30, 40, RED);
        DrawText("Respawning at your ship...", GetScreenWidth() / 2 - MeasureText("Respawning at your ship...", 20) / 2, GetScreenHeight() / 2 + 20, 20, RAYWHITE);
    }
}

// ---------------------------------------------------------------------------------------------
// Objectives
// ---------------------------------------------------------------------------------------------

namespace {
    struct Objective { const char* text; bool done; };
    std::vector<Objective> objectives = {
        { "Mine some ice", false },
        { "Craft an oxygen canister (Tab)", false },
        { "Board your ship (E)", false },
        { "Fly to another asteroid", false },
        { "Craft a power cell from iron ore", false },
    };
    std::set<int> visitedBodies;
}

void Objectives_MarkVisitedBody(int i) { visitedBodies.insert(i); }

void Objectives_Update(Scene& world) {
    if (inventory.count("ice") > 0 || inventory.count("water") > 0) objectives[0].done = true;
    if (inventory.count("oxygen_canister") > 0) objectives[1].done = true;
    if (ship.piloted) objectives[2].done = true;
    if (inventory.count("power_cell") > 0) objectives[4].done = true;
    // Visiting: get within 2 radii of any body other than the start one.
    const Vector3 pr = ship.piloted ? ship.renderPos(world) : (player.position - world.renderOrigin).toFloat();
    int index = 0;
    for (auto& obj : world.objects) {
        if (auto* p = dynamic_cast<Planetoid*>(obj.get())) {
            if (index > 0 && Vector3Distance(pr, p->position) < p->radius * p->scale * 2.0f) visitedBodies.insert(index);
            index++;
        }
    }
    if (!visitedBodies.empty()) objectives[3].done = true;
}

void Objectives_DrawHUD() {
    int y = 64;
    DrawText("OBJECTIVES", 12, y, 16, Fade(RAYWHITE, 0.8f));
    for (const Objective& o : objectives) {
        y += 20;
        DrawText(TextFormat("%s %s", o.done ? "[x]" : "[ ]", o.text), 12, y, 16, o.done ? Fade(GREEN, 0.7f) : RAYWHITE);
    }
}
