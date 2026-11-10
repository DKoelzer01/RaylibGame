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
    survival.breathable = false;
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
        // A planet's air is breathable up to ~15 km: the suit tops up its oxygen from it instead.
        const auto& atm = world.atmosphere;
        const Vector3 head = Vector3Add(pr, Vector3Scale(player.up, player.eyeHeight));
        survival.breathable = atm.enabled && Vector3Distance(pr, atm.centre) < atm.radius + survival.breathableAltitude && !IsUnderwater(world, head);
        if (survival.breathable) survival.oxygen = std::min(100.0f, survival.oxygen + survival.airRefillPerSec * dt);
        else survival.oxygen = std::max(0.0f, survival.oxygen - survival.oxygenDrainPerSec * dt);
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
    else if (survival.breathable) DrawText("BREATHABLE AIR", x + 280, y + 2, 16, Fade(SKYBLUE, 0.7f));
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
        { "Dig some stone", false },
        { "Find iron ore underground", false },
        { "Craft a power cell (Tab)", false },
        { "Board your ship (E)", false },
        { "Fly to the moon", false },
        { "Visit the asteroid ring", false },
    };
    std::set<int> visitedBodies;
}

void Objectives_MarkVisitedBody(int i) { visitedBodies.insert(i); }

void Objectives_Update(Scene& world) {
    QuestBook_Update();
    if (inventory.count("stone") > 0 || inventory.count("regolith") > 0) objectives[0].done = true;
    if (inventory.count("iron_ore") > 0) objectives[1].done = true;
    if (inventory.count("power_cell") > 0) objectives[2].done = true;
    if (ship.piloted) objectives[3].done = true;
    // Visiting: get within 2 radii of any body other than the start one.
    const Vector3 pr = ship.piloted ? ship.renderPos(world) : (player.position - world.renderOrigin).toFloat();
    int index = 0;
    for (auto& obj : world.objects) {
        if (auto* p = dynamic_cast<Planetoid*>(obj.get())) {
            if (index > 0 && Vector3Distance(pr, p->position) < p->radius * p->scale * 2.0f) visitedBodies.insert(index);
            if (index > 0 && visitedBodies.count(index)) {
                if (p->kind == BodyKind::Moon) objectives[4].done = true;
                else if (p->kind == BodyKind::Asteroid) objectives[5].done = true;
            }
            index++;
        }
    }
}


// ---------------------------------------------------------------------------------------------
// Quest book (resources/data/tiers.txt): the tier ladder from bare hands to orbit. Making (or
// placing) a tier's gate item opens the next one.
// ---------------------------------------------------------------------------------------------

namespace {
    struct Tier { std::string name, gate, text; };
    std::vector<Tier> tiers;
    int questTier = 0;
    std::set<std::string> placedItems;

    void LoadTiers() {
        if (!tiers.empty()) return;
        FILE* f = fopen("resources/data/tiers.txt", "r");
        if (!f) return;
        char buf[1024];
        while (fgets(buf, sizeof buf, f)) {
            std::string line = buf;
            if (line.empty() || line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
            std::vector<std::string> cols;
            size_t start = 0, bar;
            while ((bar = line.find('|', start)) != std::string::npos) { cols.push_back(line.substr(start, bar - start)); start = bar + 1; }
            cols.push_back(line.substr(start));
            for (auto& c : cols) { const size_t a = c.find_first_not_of(" \t"), b = c.find_last_not_of(" \t\r\n"); c = a == std::string::npos ? "" : c.substr(a, b - a + 1); }
            if (cols.size() >= 4) tiers.push_back({ cols[1], cols[2], cols[3] });
        }
        fclose(f);
    }
}

void QuestBook_NotePlaced(const std::string& item) { placedItems.insert(item); }
int QuestBook_Tier() { return questTier; }
void QuestBook_SetTier(int t) { questTier = std::max(0, t); }

void QuestBook_Update() {
    LoadTiers();
    while (questTier < (int)tiers.size() && tiers[questTier].gate != "-" &&
           (inventory.count(tiers[questTier].gate) > 0 || placedItems.count(tiers[questTier].gate))) {
        ++questTier;
        if (questTier < (int)tiers.size()) GameUI_Toast(TextFormat("Tier %d: %s", questTier, tiers[questTier].name.c_str()));
    }
}

void Objectives_DrawHUD() {
    LoadTiers();
    int y = 64;
    if (questTier < (int)tiers.size()) {
        const Tier& t = tiers[questTier];
        DrawText(TextFormat("QUEST BOOK   tier %d: %s", questTier, t.name.c_str()), 12, y, 16, Fade(GOLD, 0.9f));
        y += 20;
        if (t.gate != "-") {
            const ItemDef* d = itemDB.find(t.gate);
            DrawText(TextFormat("Next: make a %s", d ? d->name.c_str() : t.gate.c_str()), 12, y, 16, RAYWHITE);
            y += 20;
        }
        // Wrap the hint at ~64 characters.
        std::string words = t.text, line;
        size_t pos = 0;
        while (pos < words.size()) {
            size_t sp = words.find(' ', pos);
            if (sp == std::string::npos) sp = words.size();
            const std::string w = words.substr(pos, sp - pos);
            if (line.size() + w.size() + 1 > 64) { DrawText(line.c_str(), 12, y, 14, Fade(RAYWHITE, 0.75f)); y += 17; line.clear(); }
            line += (line.empty() ? "" : " ") + w;
            pos = sp + 1;
        }
        if (!line.empty()) { DrawText(line.c_str(), 12, y, 14, Fade(RAYWHITE, 0.75f)); y += 17; }
        y += 8;
    }
    // Space goals, compactly.
    for (const Objective& o : objectives) {
        if (o.done) continue;
        DrawText(TextFormat("[ ] %s", o.text), 12, y, 14, Fade(RAYWHITE, 0.6f));
        y += 16;
    }
}
