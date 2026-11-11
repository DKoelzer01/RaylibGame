#include "survival.h"
#include "suit.h"
#include "questtree.h"
#include "rocket.h"
#include "inventory.h"
#include "player.h"
#include "ship.h"
#include "audio.h"
#include "gameui.h"
#include "world.h"
#include "worldtypes.h"
#include "voxel.h"
#include "scene.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

SurvivalState survival;
int worldDifficulty = 1;
float Difficulty_HungerScale() { return worldDifficulty == 0 ? 0.6f : worldDifficulty == 2 ? 1.25f : 1.0f; }
float Difficulty_ThreatScale() { return worldDifficulty == 0 ? 0.0f : worldDifficulty == 2 ? 1.5f : 1.0f; }
const char* Difficulty_Name(int d) { return d == 0 ? "Peaceful" : d == 2 ? "Hard" : "Normal"; }

// ---------------------------------------------------------------------------------------------
// Food (resources/data/foods.txt: item | food | note) and the air temperature
// ---------------------------------------------------------------------------------------------

namespace {
    struct Food { float value = 0.0f; bool raw = false; };
    std::unordered_map<std::string, Food> foods;
    bool foodsLoaded = false;

    void LoadFoods() {
        if (foodsLoaded) return;
        foodsLoaded = true;
        FILE* f = fopen("resources/data/foods.txt", "r");
        if (!f) return;
        char line[512];
        while (fgets(line, sizeof line, f)) {
            if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
            char id[64] = {}, flag[32] = {};
            float v = 0.0f;
            // "id | value | raw" (the third column optional)
            if (sscanf(line, " %63[^ |] | %f | %31[^ |\r\n]", id, &v, flag) >= 2) foods[id] = { v, strcmp(flag, "raw") == 0 };
        }
        fclose(f);
    }

    // The air temperature where the player stands (deg C), or NaN in vacuum (the suit regulates):
    // 32 at the equator falling to about -25 at the poles, 6.5 colder per km of height, +-7 between
    // noon and midnight, and underground drifting toward a cave's steady 12.
    float AirTemperature(Scene& world) {
        Planetoid* p = player.nearest;
        if (!p || !world.atmosphere.enabled) return NAN;
        const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
        if (!pf || pf->moon) return NAN;
        const Vector3 local = p->toLocal((player.position - world.renderOrigin).toFloat());
        const float d = Vector3Length(local);
        if (d < 1.0f) return NAN;
        const Vector3 up = Vector3Scale(local, 1.0f / d);
        const float height = d - pf->radius;                  // above sea level
        if (height > survival.breathableAltitude) return NAN; // above the air: the suit regulates
        const float lat = fabsf(up.y);
        const float sun = Vector3DotProduct(Vector3Normalize(world.sunDirection), player.up); // -1 midnight .. 1 noon
        float t = 32.0f - 57.0f * powf(lat, 1.5f) - 6.5f * std::max(height, 0.0f) / 1000.0f + 7.0f * sun;
        const float depth = pf->elevation(up) - height;       // below the surface here
        if (depth > 3.0f) t += (12.0f - t) * std::min(1.0f, (depth - 3.0f) / 12.0f);
        return t;
    }
}

bool Survival_IsFood(const std::string& id) { LoadFoods(); return foods.count(id) > 0; }
float Survival_FoodValue(const std::string& id) { LoadFoods(); auto it = foods.find(id); return it == foods.end() ? 0.0f : it->second.value; }

void Survival_OnLanding(float impactSpeed) {
    const float safe = Suit_Stats().servo ? 12.0f : 9.0f; // servo legs soften landings
    if (impactSpeed <= safe || survival.dead) return;
    survival.health -= (impactSpeed - safe) * 8.0f;
    Audio_Play(Sfx::Error, 0.6f, 0.6f);
}

void Survival_Respawn(Scene& world) {
    survival.oxygen = survival.power = survival.health = 100.0f;
    survival.food = std::max(survival.food, 60.0f);
    survival.sickUntil = 0.0f;
    survival.dead = false;
    player.flyMode = false; // otherwise the free-fly camera would put the player back where they died
    // Wake up at the start, at rest with the ground there (spawnOn sets the velocity to the surface's:
    // nothing of the fall or flight that killed you carries over).
    Planetoid* start = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get());
    if (start) {
        const auto* pf = dynamic_cast<const PlanetField*>(start->field.get());
        player.spawnOn(world, *start, start->dirToWorld(pf ? pf->findSpawnDir(start->seed) : Vector3{ 0, 1, 0 }));
    } else {
        player.velocity = { 0, 0, 0 };
    }
    player.impactSpeed = 0.0f;
    GameUI_Toast("Respawned at the start");
}

void Survival_Update(Scene& world, float dt, bool inputEnabled) {
    if (survival.dead) {
        if (GetTime() - survival.deathTime > 3.0) Survival_Respawn(world);
        return;
    }
    const Vector3 pr = (player.position - world.renderOrigin).toFloat();
    survival.breathable = false;
    survival.underwater = false;
    // Life support: a rocket's pressurised cockpit (the placeholder ship that also gave it is gone).
    survival.inLifeSupport = Rocket_Boarded();
    if (survival.inLifeSupport) {
        survival.oxygen = std::min(100.0f, survival.oxygen + survival.lifeSupportRefillPerSec * dt);
        survival.power = std::min(100.0f, survival.power + survival.lifeSupportRefillPerSec * dt);
        survival.health = std::min(100.0f, survival.health + 2.0f * dt);
        player.jetpackFuel = std::min(1.0f, player.jetpackFuel + 0.5f * dt); // refill the suit's thruster tank
    } else {
        // A planet's air is breathable up to ~15 km: the suit tops up its oxygen from it instead.
        const auto& atm = world.atmosphere;
        const Vector3 head = Vector3Add(pr, Vector3Scale(player.up, player.eyeHeight));
        survival.breathable = atm.enabled && Vector3Distance(pr, atm.centre) < atm.radius + survival.breathableAltitude && !IsUnderwater(world, head);
        if (survival.breathable) survival.oxygen = std::min(100.0f, survival.oxygen + survival.airRefillPerSec * dt);
        else {
            // Underwater drains it faster (a rebreather slows that); bigger tanks last longer.
            const bool underwater = atm.enabled && IsUnderwater(world, head);
            survival.underwater = underwater;
            survival.oxygen = std::max(0.0f, survival.oxygen - survival.oxygenDrainPerSec * (underwater ? Suit_Stats().underwater : 1.0f) / Suit_Stats().oxygen * dt);
        }
        survival.power = std::max(0.0f, survival.power - (1.0f / Suit_Stats().power) * (survival.powerDrainPerSec
                                                          + (survival.headlamp && !ship.piloted ? survival.headlampPowerPerSec : 0.0f)) * dt);
        // Firing the jetpack only burns fuel; recharging the tank (on the ground) is what costs suit
        // power, and only as far as the power left can pay for.
        const float refilled = player.jetpackFuel - survival.lastJetpackFuel;
        if (refilled > 0.0f) {
            const float perTank = survival.jetpackPowerPerTank / Suit_Stats().power;
            const float affordable = survival.power / perTank;
            const float paid = std::min(refilled, affordable);
            player.jetpackFuel = survival.lastJetpackFuel + paid;
            survival.power = std::max(0.0f, survival.power - paid * perTank);
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
    // Temperature and hunger. Food drains faster running or on the jetpack, in the cold, and
    // while sick; the cockpit's life support feeds you. Heat does slow damage.
    survival.temperature = AirTemperature(world);
    {
        const Vector3 flat = Vector3Subtract(player.velocity, Vector3Scale(player.up, Vector3DotProduct(player.velocity, player.up)));
        const bool active = player.jetpackFiring || (player.grounded && Vector3Length(flat) > player.walkSpeed * 1.1f);
        float drain = survival.foodDrainPerSec * (active ? 1.6f : 1.0f) * Suit_Stats().food * Difficulty_HungerScale();
        const float T = survival.temperature;
        if (!std::isnan(T) && T < 5.0f) drain *= 1.0f + Suit_Stats().cold * (5.0f - T) / 15.0f; // insulation, the regulator
        if (GetTime() < survival.sickUntil) drain *= 2.0f;
        if (!survival.inLifeSupport) survival.food = std::max(0.0f, survival.food - drain * dt);
        if (survival.food <= 0.0f) survival.health -= survival.starvationDamagePerSec * dt;
        else if (survival.food > 80.0f && survival.oxygen > 0.0f) survival.health = std::min(100.0f, survival.health + survival.wellFedHealPerSec * dt);
        if (!std::isnan(T) && T > 35.0f && !Suit_Stats().heatProof) survival.health -= 0.05f * (T - 35.0f) * dt; // 40 C: 0.25 a second
    }

    if (inputEnabled) {
    }

    {   // A night survived (the quest tree): dark fell and day came back without a death.
        static bool inNight = false;
        static int deathsAtDusk = 0;
        const bool dark = Vector3DotProduct(Vector3Normalize(world.sunDirection), player.up) < -0.1f && player.nearest != nullptr;
        if (dark && !inNight) { inNight = true; deathsAtDusk = survival.deaths; }
        else if (!dark && inNight) { inNight = false; if (survival.deaths == deathsAtDusk) QuestBook_Count("nights"); }
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
    const int x = 20, y = GetScreenHeight() - 192;
    // The oxygen and power bars grow with the suit's capacity (basic suit: 260 px).
    const int wO2 = (int)(180.0f + 80.0f * std::min(Suit_Stats().oxygen, 2.5f)), wPwr = (int)(180.0f + 80.0f * std::min(Suit_Stats().power, 2.5f));
    Bar(x, y, wO2, "O2", survival.oxygen, { 90, 200, 255, 255 }, survival.oxygen < 20.0f);
    Bar(x, y + 22, wPwr, "PWR", survival.power, { 250, 210, 70, 255 }, survival.power < 15.0f);
    Bar(x, y + 44, 260, "HP", survival.health, { 230, 70, 70, 255 }, survival.health < 30.0f);
    Bar(x, y + 66, 260, "FOOD", survival.food, { 120, 200, 90, 255 }, survival.food < 15.0f);
    // The air temperature, coloured when it costs you (cold: hunger; heat: health).
    const float T = survival.temperature;
    if (!std::isnan(T)) {
        const Color tc = T < 5.0f ? Color{ 130, 190, 255, 255 } : T > 35.0f ? Color{ 255, 140, 70, 255 } : Fade(RAYWHITE, 0.8f);
        const SuitStats& st = Suit_Stats();
        const char* note = T < 5.0f ? (st.cold <= 0.0f ? "  cold: regulated" : st.cold < 1.0f ? "  cold: insulated" : "  cold: hungrier") : T > 35.0f ? (st.heatProof ? "  heat: regulated" : "  heat") : "";
        DrawText(TextFormat("%.0f C%s", T, note), x + std::max(wO2, wPwr) + 20, y + 68, 16, tc);
    }
    const int tx = x + std::max(wO2, wPwr) + 20;
    if (GetTime() < survival.sickUntil) DrawText("SICK", tx, y + 46, 16, ORANGE);
    if (survival.inLifeSupport) DrawText("LIFE SUPPORT", tx, y + 2, 16, SKYBLUE);
    else if (survival.breathable) DrawText("BREATHABLE AIR", tx, y + 2, 16, Fade(SKYBLUE, 0.7f));
    DrawText(TextFormat("O2 canister x%d   power cell x%d   (use from the hotbar)", inventory.count("oxygen_canister"), inventory.count("power_cell")), x, y - 22, 16, Fade(RAYWHITE, 0.8f));
    if (survival.dead) {
        DrawRectangle(0, 0, GetScreenWidth(), GetScreenHeight(), Fade(BLACK, 0.6f));
        const char* msg = survival.oxygen <= 0.0f ? "You suffocated" : "You died";
        DrawText(msg, GetScreenWidth() / 2 - MeasureText(msg, 40) / 2, GetScreenHeight() / 2 - 30, 40, RED);
        DrawText("Respawning...", GetScreenWidth() / 2 - MeasureText("Respawning...", 20) / 2, GetScreenHeight() / 2 + 20, 20, RAYWHITE);
    }
}

// Current Objectives (resources/data/tiers.txt): the tier ladder from bare hands to orbit, each
// tier a short list of tasks that tick themselves off. A tier opens the next once every task is
// done (a tier with no tasks: once its gate item is made or placed).
//   get <item> <n>   n of the item obtained, ever (crafted, mined, taken from a machine)
//   place <item>     a machine of that item placed
//   orbit            a rocket in orbit
//   visit moon|ring  near the moon / a ring rock or one of the huge rocks
// ---------------------------------------------------------------------------------------------

namespace {
    struct Task { std::string kind, target; double n = 1.0; std::string text; bool done = false; };
    struct Tier { std::string name, gate, text; std::vector<Task> tasks; };
    std::vector<Tier> tiers;
    int questTier = 0;
    std::set<std::string> placedItems;
    std::set<int> visitedBodies;
    bool reachedOrbit = false, visitedMoon = false, visitedRing = false;
    std::map<std::string, double> questCounts; // events for "do" tasks: eat, feed, penned, till, plant, harvest, waypoint, modules, explored

    std::vector<std::string> Columns(const std::string& line) {
        std::vector<std::string> cols;
        size_t start = 0, bar;
        while ((bar = line.find('|', start)) != std::string::npos) { cols.push_back(line.substr(start, bar - start)); start = bar + 1; }
        cols.push_back(line.substr(start));
        for (auto& c : cols) { const size_t a = c.find_first_not_of(" \t"), b = c.find_last_not_of(" \t\r\n"); c = a == std::string::npos ? "" : c.substr(a, b - a + 1); }
        return cols;
    }

    void LoadTiers() {
        if (!tiers.empty()) return;
        FILE* f = fopen("resources/data/tiers.txt", "r");
        if (!f) return;
        char buf[1024];
        while (fgets(buf, sizeof buf, f)) {
            std::string line = buf;
            const size_t a = line.find_first_not_of(" \t");
            if (a == std::string::npos || line[a] == '#' || line[a] == '\n' || line[a] == '\r') continue;
            const std::vector<std::string> cols = Columns(line);
            if (cols[0] == "task") {
                if (tiers.empty() || cols.size() < 3) continue;
                Task t;
                char kind[32] = {}, target[64] = {};
                double n = 1.0;
                const int got = sscanf(cols[1].c_str(), "%31s %63s %lf", kind, target, &n);
                if (got < 1) continue;
                t.kind = kind;
                t.target = got >= 2 ? target : "";
                t.n = got >= 3 ? n : 1.0;
                t.text = cols[2];
                tiers.back().tasks.push_back(t);
            } else if (cols.size() >= 4) {
                tiers.push_back({ cols[1], cols[2], cols[3], {} });
            }
        }
        fclose(f);
    }

    bool TaskMet(const Task& t) {
        if (t.kind == "get") return inventory.obtainedTotal(t.target) >= t.n - 1e-6;
        if (t.kind == "place") return placedItems.count(t.target) > 0;
        if (t.kind == "orbit") return reachedOrbit;
        if (t.kind == "visit") return t.target == "moon" ? visitedMoon : visitedRing;
        if (t.kind == "do") { auto it = questCounts.find(t.target); return it != questCounts.end() && it->second >= t.n - 1e-6; }
        return false;
    }
}

void Objectives_MarkVisitedBody(int i) { visitedBodies.insert(i); }

void Objectives_Update(Scene& world) {
    if (Rocket_Telemetry(world).periapsis > 35000.0 && Rocket_Boarded()) reachedOrbit = true;
    // Visiting: within 2 radii of any body other than the start one.
    const Vector3 pr = (player.position - world.renderOrigin).toFloat();
    int index = 0;
    for (auto& obj : world.objects) {
        if (auto* p = dynamic_cast<Planetoid*>(obj.get())) {
            if (index > 0 && Vector3Distance(pr, p->position) < p->radius * p->scale * 2.0f) {
                visitedBodies.insert(index);
                if (p->kind == BodyKind::Moon) visitedMoon = true;
                else if (p->kind == BodyKind::Asteroid) visitedRing = true;
            }
            index++;
        }
    }
    QuestBook_Update();
}

void QuestBook_NotePlaced(const std::string& item) { placedItems.insert(item); }
bool QuestBook_Met(const std::string& kind, const std::string& target, double n) { Task t; t.kind = kind; t.target = target; t.n = n; return TaskMet(t); }
void QuestBook_Count(const std::string& event, double n) { questCounts[event] += n; }
void QuestBook_AtLeast(const std::string& event, double v) { double& c = questCounts[event]; c = std::max(c, v); }
double QuestBook_CountOf(const std::string& event) { auto it = questCounts.find(event); return it == questCounts.end() ? 0.0 : it->second; }
int QuestBook_Tier() { return questTier; }
void QuestBook_SetTier(int t) { questTier = std::max(0, t); }

void QuestBook_Reset() {
    LoadTiers();
    questTier = 0;
    placedItems.clear();
    visitedBodies.clear();
    reachedOrbit = visitedMoon = visitedRing = false;
    questCounts.clear();
    QuestTree_Reset();
    for (Tier& t : tiers) for (Task& k : t.tasks) k.done = false;
}

void QuestBook_Update() {
    LoadTiers();
    while (questTier < (int)tiers.size()) {
        Tier& t = tiers[questTier];
        bool open;
        if (t.tasks.empty()) {
            open = t.gate != "-" && (inventory.obtainedTotal(t.gate) > 0.0 || inventory.count(t.gate) > 0 || placedItems.count(t.gate));
        } else {
            open = true;
            for (Task& k : t.tasks) {
                if (!k.done && TaskMet(k)) k.done = true;
                open = open && k.done;
            }
        }
        if (!open) break;
        ++questTier;
        if (questTier < (int)tiers.size()) GameUI_Toast(TextFormat("Tier %d: %s", questTier, tiers[questTier].name.c_str()));
    }
    QuestTree_Update(); // the side lines
}

int QuestBook_TierCount() { LoadTiers(); return (int)tiers.size(); }
const char* QuestBook_TierName(int i) { LoadTiers(); return i >= 0 && i < (int)tiers.size() ? tiers[i].name.c_str() : ""; }
const char* QuestBook_TierText(int i) { LoadTiers(); return i >= 0 && i < (int)tiers.size() ? tiers[i].text.c_str() : ""; }
int QuestBook_TaskCount(int i) { LoadTiers(); return i >= 0 && i < (int)tiers.size() ? (int)tiers[i].tasks.size() : 0; }
const char* QuestBook_TaskText(int i, int k) { return k >= 0 && k < QuestBook_TaskCount(i) ? tiers[i].tasks[k].text.c_str() : ""; }
bool QuestBook_TaskDone(int i, int k) { return k >= 0 && k < QuestBook_TaskCount(i) && (tiers[i].tasks[k].done || i < questTier); }

void QuestBook_Save(FILE* f) {
    fprintf(f, "quest %d\n", questTier);
    for (size_t i = 0; i < tiers.size(); ++i)
        for (size_t k = 0; k < tiers[i].tasks.size(); ++k)
            if (tiers[i].tasks[k].done) fprintf(f, "quest_done %zu %zu\n", i, k);
    for (const std::string& p : placedItems) fprintf(f, "placed %s\n", p.c_str());
    for (const auto& kv : inventory.obtained) if (kv.second > 0.0) fprintf(f, "obtained %s %.2f\n", kv.first.c_str(), kv.second);
    fprintf(f, "quest_flags %d %d %d\n", reachedOrbit ? 1 : 0, visitedMoon ? 1 : 0, visitedRing ? 1 : 0);
    for (const auto& kv : questCounts) fprintf(f, "quest_count %s %.2f\n", kv.first.c_str(), kv.second);
    QuestTree_Save(f);
}

bool QuestBook_LoadLine(const char* tok, FILE* f) {
    LoadTiers();
    char id[64] = {};
    if (strcmp(tok, "quest") == 0) { int t = 0; if (fscanf(f, " %d", &t) == 1) questTier = std::max(0, t); return true; }
    if (strcmp(tok, "quest_done") == 0) {
        int i = 0, k = 0;
        if (fscanf(f, " %d %d", &i, &k) == 2 && i >= 0 && i < (int)tiers.size() && k >= 0 && k < (int)tiers[i].tasks.size()) tiers[i].tasks[k].done = true;
        return true;
    }
    if (strcmp(tok, "placed") == 0) { if (fscanf(f, " %63s", id) == 1) placedItems.insert(id); return true; }
    if (QuestTree_LoadLine(tok, f)) return true;
    if (strcmp(tok, "quest_count") == 0) { double n = 0; if (fscanf(f, " %63s %lf", id, &n) == 2) questCounts[id] = n; return true; }
    if (strcmp(tok, "obtained") == 0) { double n = 0; if (fscanf(f, " %63s %lf", id, &n) == 2) inventory.obtained[id] = n; return true; }
    if (strcmp(tok, "quest_flags") == 0) {
        int a = 0, b = 0, c = 0;
        if (fscanf(f, " %d %d %d", &a, &b, &c) == 3) { reachedOrbit = a; visitedMoon = b; visitedRing = c; }
        return true;
    }
    return false;
}

// Wraps text at `chars` characters, drawing each line; returns the y after it.
static int DrawWrapped(const std::string& text, int x, int y, int chars, int size, Color col) {
    std::string line;
    size_t pos = 0;
    while (pos < text.size()) {
        size_t sp = text.find(' ', pos);
        if (sp == std::string::npos) sp = text.size();
        const std::string w = text.substr(pos, sp - pos);
        if (!line.empty() && line.size() + w.size() + 1 > (size_t)chars) { DrawText(line.c_str(), x, y, size, col); y += size + 3; line.clear(); }
        line += (line.empty() ? "" : " ") + w;
        pos = sp + 1;
    }
    if (!line.empty()) { DrawText(line.c_str(), x, y, size, col); y += size + 3; }
    return y;
}

void Objectives_DrawHUD() {
    LoadTiers();
    int y = 64;
    if (questTier >= (int)tiers.size()) {
        DrawText("CURRENT OBJECTIVES   all done", 12, y, 16, Fade(GOLD, 0.9f));
        return;
    }
    const Tier& t = tiers[questTier];
    DrawText(TextFormat("CURRENT OBJECTIVES   tier %d: %s", questTier, t.name.c_str()), 12, y, 16, Fade(GOLD, 0.9f));
    y += 22;
    if (t.tasks.empty() && t.gate != "-") {
        const ItemDef* d = itemDB.find(t.gate);
        DrawText(TextFormat("[ ] Make a %s", d ? d->name.c_str() : t.gate.c_str()), 12, y, 16, RAYWHITE);
        y += 20;
    }
    for (const Task& k : t.tasks) {
        y = DrawWrapped(TextFormat("%s %s", k.done ? "[x]" : "[ ]", k.text.c_str()), 12, y, 58, 15, k.done ? Fade(GREEN, 0.75f) : RAYWHITE);
        y += 1;
    }
    y += 4;
    DrawWrapped(t.text, 12, y, 64, 13, Fade(RAYWHITE, 0.6f)); // the tier's hint
}

bool Survival_IsUsable(const std::string& id) { return id == "oxygen_canister" || id == "power_cell" || Survival_IsFood(id); }

bool Survival_UseItem(const std::string& id) {
    if (Survival_IsFood(id)) {
        // Eat: not when full (no waste). Raw meat may make you sick: food drains twice as fast a minute.
        const Food& food = foods[id];
        if (survival.food < 98.0f && inventory.remove(id, 1)) {
            survival.food = std::min(100.0f, survival.food + food.value);
            QuestBook_Count("eat");
            if (food.raw && GetRandomValue(0, 99) < 30) { survival.sickUntil = (float)GetTime() + 60.0f; GameUI_Toast("You feel sick: cook your meat"); }
            Audio_Play(Sfx::Pickup, 0.5f, 0.8f);
            return true;
        }
        if (survival.food >= 98.0f) GameUI_Toast("You are full");
        Audio_Play(Sfx::Error, 0.4f);
        return false;
    }
    if (id == "oxygen_canister") {
        if (survival.oxygen < 99.0f && inventory.remove("oxygen_canister", 1)) { survival.oxygen = std::min(100.0f, survival.oxygen + 50.0f); Audio_Play(Sfx::Pickup); return true; }
    } else if (id == "power_cell") {
        // Don't waste a cell on a full store (ship while piloting, otherwise the suit).
        const bool full = survival.power >= 99.0f;
        if (!full && inventory.remove("power_cell", 1)) {
            survival.power = std::min(100.0f, survival.power + 50.0f);
            Audio_Play(Sfx::Pickup);
            return true;
        }
    }
    Audio_Play(Sfx::Error, 0.4f);
    return false;
}
