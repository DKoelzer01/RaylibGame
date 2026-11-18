#include "survival.h"
#include "net.h"
#include "gamerng.h"
#include "playerinput.h"
#include "suit.h"
#include "questtree.h"
#include "weather.h"
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
#include "combat.h"
#include "playersession.h"
#include "plots.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include "dungeontraps.h"

// (survival: the local PlayerSession's, src/space/playersession.cpp, #313)
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
    Biome hereBiome = Biome::Count; // where the player stands, as of the last update (Count: none)
    // (#337: for a given player - the local one, or a remote one a server simulates; only the local one sets the HUD's biome)
    float AirTemperature(Scene& world, const Player& player, const SurvivalState& survival, bool isLocal) {
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
        const float sun = Vector3DotProduct(Vector3Normalize(world.sunDirection), player.up); // -1 midnight .. 1 noon
        // The biome climate (Plan 5 B) at the ground here, colder with height above it; dry places swing
        // more between day and night (deserts: +-14 C, jungles and swamps: +-5).
        const float ground = pf->elevation(up);
        const PlanetField::Climate c = pf->climateAt(up, ground);
        if (isLocal) hereBiome = pf->biomeAt(up, ground);
        float t = c.temp - 6.5f * std::max(height - std::max(ground, 0.0f), 0.0f) / 1000.0f + (5.0f + 9.0f * (1.0f - c.moist)) * sun;
        const float depth = ground - height;                  // below the surface here
        if (depth > 3.0f) t += (12.0f - t) * std::min(1.0f, (depth - 3.0f) / 12.0f);
        return t;
    }
}

const char* Survival_BiomeName() { return hereBiome == Biome::Count ? "" : BiomeName(hereBiome); }
bool Survival_IsFood(const std::string& id) { LoadFoods(); return foods.count(id) > 0; }
float Survival_FoodValue(const std::string& id) { LoadFoods(); auto it = foods.find(id); return it == foods.end() ? 0.0f : it->second.value; }

// Fall damage from the speed into the floor at touchdown (#263). On Terra (g 5.9): a 3 m drop is safe, 20 m takes
// ~79, 30 m ~108 (fatal); servo legs make ~6 m safe. A landing that would take 15 or less takes nothing (#536: minor
// falls are ignored), so on Terra about 5 m is safe.
void Survival_OnLandingFor(SurvivalState& survival, const SuitStats& stats, float impactSpeed, bool local) {
    const float safe = stats.servo ? 8.5f : 6.0f; // servo legs soften landings
    constexpr float kGrace = 15.0f; // health points
    const float damage = (impactSpeed - safe) * 8.5f;
    if (damage <= kGrace || survival.dead) return;
    survival.health -= damage;
    survival.fallDamageTaken += damage;
    if (local) Audio_Play(Sfx::Error, 0.6f, 0.6f);
}
void Survival_OnLanding(float impactSpeed) { Survival_OnLandingFor(survival, Suit_Stats(), impactSpeed, true); }

void Survival_Respawn(Scene& world) { Survival_RespawnFor(world, survival, player, true); }
void Survival_RespawnFor(Scene& world, SurvivalState& survival, Player& player, bool local) {
    survival.oxygen = survival.health = 100.0f;
    // Power is kept as it was (#256): respawning used to refill it, and with it the jetpack (which recharges from it on
    // the ground). Whether the starting suit has power or a jetpack at all is #217.
    survival.food = std::max(survival.food, 60.0f);
    survival.sickUntil = 0.0f;
    survival.vacuumTime = 0.0f;
    survival.dead = false;
    player.flyMode = false; // otherwise the free-fly camera would put the player back where they died
    // Wake on your bed if it still stands (#800), else at the start, at rest with the ground there (spawnOn / placeAt set
    // the velocity to the surface's: nothing of the fall or flight that killed you carries over).
    Planetoid* bedBody = nullptr; Vector3 bedFeet{};
    if (survival.hasBed && Survival_BedSpawnFor(survival, bedBody, bedFeet)) {
        player.placeAt(world, *bedBody, bedFeet);
        player.impactSpeed = player.landingSpeed = 0.0f;
        if (local) GameUI_Toast("Respawned at your bed");
        return;
    }
    const bool lostBed = survival.hasBed;
    survival.hasBed = false; // gone or blocked: forgotten (set it again on a new bed)
    Planetoid* start = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get());
    if (start) {
        const auto* pf = dynamic_cast<const PlanetField*>(start->field.get());
        player.spawnOn(world, *start, start->dirToWorld(pf ? pf->findSpawnDir(start->seed) : Vector3{ 0, 1, 0 }));
    } else {
        player.velocity = { 0, 0, 0 };
    }
    player.impactSpeed = player.landingSpeed = 0.0f; // a landing from before the death must not hurt the respawn
    if (local) GameUI_Toast(lostBed ? "Your bed is gone or blocked: respawned at the start" : "Respawned at the start");
}

// Beds (#800) ----------------------------------------------------------------------------------------------------------
bool Survival_SetBedFor(SurvivalState& survival, const Planetoid& body, const Plot& plot, int part, bool local) {
    if (part < 0 || part >= (int)plot.parts.size() || plot.parts[part].type != PartType::Bed) return false;
    const Int3 o = plot.parts[part].origin;
    const bool same = survival.hasBed && survival.bedBody == body.name && survival.bedPlot == plot.id
                   && survival.bedCell[0] == o.x && survival.bedCell[1] == o.y && survival.bedCell[2] == o.z;
    survival.hasBed = true;
    survival.bedBody = body.name;
    survival.bedPlot = plot.id;
    survival.bedCell[0] = o.x; survival.bedCell[1] = o.y; survival.bedCell[2] = o.z;
    if (local) {
        GameUI_Toast(same ? "This is already your bed" : "Spawn set: you will wake at this bed");
        if (!same) Audio_Play(Sfx::Pickup, 0.4f, 0.8f);
    }
    return true;
}

bool Survival_BedSpawnFor(const SurvivalState& survival, Planetoid*& body, Vector3& renderFeet) {
    if (!survival.hasBed) return false;
    body = nullptr;
    for (Planetoid* p : Planetoid::All()) if (p->name == survival.bedBody) { body = p; break; }
    if (!body || !body->plots) return false;
    const Plot* plot = body->plots->byId(survival.bedPlot);
    if (!plot) return false;
    const int part = plot->partAt({ survival.bedCell[0], survival.bedCell[1], survival.bedCell[2] });
    if (part < 0 || plot->parts[part].type != PartType::Bed || plot->parts[part].origin.x != survival.bedCell[0]
        || plot->parts[part].origin.y != survival.bedCell[1] || plot->parts[part].origin.z != survival.bedCell[2]) return false;
    // On top of it, at the middle of its cells, with room to stand (2 m of free cells over the middle).
    std::vector<std::pair<Int3, CellShape>> fp;
    Plot::Footprint(plot->parts[part], fp);
    if (fp.empty()) return false;
    float cx = 0.0f, cz = 0.0f; int top = fp[0].first.y;
    for (const auto& c : fp) { cx += c.first.x + 0.5f; cz += c.first.z + 0.5f; top = std::max(top, c.first.y); }
    cx /= (float)fp.size(); cz /= (float)fp.size();
    const Int3 mid = { (int)floorf(cx), top + 1, (int)floorf(cz) };
    for (int k = 0; k < (int)(2.0f / kCellSize); ++k)
        if (plot->partAt({ mid.x, mid.y + k, mid.z }) >= 0) return false;
    const Vector3 feetPlot = { cx * kCellSize, (top + 1) * kCellSize + 0.02f, cz * kCellSize };
    renderFeet = body->toWorld(plot->toBody(feetPlot));
    return true;
}

// #337: one player's survival step - the local player (Survival_Update), or a remote one a server simulates (local false: no
// UI, sounds, quests, rocket or ship; the weather and temperature at that player). Single player runs exactly as before.
void Survival_Update(Scene& world, float dt, const PlayerInput& in) { Survival_StepFor(world, dt, in, survival, player, inventory, Suit_Stats(), true); }
void Survival_StepFor(Scene& world, float dt, const PlayerInput& in, SurvivalState& survival, Player& player, Inventory& inventory, const SuitStats& stats, bool local) {
    (void)inventory;
    if (survival.dead) {
        if (SimTime() - survival.deathTime > 3.0) Survival_RespawnFor(world, survival, player, local);
        return;
    }
    const Vector3 pr = (player.position - world.renderOrigin).toFloat();
    survival.breathable = false;
    survival.underwater = false;
    // Life support: a rocket's pressurised cockpit (the placeholder ship that also gave it is gone).
    // Life support: a cockpit in breathable air uses the air outside; elsewhere it spends the air loaded
    // into the rocket (oxygen canisters), and when that is gone the suit's own oxygen takes over (the
    // hull still holds pressure, so no vacuum clock).
    {
        const auto& atm0 = world.atmosphere;
        const bool airOutside = atm0.enabled && Vector3Distance(pr, atm0.centre) < atm0.radius + survival.breathableAltitude;
        survival.inLifeSupport = (local && Rocket_Boarded()) && (airOutside || Rocket_BreatheAboard(dt));
    }
    survival.inVacuum = false;
    if (survival.inLifeSupport) survival.vacuumTime = 0.0f;
    if (survival.inLifeSupport) {
        survival.oxygen = std::min(100.0f, survival.oxygen + survival.lifeSupportRefillPerSec * dt);
        survival.power = std::min(100.0f, survival.power + survival.lifeSupportRefillPerSec * dt);
        survival.health = std::min(100.0f, survival.health + 2.0f * dt);
        player.jetpackFuel = std::min(1.0f, player.jetpackFuel + 0.5f * dt); // refill the suit's thruster tank
    } else {
        // A planet's air is breathable up to ~15 km: the suit tops up its oxygen from it instead.
        const auto& atm = world.atmosphere;
        const Vector3 head = Vector3Add(pr, Vector3Scale(player.up, player.eyeHeight));
        survival.inGas = DungeonTraps_InGas(world, head); // (#719) a dungeon's gas pocket
        survival.breathable = atm.enabled && Vector3Distance(pr, atm.centre) < atm.radius + survival.breathableAltitude && !IsUnderwater(world, head) && !survival.inGas;
        // Vacuum (outside the air): a suit not rated for it holds 60 s, then it hurts (4 a second).
        survival.inVacuum = !(atm.enabled && Vector3Distance(pr, atm.centre) < atm.radius + survival.breathableAltitude);
        if (survival.inVacuum && !stats.vacuum) { survival.vacuumTime += dt; if (survival.vacuumTime > 60.0f) survival.health -= 4.0f * std::min(dt, survival.vacuumTime - 60.0f); }
        else survival.vacuumTime = 0.0f;
        if ((local && Rocket_Boarded())) { survival.vacuumTime = 0.0f; survival.inVacuum = false; } // inside a pressurised hull, out of air
        if (survival.breathable) survival.oxygen = std::min(100.0f, survival.oxygen + survival.airRefillPerSec * dt);
        else {
            // Underwater drains it faster (a rebreather slows that); bigger tanks last longer.
            const bool underwater = atm.enabled && IsUnderwater(world, head);
            survival.underwater = underwater;
            survival.oxygen = std::max(0.0f, survival.oxygen - survival.oxygenDrainPerSec * (underwater ? stats.underwater : 1.0f) * (survival.inGas ? 3.0f : 1.0f) / stats.oxygen * dt);
        }
        survival.power = std::max(0.0f, survival.power - (1.0f / stats.power) * (survival.powerDrainPerSec
                                                          + (survival.headlamp && !(local && ship.piloted) ? survival.headlampPowerPerSec : 0.0f)) * dt);
        // Flight pack hover (#163): a slow, steady drain in the air (~7 min of the full bar), three times that
        // boosting (Shift); at zero it cuts out.
        if (player.hovering && !player.grounded)
            survival.power = std::max(0.0f, survival.power - (in.down(Act::Sprint) ? 0.75f : 0.25f) / stats.power * dt);
        player.hoverPowered = survival.power > 0.0f;
        // Firing the jetpack only burns fuel; recharging the tank (on the ground) is what costs suit
        // power, and only as far as the power left can pay for.
        const float refilled = player.jetpackFuel - survival.lastJetpackFuel;
        if (refilled > 0.0f) {
            const float perTank = survival.jetpackPowerPerTank / stats.power;
            const float affordable = survival.power / perTank;
            const float paid = std::min(refilled, affordable);
            player.jetpackFuel = survival.lastJetpackFuel + paid;
            survival.power = std::max(0.0f, survival.power - paid * perTank);
        }
    }
    if (!stats.powered) { survival.power = 0.0f; player.hoverPowered = false; } // #217: no battery pack, no suit power at all
    survival.lastJetpackFuel = player.jetpackFuel;
    if (survival.power <= 0.0f && survival.headlamp && !(local && ship.piloted)) { survival.headlamp = false; if (local) GameUI_Toast("Headlamp out: no suit power"); }
    if (in.hit(Act::Headlamp)) {
        if (survival.power > 0.0f || (local && ship.piloted)) {
            survival.headlamp = !survival.headlamp;
            if (local) Audio_Play(Sfx::Build, 0.25f, survival.headlamp ? 1.8f : 1.5f);
        } else if (local) GameUI_Toast(stats.powered ? "No suit power for the headlamp" : "The headlamp needs suit power: fit a battery pack (LV)");
    }
    if (survival.oxygen <= 0.0f) survival.health -= survival.suffocationDamagePerSec * dt;
    // Temperature and hunger. Food drains faster running or on the jetpack, in the cold, and
    // while sick; the cockpit's life support feeds you. Heat does slow damage.
    survival.temperature = AirTemperature(world, player, survival, local);
    if (!std::isnan(survival.temperature)) { // rain and snow chill the air
        const WeatherState w = local ? Weather_Here() : (player.nearest ? Weather_At(world, Vector3Normalize(player.nearest->toLocal((player.position - world.renderOrigin).toFloat())), WorldTime()) : WeatherState{});
        if (w.kind == WeatherKind::Rain || w.kind == WeatherKind::Storm) survival.temperature -= 4.0f * w.precipitation;
        else if (w.kind == WeatherKind::Snow) survival.temperature -= 3.0f * w.precipitation;
    }
    {
        const Vector3 flat = Vector3Subtract(player.velocity, Vector3Scale(player.up, Vector3DotProduct(player.velocity, player.up)));
        const bool active = player.jetpackFiring || (player.grounded && Vector3Length(flat) > player.walkSpeed * 1.1f);
        float drain = survival.foodDrainPerSec * (active ? 1.6f : 1.0f) * stats.food * Difficulty_HungerScale();
        const float T = survival.temperature;
        if (!std::isnan(T) && T < 5.0f) drain *= 1.0f + stats.cold * (5.0f - T) / 15.0f; // insulation, the regulator
        if (SimTime() < survival.sickUntil) drain *= 2.0f;
        if (!survival.inLifeSupport) survival.food = std::max(0.0f, survival.food - drain * dt);
        if (survival.food <= 0.0f) survival.health -= survival.starvationDamagePerSec * dt;
        else if (survival.food > 80.0f && survival.oxygen > 0.0f) survival.health = std::min(100.0f, survival.health + survival.wellFedHealPerSec * dt);
        if (!std::isnan(T) && T > 35.0f && !stats.heatProof) survival.health -= 0.05f * (T - 35.0f) * dt; // 40 C: 0.25 a second
    }

    if (in.enabled) {
    }

    if (local) { // A night survived (the quest tree): dark fell and day came back without a death.
        static bool inNight = false;
        static int deathsAtDusk = 0;
        const bool dark = Vector3DotProduct(Vector3Normalize(world.sunDirection), player.up) < -0.1f && player.nearest != nullptr;
        if (dark && !inNight) { inNight = true; deathsAtDusk = survival.deaths; }
        else if (!dark && inNight) { inNight = false; if (survival.deaths == deathsAtDusk) QuestBook_Count("nights"); }
    }
    // (#337) a client process's player dies only when the server says (Vitals): its own drains here just keep the bars moving
    // between the server's numbers, and never cost the inventory it mirrors.
    if (local && gNetRole == NetRole::Client) survival.health = std::max(survival.health, 0.0f);
    else if (survival.health <= 0.0f) {
        survival.health = 0.0f;
        survival.dead = true;
        survival.deathTime = SimTime();
        survival.deaths++;
        player.autoWalk = false;
        if (local) Combat_OnDeath(LocalSession()); // what is lost, by difficulty (#165; a server's remote player: #337 slice 2)
    }
}

static float hudS = 1.0f; // (#827) the interface size, for the vitals
static void Bar(int x, int y, int w, const char* label, float value, Color c, bool warn) {
    const int bh = (int)(16 * hudS);
    DrawRectangle(x, y, w, bh, Fade(BLACK, 0.55f));
    const bool flash = warn && fmodf((float)GetTime(), 0.8f) < 0.4f;
    DrawRectangle(x, y, (int)(w * value / 100.0f), bh, flash ? RED : c);
    DrawRectangleLines(x, y, w, bh, Fade(RAYWHITE, 0.4f));
    DrawText(TextFormat("%s %3.0f", label, value), x + (int)(6 * hudS), y + (int)hudS, (int)(14 * hudS), RAYWHITE);
}

void Survival_DrawHUD() {
    hudS = GameUI_Scale();
    const float s = hudS;
    const int x = (int)(20 * s), y = GetScreenHeight() - (int)(192 * s);
    // The oxygen and power bars grow with the suit's capacity (basic suit: 260 px).
    const int wO2 = (int)((180.0f + 80.0f * std::min(Suit_Stats().oxygen, 2.5f)) * s), wPwr = (int)((180.0f + 80.0f * std::min(Suit_Stats().power, 2.5f)) * s), w260 = (int)(260 * s);
    // Full oxygen in breathable air (#277): no O2 bar and no "BREATHABLE AIR" label; the rows close up from the bottom.
    const bool showO2 = !(survival.oxygen >= 99.9f && survival.breathable && !survival.underwater);
    if (showO2) Bar(x, y, wO2, "O2", survival.oxygen, { 90, 200, 255, 255 }, survival.oxygen < 20.0f);
    if (Suit_Stats().powered) Bar(x, y + (int)(22 * s), wPwr, "PWR", survival.power, { 250, 210, 70, 255 }, survival.power < 15.0f); // hidden with no battery pack (#217)
    Bar(x, y + (int)(44 * s), w260, "HP", survival.health, { 230, 70, 70, 255 }, survival.health < 30.0f);
    Bar(x, y + (int)(66 * s), w260, "FOOD", survival.food, { 120, 200, 90, 255 }, survival.food < 15.0f);
    // Stamina (#165): spent by blocks; only while it is below full, under the others.
    if (survival.stamina < 99.9f) Bar(x, y + (int)(88 * s), w260, "STA", survival.stamina, { 235, 235, 235, 255 }, survival.stamina < 20.0f);
    // The air temperature, coloured when it costs you (cold: hunger; heat: health).
    const float T = survival.temperature;
    if (!std::isnan(T)) {
        const Color tc = T < 5.0f ? Color{ 130, 190, 255, 255 } : T > 35.0f ? Color{ 255, 140, 70, 255 } : Fade(RAYWHITE, 0.8f);
        const SuitStats& st = Suit_Stats();
        const char* note = T < 5.0f ? (st.cold <= 0.0f ? "  cold: regulated" : st.cold < 1.0f ? "  cold: insulated" : "  cold: hungrier") : T > 35.0f ? (st.heatProof ? "  heat: regulated" : "  heat") : "";
        DrawText(TextFormat("%.0f C%s", T, note), x + std::max(wO2, wPwr) + (int)(20 * s), y + (int)(68 * s), (int)(16 * s), tc);
    }
    const int tx = x + std::max(wO2, wPwr) + (int)(20 * s), fs = (int)(16 * s);
    if (SimTime() < survival.sickUntil) DrawText("SICK", tx, y + (int)(46 * s), fs, ORANGE);
    if (survival.inLifeSupport) DrawText("LIFE SUPPORT", tx, y + (int)(2 * s), fs, SKYBLUE);
    else if (survival.inGas) DrawText("TOXIC GAS", tx, y + (int)(2 * s), fs, Fade({ 170, 220, 80, 255 }, 0.9f));
    else if (survival.breathable && showO2) DrawText("BREATHABLE AIR", tx, y + (int)(2 * s), fs, Fade(SKYBLUE, 0.7f));
    else if (survival.inVacuum) {
        if (Suit_Stats().vacuum) DrawText("VACUUM: suit rated", tx, y + (int)(2 * s), fs, Fade(SKYBLUE, 0.7f));
        else { const bool flash = survival.vacuumTime > 50.0f && fmodf((float)GetTime(), 0.6f) < 0.3f; DrawText(survival.vacuumTime < 60.0f ? TextFormat("VACUUM: suit not rated, %.0f s", 60.0f - survival.vacuumTime) : "VACUUM: suit failing - get to air or a cockpit", tx, y + (int)(2 * s), fs, flash ? RED : ORANGE); }
    }
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
    struct Task { std::string kind, target; double n = 1.0; std::string text; };
    struct Tier { std::string name, gate, text; std::vector<Task> tasks; };
    std::vector<Tier> tiers;
    // (#338) the selected book: the local game's, or (on a server) the team being counted or evaluated
    QuestBookState localBook;
    QuestBookState* book = &localBook;
    QuestBookState& Q() { return *book; }
    double ObtainedTotal(const std::string& id) {
        if (!Q().ownObtained) return inventory.obtainedTotal(id);
        const auto& m = Q().obtained;
        if (id.rfind("any_", 0) == 0) { double t = 0.0; for (const auto& kv : m) if (ItemMatches(id, kv.first)) t += kv.second; return t; }
        auto it = m.find(id); return it == m.end() ? 0.0 : it->second;
    }
    bool Done(int tier, int task) { return Q().tasksDone.count({ tier, task }) > 0; }

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
        if (t.kind == "get") return ObtainedTotal(t.target) >= t.n - 1e-6;
        if (t.kind == "place") return Q().placed.count(t.target) > 0;
        if (t.kind == "orbit") return Q().reachedOrbit;
        if (t.kind == "visit") return t.target == "moon" ? Q().visitedMoon : Q().visitedRing;
        if (t.kind == "do") { auto it = Q().counts.find(t.target); return it != Q().counts.end() && it->second >= t.n - 1e-6; }
        return false;
    }
    // A counted task's text with how far along it is (#563: "Hammer 4 bronze plates (2/4)"); tasks of one thing as written.
    std::string TaskLine(const Task& t, bool done) {
        if ((t.kind != "get" && t.kind != "do") || t.n < 2.0 - 1e-6) return t.text;
        double have = 0.0;
        if (t.kind == "get") have = ObtainedTotal(t.target);
        else if (auto it = Q().counts.find(t.target); it != Q().counts.end()) have = it->second;
        if (done) have = t.n;
        return TextFormat("%s (%d/%d)", t.text.c_str(), (int)std::min(floor(have + 1e-6), t.n), (int)t.n);
    }
}

void Objectives_MarkVisitedBody(int i) { Q().visitedBodies.insert(i); }

void Objectives_Update(Scene& world) {
    if (Rocket_Telemetry(world).periapsis > 35000.0 && Rocket_Boarded()) Q().reachedOrbit = true;
    // Visiting: within 2 radii of any body other than the start one.
    const Vector3 pr = (player.position - world.renderOrigin).toFloat();
    int index = 0;
    for (auto& obj : world.objects) {
        if (auto* p = dynamic_cast<Planetoid*>(obj.get())) {
            if (index > 0 && Vector3Distance(pr, p->position) < p->radius * p->scale * 2.0f) {
                Q().visitedBodies.insert(index);
                if (p->kind == BodyKind::Moon) Q().visitedMoon = true;
                else if (p->kind == BodyKind::Asteroid) Q().visitedRing = true;
            }
            index++;
        }
    }
    QuestBook_Update();
}

QuestBookState& QuestBook_LocalState() { return localBook; }
void QuestBook_Select(QuestBookState* b) { book = b ? b : &localBook; }
QuestBookState& QuestBook_Selected() { return Q(); }
void QuestBook_NotePlaced(const std::string& item) { Q().placed.insert(item); }
bool QuestBook_Met(const std::string& kind, const std::string& target, double n) { Task t; t.kind = kind; t.target = target; t.n = n; return TaskMet(t); }
void QuestBook_Count(const std::string& event, double n) { Q().counts[event] += n; }
void QuestBook_AtLeast(const std::string& event, double v) { double& c = Q().counts[event]; c = std::max(c, v); }
double QuestBook_CountOf(const std::string& event) { auto it = Q().counts.find(event); return it == Q().counts.end() ? 0.0 : it->second; }
int QuestBook_Tier() { return Q().tier; }
void QuestBook_SetTier(int t) { Q().tier = std::max(0, t); }

void QuestBook_Reset() {
    LoadTiers();
    Q().tier = 0;
    Q().placed.clear();
    Q().visitedBodies.clear();
    Q().reachedOrbit = Q().visitedMoon = Q().visitedRing = false;
    Q().counts.clear();
    QuestTree_Reset();
    Q().tasksDone.clear();
    if (Q().ownObtained) Q().obtained.clear(); // (#338) a team book's own totals (the local book's live in the inventory)
}

void QuestBook_Update() {
    LoadTiers();
    // (#338, #816) a joined client's book is its team's, decided and paid on the server (NetMsg QuestBook): nothing here
    if (gNetRole == NetRole::Client && book == &localBook) return;
    while (Q().tier < (int)tiers.size()) {
        Tier& t = tiers[Q().tier];
        bool open;
        if (t.tasks.empty()) {
            open = t.gate != "-" && (ObtainedTotal(t.gate) > 0.0 || (!Q().ownObtained && inventory.count(t.gate) > 0) || Q().placed.count(t.gate));
        } else {
            open = true;
            for (int k = 0; k < (int)t.tasks.size(); ++k) {
                if (!Done(Q().tier, k) && TaskMet(t.tasks[k])) Q().tasksDone.insert({ Q().tier, k });
                open = open && Done(Q().tier, k);
            }
        }
        if (!open) break;
        ++Q().tier;
        if (Q().tier < (int)tiers.size()) {
            const char* news = TextFormat("Tier %d: %s", Q().tier, tiers[Q().tier].name.c_str());
            if (gQuestNews) gQuestNews(news); else GameUI_Toast(news);
        }
    }
    QuestTree_Update(); // the side lines
}

int QuestBook_TierCount() { LoadTiers(); return (int)tiers.size(); }
const char* QuestBook_TierName(int i) { LoadTiers(); return i >= 0 && i < (int)tiers.size() ? tiers[i].name.c_str() : ""; }
const char* QuestBook_TierText(int i) { LoadTiers(); return i >= 0 && i < (int)tiers.size() ? tiers[i].text.c_str() : ""; }
int QuestBook_TaskCount(int i) { LoadTiers(); return i >= 0 && i < (int)tiers.size() ? (int)tiers[i].tasks.size() : 0; }
const char* QuestBook_TaskText(int i, int k) { // with its count (#563)
    static std::string s;
    if (k < 0 || k >= QuestBook_TaskCount(i)) return "";
    s = TaskLine(tiers[i].tasks[k], QuestBook_TaskDone(i, k));
    return s.c_str();
}
bool QuestBook_TaskDone(int i, int k) { return k >= 0 && k < QuestBook_TaskCount(i) && (Done(i, k) || i < Q().tier); }

void QuestBook_Save(FILE* f) {
    fprintf(f, "quest %d\n", Q().tier);
    for (size_t i = 0; i < tiers.size(); ++i)
        for (size_t k = 0; k < tiers[i].tasks.size(); ++k)
            if (Done((int)i, (int)k)) fprintf(f, "quest_done %zu %zu\n", i, k);
    for (const std::string& p : Q().placed) fprintf(f, "placed %s\n", p.c_str());
    for (const auto& kv : (Q().ownObtained ? Q().obtained : inventory.obtained)) if (kv.second > 0.0) fprintf(f, "obtained %s %.2f\n", kv.first.c_str(), kv.second);
    fprintf(f, "quest_flags %d %d %d\n", Q().reachedOrbit ? 1 : 0, Q().visitedMoon ? 1 : 0, Q().visitedRing ? 1 : 0);
    for (const auto& kv : Q().counts) fprintf(f, "quest_count %s %.2f\n", kv.first.c_str(), kv.second);
    QuestTree_Save(f);
}

bool QuestBook_LoadLine(const char* tok, FILE* f) {
    LoadTiers();
    char id[64] = {};
    if (strcmp(tok, "quest") == 0) { int t = 0; if (fscanf(f, " %d", &t) == 1) Q().tier = std::max(0, t); return true; }
    if (strcmp(tok, "quest_done") == 0) {
        int i = 0, k = 0;
        if (fscanf(f, " %d %d", &i, &k) == 2 && i >= 0 && i < (int)tiers.size() && k >= 0 && k < (int)tiers[i].tasks.size()) Q().tasksDone.insert({ i, k });
        return true;
    }
    if (strcmp(tok, "placed") == 0) { if (fscanf(f, " %63s", id) == 1) Q().placed.insert(id); return true; }
    if (QuestTree_LoadLine(tok, f)) return true;
    if (strcmp(tok, "quest_count") == 0) { double n = 0; if (fscanf(f, " %63s %lf", id, &n) == 2) Q().counts[id] = n; return true; }
    if (strcmp(tok, "obtained") == 0) { double n = 0; if (fscanf(f, " %63s %lf", id, &n) == 2) (Q().ownObtained ? Q().obtained : inventory.obtained)[id] = n; return true; }
    if (strcmp(tok, "quest_flags") == 0) {
        int a = 0, b = 0, c = 0;
        if (fscanf(f, " %d %d %d", &a, &b, &c) == 3) { Q().reachedOrbit = a; Q().visitedMoon = b; Q().visitedRing = c; }
        return true;
    }
    return false;
}

// Wraps text at `chars` characters, drawing each line; returns the y after it.
// draw false: only measure (wrapMaxW: the widest line so far, for a backing)
static int wrapMaxW = 0;
static int DrawWrapped(const std::string& text, int x, int y, int chars, int size, Color col, bool draw = true) {
    std::string line;
    auto put = [&](const std::string& l) { if (draw) DrawText(l.c_str(), x, y, size, col); wrapMaxW = std::max(wrapMaxW, MeasureText(l.c_str(), size)); };
    size_t pos = 0;
    while (pos < text.size()) {
        size_t sp = text.find(' ', pos);
        if (sp == std::string::npos) sp = text.size();
        const std::string w = text.substr(pos, sp - pos);
        if (!line.empty() && line.size() + w.size() + 1 > (size_t)chars) { put(line); y += size + 3; line.clear(); }
        line += (line.empty() ? "" : " ") + w;
        pos = sp + 1;
    }
    if (!line.empty()) { put(line); y += size + 3; }
    return y;
}

// Pinned quests (#825): up to 3 tasks (tier, task) the player picked on the Quests screen; the HUD shows them. A player's
// own (saved in their record), the quest book itself the team's. A finished task drops out of the pins.
namespace { std::vector<std::pair<int, int>> pins; }
const std::vector<std::pair<int, int>>& QuestBook_Pins() { return pins; }
bool QuestBook_Pinned(int tier, int task) { return std::find(pins.begin(), pins.end(), std::make_pair(tier, task)) != pins.end(); }
bool QuestBook_Pin(int tier, int task, bool on) {
    const auto it = std::find(pins.begin(), pins.end(), std::make_pair(tier, task));
    if (!on) { if (it != pins.end()) pins.erase(it); return true; }
    if (it != pins.end()) return true;
    if ((int)pins.size() >= kMaxPins || task < 0 || task >= QuestBook_TaskCount(tier)) return false;
    pins.push_back({ tier, task });
    return true;
}
void QuestBook_ClearPins() { pins.clear(); }
void QuestBook_SavePins(FILE* f) { for (const auto& [t, k] : pins) fprintf(f, "quest_pin %d %d\n", t, k); }
bool QuestBook_LoadPinLine(const char* tok, FILE* f) {
    if (strcmp(tok, "quest_pin") != 0) return false;
    int t = 0, k = 0;
    if (fscanf(f, " %d %d", &t, &k) == 2) QuestBook_Pin(t, k, true);
    return true;
}

void Objectives_DrawHUD() {
    LoadTiers();
    // Finished pins drop out (the next ones are the player's to pick).
    pins.erase(std::remove_if(pins.begin(), pins.end(), [](const std::pair<int, int>& p) { return QuestBook_TaskDone(p.first, p.second); }), pins.end());
    // (#827) Laid out twice: measured, then drawn over a dark backing (white text over bright sky was hard to read); the
    // interface size scales it.
    const float s = GameUI_Scale();
    const int x = (int)(12 * s), top = (int)(64 * s), big = (int)(16 * s), mid = (int)(15 * s), small = (int)(13 * s);
    auto layout = [&](bool draw) {
        int y = top;
        wrapMaxW = 0;
        auto text = [&](const char* t, int size, Color c) { if (draw) DrawText(t, x, y, size, c); wrapMaxW = std::max(wrapMaxW, MeasureText(t, size)); y += size + (int)(6 * s); };
        if (Q().tier >= (int)tiers.size() && pins.empty()) { text("QUESTS   all done", big, Fade(GOLD, 0.9f)); return y; }
        if (!pins.empty()) { // the pinned ones, each with its tier
            text("PINNED QUESTS", big, Fade(GOLD, 0.9f));
            for (const auto& [ti, ki] : pins) {
                if (ti < 0 || ti >= (int)tiers.size() || ki < 0 || ki >= (int)tiers[ti].tasks.size()) continue;
                y = DrawWrapped(TextFormat("[ ] %s", TaskLine(tiers[ti].tasks[ki], false).c_str()), x, y, 58, mid, RAYWHITE, draw);
                y = DrawWrapped(TextFormat("    tier %d: %s", ti, tiers[ti].name.c_str()), x, y, 64, small, Fade(RAYWHITE, 0.55f), draw);
                y += (int)(2 * s);
            }
            return y;
        }
        // Nothing pinned: the current tier's first three open tasks (pin others on the Quests screen).
        const Tier& t = tiers[Q().tier];
        text(TextFormat("QUESTS   tier %d: %s", Q().tier, t.name.c_str()), big, Fade(GOLD, 0.9f));
        if (t.tasks.empty() && t.gate != "-") {
            const ItemDef* d = itemDB.find(t.gate);
            text(TextFormat("[ ] Make a %s", d ? d->name.c_str() : t.gate.c_str()), big, RAYWHITE);
        }
        int shown = 0, more = 0;
        for (int ki = 0; ki < (int)t.tasks.size(); ++ki) {
            if (Done(Q().tier, ki)) continue;
            if (shown >= kMaxPins) { ++more; continue; }
            y = DrawWrapped(TextFormat("[ ] %s", TaskLine(t.tasks[ki], false).c_str()), x, y, 58, mid, RAYWHITE, draw);
            y += (int)(1 * s);
            ++shown;
        }
        if (more > 0) text(TextFormat("+%d more (pin the ones you want: Quests)", more), small, Fade(RAYWHITE, 0.55f));
        return y;
    };
    const int bottom = layout(false);
    const int pad = (int)(6 * s);
    DrawRectangleRounded({ (float)(x - pad), (float)(top - pad), (float)(wrapMaxW + 2 * pad), (float)(bottom - top + pad) }, 0.12f, 6, Fade(BLACK, 0.42f));
    layout(true);
}

bool Survival_IsUsable(const std::string& id) { return id == "oxygen_canister" || id == "power_cell" || id == "repair_kit" || Survival_IsFood(id); }

// #337: using an item - the local player's (a client process asks the server instead: its inventory and vitals are the
// server's; nothing changes here until the server's deltas come back), or a remote player's on a server (local false: no
// UI, sounds or quests).
bool Survival_UseItem(const std::string& id) {
    if (gNetRole == NetRole::Client && gSurvivalUseRequest && id != "repair_kit") { gSurvivalUseRequest(id); return true; }
    return Survival_UseItemFor(id, survival, inventory, Suit_Stats(), true);
}
bool Survival_UseItemFor(const std::string& id, SurvivalState& survival, Inventory& inventory, const SuitStats& stats, bool local) {
    if (id == "repair_kit") { // mends the most damaged part of the rocket in reach (rocket.h)
        if (local && inventory.count(id) > 0 && Rocket_UseRepairKit()) { inventory.remove(id, 1); return true; }
        return false;
    }
    if (Survival_IsFood(id)) {
        // Eat: not when full (no waste). Raw meat may make you sick: food drains twice as fast a minute.
        const Food& food = foods[id];
        if (survival.food < 98.0f && inventory.remove(id, 1)) {
            survival.food = std::min(100.0f, survival.food + food.value);
            if (local) QuestBook_Count("eat");
            if (food.raw && GameRng_Int(0, 99) < 30) { survival.sickUntil = (float)SimTime() + 60.0f; if (local) GameUI_Toast("You feel sick: cook your meat"); }
            if (local) Audio_Play(Sfx::Pickup, 0.5f, 0.8f);
            return true;
        }
        if (survival.food >= 98.0f) if (local) GameUI_Toast("You are full");
        if (local) Audio_Play(Sfx::Error, 0.4f);
        return false;
    }
    if (id == "oxygen_canister") {
        if (survival.oxygen < 99.0f && inventory.remove("oxygen_canister", 1)) { survival.oxygen = std::min(100.0f, survival.oxygen + 50.0f); if (local) Audio_Play(Sfx::Pickup); return true; }
    } else if (id == "power_cell") {
        if (!stats.powered) { if (local) GameUI_Toast("A power cell needs a battery pack in the suit (#217)"); if (local) Audio_Play(Sfx::Error, 0.4f); return false; }
        // Don't waste a cell on a full store (ship while piloting, otherwise the suit).
        const bool full = survival.power >= 99.0f;
        if (!full && inventory.remove("power_cell", 1)) {
            survival.power = std::min(100.0f, survival.power + 50.0f);
            if (local) Audio_Play(Sfx::Pickup);
            return true;
        }
    }
    if (local) Audio_Play(Sfx::Error, 0.4f);
    return false;
}
