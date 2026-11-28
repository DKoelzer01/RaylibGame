// hazards/heat-depth (#1472): heat with depth against the suit's heat rating (the user's numbers, 2026-11-27). Two rooms are
// carved, 300 m down (heat 1) and 600 m down (heat 2), and the player waits in each (hovering: the test is the heat, not the
// deep ground's collision). 300 m in a suit of rating 0: the bar is half full after 90 s (the HUD warns) and no health is
// lost; a cooling vest (rating 1) cools it; 600 m with nothing: two short, full within 30 s, then it hurts; the thermal
// regulator with a battery pack (rating 2) cools it there.
#include "scenario.h"
#include "inventory.h"
#include "player.h"
#include "suit.h"
#include "survival.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <vector>

namespace {
struct State {
    Vector3 home{}, ax{}, ay{};
    float r0 = 0.0f;
    float heatAt = 0.0f, healthAt = 100.0f;
    Vector3 pt(float x, float y, float h) const {
        const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(home, r0), Vector3Add(Vector3Scale(ax, x), Vector3Scale(ay, y))));
        return Vector3Scale(d, r0 + h);
    }
    TerrainEdit room(float x, float depth) const { // 8 x 8 m, 4 m high, its floor `depth` m down
        TerrainEdit e{};
        e.center = pt(x, 0.0f, -depth + 2.0f);
        e.carve = true; e.box = true; e.half = { 4.0f, 4.0f, 2.0f };
        e.ax = ax; e.ay = ay; e.az = Vector3Normalize(e.center);
        e.radius = Vector3Length(e.half);
        return e;
    }
} S;
constexpr float kRoomA = 300.0f, kRoomB = 600.0f, kRoomBx = 40.0f;

void Unfit() { for (int i = 0; i < (int)Suit_Modules().size(); ++i) if (!Suit_Modules()[i].empty()) Suit_Unfit(i); }
bool Fit(const char* id) { if (inventory.count(id) == 0) inventory.add(id, 1); return Suit_Fit(id); }
// The player hovers in the room (placed once: a placement starts the world loading around it again).
void Place(ScenarioCtx& c, float x, float depth) {
    player.flyMode = true;
    player.placeAt(c.world, *c.body, c.body->toWorld(S.pt(x, -2.5f, -depth + 0.2f)));
}
// Each frame: the camera at eye height in the room, the room in front (the HUD's bar on screen).
void Hold(ScenarioCtx& c, float x, float depth) { c.camera(S.pt(x, -3.0f, -depth + 1.8f), S.pt(x, 3.5f, -depth + 0.8f), S.home); }
}

SCENARIO("hazards/heat-depth") {
    s.tags = { "hazards", "survival", "world" }; // (world: the runner's part 2)
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 420.0;
    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("carve a room 300 m down and one 600 m down; nothing fitted", 1.0, [](ScenarioCtx& c) {
        S.home = c.spawnDirLocal();
        S.ax = Vector3Normalize(Vector3CrossProduct(S.home, { 0, 1, 0 })); S.ay = Vector3CrossProduct(S.ax, S.home);
        S.r0 = c.field->radius + c.field->elevation(S.home);
        const TerrainEdit rooms[2] = { S.room(0.0f, kRoomA), S.room(kRoomBx, kRoomB) };
        ApplyTerrainEdits(*c.body, rooms, 2);
        Unfit();
        c.eq("the suit's heat rating, nothing fitted", Suit_Stats().heat, 0);
        survival.heat = 0.0f; survival.health = 100.0f;
        Place(c, 0.0f, kRoomA);
        Hold(c, 0.0f, kRoomA);
        return true;
    });
    s.step("remesh", 60.0, [](ScenarioCtx& c) {
        Hold(c, 0.0f, kRoomA);
        survival.heat = 0.0f; survival.heatWarned = false; // (the clock starts once the room is there)
        c.waitingFor("the rooms meshed (an idle world)");
        return c.t > 1.0 && WorldIdle();
    });
    s.step("300 m, rating 0: half full in 90 s, the warning, no harm yet", 100.0, [](ScenarioCtx& c) {
        Hold(c, 0.0f, kRoomA);
        if (c.t < 91.0) return false; // (half at 90 s; a second past it, the warning has fired)
        c.eq("heat level 300 m down", survival.heatLevel, 1);
        c.check("the bar about half full", survival.heat > 48.0f && survival.heat < 53.0f, TextFormat("%.1f", survival.heat));
        c.check("the half-full warning shown", survival.heatWarned);
        c.check("no health lost before it is full", survival.health >= 99.9f, TextFormat("%.1f", survival.health));
        c.metric("heat after 90 s at 300 m", survival.heat, "%");
        c.keepCapture("the HEAT bar at half, its warning");
        return true;
    });
    s.step("a cooling vest (rating 1) cools it", 40.0, [](ScenarioCtx& c) {
        static bool fitted = false;
        if (!fitted) { fitted = true; c.check("a cooling vest fitted", Fit("cooling_vest")); S.heatAt = survival.heat; }
        Hold(c, 0.0f, kRoomA);
        if (c.t < 30.0) return false;
        c.eq("rating 1", Suit_Stats().heat, 1);
        c.check("cooled by about a third in 30 s", S.heatAt - survival.heat > 30.0f, TextFormat("%.1f -> %.1f", S.heatAt, survival.heat));
        return true;
    });
    s.step("600 m with nothing fitted: two short, full within 30 s, then it hurts", 100.0, [](ScenarioCtx& c) {
        static bool moved = false;
        if (!moved) { moved = true; Unfit(); S.heatAt = survival.heat; Place(c, kRoomBx, kRoomB); }
        c.waitingFor("the deep room loaded");
        if (!WorldIdle()) return false;
        static double since = -1.0; // (the 40 s count from the world being there, not the move)
        if (since < 0.0) since = c.t;
        Hold(c, kRoomBx, kRoomB);
        if (c.t - since < 40.0) return false;
        c.eq("heat level 600 m down", survival.heatLevel, 2);
        c.eq("rating 0 again", Suit_Stats().heat, 0);
        c.check("full", survival.heat >= 100.0f, TextFormat("%.1f (from %.1f)", survival.heat, S.heatAt));
        c.check("it hurts once full", survival.health < 99.0f, TextFormat("%.1f", survival.health));
        c.check("the harm is the heat", std::string(survival.lastHarm) == "heat", survival.lastHarm);
        S.healthAt = survival.health;
        return true;
    });
    s.step("the thermal regulator with suit power (rating 2) cools it at 600 m", 40.0, [](ScenarioCtx& c) {
        static bool fitted = false;
        if (!fitted) { fitted = true; c.check("a regulator fitted", Fit("thermal_regulator")); c.check("a battery pack fitted", Fit("battery_pack")); }
        Hold(c, kRoomBx, kRoomB);
        if (c.t < 30.0) return false;
        c.eq("rating 2", Suit_Stats().heat, 2);
        c.check("cooling", survival.heat < 70.0f, TextFormat("%.1f", survival.heat));
        c.check("no more harm", survival.health >= S.healthAt - 1.0f, TextFormat("%.1f (was %.1f)", survival.health, S.healthAt));
        return true;
    });
}
