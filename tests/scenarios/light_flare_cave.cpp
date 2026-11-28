// light/flare-cave (#2062, underground plan #1465): a thrown flare lights a dark cave through the light grid, then burns out.
// A sealed 16 x 16 x 4 m room is carved 10 m down at night: dark. The player, standing at one end, throws a flare (the RMB
// path, Flares_ThrowFromPlayer, with a short burn for the test). Checks: one flare left the inventory; it lands on the floor;
// once at rest the grid at the flare and 4 m from it is lit (a source, no point light); the room's far corner stays dark
// past its reach; after the burn the flare is gone and the grid there is dark again. The capture shows dark, lit, dark.
#include "scenario.h"
#include "flares.h"
#include "inventory.h"
#include "lightlevel.h"
#include "player.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <string>

namespace {
struct State {
    Vector3 home{}, ax{}, ay{};
    float r0 = 0.0f;
    Vector3 pt(float x, float y, float h) const {
        const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(home, r0), Vector3Add(Vector3Scale(ax, x), Vector3Scale(ay, y))));
        return Vector3Scale(d, r0 + h);
    }
    TerrainEdit box(float x0, float x1, float y0, float y1, float h0, float h1) const {
        TerrainEdit e{};
        e.center = pt(0.5f * (x0 + x1), 0.5f * (y0 + y1), 0.5f * (h0 + h1));
        e.carve = true; e.box = true; e.half = { 0.5f * (x1 - x0), 0.5f * (y1 - y0), 0.5f * (h1 - h0) };
        e.ax = ax; e.ay = ay; e.az = Vector3Normalize(e.center);
        e.radius = Vector3Length(e.half);
        return e;
    }
    Vector3 rest{};      // where the flare came to rest (body-local)
    int peak = 0;        // the grid's level at the flare while it burned
} S;
constexpr float kFloor = -10.0f, kCeil = -6.0f; // m under the surface
constexpr float kBurn = 8.0f;                    // s: a short burn for the test (the item burns kFlareBurnSeconds)
int Level(Vector3 local) { return LightLevel_At(*player.nearest, local).level; }
int LevelAtFlare() { return Level(Vector3Add(S.rest, Vector3Scale(Vector3Normalize(S.rest), 0.25f))); }
void View(ScenarioCtx& c) { c.camera(S.pt(-7.0f, -1.5f, kFloor + 1.7f), S.pt(3.0f, 1.0f, kFloor + 0.3f), S.home); }
}

SCENARIO("light/flare-cave") {
    s.tags = { "light", "world" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--sun-elev", "-30" }; // night: only the flare lights the room
    s.timeout = 150.0;

    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("carve a dark room", 1.0, [](ScenarioCtx& c) {
        S.home = c.field->findSpawnDir(c.body->seed);
        S.ax = Vector3Normalize(Vector3CrossProduct(S.home, { 0, 1, 0 })); S.ay = Vector3CrossProduct(S.ax, S.home);
        S.r0 = c.field->radius + c.field->elevation(S.home);
        TerrainEdit room = S.box(-8.0f, 8.0f, -8.0f, 8.0f, kFloor, kCeil);
        ApplyTerrainEdits(*c.body, &room, 1);
        Flares_Clear();
        View(c);
        return true;
    });
    s.step("meshed, the room dark", 30.0, [](ScenarioCtx& c) {
        View(c);
        c.waitingFor("an idle world and grid");
        if (!(c.t > 1.0 && WorldIdle() && !LightLevel_Stats().busy)) return false;
        c.eq("the room is dark before the flare", Level(S.pt(0.0f, 0.0f, kFloor + 1.0f)), 0);
        return true;
    });
    s.step("dark: before the throw", 3.0, [](ScenarioCtx& c) {
        static bool cued = false;
        if (!cued) { cued = true; c.cue(); }
        View(c);
        return c.t > 1.0;
    });
    s.step("throw a flare", 1.0, [](ScenarioCtx& c) {
        View(c);
        inventory.add("flare", 2);
        const bool thrown = Flares_ThrowFromPlayer(c.world, kBurn);
        c.check("the throw went", thrown);
        c.eq("one flare left the inventory", inventory.count("flare"), 1);
        c.eq("one flare in flight", Flares_Count(), 1);
        return true;
    });
    s.step("it lands", 5.0, [](ScenarioCtx& c) {
        View(c);
        Vector3 at{}; bool resting = false; float reach = 0.0f;
        if (!Flares_First(at, resting, reach)) { c.fail("the flare vanished in flight"); return true; }
        c.waitingFor("the flare on the floor");
        if (!resting) return false;
        S.rest = at;
        const float h = Vector3Length(at) - Vector3Length(S.pt(0.0f, 0.0f, kFloor));
        c.between("it rests on the floor (m above the carved floor)", h, -0.6, 0.6);
        c.gt("it flew from the hand (m along the room)", Vector3DotProduct(Vector3Subtract(at, S.pt(-7.0f, -1.5f, kFloor)), S.ax), 2.0);
        return true;
    });
    s.step("the grid lights up", 4.0, [](ScenarioCtx& c) {
        View(c);
        Vector3 at{}; bool resting = false; float reach = 0.0f;
        if (Flares_First(at, resting, reach)) S.rest = at;
        S.peak = std::max(S.peak, LevelAtFlare());
        c.waitingFor("the grid lit round the flare: level " + std::to_string(LevelAtFlare()));
        return S.peak >= 100 && !LightLevel_Stats().busy;
    });
    s.step("lit: while it burns", 2.5, [](ScenarioCtx& c) {
        View(c);
        static bool checked = false;
        if (!checked) {
            checked = true;
            c.gt("the grid at the flare is lit (level)", LevelAtFlare(), 100);
            const Vector3 up = Vector3Normalize(S.rest);
            const Vector3 near4 = Vector3Add(Vector3Add(S.rest, Vector3Scale(S.ay, 4.0f)), Vector3Scale(up, 1.0f));
            c.gt("4 m from the flare is lit (level)", Level(near4), 40);
            c.gt("the flare is a grid source", LightLevel_Stats().sources, 0);
        }
        return c.t > 2.0;
    });
    s.step("it burns out", kBurn + 4.0, [](ScenarioCtx& c) {
        View(c);
        c.waitingFor("the flare gone and the grid dark again: level " + std::to_string(LevelAtFlare()));
        return Flares_Count() == 0 && !LightLevel_Stats().busy && LevelAtFlare() == 0;
    });
    s.step("dark again", 2.0, [](ScenarioCtx& c) {
        View(c);
        static bool checked = false;
        if (!checked) {
            checked = true;
            c.eq("the grid at the flare is dark again", LevelAtFlare(), 0);
            c.eq("the room is dark again", Level(S.pt(0.0f, 0.0f, kFloor + 1.0f)), 0);
            c.eq("no grid sources left", LightLevel_Stats().sources, 0);
        }
        return c.t > 1.5;
    });
    s.step("done", 1.0, [](ScenarioCtx& c) {
        c.keepCapture("a flare lights a dark cave and burns out (#2062): for a person to judge");
        return true;
    });
}
