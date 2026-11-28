// light/torch-room (#1739 PR 2): many torches light a room through the light grid, not only the nearest few.
// User 2026-11-23: "placing a lot of torches didnt actually fully light up the area, just what was around the player
// (closest lights)". A sealed 24 x 24 x 4 m room is carved 10 m down and 36 torches stand on its floor, 4 m apart. Checks:
// every torch is a grid source; the room's far corner and the middle between torches are lit; the rock over the ceiling
// and the open air on the surface above stay dark (no light through rock); a tunnel dug out of the room lights up along
// its first metres (terrain edits relight). Views from a corner and the middle are kept for a person to judge.
#include "scenario.h"
#include "animals.h"
#include "creatures.h"
#include "lightlevel.h"
#include "player.h"
#include "plots.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
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
    int torches = 0;
    float worstMs = 0.0f;
} S;
constexpr float kFloor = -10.0f, kCeil = -6.0f; // m under the surface
int Level(Vector3 local) { return LightLevel_At(*player.nearest, local).level; }
}

SCENARIO("light/torch-room") {
    s.tags = { "light", "world" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--sun-elev", "-30" }; // night: only the torches light the views
    s.timeout = 150.0;

    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("carve the room, 36 torches", 1.0, [](ScenarioCtx& c) {
        S.home = c.field->findSpawnDir(c.body->seed);
        S.ax = Vector3Normalize(Vector3CrossProduct(S.home, { 0, 1, 0 })); S.ay = Vector3CrossProduct(S.ax, S.home);
        S.r0 = c.field->radius + c.field->elevation(S.home);
        TerrainEdit room = S.box(-12.0f, 12.0f, -12.0f, 12.0f, kFloor, kCeil);
        ApplyTerrainEdits(*c.body, &room, 1);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        Plot* plot = c.body->plots->plotFor(S.pt(0.0f, 0.0f, kFloor + 0.25f), 0.0f, true);
        for (int i = -10; i <= 10; i += 4)          // plot cells are 0.5 m: torches every 4 m (8 cells), x and z from -10 m to 10 m
            for (int j = -10; j <= 10; j += 4)
                if (plot->place({ PartType::Torch, { i * 2, 0, j * 2 }, 0 }) >= 0) ++S.torches;
        c.eq("torches placed", S.torches, 36);
        c.camera(S.pt(-11.0f, -11.0f, kFloor + 1.7f), S.pt(4.0f, 4.0f, kFloor + 0.5f), S.home);
        return true;
    });
    s.step("meshed and relit", 30.0, [](ScenarioCtx& c) {
        c.camera(S.pt(-11.0f, -11.0f, kFloor + 1.7f), S.pt(4.0f, 4.0f, kFloor + 0.5f), S.home);
        const LightLevelStats st = LightLevel_Stats();
        S.worstMs = std::max(S.worstMs, st.ms); // (#1739 profile: one relight took 42-45 ms in a frame) spread over frames now
        c.waitingFor("an idle world and the grid relit: " + std::to_string(st.sources) + " sources, busy " + std::to_string(st.busy));
        if (!(c.t > 1.0 && WorldIdle() && !st.busy && st.sources > 0)) return false;
        c.metric("light grid: worst frame while 36 torches lit up (ms)", S.worstMs, "ms", false);
        c.metric("light grid: slowest brick read (ms)", st.brickMsMax, "ms", false);
        // A frame does 1 ms of grid work plus at most one brick read past it (GUESS 8 ms on a shared machine: a brick read
        // alone was ~3 ms on the sandbox); before #1739's budget the whole relight ran in one frame (42-45 ms in play).
        // A timing gate (#1034): on the shared runner one brick read alone took 7.5-52 ms (2/20 passes, #1829), so it is
        // judged on the quiet PC only.
        c.timingLt("the relight is spread over frames (worst frame, ms)", S.worstMs, 8.0);
        return true;
    });
    s.step("the room is lit, the rock is not", 1.0, [](ScenarioCtx& c) {
        const LightLevelStats st = LightLevel_Stats();
        c.eq("every torch is a grid source", st.sources, 36);
        c.note(TextFormat("grid: %d bricks, %.2f ms last frame", st.bricks, st.ms));
        const int middle = Level(S.pt(2.0f, 2.0f, kFloor + 1.0f)), corner = Level(S.pt(-11.5f, 11.5f, kFloor + 1.0f));
        c.check("lit between torches", middle >= 60, TextFormat("level %d", middle));
        c.check("the far corner is lit", corner >= 20, TextFormat("level %d", corner));
        c.eq("the rock 2 m over the ceiling is dark", Level(S.pt(0.0f, 0.0f, kCeil + 2.0f)), 0);
        c.eq("the surface over the room is dark", Level(S.pt(0.0f, 0.0f, 1.5f)), 0);
        return true;
    });
    s.step("from the corner", 3.0, [](ScenarioCtx& c) {
        static bool cued = false;
        if (!cued) { cued = true; c.cue(); }
        c.camera(S.pt(-11.0f, -11.0f, kFloor + 1.7f), S.pt(4.0f, 4.0f, kFloor + 0.5f), S.home);
        return c.t > 2.0;
    });
    s.step("from the middle", 3.0, [](ScenarioCtx& c) {
        c.camera(S.pt(0.0f, -1.0f, kFloor + 1.7f), S.pt(0.0f, 11.0f, kFloor + 1.0f), S.home);
        return c.t > 2.0;
    });
    s.step("dig a tunnel out", 30.0, [](ScenarioCtx& c) {
        static bool dug = false;
        if (!dug) {
            dug = true;
            c.check("the tunnel's rock is dark first", Level(S.pt(14.0f, 2.0f, kFloor + 1.0f)) == 0);
            TerrainEdit t = S.box(11.0f, 34.0f, 1.0f, 3.0f, kFloor, kFloor + 2.5f); // 22 m past the room: its far end is out of every torch's reach
            ApplyTerrainEdits(*c.body, &t, 1);
        }
        c.camera(S.pt(6.0f, 2.0f, kFloor + 1.7f), S.pt(20.0f, 2.0f, kFloor + 1.0f), S.home);
        c.waitingFor("the grid relit after the dig");
        return c.t > 0.5 && !LightLevel_Stats().busy && Level(S.pt(14.0f, 2.0f, kFloor + 1.0f)) > 0;
    });
    s.step("the tunnel is lit near the room", 3.0, [](ScenarioCtx& c) {
        static bool checked = false;
        if (!checked) {
            checked = true;
            const int in = Level(S.pt(14.0f, 2.0f, kFloor + 1.0f));
            c.check("3 m into the tunnel is lit", in >= 20, TextFormat("level %d", in));
            // (#1739 PR 4) the cave spawn rule reads it: no spawn in the lit room, one possible at the dark far end
            const Vector3 up = Vector3Normalize(S.home);
            c.check("a cave spawn site in the room is not dark", !Animals_CaveSiteDark(*c.body, S.pt(2.0f, 2.0f, kFloor + 0.1f), up));
            c.check("a site at the tunnel's far end is dark", Animals_CaveSiteDark(*c.body, S.pt(32.0f, 2.0f, kFloor + 0.1f), up),
                    TextFormat("level %d", Level(S.pt(32.0f, 2.0f, kFloor + 1.1f))));
            // (#1740) the site's light (SpawnSite::light, what Spawn_Fits compares with dark=) is the grid's, not a guess:
            // the lost miner (dark=0.1) fits the dark end, not the lit room.
            const float litRoom = Animals_CaveSiteLight(*c.body, S.pt(2.0f, 2.0f, kFloor + 0.1f), up, false, 10.0f, false);
            const float darkEnd = Animals_CaveSiteLight(*c.body, S.pt(32.0f, 2.0f, kFloor + 0.1f), up, false, 10.0f, false);
            c.metric("spawn site light in the lit room", litRoom, "", false);
            c.check("a spawn site in the lit room reads the grid (light > 0.2)", litRoom > 0.2f, TextFormat("%.3f", litRoom));
            c.check("a spawn site at the dark far end reads dark (light <= 0.07)", darkEnd <= 0.07f, TextFormat("%.3f", darkEnd));
            const CreatureSpec* m = Creature_Spec("lost_miner");
            if (c.check("the lost miner's row is loaded", m != nullptr)) {
                SpawnSite s; s.kind = 0; s.depth = 50.0f; // the miner's depth range; only the light differs
                s.light = litRoom; c.check("the lost miner does not fit the lit site", !Spawn_Fits(*m, s));
                s.light = darkEnd; c.check("the lost miner fits the dark site", Spawn_Fits(*m, s));
            }
        }
        c.camera(S.pt(6.0f, 2.0f, kFloor + 1.7f), S.pt(20.0f, 2.0f, kFloor + 1.0f), S.home);
        return c.t > 2.0;
    });
    s.step("done", 1.0, [](ScenarioCtx& c) {
        c.keepCapture("a torch-lit room through the light grid (#1739): for a person to judge");
        return true;
    });
}
