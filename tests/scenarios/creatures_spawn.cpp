// creatures/cave-spawn (#1223, Plan 12 spawn statistics): cave creatures spawn in the cave graph near the player, by its
// rules. A test row makes the hornhide a cave dweller (habitat=cave, 6 m+ deep, dark, 2 a chunk, nests); the player stands
// on the surface over seed 42's nearest chamber (#1740: none spawn there), then in it. Within a few passes: sites were found and creatures placed, every one of
// them underground (6 m+ under the surface), standing in open space, never within 24 m of the player, and no chunk over its
// cap. On peaceful (forgotten chunks, a fresh pass), none spawn. The counts are reported (deterministic from the seed).
#include "scenario.h"
#include "animals.h"
#include "creatures.h"
#include "player.h"
#include "survival.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <map>
#include <string>
#include <tuple>

namespace {
struct St { Vector3 chamber{}, out{}; int count0 = 0, spawned0 = 0; } S;
Vector3 FeetLocal(ScenarioCtx& c) { return c.body->toLocal((player.position - c.world.renderOrigin).toFloat()); }
}

SCENARIO("creatures/cave-spawn") {
    s.tags = { "creatures", "caves" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 240.0;
    s.record = false;

    s.step("settle", 120.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("a cave dweller, over the nearest chamber", 5.0, [](ScenarioCtx& c) {
        std::string err;
        if (!c.check("the test row parses", Creatures_TestAdd("hornhide | livestock | aggro=passive habitat=cave depth=6-500 dark=0.1 cap=2 respawn=30 nests=1", &err), err)) { c.fail(err); return true; }
        const Vector3 home = c.field->findSpawnDir(c.body->seed);
        if (!c.check("a chamber within 3 km", c.field->nearestCaveMouth(home, 3000.0f, S.chamber, S.out, 3))) { c.fail("no chamber"); return true; }
        c.metric("chamber from the start (m)", Vector3Distance(Vector3Scale(Vector3Normalize(S.chamber), c.field->radius), Vector3Scale(home, c.field->radius)), "m");
        Animals_SetNaturalSpawns(false);
        Animals_SetThreats(true);
        worldDifficulty = 1;
        Animals_ResetCaveSpawns();
        player.flyMode = false;
        player.spawnOn(c.world, *c.body, c.body->dirToWorld(Vector3Normalize(S.chamber)));
        S.count0 = Animals_Count();
        return true;
    });
    // (#1740) On the surface over the chamber: nothing spawns (the player is not underground).
    s.step("on the surface over it: none spawn", 60.0, [](ScenarioCtx& c) {
        c.waitingFor("the terrain over the chamber built, then 5 s (two passes)");
        static double idleAt = -1.0;
        if (idleAt < 0.0 && c.t > 1.0 && WorldIdle()) { idleAt = c.t; S.spawned0 = Animals_CaveSpawned(); Animals_CaveSpawnNow(); }
        if (idleAt < 0.0 || c.t < idleAt + 5.0) return false;
        c.eq("no cave creatures while the player is on the surface", Animals_CaveSpawned() - S.spawned0, 0);
        return true;
    });
    // (#1740) Down in the chamber (its floor under the centre): now they spawn, 24-64 m away, unseen, underground.
    s.step("in the chamber: a few spawn passes", 60.0, [](ScenarioCtx& c) {
        static double idleAt = -1.0;
        if (c.t < 0.05) {
            const Vector3 up = Vector3Normalize(S.chamber);
            Vector3 feet = S.chamber;
            for (float d = 0.0f; d < 40.0f; d += 0.25f) { const Vector3 q = Vector3Subtract(S.chamber, Vector3Scale(up, d)); if (c.body->densityAt(q) >= ISO_LEVEL) { feet = Vector3Add(q, Vector3Scale(up, 0.3f)); break; } }
            player.placeAt(c.world, *c.body, c.body->toWorld(feet));
            idleAt = -1.0;
            return false;
        }
        c.waitingFor("the chamber built, then 5 s (two passes)");
        if (idleAt < 0.0 && c.t > 1.0 && WorldIdle()) { idleAt = c.t; Animals_CaveSpawnNow(); }
        if (idleAt < 0.0 || c.t < idleAt + 5.0) return false;
        c.metric("spawn sites in the last pass", Animals_CaveSites(), "", false);
        c.metric("cave creatures spawned", Animals_CaveSpawned(), "", false);
        c.gt("spawn sites found around the chamber", Animals_CaveSites(), 0.5);
        c.gt("cave creatures spawned", Animals_CaveSpawned(), 0.5);
        const Vector3 me = FeetLocal(c);
        int under = 0, open = 0, far = 0, within = 0, n = 0;
        std::map<std::tuple<int, int, int>, int> perChunk;
        for (int i = S.count0; i < Animals_Count(); ++i) {
            if (Animals_Body(i) != c.body->name || (Animals_SpeciesOf(i).first != "hornhide" && Animals_SpeciesOf(i).first != "creature_nest")) continue; // (the cave dwellers and nests only)
            const Vector3 p = Animals_Position(i);
            const Vector3 dir = Vector3Normalize(p);
            const float depth = c.field->radius + c.field->elevation(dir) - Vector3Length(p);
            ++n;
            under += depth >= 5.5f;
            open += c.body->densityAt(Vector3Add(p, Vector3Scale(dir, 0.6f))) < ISO_LEVEL;
            far += Vector3Distance(p, me) >= 23.5f;
            within += Vector3Distance(p, me) <= 64.5f || Animals_SpeciesOf(i).first == "creature_nest";
            ++perChunk[{ (int)floorf(p.x / 64.0f), (int)floorf(p.y / 64.0f), (int)floorf(p.z / 64.0f) }];
        }
        c.metric("new animals", n, "", false);
        c.eq("all underground (6 m+)", under, n);
        c.eq("all in open space", open, n);
        c.eq("none within 24 m of the player", far, n);
        c.eq("all within 64 m of the player (or nests)", within, n);
        int worst = 0; std::tuple<int, int, int> worstKey{};
        for (const auto& kv : perChunk) if (kv.second > worst) { worst = kv.second; worstKey = kv.first; }
        if (worst > 3) // (what is in the fullest chunk, and how near each stands to its faces: a wanderer in from the next one?)
            for (int i = S.count0; i < Animals_Count(); ++i) {
                const Vector3 p = Animals_Position(i);
                if (std::tuple<int, int, int>{ (int)floorf(p.x / 64.0f), (int)floorf(p.y / 64.0f), (int)floorf(p.z / 64.0f) } != worstKey) continue;
                const float ex = std::min(p.x - 64.0f * floorf(p.x / 64.0f), 64.0f - (p.x - 64.0f * floorf(p.x / 64.0f)));
                const float ey = std::min(p.y - 64.0f * floorf(p.y / 64.0f), 64.0f - (p.y - 64.0f * floorf(p.y / 64.0f)));
                const float ez = std::min(p.z - 64.0f * floorf(p.z / 64.0f), 64.0f - (p.z - 64.0f * floorf(p.z / 64.0f)));
                c.note(TextFormat("fullest chunk: %s, %.1f m from its nearest face", Animals_SpeciesOf(i).first.c_str(), std::min(ex, std::min(ey, ez))));
            }
        c.lt("most in one chunk (cap 2 plus a nest)", worst, 3.5);
        return true;
    });
    // (#1717 / #1227, playtest: "a lost miner wandering around near dawn" on open grass) 20 s of living there: none ends up on
    // the surface (a creature off its feet used to take the surface fallback and stand over its cave by day).
    s.step("20 s later: all still underground", 30.0, [](ScenarioCtx& c) {
        if (c.t < 20.0) { c.waitingFor("20 s of the cave creatures living"); return false; }
        int under = 0, n = 0;
        for (int i = S.count0; i < Animals_Count(); ++i) {
            if (Animals_Body(i) != c.body->name || Animals_SpeciesOf(i).first != "hornhide") continue;
            const Vector3 p = Animals_Position(i), dir = Vector3Normalize(p);
            ++n;
            under += c.field->radius + c.field->elevation(dir) - Vector3Length(p) >= 3.0f;
        }
        c.metric("cave creatures after 20 s", n, "", false);
        c.eq("all still underground (3 m+) after 20 s", under, n);
        return true;
    });
    // (#1740) The player leaves (150 m away, on the surface): the cave dwellers go home (the nests stay).
    s.step("far away: they go home", 30.0, [](ScenarioCtx& c) {
        static bool moved = false;
        if (!moved) {
            moved = true;
            const Vector3 up = Vector3Normalize(S.chamber);
            const Vector3 side = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
            player.spawnOn(c.world, *c.body, c.body->dirToWorld(Vector3Normalize(Vector3Add(S.chamber, Vector3Scale(side, 150.0f)))));
            return false;
        }
        // (#1756) The move brings up a loading screen and animals do not update under it (main.cpp: !loading): wait for the
        // world to be idle, then 2 s of updates (past 80 m they run every 4th frame), not 3 s of sim time from the move.
        static double idleAt = -1.0;
        c.waitingFor("the terrain at the new spot built (an idle world), then 2 s");
        if (idleAt < 0.0 && c.t > 0.5 && WorldIdle()) idleAt = c.t;
        if (idleAt < 0.0 || c.t < idleAt + 2.0) return false;
        int left = 0;
        for (int i = S.count0; i < Animals_Count(); ++i) left += Animals_Body(i) == c.body->name && Animals_SpeciesOf(i).first == "hornhide";
        c.eq("cave dwellers left after the player went 150 m away", left, 0);
        return true;
    });
    s.step("peaceful: none", 15.0, [](ScenarioCtx& c) {
        static bool armed = false;
        if (!armed) { armed = true; worldDifficulty = 0; Animals_ResetCaveSpawns(); Animals_CaveSpawnNow(); S.spawned0 = Animals_CaveSpawned(); return false; }
        if (c.t < 5.0) return false;
        c.eq("no cave creatures on peaceful", Animals_CaveSpawned() - S.spawned0, 0);
        worldDifficulty = 1;
        return true;
    });
}
