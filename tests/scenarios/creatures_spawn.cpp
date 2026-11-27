// creatures/cave-spawn (#1223, Plan 12 spawn statistics): cave creatures spawn in the cave graph near the player, by its
// rules. A test row makes the hornhide a cave dweller (habitat=cave, 6 m+ deep, dark, 2 a chunk, nests); the player stands
// on the surface over seed 42's nearest chamber. Within a few passes: sites were found and creatures placed, every one of
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
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <tuple>

namespace {
using ChunkKey = std::tuple<int, int, int>;
struct Spawned { ChunkKey chunk; bool open; };
struct St { Vector3 chamber{}, out{}; int count0 = 0, spawned0 = 0; std::map<uint64_t, Spawned> spawnedIn; } S;
// Open space: air 0.6 m over the feet (radially).
bool OpenAt(ScenarioCtx& c, Vector3 p) { return c.body->densityAt(Vector3Add(p, Vector3Scale(Vector3Normalize(p), 0.6f))) < ISO_LEVEL; }
Vector3 FeetLocal(ScenarioCtx& c) { return c.body->toLocal((player.position - c.world.renderOrigin).toFloat()); }
ChunkKey ChunkOfLocal(Vector3 p) { return { (int)floorf(p.x / 64.0f), (int)floorf(p.y / 64.0f), (int)floorf(p.z / 64.0f) }; }
bool CaveDweller(ScenarioCtx& c, int i) { return Animals_Body(i) == c.body->name && (Animals_SpeciesOf(i).first == "hornhide" || Animals_SpeciesOf(i).first == "creature_nest"); }
// Each new cave dweller's (or nest's) 64 m chunk on the first frame it is seen: where the spawner placed it, before it
// wanders (the pass runs on the next update after Animals_CaveSpawnNow, so this runs every frame of the step).
void NoteSpawns(ScenarioCtx& c) {
    for (int i = S.count0; i < Animals_Count(); ++i)
        if (CaveDweller(c, i) && !S.spawnedIn.count(Animals_NetId(i))) {
            const Vector3 p = Animals_Position(i);
            S.spawnedIn[Animals_NetId(i)] = { ChunkOfLocal(p), OpenAt(c, p) };
        }
}
int FullestAsSpawned() {
    std::map<ChunkKey, int> perChunk;
    int worst = 0;
    for (const auto& kv : S.spawnedIn) worst = std::max(worst, ++perChunk[kv.second.chunk]);
    return worst;
}
int OpenAsSpawned() { int n = 0; for (const auto& kv : S.spawnedIn) n += kv.second.open; return n; }
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
    s.step("a few spawn passes", 60.0, [](ScenarioCtx& c) {
        c.waitingFor("the terrain over the chamber built, then 5 s (two passes)");
        static double idleAt = -1.0;
        if (idleAt < 0.0 && c.t > 1.0 && WorldIdle()) { idleAt = c.t; Animals_CaveSpawnNow(); }
        NoteSpawns(c);
        if (idleAt < 0.0 || c.t < idleAt + 5.0) return false;
        c.metric("spawn sites in the last pass", Animals_CaveSites(), "", false);
        c.metric("cave creatures spawned", Animals_CaveSpawned(), "", false);
        c.gt("spawn sites found around the chamber", Animals_CaveSites(), 0.5);
        c.gt("cave creatures spawned", Animals_CaveSpawned(), 0.5);
        const Vector3 me = FeetLocal(c);
        int under = 0, open = 0, far = 0, n = 0;
        std::map<std::tuple<int, int, int>, int> perChunk;
        for (int i = S.count0; i < Animals_Count(); ++i) {
            if (Animals_Body(i) != c.body->name || (Animals_SpeciesOf(i).first != "hornhide" && Animals_SpeciesOf(i).first != "creature_nest")) continue; // (the cave dwellers and nests only)
            const Vector3 p = Animals_Position(i);
            const Vector3 dir = Vector3Normalize(p);
            const float depth = c.field->radius + c.field->elevation(dir) - Vector3Length(p);
            ++n;
            under += depth >= 5.5f;
            const bool isOpen = OpenAt(c, p);
            open += isOpen;
            if (!isOpen) c.note(TextFormat("not in open space after 5 s: %s, %.1f m deep, spawned %s", Animals_SpeciesOf(i).first.c_str(), depth,
                                           S.spawnedIn.count(Animals_NetId(i)) && S.spawnedIn[Animals_NetId(i)].open ? "in open space" : "in rock"));
            far += Vector3Distance(p, me) >= 23.5f;
            ++perChunk[{ (int)floorf(p.x / 64.0f), (int)floorf(p.y / 64.0f), (int)floorf(p.z / 64.0f) }];
        }
        c.metric("new animals", n, "", false);
        c.eq("all underground (6 m+)", under, n);
        c.eq("all in open space as spawned", OpenAsSpawned(), (int)S.spawnedIn.size());
        c.metric("in open space after 5 s of wandering", open, "", false);
        c.eq("none within 24 m of the player", far, n);
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
        // The cap holds where the pass placed them (the first visit fills each chunk; respawn is 30 s, past this step).
        // 5 s later some have wandered over a chunk face (seed 42's #1883 start: 4 hornhides, two within 1.2 m of a
        // face), which the spawner counts as the next chunk's: reported, not checked.
        c.metric("most in one chunk after 5 s of wandering", worst, "", false);
        c.metric("cave dwellers and nests seen spawning", (double)S.spawnedIn.size(), "", false);
        c.lt("most in one chunk as spawned (cap 2 plus a nest)", FullestAsSpawned(), 3.5);
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
