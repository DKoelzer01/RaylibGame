// --test-ring (milestone 6 gate): the ring at scale and a round trip to one of its rocks.
//   1. the ring: rock count, size mix, lanes cleared by the huge rocks
//   2. from the ground: the band in the sky (bench_ring_ground.png), frame time
//   3. inside the ring next to a large rock: instances drawn, frame time (bench_ring_inside.png);
//      the rock becomes a voxel body
//   4. dig it, fly 20 km away (it turns back into an instance), come back: the hole is still there;
//      the same after a save and a load
//   5. stand on a huge rock (its own low gravity)

#include "ring_test.h"
#include "player.h"
#include "survival.h"
#include "ring.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

namespace {
    int phase = 0;
    double phaseStart = 0.0;
    int rock = -1;
    Vector3 digLocal{};
    double frameSum = 0.0, frameMax = 0.0;
    int frames = 0, instances = 0;
    bool ok = true;
    Planetoid* huge = nullptr;
    int groundedFrames = 0, standFrames = 0;
    DVec3 cameraOffset{};   // from the rock, world space

    Planetoid* Planet(Scene& w) { return dynamic_cast<Planetoid*>(w.objects[0].get()); }
    void Next() { ++phase; phaseStart = GetTime(); frameSum = frameMax = 0.0; frames = 0; }
    double T() { return GetTime() - phaseStart; }
    void Measure() { const double dt = GetFrameTime() * 1000.0; frameSum += dt; frameMax = std::max(frameMax, dt); ++frames; }
    void PlaceCamera(Scene& world, RingSystem& ring) {
        const DVec3 at = ring.rockPosition(rock);
        world.setCameraWorld(at + cameraOffset, at);
        player.position = world.renderOrigin + world.camera.position;
    }
    float DugDensity(RingSystem& ring) {
        Planetoid* b = ring.promotedBody(rock);
        return b ? b->densityAt(digLocal) : -99.0f;
    }
}

int RingTest_Update(Scene& world) {
    if (Rings().empty()) { printf("bench: --test-ring needs a ring (the default bodies.txt has one)\n"); return -1; }
    RingSystem& ring = *Rings()[0];
    Planetoid* planet = Planet(world);
    player.flyMode = true;
    player.flyFollow = false;
    switch (phase) {
    case 0: {
        int small = 0, big = 0;
        float maxR = 0.0f;
        for (const RingRock& r : ring.rocks) { if (r.radius <= 0.0f) continue; (r.radius * 0.72f > 30.0f ? big : small)++; maxR = std::max(maxR, r.radius * 0.72f); }
        printf("bench: ring: %zu rocks (%d small, %d large up to %.0f m across, %d cleared from lanes) in a %.1f km band, %.0f m thick; %zu KB of rock data\n",
               ring.rocks.size(), small, big, 2.0f * maxR, ring.clearedCount, (ring.outerRadius - ring.innerRadius) / 1000.0, ring.thickness, ring.rocks.size() * sizeof(RingRock) / 1024);
        ok &= ring.rocks.size() >= 50000;
        // Look at the band from the ground: the nearest point of the ring's centre line.
        const DVec3 feet = player.position;
        const DVec3 rel = feet - planet->worldPos;
        const double a = atan2(rel.z, rel.x), rr = (ring.innerRadius + ring.outerRadius) * 0.5;
        const DVec3 target = planet->worldPos + DVec3(cos(a) * rr, 0.0, sin(a) * rr);
        const DVec3 up = rel * (1.0 / rel.length());
        world.setCameraWorld(feet + up * 2.0, target);
        Next();
        break;
    }
    case 1:
        if (T() > 2.0) Measure();
        if (T() > 4.0) {
            TakeScreenshot("bench_ring_ground.png");
            printf("bench: ring from the ground: %.2f ms avg, %.2f ms worst over %d frames\n", frameSum / std::max(1, frames), frameMax, frames);
            // A large rock to visit, 1 km out from it (radially, so the band is in view).
            int best = -1;
            float bestR = 0.0f;
            for (size_t k = 0; k < ring.rocks.size(); ++k) if (ring.rocks[k].radius > bestR && ring.rocks[k].radius * 0.72f < 90.0f) { bestR = ring.rocks[k].radius; best = (int)k; }
            rock = best;
            const DVec3 at = ring.rockPosition(rock);
            const DVec3 out = (at - planet->worldPos) * (1.0 / (at - planet->worldPos).length());
            cameraOffset = out * 900.0 + DVec3(0.0, 120.0, 0.0);
            PlaceCamera(world, ring);
            Next();
        }
        break;
    case 2:
        PlaceCamera(world, ring);
        if (T() > 3.0 && WorldIdle()) { Measure(); instances = std::max(instances, ring.drawnInstances); }
        if (T() > 6.0 && frames > 30) {
            TakeScreenshot("bench_ring_inside.png");
            Planetoid* b = ring.promotedBody(rock);
            printf("bench: ring inside: %.2f ms avg, %.2f ms worst over %d frames; %d instances drawn, %d rocks are voxel bodies; rock %d (%.0f m across) promoted: %s\n",
                   frameSum / std::max(1, frames), frameMax, frames, instances, ring.promotedCount(), rock, ring.rocks[rock].radius * 1.44f, b ? "yes" : "NO");
            ok &= b != nullptr && frameSum / std::max(1, frames) < 8.0;
            if (!b) return -1;
            // Dig into its surface facing the camera.
            const Vector3 toCam = Vector3Normalize(b->toLocal(world.camera.position));
            float rS = b->radius;
            while (rS > 1.0f && b->densityAt(Vector3Scale(toCam, rS)) < ISO_LEVEL) rS -= 0.5f;
            digLocal = Vector3Scale(toCam, rS - 1.0f);
            for (int k = 0; k < 6; ++k) ApplyTerrainEdit(*b, TerrainEdit{ Vector3Scale(toCam, rS - 0.5f * k), 2.5f, -3.0f });
            printf("bench: ring: dug into rock %d at %.1f m from its centre: density there %.2f (solid >= %.2f)\n", rock, rS, DugDensity(ring), ISO_LEVEL);
            ok &= DugDensity(ring) < ISO_LEVEL;
            Next();
        }
        break;
    case 3: // leave
        cameraOffset = DVec3(0.0, 20000.0, 0.0);
        PlaceCamera(world, ring);
        if (T() > 1.5) {
            const bool gone = ring.promotedBody(rock) == nullptr;
            printf("bench: ring: 20 km away the rock is an instance again: %s (%d voxel rocks)\n", gone ? "yes" : "NO", ring.promotedCount());
            ok &= gone;
            Next();
        }
        break;
    case 4: // return
        cameraOffset = DVec3(0.0, 120.0, 0.0) + (ring.rockPosition(rock) - planet->worldPos) * (900.0 / (ring.rockPosition(rock) - planet->worldPos).length());
        PlaceCamera(world, ring);
        if (T() > 1.0) {
            const float d = DugDensity(ring);
            printf("bench: ring: back at the rock: voxel body %s, the hole %s (density %.2f)\n", ring.promotedBody(rock) ? "yes" : "NO", d < ISO_LEVEL && d > -90.0f ? "kept" : "LOST", d);
            ok &= d < ISO_LEVEL && d > -90.0f;
            const bool saved = SaveWorld(world), loaded = saved && LoadWorldSave(world);
            printf("bench: ring: save %s, load %s\n", saved ? "ok" : "FAILED", loaded ? "ok" : "FAILED");
            ok &= loaded;
            Next();
        }
        break;
    case 5: { // after the load (the world was rebuilt): back to the rock
        RingSystem& r2 = *Rings()[0];
        PlaceCamera(world, r2);
        if (T() > 1.0) {
            const float d = DugDensity(r2);
            printf("bench: ring: after the load the hole %s (density %.2f)\n", d < ISO_LEVEL && d > -90.0f ? "is still there" : "is GONE", d);
            ok &= d < ISO_LEVEL && d > -90.0f;
            for (auto& obj : world.objects)
                if (auto* p = dynamic_cast<Planetoid*>(obj.get()); p && p->name == "bastion") huge = p;
            if (!huge) { printf("bench: ring: no huge rock 'bastion'\n"); return -1; }
            player.flyMode = false;
            player.spawnOn(world, *huge, huge->dirToWorld(Vector3Normalize({ 0.2f, 1.0f, 0.3f })));
            Next();
        }
        break;
    }
    case 6:
        player.flyMode = false;
        if (T() > 1.5) { ++standFrames; if (player.grounded) ++groundedFrames; }
        if (T() > 4.5) {
            const Vector3 feet = (player.position - world.renderOrigin).toFloat();
            printf("bench: ring: standing on %s (%.0f m across): grounded %d%% of 3 s, gravity %.2f m/s^2, %.0f m from its centre\n",
                   huge->name.c_str(), huge->radius * 1.44f, 100 * groundedFrames / std::max(1, standFrames), player.gravity, Vector3Distance(feet, huge->position));
            ok &= groundedFrames > standFrames * 9 / 10 && player.gravity > 0.5f && player.gravity < 3.0f && survival.deaths == 0;
            printf("bench: ring gate: %s\n", ok ? "PASS" : "FAIL");
            return ok ? 1 : -1;
        }
        break;
    }
    return 0;
}
