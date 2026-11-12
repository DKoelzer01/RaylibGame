#ifndef FLUID_H
#define FLUID_H

// Flowing water (Plan 4): a sparse 3D cellular simulation, Terraria-style, with 1 m cells.
//
// Static water (the sea, rivers, lakes, ponds: PlanetField::waterElevation) stays the default and is
// drawn by the terrain chunks. The simulation starts where a terrain edit (or a bucket, a pump) opens
// space next to water, and keeps every cell it has touched afterwards.
//
// Sites: each disturbed area gets a site, a grid in a tangent frame (east, north, up) anchored at the
// first disturbance; over a site's 512 m the planet curves by ~0.5 m, so "down" is the frame's -up.
// A site holds bricks of 16^3 cells allocated where water goes. A cell is solid (terrain density at
// its centre, edits included) or open, and open cells hold a fill of 0-255 (255 = 1 m^3).
//
// Sources: an open cell inside static water is a source. The sea, rivers and lakes are infinite:
// their cells are topped up to their level every step. Ponds are finite: the first time a site
// reaches one, the pond is taken over (its static surface disappears, PlanetField::claimPond) and
// its water becomes ordinary cells that can drain away.
//
// Step (3 a tick at 20 Hz, on a worker): bottom to top, water falls into the cell below, then what is left evens
// out with the four side neighbours (a fifth of the difference each; films under ~3 cm stay put). Twice
// a second a levelling pass makes each connected body of water communicate (U-bends, a flooded cave
// rising to its source's level), easing toward the level a little at a time. Bricks with no change
// for 2 s sleep; edits, buckets and moving neighbours wake them.
//
// Drawing: marching cubes on the fill field of each changed brick (on the worker), the surface
// carried into solid cells so it meets the terrain without a gap; cells of static water are left to
// the chunks' own surface. Uploaded and drawn on the main thread with the water material.

#include <raylib.h>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

class Scene;
class Planetoid;
class ChunkBatch;
struct RenderContext;

namespace Fluid {
    // A new world or a load: forget everything (waits for a running step).
    void Reset();
    // Main thread, once a frame: hands queued commands to the next step, starts a step job when one is
    // due and the last has finished, uploads finished brick meshes.
    void Update(Scene& world, float dt);
    // A terrain edit on a body (body-local centre, radius): wakes water there, refreshes solidity, and
    // starts a site if the edit opens space next to water.
    void OnTerrainEdit(Planetoid& body, Vector3 localCentre, float radius);
    // In Planetoid::drawTransparent, with the water material bound.
    void Draw(const Planetoid& body, const RenderContext& ctx, ChunkBatch& batch);

    // Queries (main thread, from the last published step).
    bool Wet(const Planetoid& body, Vector3 local);                       // under moving water
    float SurfaceAbove(const Planetoid& body, Vector3 local, float reach); // height of the water surface over this point (m, < 0: none within reach)
    bool StandingWaterNear(const Planetoid& body, Vector3 local, float range); // irrigation
    float WaterAt(const Planetoid& body, Vector3 local, float radius);    // m^3 of water (moving or static) within radius

    // Fluid kinds (resources/data/fluids.txt; Plan 5): kind 0 is water, the rest by the table (crude oil,
    // later lava...). Every cell holds one kind; kinds don't mix (a fluid flows only into empty cells or
    // its own kind, so oil rests on water). Up to 8 kinds. Per kind: how readily it spreads sideways
    // (viscosity), the thinnest film that still spreads, its colour, opacity and glow, and what a bucket
    // scooping it yields.
    struct Kind {
        std::string id = "water", name = "Water";
        Color color = { 40, 90, 160, 255 };
        float opacity = 0.6f, glow = 0.0f;   // glow: emissive strength (lava)
        int viscosity = 5;                   // sideways flow per step = difference / viscosity (water 5)
        int minSpread = 10;                  // films thinner than this (/255 m) stay put
        std::string bucketItem;              // what a bucket scooping 1 m^3 yields ("" = the water bucket)
        int bucketCount = 0;
    };
    void LoadKinds(const char* path);       // at startup; water alone if missing
    int KindCount();
    const Kind& KindInfo(int kind);
    int KindByName(const std::string& id);   // -1 if unknown

    // Moving fluid by hand or machine (applied at the next step). Static water counts too: taking from
    // the sea or a river is free (they are infinite); taking from a pond takes it over first.
    float Take(Planetoid& body, Vector3 local, float m3, float radius = 1.5f, int kind = 0); // returns what will be taken
    float Pour(Planetoid& body, Vector3 local, float m3, int kind = 0);                     // returns what will be poured
    float FluidAt(const Planetoid& body, Vector3 local, float radius, int kind);          // m^3 of moving fluid of a kind within radius
    int KindAt(const Planetoid& body, Vector3 local, float radius);                      // the most plentiful non-water kind there, or -1

    // Sources (springs, seeps): add fluid at up to `rate` m^3/s while the fluid surface over the source
    // is below its head (`head` m above the source point). A basin fills to the head and the source
    // rests (equilibrium); where the rim is lower than the head the overflow runs off as a stream.
    // Saved with the site. A second source on the same cell replaces the first.
    void AddSource(Planetoid& body, Vector3 local, int kind, float rate, float head);
    int SourceCount();

    // Saves: per body, after the factory block.
    void Save(FILE* f, uint64_t bodyId);
    void Load(FILE* f, uint64_t bodyId);

    struct Stats {
        int sites = 0, bricks = 0, activeBricks = 0, meshes = 0, claimedPonds = 0;
        double wetCells = 0.0, volume = 0.0;  // cells with water, m^3 (non-source)
        double otherVolume = 0.0;             // m^3 of other fluids (oil, ...)
        double stepMs = 0.0, maxStepMs = 0.0, meshMs = 0.0;
        long long steps = 0;
        long long dbg[4] = {};     // tests: cell changes in the last tick by kind (down, side, source refill, levelling)
    };
    Stats GetStats();
    void WaitIdle();            // tests: until no step job is running
    void SetSpeed(int steps);   // tests: steps per tick (1 = 20 Hz real time)
    bool Asleep();              // every brick asleep and nothing queued
    std::string DebugCell(const Planetoid& body, Vector3 local); // tests: "fill N flags SRTD"
}

#endif
