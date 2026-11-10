#ifndef RING_H
#define RING_H

// Planetary rings (milestone 6, design doc section 10): tens of thousands of small asteroids on
// circular orbits in a thin band around a planet.
//
// Every rock is derived from the world seed and its index (position in the band, size, spin,
// shape seed), so untouched rocks cost nothing to store. Its orbit is analytic, like the moon's.
//   Far: rocks draw as GPU instances of a few dozen pre-made rock meshes; rocks too small to see
//        are skipped, and beyond ~40 km the ring is a translucent band, shaded by the planet's
//        shadow, that fades out as you approach the ring plane.
//   Near: within a reach that grows with its size a rock becomes a full voxel body (diggable,
//        walkable, with gravity), and turns back into an instance when you leave, keeping its
//        edits (and any build plot) for the next visit and for saves.
// Rocks never collide with each other (non-crossing circular orbits).

#include <raylib.h>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "dvec3.h"
#include "terraindeltas.h"

class Scene;
class Planetoid;
class PlotSet;
struct RenderContext;
struct BodyDef;

struct RingRock {
    double orbitRadius;   // from the parent's centre, in its equatorial plane
    double height;        // along its axis
    double phase;         // rad at world time 0
    double rate;          // rad/s
    float radius;         // the voxel body's radius (its mean surface is ~0.72 of it)
    float spin;           // deg/s about its own axis
    uint32_t seed;
    uint8_t variant;      // instance mesh
    Color tint;
};

// What a rock keeps while it is an instance (after a visit, or from a save).
struct RingRockState {
    std::shared_ptr<const TerrainDeltas> edits;
    std::unordered_set<TerrainDeltas::Key> dirtyRegions;
    bool allRegionsDirty = true;
    std::string regionsSavedTo;
    std::unique_ptr<PlotSet> plots;
};

class RingSystem {
public:
    RingSystem(const BodyDef& def, Planetoid& parent, uint32_t seed, int ringNo);
    ~RingSystem();

    int ringNo;
    std::string name;
    std::vector<RingRock> rocks;
    double innerRadius = 0, outerRadius = 0, thickness = 0;
    int largeCount = 0;

    // Main thread, once a frame after the bodies moved: promote rocks near `focus` (world space) to
    // voxel bodies, demote those left behind.
    void update(Scene& world, DVec3 focus);
    // Draw the instances (mid and near passes) and the band (far pass).
    void drawRocks(Scene& world, const RenderContext& ctx, bool rebuild); // rebuild: first pass of the frame
    void drawBand(Scene& world, const RenderContext& ctx);
    // Distance range of the band from the camera (for the far pass's clip planes); false if hidden.
    bool bandRange(Scene& world, double& dmin, double& dmax) const;

    DVec3 rockPosition(int index) const; // world space, now
    int promotedCount() const { return (int)promoted.size(); }
    Planetoid* promotedBody(int index) const; // the voxel body of a promoted rock, or null
    int nearestRock(DVec3 world, double* distance = nullptr) const;
    // A big body on a circular orbit in the ring sweeps a lane clear: rocks whose orbits it crosses are
    // left out (their indices stay, so saves keep their meaning). Returns how many.
    int clearLane(double orbitRadius, double height, double radius);
    int clearedCount = 0;
    int drawnInstances = 0;

    // Saves: rocks with edits or plots (promoted or not). Bodies use region ids 1000000 + ...
    void save(FILE* f, const std::string& regionDir);
    bool loadRock(FILE* f, int rock, size_t strokes, const std::string& regionDir); // after "ringrock R I strokes N"
    int regionBodyId(int rock) const { return 1000000 + ringNo * 100000 + rock; }

private:
    uint64_t parentId;
    Planetoid* parent;
    std::unordered_map<int, uint64_t> promoted;      // rock -> body id
    std::unordered_map<int, RingRockState> dormant;  // rocks with state, as instances
    // Angular bins (rebuilt every few seconds): rocks sorted by their angle then.
    std::vector<std::vector<int>> bins;
    std::vector<int> large;                          // large rocks: drawn from much further away
    double binnedAt = -1e9;
    void rebin();
    double angleOf(const RingRock& r) const;
    void promote(Scene& world, int rock);
    void demote(Scene& world, int rock, Planetoid* body);
    // GPU
    std::vector<unsigned int> vao, vbo;              // per mesh variant and detail level
    std::vector<int> instCount;
};

// Gravity of ring rocks (and hand-placed huge ones): low, growing slowly with size (the asteroid
// rule would give a 1 km rock nearly Earth gravity).
float RingRockGravity(float radius);

#endif
