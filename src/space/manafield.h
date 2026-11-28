#pragma once

// The mana field (#1248, magic step 1; the magic doc section 3): ley lines and nodes of the four base kinds, and a density
// query for the first slice's nine kinds. Pure and seeded: the lines and nodes come only from the world seed and the planet
// radius (and the seed's climate for a node's kind), so the server and every client draw the same field without messages.
//
// Lines: each base kind has its own seeded noise, and its ley lines are that noise's zero contours (long, curving, closed
// lines a few km apart), shown where the kind's regional strength allows. The strength rises with the kind's geology: Earth
// near the fault regions (and Earth also runs along the real fault lines, #590), Water in karst country (#592), Fire near
// volcanoes (#593); Air has its own regions and its flow rises with height, wind and storms at query time.
// Nodes: where lines of two different kinds cross. The algebra (mana.txt) names the node's kind: Fire x Air = Storm,
// Earth x Air = Shard, ...; a Fire line crossing frozen ground cancels into Aether.
//
// Everything magic is behind RAYLIBGAME_MAGIC until the first playable slice (magic.h).

#include "mana.h"

#include <raylib.h>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <vector>

enum ManaBase : int { kManaEarth = 0, kManaWater, kManaFire, kManaAir, kManaBases };

struct ManaFieldParams {
    unsigned int seed = 0;   // the planet's terrain seed (PlanetField::fieldSeed)
    float radius = 47000.0f; // sea level (m)
    bool volcanoes = true, karst = true; // the world's generation switches (WorldGenParams)
    // The annual mean surface temperature (C) at a direction, seed-pure (PlanetField::climateAt). Optional: without it no
    // node is frozen ground, so no Fire line makes an Aether node.
    std::function<float(Vector3 dir)> meanTempC;
};

struct ManaNode {
    Vector3 dir{};      // unit, body-local
    int a = -1, b = -1; // the crossing lines' base kinds (ManaBase)
    int kind = -1;      // its kind (mana.txt index): the compound of a and b, or Aether; -1 if the table lacks it
    float strength = 0; // 0..1: the weaker line's strength there
    uint32_t id = 0;    // stable for the seed
};

// What the place adds at query time (the game fills it from the planet and the weather; tests set it by hand).
struct ManaEnv {
    float elevation = 0.0f; // m above sea level
    float depth = 0.0f;     // m below the surface (0 in the open)
    float tempC = 15.0f;    // the air temperature now
    float night = 0.0f;     // 0 day .. 1 deep night
    float wind = 0.0f;      // m/s
    float storm = 0.0f;     // 0..1: a thunderstorm here
    int biome = -1;         // Biome (voxel.h), -1 unknown
    bool water = false;     // in or beside open water
};

struct ManaSample {
    float v[kManaMaxKinds] = {}; // flow per kind (mana.txt index), 0..2; only the slice's kinds are filled
    float of(int kind) const { return kind >= 0 && kind < kManaMaxKinds ? v[kind] : 0.0f; }
};

struct ManaSegment { Vector3 a, b; float strength; int base; }; // a ley line piece: unit directions (body-local)

class ManaField {
public:
    explicit ManaField(const ManaFieldParams& p);

    // A base kind's ley lines at a unit direction.
    float strength(int base, Vector3 dir) const;     // 0..1: the regional strength (seed-pure)
    float lineDistance(int base, Vector3 dir) const; // m to its nearest line (an estimate: |value| / |gradient|)
    float lineFlow(int base, Vector3 dir) const;     // strength x the line profile (1 on a full-strength line)

    // Nodes within `range` metres (along the surface) of a direction; built per 4 km tile on first use and cached.
    void nodesNear(Vector3 dir, float range, std::vector<ManaNode>& out) const;

    // The flow of the slice's kinds at a direction (the magic doc 3: lines, nodes, and what the place adds).
    ManaSample density(Vector3 dir, const ManaEnv& env) const;

    // The ley line pieces within `range` m of a direction (marching squares on a `step` m tangent grid), for the debug
    // overlay and magic sight. Pieces where the kind is weaker than minStrength are left out.
    void lineSegments(Vector3 dir, float range, float step, float minStrength, std::vector<ManaSegment>& out) const;

    const ManaFieldParams& params() const { return p; }

    static constexpr float kLineWidth = 30.0f;  // m: a line's core half-width (guess)
    static constexpr float kTileSize = 4000.0f; // m: node tiles
    static constexpr float kNodeMinStrength = 0.3f;
    static constexpr float kNodeReach = 60.0f;  // m: a node's glow in the density (guess)

private:
    struct Family { int base; int noise; }; // a set of lines: a noise's zero contours (noise -1: the fault field)
    float familyValue(int fam, Vector3 s) const; // the signed contour value at a surface point (m)
    float familyDistance(int fam, Vector3 dir) const;
    int tilesPerFace() const;
    const std::vector<ManaNode>& tile(int face, int i, int j) const;

    ManaFieldParams p;
    std::vector<Family> families;
    struct Noise; std::shared_ptr<Noise> noise;
    mutable std::mutex mu;
    mutable std::unordered_map<uint64_t, std::vector<ManaNode>> tiles;
    int kindOfPair[kManaBases][kManaBases];
    int kBase[kManaBases], kAether = -1, kFrost = -1, kStorm = -1, kShard = -1, kRot = -1; // mana.txt indices
};
