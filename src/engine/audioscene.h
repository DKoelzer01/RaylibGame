#ifndef AUDIOSCENE_H
#define AUDIOSCENE_H

// The world as sound sees it (#1293, audio step 3b): rays against the terrain's density field (with the player's edits)
// and the build plots' part cells, for Steam Audio's occlusion and transmission (a custom scene: its ray callbacks call
// these). No meshes are built and no Embree is needed. A snapshot is immutable: the body's field and edits are the
// shared read-only snapshots the chunk workers use, the plots the grass's part-cell snapshot (Grass_BlockersOf), so it
// is read on the simulation thread while the main thread goes on.
//
// "Audio space": global axes, metres, origin at `origin` (a global point near the listener, the snapshot's), in float.

#include "dvec3.h"
#include "foliage.h"     // GrassBlockers
#include "voxel.h"       // DensityField, VoxelMaterial

#include <memory>
#include <vector>

class Scene;

// Acoustic material classes (Steam Audio's reference values, steamaudio.cpp).
enum class AudioMat : uint8_t { Generic, Rock, Soil, Gravel, Snow, Wood, Metal, Glass, Brick, Count };
AudioMat AudioMatOfVoxel(VoxelMaterial m);

struct AudioScene {
    struct Body {
        std::shared_ptr<const DensityField> field;
        std::shared_ptr<const TerrainDeltas> edits;
        std::shared_ptr<const GrassBlockers> plots;   // may be null
        DVec3 b{};                                    // body-local position of the audio-space origin
        double A[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } }; // body-local = b + A * audio
        float localPerMetre = 1.0f;                   // body-local units in one metre (1 / the body's scale)
    };
    DVec3 origin{};                                   // global
    Vector3 listener{ 0, 0, 0 };                      // the listener in audio space (the ray ends near it get no margin)
    std::vector<Body> bodies;

    static constexpr float kSourceMargin = 0.75f;     // solid this close to a sound's end of a ray is ignored (a hum
                                                      // inside its machine, a source half in the ground)
    // The first solid along origin + dir t, t in [tMin, tMax] (dir unit, audio space): its t, or a negative number.
    // listenerAt: where the listener is now (audio space), if it moved since the snapshot (default: `listener`).
    float closestHit(Vector3 o, Vector3 dir, float tMin, float tMax, AudioMat* mat = nullptr, const Vector3* listenerAt = nullptr) const;
    bool anyHit(Vector3 o, Vector3 dir, float tMin, float tMax, const Vector3* listenerAt = nullptr) const { return closestHit(o, dir, tMin, tMax, nullptr, listenerAt) >= 0.0f; }

private:
    float terrainHit(const Body& b, Vector3 o, Vector3 dir, float lo, float hi, AudioMat* mat) const;
    float plotHit(const Body& b, const GrassBlockPlot& p, Vector3 o, Vector3 dir, float lo, float hi, AudioMat* mat) const;
};

// Main thread: the bodies whose ground can be within `range` metres of the listener, around it.
std::shared_ptr<const AudioScene> AudioScene_Snapshot(Scene& world, const DVec3& listenerGlobal, float range = 80.0f);

#endif
