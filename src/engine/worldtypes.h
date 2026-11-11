#ifndef worldtypes_h
#define worldtypes_h

#include <cstdint>
#include <atomic>
#include <memory>
#include <unordered_map>
#include <vector>
#include <raylib.h>

#include "object.h"
#include "voxel.h"
#include "bodies.h"
#include "foliage.h"
#include "trees.h"
#include "farsphere.h"

// One octree LOD node of a planetoid's surface. Its mesh is built from the planetoid's density
// field on a worker thread and only uploaded here.
class Chunk : public Object {
public:
    Mesh mesh{};          // GPU mesh (empty if the node contains no surface)
    Mesh waterMesh{};     // sea surface in this node (empty if none)
    FoliageGpu grass{};   // grass tufts (finest level only)
    // The build version this node still wants (0 once freed). Queued jobs check it and skip work
    // nobody wants any more (flying fast used to fill the workers with builds for nodes long passed).
    std::shared_ptr<std::atomic<uint32_t>> wanted = std::make_shared<std::atomic<uint32_t>>(0);
    NodeKey key{};
    bool ready = false;   // a build result has arrived at least once
    // Body-local upper bound for fragments. Nodes overlap their +x/+y/+z neighbours by one cell to
    // hide cracks between LOD levels; where the neighbour on an axis is drawn at the same level the
    // overlap would just duplicate (and z-fight with) its surface, so it is clipped on that axis.
    Vector3 clipMax = { 1e30f, 1e30f, 1e30f };

    Chunk() : Object("chunk", "chunk", {0, 0, 0}, {0, 0, 0}, WHITE, 1.0f) {}
    virtual ~Chunk();

    // Replace the mesh with freshly built data (main thread: uploads to the GPU).
    void setMesh(const ChunkMeshData& data);
    bool hasMesh() const { return mesh.vaoId > 0; }

    // `model` = the chunk's full model matrix (Planetoid::chunkModel: translation in double).
    void drawChunk(const RenderContext& ctx, const Matrix& model);
    void drawChunkDepth(const RenderContext& ctx, const Matrix& model);
    void drawWater(const RenderContext& ctx, const Matrix& model);  // with ctx.batch bound to the water material
    void drawGrass(const RenderContext& ctx, const Matrix& model);  // grass shader bound (Planetoid::drawFoliage)

private:
    void releaseMesh();
    bool visible(const RenderContext& ctx, const Matrix& model) const;
};

class Planetoid : public Object {
public:
    const uint64_t id;                                   // stable handle for async results
    unsigned int seed;
    BodyKind kind = BodyKind::Asteroid;
    float surfaceGravity = 9.81f;                        // m/s^2 at gravityRadius
    float gravityRadius = 0.0f;                          // mean surface: inverse-square outside, linear inside
    float radius;
    Vector3 spin = { 0, 0, 0 };                          // degrees per second about X/Y/Z
    // Circular orbit around another body (bodies.txt orbit=): advanced analytically from the world
    // clock in advanceBodies(), in the parent's equatorial (XZ) plane. orbitParent is validated by id.
    Planetoid* orbitParent = nullptr;
    uint64_t orbitParentId = 0;
    double orbitRadius = 0.0;   // from the parent's centre (m)
    double orbitRate = 0.0;     // rad/s, from the parent's gravity: sqrt(GM / r^3)
    double orbitPhase = 0.0;    // rad at world time 0
    double orbitHeight = 0.0;   // offset along the parent's axis (ring members spread above/below the equator)
    Vector3 linearVelocity = { 0, 0, 0 }; // m/s of the body itself (orbit), part of its surface velocity
    std::shared_ptr<const DensityField> field;           // immutable, shared with workers
    // Flat (open-world) bodies: +Y up, constant gravity, LOD tiled around the camera within the
    // height band [minHeight, maxHeight] instead of a sphere of `radius`.
    bool flat = false;
    float minHeight = 0.0f, maxHeight = 0.0f;
    std::shared_ptr<const TerrainDeltas> edits;          // player edits as density deltas; immutable snapshot, shared with workers
    std::unordered_set<TerrainDeltas::Key> dirtyRegions; // regions changed since the last save (incremental saves)
    bool allRegionsDirty = true;                         // a fresh body: the first save writes every region
    std::string regionsSavedTo;                          // the region directory it was last saved to / loaded from

    // LOD state (main thread only).
    std::unordered_map<NodeKey, std::unique_ptr<Chunk>> nodes;   // built or being built
    std::unordered_map<NodeKey, uint32_t> requestedVersion;      // latest request per node
    std::vector<Chunk*> drawList;                                // what is rendered this frame
    int rootLevel = 0;                                           // level whose 2x2x2 nodes cover the body
    // UpdateLod skips work while nothing changed: the camera (in body space) stayed within 1 m of
    // where the node selection was last made and no chunk results or edits arrived.
    bool lodDirty = true;
    Vector3 lodCamLocal = { 1e30f, 1e30f, 1e30f };

    Planetoid(std::string name, DVec3 worldPos, Vector3 rotation, float scale, unsigned int seed, float radius);
    // Open-world terrain body.
    Planetoid(std::string name, DVec3 worldPos, std::shared_ptr<const TerrainField> terrain);
    virtual ~Planetoid();

    void draw(const RenderContext& ctx) override;
    void drawDepthOnly(const RenderContext& ctx) override;
    void drawTransparent(const RenderContext& ctx) override; // sea surface
    void drawFoliage(const RenderContext& ctx);             // grass (ctx.grassShader bound)

    // World (render-relative) <-> planetoid-local (voxel) space.
    Vector3 toLocal(Vector3 world) const;
    Vector3 toWorld(Vector3 local) const;
    // Double-precision versions in absolute world space (not render space): exact far from the
    // body centre, for anchors and chunk placement on large bodies (floats step ~4-8 mm at 60 km).
    DVec3 toWorldD(const DVec3& local) const;
    DVec3 toLocalD(const DVec3& worldAbs) const;
    // Model matrix of one of this body's chunks in render space, its translation computed in double.
    Matrix chunkModel(const Chunk& chunk, const DVec3& renderOrigin) const;
    // Model matrix for geometry positioned relative to a body-local anchor (foliage batches).
    Matrix localModel(const DVec3& anchorLocal, const DVec3& renderOrigin) const;
    std::unique_ptr<TreeSystem> trees;                   // planets with life: trees and bushes
    std::unique_ptr<FarSphere> farSphere;                // planets and moons: drawn as a sphere from far away
    bool useFarSphere(const RenderContext& ctx) const;  // far enough (zoom included) and baked

    float densityAt(Vector3 local) const; // includes edits; main-thread helper (raycasts)
    // Velocity (world units/s) of the body-fixed point currently at renderPoint, from the spin.
    Vector3 surfaceVelocity(Vector3 renderPoint) const;
    // Body-local direction <-> world direction (rotation only).
    Vector3 dirToLocal(Vector3 worldDir) const;
    Vector3 dirToWorld(Vector3 localDir) const;
    Quaternion rotationQuat() const; // body orientation (local -> world), for riding along with spin

    static Planetoid* find(uint64_t id); // nullptr once destroyed
};

#endif
