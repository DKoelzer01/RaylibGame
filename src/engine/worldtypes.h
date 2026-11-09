#ifndef worldtypes_h
#define worldtypes_h

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>
#include <raylib.h>

#include "object.h"
#include "voxel.h"

// One octree LOD node of a planetoid's surface. Its mesh is built from the planetoid's density
// field on a worker thread and only uploaded here.
class Chunk : public Object {
public:
    Mesh mesh{};          // GPU mesh (empty if the node contains no surface)
    NodeKey key{};
    bool ready = false;   // a build result has arrived at least once

    Chunk() : Object("chunk", "chunk", {0, 0, 0}, {0, 0, 0}, WHITE, 1.0f) {}
    virtual ~Chunk();

    // Replace the mesh with freshly built data (main thread: uploads to the GPU).
    void setMesh(const ChunkMeshData& data);
    bool hasMesh() const { return mesh.vaoId > 0; }

    void drawChunk(const RenderContext& ctx, const Matrix& parentTransform);
    void drawChunkDepth(const RenderContext& ctx, const Matrix& parentTransform);

private:
    void releaseMesh();
    bool visible(const RenderContext& ctx, const Matrix& model) const;
    Matrix modelMatrix(const Matrix& parentTransform) const;
};

class Planetoid : public Object {
public:
    const uint64_t id;                                   // stable handle for async results
    unsigned int seed;
    float radius;
    Vector3 spin = { 0, 0, 0 };                          // degrees per second about X/Y/Z
    std::shared_ptr<const PlanetoidField> field;         // immutable, shared with workers
    std::shared_ptr<const EditList> edits;               // copy-on-write, shared with workers

    // LOD state (main thread only).
    std::unordered_map<NodeKey, std::unique_ptr<Chunk>> nodes;   // built or being built
    std::unordered_map<NodeKey, uint32_t> requestedVersion;      // latest request per node
    std::vector<Chunk*> drawList;                                // what is rendered this frame
    int rootLevel = 0;                                           // level whose 2x2x2 nodes cover the body

    Planetoid(std::string name, DVec3 worldPos, Vector3 rotation, float scale, unsigned int seed, float radius);
    virtual ~Planetoid();

    void draw(const RenderContext& ctx) override;
    void drawDepthOnly(const RenderContext& ctx) override;

    // World (render-relative) <-> planetoid-local (voxel) space.
    Vector3 toLocal(Vector3 world) const;
    Vector3 toWorld(Vector3 local) const;

    float densityAt(Vector3 local) const; // includes edits; main-thread helper (raycasts)

    static Planetoid* find(uint64_t id); // nullptr once destroyed
};

#endif
