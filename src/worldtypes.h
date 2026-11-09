#ifndef worldtypes_h
#define worldtypes_h

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>
#include <raylib.h>

#include "object.h"
#include "voxel.h"

// One CHUNK_SIZE^3 block of a planetoid's surface mesh. Chunks are independent: their mesh is
// built from the planetoid's density field on a worker thread and only uploaded here.
class Chunk : public Object {
public:
    Mesh mesh{};      // GPU mesh; CPU arrays allocated with MemAlloc so UnloadMesh frees them
    Int3 coord{};     // chunk coordinate; voxel offset = coord * CHUNK_SIZE

    Chunk() : Object("chunk", "chunk", {0, 0, 0}, {0, 0, 0}, WHITE, 1.0f) {}
    virtual ~Chunk();

    // Replace the mesh with freshly built data (main thread: uploads to the GPU).
    void setMesh(const ChunkMeshData& data);

    // parentTransform: the owning planetoid's world transform.
    void drawChunk(const RenderContext& ctx, const Matrix& parentTransform);
    void drawChunkDepth(const RenderContext& ctx, const Matrix& parentTransform);

    Vector3 voxelOffset() const { return { (float)(coord.x * CHUNK_SIZE), (float)(coord.y * CHUNK_SIZE), (float)(coord.z * CHUNK_SIZE) }; }

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
    std::shared_ptr<const PlanetoidField> field;         // immutable, shared with workers
    std::shared_ptr<const EditList> edits;               // copy-on-write, shared with workers
    std::unordered_map<Int3, std::unique_ptr<Chunk>> chunks;      // only chunks with a mesh
    std::unordered_map<Int3, uint32_t> requestedVersion;          // latest mesh request per chunk
    bool initialRequestDone = false;

    Planetoid(std::string name, Vector3 position, Vector3 rotation, float scale, unsigned int seed, float radius);
    virtual ~Planetoid();

    void draw(const RenderContext& ctx) override;
    void drawDepthOnly(const RenderContext& ctx) override;

    // World <-> planetoid-local (voxel) space.
    Vector3 toLocal(Vector3 world) const;
    Vector3 toWorld(Vector3 local) const;

    float densityAt(Vector3 local) const; // includes edits; main-thread helper (raycasts)

    static Planetoid* find(uint64_t id); // nullptr once destroyed
};

#endif
