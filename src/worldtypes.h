#ifndef worldtypes_h
#define worldtypes_h

#include <tuple>
#include <unordered_map>
#include <vector>
#include "object.h"
#include "SimplexNoise.h"
#include <raylib.h>


// Hash for tuple<int, int, int, int>
struct Tuple4Hash {
    std::size_t operator()(const std::tuple<int, int, int, int>& t) const {
        std::size_t h1 = std::hash<int>()(std::get<0>(t));
        std::size_t h2 = std::hash<int>()(std::get<1>(t));
        std::size_t h3 = std::hash<int>()(std::get<2>(t));
        std::size_t h4 = std::hash<int>()(std::get<3>(t));
        // Combine hashes
        return (((h1 ^ (h2 << 1)) >> 1) ^ (h3 << 1)) ^ (h4 << 2);
    }
};

struct Int3 {
    int x, y, z;
    bool operator==(const Int3& other) const { return x == other.x && y == other.y && z == other.z; }
};

namespace std {
    template<>
    struct hash<Int3> {
        size_t operator()(const Int3& k) const {
            return ((hash<int>()(k.x) ^ (hash<int>()(k.y) << 1)) >> 1) ^ (hash<int>()(k.z) << 2);
        }
    };
}

struct Vector3Hash {
    std::size_t operator()(const Vector3& v) const {
        std::size_t hx = std::hash<float>()(v.x);
        std::size_t hy = std::hash<float>()(v.y);
        std::size_t hz = std::hash<float>()(v.z);
        return ((hx ^ (hy << 1)) >> 1) ^ (hz << 1);
    }
};


class Chunk : public Object {
public:
    std::vector<float> noiseValues;   // (CHUNK_SIZE+1)^3 density samples, chunk-local
    std::vector<Vector3> vertices;    // CPU copy of mesh positions, chunk-local
    std::vector<int> indices;
    Mesh mesh{};                      // CPU arrays allocated with MemAlloc so UnloadMesh can free them

    Chunk* neighbors[26] = {nullptr};
    bool normalsPending = false;
    uint32_t neighborMask = 0;

    Chunk() : Object("chunk", "chunk", {0, 0, 0}, {0, 0, 0}, WHITE, 1.0f), offset({0, 0, 0}) {}
    virtual ~Chunk();

    Int3 offset;   // Chunk offset inside the planetoid, in voxels (multiple of CHUNK_SIZE)

    void setNeighbor(int idx, Chunk* neighbor) {
        neighbors[idx] = neighbor;
        if (neighbor) neighborMask |= (1u << idx);
        else neighborMask &= ~(1u << idx);
    }

    bool allNeighborsPresent() const { return neighborMask == 0x3FFFFFF; } // 26 bits set

    void onNeighborAdded(int idx, Chunk* neighbor) {
        setNeighbor(idx, neighbor);
        if (normalsPending && allNeighborsPresent()) {
            calculateNormals();
            normalsPending = false;
        }
    }

    void tryCalculateNormals();
    // Computes normals and uploads (or updates) the mesh on the GPU.
    void calculateNormals();

    static const int neighborOffsets[26][3];
    void assignNeighborsAndNotify(std::unordered_map<Int3, std::unique_ptr<Chunk>>& chunkChildren);

    // parentTransform: the owning planetoid's world transform.
    void drawChunk(const RenderContext& ctx, const Matrix& parentTransform);
    void drawChunkDepth(const RenderContext& ctx, const Matrix& parentTransform);

private:
    bool visible(const RenderContext& ctx, const Matrix& model) const;
    Matrix modelMatrix(const Matrix& parentTransform) const;
};

class Planetoid: public Object {
public:
    int size;
    Vector3 seed;
    std::unordered_map<Int3,bool> generatedChunks; // Store generated chunk positions
    std::unordered_map<Int3, std::unique_ptr<Chunk>> chunkChildren;
    Planetoid(std::string name, Vector3 position, Vector3 rotation, Color color, float scale, size_t size);
    virtual ~Planetoid();
    void draw(const RenderContext& ctx) override;
    void drawDepthOnly(const RenderContext& ctx) override;
};

extern const int CHUNK_SIZE;

#endif
