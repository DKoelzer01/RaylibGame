#ifndef worldtypes_h
#define worldtypes_h

#include <tuple>
#include <unordered_map>
#include <vector>
#include "object.h"
#include "SimplexNoise.h"
#include "edgecache.h" // EdgeCacheEntry is now defined here
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
    std::vector<float> noiseValues;
    std::vector<Vector3> vertices;
    std::vector<int> indices;
    Mesh mesh;
    Model model;
    Chunk* neighbors[26] = {nullptr};
    bool normalsPending = false;
    uint32_t neighborMask = 0;

    Chunk()
        : Object("chunk", "chunk", {0, 0, 0}, {0, 0, 0}, WHITE, 1.0f), position({0, 0, 0}) {
        mesh = { 0 };
        mesh.vertices = nullptr;
        mesh.indices = nullptr;
        mesh.normals = nullptr;
    }

    Chunk(const Int3& pos, Vector3 worldPos, Vector3 rotation, Color color, float scale)
        : Object("chunk", "chunk", worldPos, rotation, color, scale), position(pos) {}
    virtual ~Chunk();

    Int3 position; // Position in chunk grid

    // Assign a neighbor at a given index (0-25)
    void setNeighbor(int idx, Chunk* neighbor) {
        neighbors[idx] = neighbor;
        if (neighbor) neighborMask |= (1u << idx);
        else neighborMask &= ~(1u << idx);
    }

    // Check if all 26 neighbors are present
    bool allNeighborsPresent() const {
        return neighborMask == 0x3FFFFFF; // 26 bits set
    }

    // Called when a neighbor is added
    void onNeighborAdded(int idx, Chunk* neighbor) {
        setNeighbor(idx, neighbor);
        if (normalsPending && allNeighborsPresent()) {
            calculateNormals();
            normalsPending = false;
        }
    }

    // Attempt to calculate normals if ready
    void tryCalculateNormals();

    // Placeholder for your normal calculation logic
    void calculateNormals();

    // Helper: 26 neighbor offsets (faces, edges, corners)
    static const int neighborOffsets[26][3];

    // Method for Chunk: assign all 26 neighbors and notify them
    void assignNeighborsAndNotify(std::unordered_map<Int3, std::unique_ptr<Chunk>>& chunkChildren);

    void draw(Shader* lightingShader) override;
    void drawDepthOnly(const Matrix& lightSpaceMatrix, Shader* depthShader) override;
};

class Planetoid: public Object {
public:
    int size;
    Vector3 seed;
    std::unordered_map<std::tuple<int, int, int, int>, std::vector<EdgeCacheEntry>, Tuple4Hash> sharedEdgeCaches; // Shared edge caches for chunks
    std::unordered_map<Int3,bool> generatedChunks; // Store generated chunk positions
    std::unordered_map<Int3, std::unique_ptr<Chunk>> chunkChildren;
    Planetoid(std::string name, Vector3 position, Vector3 rotation, Color color, float scale, size_t size);
    virtual ~Planetoid();
    float GetNoise(float wx, float wy, float wz); // Get noise value at world coordinates
    void draw(Shader* lightingShader) override;
    void drawDepthOnly(const Matrix& lightSpaceMatrix, Shader* depthShader) override;
    // Update: Now takes chunkWorldPos (relative to planetoid) and origin (planetoid world position)
    Chunk generateChunk(const Vector3& chunkWorldPos, const Vector3& origin, SimplexNoise* noise);
};

#include "cubemarch.h" // Only include after Planetoid if needed

extern const int CHUNK_SIZE;

#endif