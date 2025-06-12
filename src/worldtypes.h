#ifndef worldtypes_h
#define worldtypes_h

#include <tuple>
#include <unordered_map>
#include <vector>
#include "object.h"
#include "SimplexNoise.h"
#include "edgecache.h" // EdgeCacheEntry is now defined here


// Forward declaration for ChunkObject
class ChunkObject;

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


struct Chunk {
    std::vector<float> noiseValues;
    std::vector<Vector3> vertices;
    std::vector<int> indices;
    Int3 position; // Position of the chunk in the world
    Mesh mesh;     // Store mesh for chunk lifetime
    Model model;   // Store model for chunk lifetime
    // Add more as needed (e.g., mesh, cache)
};

class Planetoid: public Object {
public:
    int size;
    Vector3 seed;
    std::unordered_map<std::tuple<int, int, int, int>, std::vector<EdgeCacheEntry>, Tuple4Hash> sharedEdgeCaches; // Shared edge caches for chunks
    std::unordered_map<Int3,bool> generatedChunks; // Store generated chunk positions
    std::unordered_map<Int3, std::unique_ptr<ChunkObject>> chunkChildren;
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