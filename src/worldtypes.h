#ifndef worldtypes_h
#define worldtypes_h

#include <tuple>
#include <unordered_map>
#include <vector>

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
    std::vector<size_t> indices;
    Int3 position; // Position of the chunk in the world
    Mesh mesh;     // Store mesh for chunk lifetime
    Model model;   // Store model for chunk lifetime
    // Add more as needed (e.g., mesh, cache)
};

class Planetoid: public Object {
public:
    int size;
    std::unordered_map<std::tuple<int, int, int, int>, std::vector<int>, Tuple4Hash> sharedEdgeCaches; // Shared edge caches for chunks
    std::unordered_map<Int3,bool> generatedChunks; // Store generated chunk positions
    Planetoid(std::string name, Vector3 position, Vector3 rotation, Color color, float scale, size_t size);
    virtual ~Planetoid();
    void draw() override;
};

#endif