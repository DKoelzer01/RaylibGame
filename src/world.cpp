#include "world.h"
#include "utils/logger.h"


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

// Helper for 3D integer coordinates as chunk keys
namespace std {
    template<>
    struct hash<Int3> {
        size_t operator()(const Int3& k) const {
            return ((hash<int>()(k.x) ^ (hash<int>()(k.y) << 1)) >> 1) ^ (hash<int>()(k.z) << 2);
        }
    };
}

constexpr size_t CHUNK_SIZE = 32; // You can adjust this for performance/memory

void worldHandler(Scene& world) {
    // Initialize the camera for the world scene
    world.camera.position = { 0.0f, 10.0f, 10.0f };
    world.camera.target = { 0.0f, 0.0f, 0.0f };
    world.camera.up = { 0.0f, 1.0f, 0.0f };
    world.camera.fovy = 45.0f;
    world.camera.projection = CAMERA_PERSPECTIVE;

    generateWorld(world);
    loadWorld(world);
} 

void writeNoiseValuesToFile(const std::vector<float>& noiseValues, int size, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open file for writing noise values: " << filename << std::endl;
        return;
    }
    // Write size as header (optional, for reference)
    out << "Size: " << size << "\n";
    // Write noise values in a human-readable 3D format
    for (int z = 0; z < size; ++z) {
        out << "z = " << z << ":\n";
        for (int y = 0; y < size; ++y) {
            for (int x = 0; x < size; ++x) {
                int idx = x + y * size + z * size * size;
                out << noiseValues[idx];
                if (x < size - 1) out << ", ";
            }
            out << "\n";
        }
        out << "\n";
    }
    out.close();
}

void checkChunkBorderPositions(
    const Vector3& chunkOffsetA, const Vector3& chunkOffsetB,
    int size, char axis)
{
    // axis: 'x', 'y', or 'z'
    int mismatches = 0;
    float epsilon = 1e-6f;
    if (axis == 'x') {
        for (int z = 0; z < size; ++z) {
            for (int y = 0; y < size; ++y) {
                float wxA = chunkOffsetA.x + (size - 1);
                float wyA = chunkOffsetA.y + y;
                float wzA = chunkOffsetA.z + z;
                float wxB = chunkOffsetB.x + 0;
                float wyB = chunkOffsetB.y + y;
                float wzB = chunkOffsetB.z + z;
                if (fabs(wxA - wxB) > epsilon || fabs(wyA - wyB) > epsilon || fabs(wzA - wzB) > epsilon) {
                    logger.logf("X border position mismatch at (y=%d, z=%d): (%.8f,%.8f,%.8f) vs (%.8f,%.8f,%.8f)\n",
                        y, z, wxA, wyA, wzA, wxB, wyB, wzB);
                    ++mismatches;
                }
            }
        }
    } else if (axis == 'y') {
        for (int z = 0; z < size; ++z) {
            for (int x = 0; x < size; ++x) {
                float wxA = chunkOffsetA.x + x;
                float wyA = chunkOffsetA.y + (size - 1);
                float wzA = chunkOffsetA.z + z;
                float wxB = chunkOffsetB.x + x;
                float wyB = chunkOffsetB.y + 0;
                float wzB = chunkOffsetB.z + z;
                if (fabs(wxA - wxB) > epsilon || fabs(wyA - wyB) > epsilon || fabs(wzA - wzB) > epsilon) {
                    logger.logf("Y border position mismatch at (x=%d, z=%d): (%.8f,%.8f,%.8f) vs (%.8f,%.8f,%.8f)\n",
                        x, z, wxA, wyA, wzA, wxB, wyB, wzB);
                    ++mismatches;
                }
            }
        }
    } else if (axis == 'z') {
        for (int y = 0; y < size; ++y) {
            for (int x = 0; x < size; ++x) {
                float wxA = chunkOffsetA.x + x;
                float wyA = chunkOffsetA.y + y;
                float wzA = chunkOffsetA.z + (size - 1);
                float wxB = chunkOffsetB.x + x;
                float wyB = chunkOffsetB.y + y;
                float wzB = chunkOffsetB.z + 0;
                if (fabs(wxA - wxB) > epsilon || fabs(wyA - wyB) > epsilon || fabs(wzA - wzB) > epsilon) {
                    logger.logf("Z border position mismatch at (x=%d, y=%d): (%.8f,%.8f,%.8f) vs (%.8f,%.8f,%.8f)\n",
                        x, y, wxA, wyA, wzA, wxB, wyB, wzB);
                    ++mismatches;
                }
            }
        }
    }
    if (mismatches == 0) {
        logger.logf("No floating-point mismatches found on %c face!\n", axis);
    } else {
        logger.logf("Total %d floating-point mismatches found on %c face.\n", mismatches, axis);
    }
}

void checkChunkBorderNoise(
    const std::vector<float>& chunkA, const std::vector<float>& chunkB,
    int size, char axis)
{
    // axis: 'x', 'y', or 'z'
    int mismatches = 0;
    if (axis == 'x') {
        for (int z = 0; z < size; ++z) {
            for (int y = 0; y < size; ++y) {
                int idxA = (size - 1) + y * size + z * size * size;
                int idxB = 0 + y * size + z * size * size;
                if (fabs(chunkA[idxA] - chunkB[idxB]) > 1e-6f) {
                    logger.logf("X mismatch at (y=%d, z=%d): %.8f vs %.8f\n", y, z, chunkA[idxA], chunkB[idxB]);
                    ++mismatches;
                }
            }
        }
    } else if (axis == 'y') {
        for (int z = 0; z < size; ++z) {
            for (int x = 0; x < size; ++x) {
                int idxA = x + (size - 1) * size + z * size * size;
                int idxB = x + 0 * size + z * size * size;
                if (fabs(chunkA[idxA] - chunkB[idxB]) > 1e-6f) {
                    logger.logf("Y mismatch at (x=%d, z=%d): %.8f vs %.8f\n", x, z, chunkA[idxA], chunkB[idxB]);
                    ++mismatches;
                }
            }
        }
    } else if (axis == 'z') {
        for (int y = 0; y < size; ++y) {
            for (int x = 0; x < size; ++x) {
                int idxA = x + y * size + (size - 1) * size * size;
                int idxB = x + y * size + 0 * size * size;
                if (fabs(chunkA[idxA] - chunkB[idxB]) > 1e-6f) {
                    logger.logf("Z mismatch at (x=%d, y=%d): %.8f vs %.8f\n", x, y, chunkA[idxA], chunkB[idxB]);
                    ++mismatches;
                }
            }
        }
    }
    if (mismatches == 0) {
        logger.logf("No mismatches found on %c face!\n", axis);
    } else {
        logger.logf("Total %d mismatches found on %c face.\n", mismatches, axis);
    }
}

void weightNoise(std::vector<float>& noiseValues, int size, Vector3* chunkOffset, Vector3* center, float weight, float maxDist) {
    if (size <= 0) return;
    for (int z = 0; z < size; ++z) {
        for (int y = 0; y < size; ++y) {
            for (int x = 0; x < size; ++x) {
                int idx = x + y * size + z * size * size;
                float wx = x + chunkOffset->x;
                float wy = y + chunkOffset->y;
                float wz = z + chunkOffset->z;
                float dx = wx - center->x;
                float dy = wy - center->y;
                float dz = wz - center->z;
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                float scale = 2.0f - (dist / maxDist) * weight;
                if (scale < 0.0f) scale = 0.0f;
                noiseValues[idx] *= scale;
            }
        }
    }
}

void generateWorld(Scene& world) {
    logger.logf("Generating world scene...\n");
    clock_t worldGenStart = clock();

    //Delete old planetoids
    //TODO: Implement a proper cleanup of old planetoids

    float randScale = 100000.0f; // Scale for random number generation
    float genRange = 500.0f; // Range for planetoid generation
    float minSize = 50.0f; // Minimum size of planetoids
    float maxSize = 250.0f; // Maximum size of planetoids
    logger.logf("Generated %d planetoids.\n", generatePlanetoids(randScale, world, genRange, minSize, maxSize)); // Generate planetoids in the world

    clock_t worldGenEnd = clock();
    double worldGenTime = static_cast<double>(worldGenEnd - worldGenStart) / CLOCKS_PER_SEC;
    logger.logf("World generation completed in %.2f seconds.\n", worldGenTime);
}

int generatePlanetoids(float randScale, Scene& world, float genRange, float minSize, float maxSize) {
    int planetoidsGenerated = 0;

    //Setup RNG
    std::random_device rd; // Obtain a seed from the system
    std::mt19937 gen(rd()); // Initialize the Mersenne Twister engine with the seed
    std::uniform_int_distribution<> distrib(1, randScale);

    //Generate planetoids
    int numPlanetoids = (distrib(gen)/randScale) * 10; // Number of planetoids to generate
    //DEBUG
    numPlanetoids = 1; // For testing purposes, generate only one planetoid
    logger.logf("Generating %d planetoids...\n", numPlanetoids);

    for(int i = 0; i < numPlanetoids; ++i) {
        // Generate random position, rotation, color, size, and scale for the planetoid

        // Position is centered around (0,0,0) with a range of genRange
        Vector3 position = {    static_cast<float>(((distrib(gen)/randScale)*genRange)-(genRange/2)), 
                                static_cast<float>(((distrib(gen)/randScale)*genRange)-(genRange/2)), 
                                static_cast<float>(((distrib(gen)/randScale)*genRange)-(genRange/2))};

        position = { 0.0f, 0.0f, 0.0f }; // Fixed position for testing

        // Rotation is random in degrees, scaled by randScale
        Vector3 rotation = {    static_cast<float>(distrib(gen)/randScale * 360), 
                                static_cast<float>(distrib(gen)/randScale * 360), 
                                static_cast<float>(distrib(gen)/randScale * 360)};

        // Using a random color generator with a range of 0-255 for RGB values
        Color color = { static_cast<unsigned char>(distrib(gen)/randScale * 255), 
                        static_cast<unsigned char>(distrib(gen)/randScale * 255), 
                        static_cast<unsigned char>(distrib(gen)/randScale * 255), 255};

        // Size is a random float between minSize and maxSize
        size_t size = static_cast<size_t>((distrib(gen)/randScale)*maxSize) + static_cast<size_t>(minSize);
        size = 250; // Fixed size for testing

        // Scale is a fixed value for now, can be adjusted later
        float scale = 1.0f;

        logger.logf("Generating planetoid %d...\n", i);
        clock_t planetoidGenStart = clock();

        float frequencyNoise = static_cast<float>(distrib(gen)/randScale)-0.5f; // Random frequency noise to add some variation from -0.5 to 0.5

        SimplexNoise* noise = new SimplexNoise((1.5f+frequencyNoise)/static_cast<float>(size), 2.0f, 2.0f, 0.5f); // Create a new instance of SimplexNoise
        if (!noise) {
            std::cerr << "Failed to create SimplexNoise instance." << std::endl;
            continue; // Skip this planetoid if noise generation fails;
        }

        logger.logf("Generating noise for planetoid %d at position (%f, %f, %f) with size %zu and frequency %f\n", i, position.x, position.y, position.z, size, 1.5f+frequencyNoise);

        // Instead of one huge mesh, generate a grid of chunks
        std::unordered_map<Int3, Chunk> chunks;
        int chunkGrid = size/CHUNK_SIZE + (size % CHUNK_SIZE != 0 ? 1 : 0); // Calculate the number of chunks needed in each dimension
        logger.logf("Chunk grid size: %d x %d x %d\n", chunkGrid, chunkGrid, chunkGrid);

        // Center the chunk grid around the planetoid position
        Vector3 minCorner = {
            position.x - chunkGrid * CHUNK_SIZE * 0.5f,
            position.y - chunkGrid * CHUNK_SIZE * 0.5f,
            position.z - chunkGrid * CHUNK_SIZE * 0.5f
        };

        // --- Shared edge cache setup ---
        // Key: (minChunkX, minChunkY, minChunkZ, faceDir), faceDir: 0=X, 1=Y, 2=Z
        std::unordered_map<std::tuple<int, int, int, int>, std::vector<int>, Tuple4Hash> sharedEdgeCaches;

        for (int cz = 0; cz < chunkGrid; ++cz) {
            for (int cy = 0; cy < chunkGrid; ++cy) {
                for (int cx = 0; cx < chunkGrid; ++cx) {
                    logger.logf("Generating chunk at (%d, %d, %d)\n", cx, cy, cz);

                    Int3 chunkPos = {cx, cy, cz};
                    Chunk& chunk = chunks[chunkPos];
                    // Allocate noise for this chunk with a 1-voxel border
                    chunk.noiseValues.resize((CHUNK_SIZE + 1) * (CHUNK_SIZE + 1) * (CHUNK_SIZE + 1));
                    // Calculate chunk's world offset
                    Vector3 chunkOffset = {
                        minCorner.x + cx * CHUNK_SIZE,
                        minCorner.y + cy * CHUNK_SIZE,
                        minCorner.z + cz * CHUNK_SIZE
                    };
                    // Generate noise for this chunk (including border)
                    for (size_t z = 0; z <= CHUNK_SIZE; ++z) {
                        for (size_t y = 0; y <= CHUNK_SIZE; ++y) {
                            for (size_t x = 0; x <= CHUNK_SIZE; ++x) {
                                size_t idx = x + y * (CHUNK_SIZE + 1) + z * (CHUNK_SIZE + 1) * (CHUNK_SIZE + 1);
                                float wx = chunkOffset.x + x;
                                float wy = chunkOffset.y + y;
                                float wz = chunkOffset.z + z;
                                chunk.noiseValues[idx] = (noise->fractal(3, wx, wy, wz) + 1.0f);
                            }
                        }
                    }

                    float minNoise = chunk.noiseValues[0];
                    float maxNoise = chunk.noiseValues[0];
                    for (float v : chunk.noiseValues) {
                        if (v < minNoise) minNoise = v;
                        if (v > maxNoise) maxNoise = v;
                    }
                    logger.logf("Chunk (%d,%d,%d) noise min: %f, max: %f\n", cx, cy, cz, minNoise, maxNoise);

                    // Weight noise values for this chunk (use CHUNK_SIZE+1)
                    // Use chunkOffset for weighting
                    float planetoidRadius = chunkGrid * CHUNK_SIZE * 0.5f;
                    weightNoise(chunk.noiseValues, CHUNK_SIZE + 1, &chunkOffset, &position, 4.0f, planetoidRadius);

                    minNoise = chunk.noiseValues[0];
                    maxNoise = chunk.noiseValues[0];
                    for (float v : chunk.noiseValues) {
                        if (v < minNoise) minNoise = v;
                        if (v > maxNoise) maxNoise = v;
                    }
                    logger.logf("Chunk (%d,%d,%d) weighted noise min: %f, max: %f\n", cx, cy, cz, minNoise, maxNoise);

                    // --- Shared edge cache pointers for this chunk ---
                    // Each chunk needs up to 3 shared caches (for +X, +Y, +Z faces)
                    // Each shared cache is sized CHUNK_SIZE*CHUNK_SIZE*12 (for the face)
                    std::vector<int> localEdgeCache(CHUNK_SIZE * CHUNK_SIZE * CHUNK_SIZE * 12, -1);
                    std::fill(localEdgeCache.begin(), localEdgeCache.end(), -1);
                    std::vector<int>* edgeCacheX = nullptr;
                    std::vector<int>* edgeCacheY = nullptr;
                    std::vector<int>* edgeCacheZ = nullptr;
                    if (cx < chunkGrid - 1) {
                        auto key = std::make_tuple(cx + 1, cy, cz, 0); // Use neighbor's index!
                        edgeCacheX = &sharedEdgeCaches[key];
                        if (edgeCacheX->empty()) {
                            edgeCacheX->resize(CHUNK_SIZE * CHUNK_SIZE * 12, -1);
                        }
                    }
                    if (cy < chunkGrid - 1) {
                        auto key = std::make_tuple(cx, cy + 1, cz, 1);
                        edgeCacheY = &sharedEdgeCaches[key];
                        if (edgeCacheY->empty()) {
                            edgeCacheY->resize(CHUNK_SIZE * CHUNK_SIZE * 12, -1);
                        }
                    }
                    if (cz < chunkGrid - 1) {
                        auto key = std::make_tuple(cx, cy, cz + 1, 2);
                        edgeCacheZ = &sharedEdgeCaches[key];
                        if (edgeCacheZ->empty()) {
                            edgeCacheZ->resize(CHUNK_SIZE * CHUNK_SIZE * 12, -1);
                        }
                    }
                    
                    chunk.vertices.clear();
                    chunk.indices.clear();

                    logger.logf("Chunk (%d,%d,%d) marching cubes...\n", cx, cy, cz);
                    size_t vertsBefore = chunk.vertices.size();
                    size_t indsBefore = chunk.indices.size();

                    // Marching cubes: use CHUNK_SIZE+1 for noise, but mesh is CHUNK_SIZE-1
                    // For each cube, select the correct edge cache (local or shared)
                    for (size_t z = 0; z < CHUNK_SIZE; ++z) {
                        for (size_t y = 0; y < CHUNK_SIZE; ++y) {
                            for (size_t x = 0; x < CHUNK_SIZE; ++x) {
                                int cacheType = 0; // 0 = local, 1 = X, 2 = Y, 3 = Z
                                std::vector<int>* edgeCache = &localEdgeCache;

                                // Only use the shared cache for the +X, +Y, or +Z face (never for x==0, y==0, z==0)
                                if (x == CHUNK_SIZE - 1 && edgeCacheX) {
                                    edgeCache = edgeCacheX;
                                    cacheType = 1;
                                } else if (y == CHUNK_SIZE - 1 && edgeCacheY) {
                                    edgeCache = edgeCacheY;
                                    cacheType = 2;
                                } else if (z == CHUNK_SIZE - 1 && edgeCacheZ) {
                                    edgeCache = edgeCacheZ;
                                    cacheType = 3;
                                }

                                marchCube(
                                    x, y, z,
                                    chunk.noiseValues, CHUNK_SIZE + 1,
                                    chunk.vertices, chunk.indices,
                                    *edgeCache, 0.5f,
                                    cx, cy, cz, chunkGrid, cacheType
                                );
                            }
                        }
                    }
                    logger.logf("Chunk (%d,%d,%d) verts: %zu -> %zu, inds: %zu -> %zu", cx, cy, cz, vertsBefore, chunk.vertices.size(), indsBefore, chunk.indices.size());

                    // --- Create mesh and model for this chunk ---
                    if (chunk.vertices.empty() || chunk.indices.empty()) {
                        //logger.logf("Chunk at (%d, %d, %d) has no vertices or indices, skipping.\n", cx, cy, cz);
                        continue; // Skip this chunk if it has no vertices or indices
                    }
                    Mesh mesh = { 0 };
                    // Allocate mesh.vertices and mesh.indices
                    mesh.vertexCount = static_cast<int>(chunk.vertices.size());
                    mesh.vertices = new float[mesh.vertexCount * 3];
                    for (size_t j = 0; j < chunk.vertices.size(); ++j) {
                        const Vector3& v = chunk.vertices[j];
                        mesh.vertices[j * 3 + 0] = v.x;
                        mesh.vertices[j * 3 + 1] = v.y;
                        mesh.vertices[j * 3 + 2] = v.z;
                    }
                    mesh.triangleCount = static_cast<int>(chunk.indices.size() / 3);
                    mesh.indices = new unsigned short[chunk.indices.size()];
                    for (size_t j = 0; j < chunk.indices.size(); ++j) {
                        mesh.indices[j] = static_cast<unsigned short>(chunk.indices[j]);
                    }
                    UploadMesh(&mesh, false);
                    Model model = LoadModelFromMesh(mesh);



                    // Create a GameObject for this chunk
                    std::string chunkName = "chunk_" + std::to_string(cx) + "_" + std::to_string(cy) + "_" + std::to_string(cz);
                    Vector3 chunkCoords = chunkOffset;
                    Color chunkColor = color;
                    float chunkScale = scale;
                    Vector3 chunkRot = rotation;
                    auto chunkObj = std::make_unique<GameObject>("gameobject", chunkName, chunkCoords, chunkRot, chunkColor, chunkScale, model);
                    chunkObj->isActive = true;
                    chunkObj->parent = &world.rootObject;
                    world.rootObject.children.push_back(std::move(chunkObj));

                    // --- DEBUG: Compare shared edge cache indices at chunk borders ---
                    logger.logf("Comparing shared edge caches for chunk (%d,%d,%d)...\n", cx, cy, cz);
                    for (int cz = 0; cz < chunkGrid; ++cz) {
                        for (int cy = 0; cy < chunkGrid; ++cy) {
                            for (int cx = 0; cx < chunkGrid; ++cx) {
                                // +X neighbor
                                if (cx < chunkGrid - 1) {
                                    auto keyA = std::make_tuple(cx + 1, cy, cz, 0);
                                    auto keyB = std::make_tuple(cx + 1, cy, cz, 0);
                                    const auto& cacheA = sharedEdgeCaches[keyA];
                                    const auto& cacheB = sharedEdgeCaches[keyB]; // Both chunks should use the same cache object!
                                    logger.logf("Comparing shared edge cache for X face at (%d,%d,%d):\n", cx+1, cy, cz);
                                    for (size_t i = 0; i < cacheA.size(); ++i) {
                                        if (cacheA[i] != cacheB[i]) {
                                            logger.logf("  MISMATCH at i=%zu: cacheA=%d, cacheB=%d\n", i, cacheA[i], cacheB[i]);
                                        }
                                    }
                                }
                                // +Y neighbor
                                if (cy < chunkGrid - 1) {
                                    auto keyA = std::make_tuple(cx, cy + 1, cz, 1);
                                    auto keyB = std::make_tuple(cx, cy + 1, cz, 1);
                                    const auto& cacheA = sharedEdgeCaches[keyA];
                                    const auto& cacheB = sharedEdgeCaches[keyB];
                                    logger.logf("Comparing shared edge cache for Y face at (%d,%d,%d):\n", cx, cy+1, cz);
                                    for (size_t i = 0; i < cacheA.size(); ++i) {
                                        if (cacheA[i] != cacheB[i]) {
                                            logger.logf("  MISMATCH at i=%zu: cacheA=%d, cacheB=%d\n", i, cacheA[i], cacheB[i]);
                                        }
                                    }
                                }
                                // +Z neighbor
                                if (cz < chunkGrid - 1) {
                                    auto keyA = std::make_tuple(cx, cy, cz + 1, 2);
                                    auto keyB = std::make_tuple(cx, cy, cz + 1, 2);
                                    const auto& cacheA = sharedEdgeCaches[keyA];
                                    const auto& cacheB = sharedEdgeCaches[keyB];
                                    logger.logf("Comparing shared edge cache for Z face at (%d,%d,%d):\n", cx, cy, cz+1);
                                    for (size_t i = 0; i < cacheA.size(); ++i) {
                                        if (cacheA[i] != cacheB[i]) {
                                            logger.logf("  MISMATCH at i=%zu: cacheA=%d, cacheB=%d\n", i, cacheA[i], cacheB[i]);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // --- DEBUG: Print vertex positions at chunk boundaries ---
                    logger.logf("Vertex boundary check for chunk (%d,%d,%d):\n", cx, cy, cz);
                    if (cx < chunkGrid - 1) {
                        logger.logf("Vertex boundary check: chunk (%d,%d,%d) +X face\n", cx, cy, cz);
                        for (const Vector3& v : chunk.vertices) {
                            if (fabs(v.x - (chunkOffset.x + CHUNK_SIZE)) < 1e-4) {
                                logger.logf("v=(%.4f,%.4f,%.4f)\n", v.x, v.y, v.z);
                            }
                        }
                    }
                    if (cy < chunkGrid - 1) {
                        logger.logf("Vertex boundary check: chunk (%d,%d,%d) +Y face\n", cx, cy, cz);
                        for (const Vector3& v : chunk.vertices) {
                            if (fabs(v.y - (chunkOffset.y + CHUNK_SIZE)) < 1e-4) {
                                logger.logf("v=(%.4f,%.4f,%.4f)\n", v.x, v.y, v.z);
                            }
                        }
                    }
                    if (cz < chunkGrid - 1) {
                        logger.logf("Vertex boundary check: chunk (%d,%d,%d) +Z face\n", cx, cy, cz);
                        for (const Vector3& v : chunk.vertices) {
                            if (fabs(v.z - (chunkOffset.z + CHUNK_SIZE)) < 1e-4) {
                                logger.logf("v=(%.4f,%.4f,%.4f)\n", v.x, v.y, v.z);
                            }
                        }
                    }
                }
            }
        }

        // --- DEBUG: Compare world-space border vertices between adjacent chunks ---
        logger.logf("Comparing world-space border vertices between adjacent chunks...\n");
        for (int cz = 0; cz < chunkGrid; ++cz) {
            for (int cy = 0; cy < chunkGrid; ++cy) {
                for (int cx = 0; cx < chunkGrid; ++cx) {
                    Int3 chunkPos = {cx, cy, cz};
                    const auto& chunkA = chunks[chunkPos];
                    Vector3 chunkOffsetA = {
                        minCorner.x + cx * CHUNK_SIZE,
                        minCorner.y + cy * CHUNK_SIZE,
                        minCorner.z + cz * CHUNK_SIZE
                    };

                    // +X neighbor
                    if (cx < chunkGrid - 1) {
                        Int3 neighborX = {cx + 1, cy, cz};
                        const auto& chunkB = chunks[neighborX];
                        Vector3 chunkOffsetB = {
                            minCorner.x + (cx + 1) * CHUNK_SIZE,
                            minCorner.y + cy * CHUNK_SIZE,
                            minCorner.z + cz * CHUNK_SIZE
                        };
                        logger.logf("Comparing +X border vertices: chunk (%d,%d,%d) <-> (%d,%d,%d)\n", cx, cy, cz, cx+1, cy, cz);
                        for (const Vector3& vA : chunkA.vertices) {
                            if (fabs(vA.x - (chunkOffsetA.x + CHUNK_SIZE)) < 1e-4) {
                                // Look for a matching vertex in chunkB at x == chunkOffsetB.x
                                bool found = false;
                                for (const Vector3& vB : chunkB.vertices) {
                                    if (fabs(vB.x - chunkOffsetB.x) < 1e-4 &&
                                        fabs(vA.y - vB.y) < 1e-4 &&
                                        fabs(vA.z - vB.z) < 1e-4) {
                                        found = true;
                                        break;
                                    }
                                }
                                if (!found) {
                                    logger.logf("  MISMATCH: vA=(%.4f,%.4f,%.4f) has no match in neighbor\n", vA.x, vA.y, vA.z);
                                }
                            }
                        }
                    }

                    // +Y neighbor
                    if (cy < chunkGrid - 1) {
                        Int3 neighborY = {cx, cy + 1, cz};
                        const auto& chunkB = chunks[neighborY];
                        Vector3 chunkOffsetB = {
                            minCorner.x + cx * CHUNK_SIZE,
                            minCorner.y + (cy + 1) * CHUNK_SIZE,
                            minCorner.z + cz * CHUNK_SIZE
                        };
                        logger.logf("Comparing +Y border vertices: chunk (%d,%d,%d) <-> (%d,%d,%d)\n", cx, cy, cz, cx, cy+1, cz);
                        for (const Vector3& vA : chunkA.vertices) {
                            if (fabs(vA.y - (chunkOffsetA.y + CHUNK_SIZE)) < 1e-4) {
                                // Look for a matching vertex in chunkB at y == chunkOffsetB.y
                                bool found = false;
                                for (const Vector3& vB : chunkB.vertices) {
                                    if (fabs(vB.y - chunkOffsetB.y) < 1e-4 &&
                                        fabs(vA.x - vB.x) < 1e-4 &&
                                        fabs(vA.z - vB.z) < 1e-4) {
                                        found = true;
                                        break;
                                    }
                                }
                                if (!found) {
                                    logger.logf("  MISMATCH: vA=(%.4f,%.4f,%.4f) has no match in neighbor\n", vA.x, vA.y, vA.z);
                                }
                            }
                        }
                    }

                    // +Z neighbor
                    if (cz < chunkGrid - 1) {
                        Int3 neighborZ = {cx, cy, cz + 1};
                        const auto& chunkB = chunks[neighborZ];
                        Vector3 chunkOffsetB = {
                            minCorner.x + cx * CHUNK_SIZE,
                            minCorner.y + cy * CHUNK_SIZE,
                            minCorner.z + (cz + 1) * CHUNK_SIZE
                        };
                        logger.logf("Comparing +Z border vertices: chunk (%d,%d,%d) <-> (%d,%d,%d)\n", cx, cy, cz, cx, cy, cz+1);
                        for (const Vector3& vA : chunkA.vertices) {
                            if (fabs(vA.z - (chunkOffsetA.z + CHUNK_SIZE)) < 1e-4) {
                                // Look for a matching vertex in chunkB at z == chunkOffsetB.z
                                bool found = false;
                                for (const Vector3& vB : chunkB.vertices) {
                                    if (fabs(vB.z - chunkOffsetB.z) < 1e-4 &&
                                        fabs(vA.x - vB.x) < 1e-4 &&
                                        fabs(vA.y - vB.y) < 1e-4) {
                                        found = true;
                                        break;
                                    }
                                }
                                if (!found) {
                                    logger.logf("  MISMATCH: vA=(%.4f,%.4f,%.4f) has no match in neighbor\n", vA.x, vA.y, vA.z);
                                }
                            }
                        }
                    }
                }
            }
        }

        //--- DEBUG: Check chunk border positions ---
        logger.logf("Checking chunk border positions...\n");
        for (int cz = 0; cz < chunkGrid; ++cz) {
            for (int cy = 0; cy < chunkGrid; ++cy) {
                for (int cx = 0; cx < chunkGrid; ++cx) {
                    Int3 chunkPos = {cx, cy, cz};
                    Vector3 chunkOffsetA = {
                        minCorner.x + cx * CHUNK_SIZE,
                        minCorner.y + cy * CHUNK_SIZE,
                        minCorner.z + cz * CHUNK_SIZE
                    };

                    // +X neighbor
                    if (cx < chunkGrid - 1) {
                        Vector3 chunkOffsetB = {
                            minCorner.x + (cx + 1) * CHUNK_SIZE,
                            minCorner.y + cy * CHUNK_SIZE,
                            minCorner.z + cz * CHUNK_SIZE
                        };
                        logger.logf("Checking X border positions: (%d,%d,%d) <-> (%d,%d,%d)\n", cx, cy, cz, cx+1, cy, cz);
                        checkChunkBorderPositions(chunkOffsetA, chunkOffsetB, CHUNK_SIZE + 1, 'x');
                    }
                    // +Y neighbor
                    if (cy < chunkGrid - 1) {
                        Vector3 chunkOffsetB = {
                            minCorner.x + cx * CHUNK_SIZE,
                            minCorner.y + (cy + 1) * CHUNK_SIZE,
                            minCorner.z + cz * CHUNK_SIZE
                        };
                        logger.logf("Checking Y border positions: (%d,%d,%d) <-> (%d,%d,%d)\n", cx, cy, cz, cx, cy+1, cz);
                        checkChunkBorderPositions(chunkOffsetA, chunkOffsetB, CHUNK_SIZE + 1, 'y');
                    }
                    // +Z neighbor
                    if (cz < chunkGrid - 1) {
                        Vector3 chunkOffsetB = {
                            minCorner.x + cx * CHUNK_SIZE,
                            minCorner.y + cy * CHUNK_SIZE,
                            minCorner.z + (cz + 1) * CHUNK_SIZE
                        };
                        logger.logf("Checking Z border positions: (%d,%d,%d) <-> (%d,%d,%d)\n", cx, cy, cz, cx, cy, cz+1);
                        checkChunkBorderPositions(chunkOffsetA, chunkOffsetB, CHUNK_SIZE + 1, 'z');
                    }
                }
            }
        }

        //--- DEBUG: Check shared edge caches for mismatches ---
        logger.logf("Checking shared edge caches for mismatches...\n");
        for (int cz = 0; cz < chunkGrid; ++cz) {
            for (int cy = 0; cy < chunkGrid; ++cy) {
                for (int cx = 0; cx < chunkGrid; ++cx) {
                    Int3 chunkPos = {cx, cy, cz};
                    const auto& chunkA = chunks[chunkPos];
                    // +X neighbor
                    if (cx < chunkGrid - 1) {
                        Int3 neighborX = {cx + 1, cy, cz};
                        const auto& chunkB = chunks[neighborX];
                        logger.logf("Checking X border: (%d,%d,%d) <-> (%d,%d,%d)\n", cx, cy, cz, cx+1, cy, cz);
                        checkChunkBorderNoise(chunkA.noiseValues, chunkB.noiseValues, CHUNK_SIZE + 1, 'x');
                    }
                    // +Y neighbor
                    if (cy < chunkGrid - 1) {
                        Int3 neighborY = {cx, cy + 1, cz};
                        const auto& chunkB = chunks[neighborY];
                        logger.logf("Checking Y border: (%d,%d,%d) <-> (%d,%d,%d)\n", cx, cy, cz, cx, cy+1, cz);
                        checkChunkBorderNoise(chunkA.noiseValues, chunkB.noiseValues, CHUNK_SIZE + 1, 'y');
                    }
                    // +Z neighbor
                    if (cz < chunkGrid - 1) {
                        Int3 neighborZ = {cx, cy, cz + 1};
                        const auto& chunkB = chunks[neighborZ];
                        logger.logf("Checking Z border: (%d,%d,%d) <-> (%d,%d,%d)\n", cx, cy, cz, cx, cy, cz+1);
                        checkChunkBorderNoise(chunkA.noiseValues, chunkB.noiseValues, CHUNK_SIZE + 1, 'z');
                    }
                }
            }
        }

        clock_t planetoidGenEnd = clock();
        double planetoidGenTime = static_cast<double>(planetoidGenEnd - planetoidGenStart) / CLOCKS_PER_SEC;
        logger.logf("Planetoid %d generated in %.2f seconds.\n\n", i, planetoidGenTime);
        delete noise; // Clean up the noise instance
        planetoidsGenerated++;
    }
    return planetoidsGenerated; // Return the number of planetoids generated
}


void loadWorld(Scene& world) {
    // Load the world scene from a file or initialize it
    logger.logf("Loading world scene...\n");
    std::ifstream file("assets/scenes/world");
    if (!file.is_open()) {
        std::cerr << "Failed to open world scene file." << std::endl;
        return;
    }

    std::string line;
    while (getline(file, line)) {
        if (line.empty()) continue; // Skip empty lines

        logger.logf("Processing line: %s\n", line.c_str());
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> result;
        while(std::getline(ss, cell, ',')) {
            result.push_back(cell);
        }
        if (result.size() == 13) {
            std::string type = result[0];
            std::string name = result[1];
            Vector3 position = { std::stof(result[2]), std::stof(result[3]), std::stof(result[4]) };
            Vector3 rotation = { std::stof(result[5]), std::stof(result[6]), std::stof(result[7]) };
            Color color = { static_cast<unsigned char>(std::stoi(result[8])),
                            static_cast<unsigned char>(std::stoi(result[9])),
                            static_cast<unsigned char>(std::stoi(result[10])),
                            static_cast<unsigned char>(std::stoi(result[11])) };
            float scale = std::stof(result[12]);

            if (type == "gameobject") {
                logger.logf("Creating GameObject: %s\n", name.c_str());
                std::string modelPath = "assets/models/" + name + "/" + name + ".obj";
                std::string texturePath = "assets/models/" + name + "/diffuse.png";
                Model model;
                // Model model = LoadModelFromMesh(GenMeshCube(1.0f, 1.0f, 1.0f)); // Default model if loading fails
                if (!FileExists(modelPath.c_str())) {
                    std::cerr << "Model file does not exist: " << modelPath << std::endl;
                    continue;
                } else {
                    model = LoadModel(modelPath.c_str());
                    logger.logf("meshCount: %d, materialCount: %d, boneCount: %d\n", model.meshCount, model.materialCount, model.boneCount);
                    logger.logf("bones: %p, bindPose: %p\n", model.bones, model.bindPose);
                    if (model.boneCount > 0 && (model.bones == nullptr || model.bindPose == nullptr)) {
                        model.boneCount = 0;
                        model.bones = nullptr;
                        model.bindPose = nullptr;
                    }
                }

                if (model.meshCount == 0 || model.materialCount == 0) {
                    std::cerr << "Failed to load model or materials: " << modelPath << std::endl;
                    continue;
                }

                if (FileExists(texturePath.c_str()) && model.materials != nullptr && model.materials[0].maps != nullptr) {
                    model.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = LoadTexture(texturePath.c_str());
                } else {
                    logger.logf("Texture file does not exist or materials/maps not allocated: %s\n", texturePath.c_str());
                }

                world.objects.push_back(std::make_unique<GameObject>(type, name, position, rotation, color, scale, model));
                logger.logf("GameObject created: %s\n", name.c_str());
            } else if (type == "text") {
                world.uiObjects.push_back(std::make_unique<TextObject>(name, position, rotation, color, scale, static_cast<int>(scale * 20))); // Assuming scale is used for font size
            } else if (type == "cube" || type == "sphere") {
                world.objects.push_back(std::make_unique<GameObject>(type, name, position, rotation, color, scale));
            } else {
                std::cerr << "Unknown object type: " << type << std::endl;
                continue; // Skip unknown types
            }
        } else {
            std::cerr << "Invalid line format: " << line << std::endl;
            continue; // Skip invalid lines
        }
    }
    file.close();
}

