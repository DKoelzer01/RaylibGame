#include "world.h"
#include "utils/logger.h"

//Setup RNG
float randScale = 100000.0f; // Scale for random number generation
std::random_device rd; // Obtain a seed from the system
std::mt19937 gen(rd()); // Initialize the Mersenne Twister engine with the seed
std::uniform_int_distribution<> distrib(1, randScale);

// Helpers for 3D integer coordinates as chunk keys


constexpr int CHUNK_SIZE = 32; // You can adjust this for performance/memory

Planetoid::Planetoid(std::string name, Vector3 position, Vector3 rotation, Color color, float scale, size_t size)
    : Object("planetoid", name, position, rotation, color, scale), size(size) {
    // Initialize the planetoid with the given size
    this->sharedEdgeCaches.clear();
    this->generatedChunks.clear();
}

float Planetoid::GetNoise(float wx, float wy, float wz) {
    // Convert world coordinates to chunk coordinates
    int chunkX = static_cast<int>((wx-position.x) / CHUNK_SIZE);
    int chunkY = static_cast<int>((wy-position.y)  / CHUNK_SIZE);
    int chunkZ = static_cast<int>((wz-position.z)  / CHUNK_SIZE);
    Int3 chunkKey = {chunkX, chunkY, chunkZ};
    // Check if the chunk has been generated
    if (chunkChildren.find(chunkKey) == chunkChildren.end()) {
       return 0.0f; // Return 0 if the chunk does not exist
    } else {
        // If the chunk exists, retrieve it
        auto& chunkObject = chunkChildren[chunkKey];
        if (chunkObject) {
            // Calculate local coordinates within the chunk
            float localX = wx - (position.x + chunkX * CHUNK_SIZE);
            float localY = wy - (position.y + chunkY * CHUNK_SIZE);
            float localZ = wz - (position.z + chunkZ * CHUNK_SIZE);
            // Get noise value from the chunk's noise values
            return chunkObject->chunk.noiseValues[static_cast<int>(localX) + static_cast<int>(localY) * size + static_cast<int>(localZ) * size * size];
        }
        return 0.0f; // Return 0 if the chunk does not exist
    }
}

Planetoid::~Planetoid() = default;

void Planetoid::draw(Shader* lightingShader) {
    if (!isActive) return; // Skip drawing if the object is not active
    for(const auto& objPtr : children) {
        objPtr->draw(lightingShader);
    }
    for(const auto& objPtr : chunkChildren) {
        objPtr.second->draw(lightingShader);
    }
}

void Planetoid::drawDepthOnly(const Matrix& lightSpaceMatrix, Shader* depthShader) {
    if (!isActive) return;
    for(const auto& objPtr : children) {
        objPtr->drawDepthOnly(lightSpaceMatrix, depthShader);
    }
    for(const auto& objPair : chunkChildren) {
        objPair.second->drawDepthOnly(lightSpaceMatrix, depthShader);
    }
}

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

bool isAllBelowThreshold(const std::vector<float>& vec, float threshold) {
    return std::all_of(vec.begin(), vec.end(), [threshold](float i){return i < threshold; });
}

void logNoiseValues(const std::vector<float>& noiseValues, int size) {
    if (size <= 0) return;
    logger.logf("Noise values (size %d):\n", size);
    float minNoise = noiseValues[0];
    float maxNoise = noiseValues[0];
    for (int z = 0; z < size; ++z) {
        logger.logf("z = %d:\n", z);
        for (int y = 0; y < size; ++y) {
            for (int x = 0; x < size; ++x) {
                int idx = x + y * size + z * size * size;
                if (noiseValues[idx] < minNoise) minNoise = noiseValues[idx];
                if (noiseValues[idx] > maxNoise) maxNoise = noiseValues[idx];
                logger.logf("%02.4f ", noiseValues[idx]);
            }
            logger.logf("\n");
        }
        logger.logf("\n");
    }
    logger.logf("Min noise value: %.4f, Max noise value: %.4f\n", minNoise, maxNoise);
}

void weightNoise(std::vector<float>& noiseValues, int size, const Vector3& chunkWorldPos, const Vector3& planetoidCenter, float falloff, float maxDist) {
    float minScale = 1.0f;
    float maxScale = 0.0f;
    for (int z = 0; z < size; ++z) {
        for (int y = 0; y < size; ++y) {
            for (int x = 0; x < size; ++x) {
                int idx = x + y * size + z * size * size;
                // Compute voxel position in world space for correct distance calculation
                Vector3 voxelPos = { chunkWorldPos.x + x + planetoidCenter.x, chunkWorldPos.y + y + planetoidCenter.y, chunkWorldPos.z + z + planetoidCenter.z };
                float dist = Vector3Distance(voxelPos, planetoidCenter);
                float scale = std::max(0.0f, 1.0f - (dist / maxDist));
                minScale = std::min(minScale, scale);
                maxScale = std::max(maxScale, scale);
                noiseValues[idx] *= scale;
            }
        }
    }
    logger.logf("[weightNoise] minScale=%.4f maxScale=%.4f\n", minScale, maxScale);
}

void normalizeNoise(std::vector<float>& noiseValues, int size, float scale = 1.0f) {
    float minValue = 0.0f; 
    float maxValue = 1.0f;

    if (size <= 0) return;
    float minNoise = noiseValues[0];
    float maxNoise = noiseValues[0];
    
    // Find min and max noise values
    for (float v : noiseValues) {
        if (v < minNoise) minNoise = v;
        if (v > maxNoise) maxNoise = v;
    }

    // Normalize the noise values to the range [minValue, maxValue]
    for (float& v : noiseValues) {
        v = ((v - minNoise) / (maxNoise - minNoise)) * (maxValue - minValue) + minValue; 
        v *= scale; // Scale the normalized value
        v = round(v); // Round to the nearest integer
    }
}

void generateWorld(Scene& world) {
    logger.logf("Generating world scene...\n");
    clock_t worldGenStart = clock();

    float genRange = 500.0f; // Radius for planetoid generation
    float minSize = 250.0f; // Minimum size of planetoids
    float maxSize = 500.0f; // Maximum size of planetoids
    logger.logf("Generated %d planetoids.\n", generatePlanetoids(randScale, world, genRange, minSize, maxSize)); // Generate planetoids in the world

    clock_t worldGenEnd = clock();
    double worldGenTime = static_cast<double>(worldGenEnd - worldGenStart) / CLOCKS_PER_SEC;
    logger.logf("World generation completed in %.2f seconds.\n", worldGenTime);
}

// Update: Now takes chunkWorldPos (relative to planetoid) and origin (planetoid world position)
Chunk generateChunk(const Vector3& chunkWorldPos, const Vector3& origin, SimplexNoise* noise) {
    Chunk chunk = {};
    logger.logf("[generateChunk] Generating chunk at chunkWorldPos (%.2f, %.2f, %.2f) and origin (%.2f, %.2f, %.2f)\n", chunkWorldPos.x, chunkWorldPos.y, chunkWorldPos.z, origin.x, origin.y, origin.z);
    // Allocate noise for this chunk with a 1-voxel border
    const int chunkNoiseSize = (CHUNK_SIZE + 1);
    const int totalNoise = chunkNoiseSize * chunkNoiseSize * chunkNoiseSize;
    chunk.noiseValues.resize(totalNoise);
    chunk.position = {(int)chunkWorldPos.x, (int)chunkWorldPos.y, (int)chunkWorldPos.z}; // Set to relative position only
    // Use local chunk coordinates for mesh, but sample noise at world position
    size_t idx = 0;
    for (int z = 0; z <= CHUNK_SIZE; ++z) {
        for (int y = 0; y <= CHUNK_SIZE; ++y) {
            for (int x = 0; x <= CHUNK_SIZE; ++x, ++idx) {
                float wx = origin.x + chunkWorldPos.x + x;
                float wy = origin.y + chunkWorldPos.y + y;
                float wz = origin.z + chunkWorldPos.z + z;
                chunk.noiseValues[idx] = noise->fractal(4, wx, wy, wz) + 1.0f;
            }
        }
    }
    return chunk;
}

// Iterative chunk generation using BFS
void iterativeChunk(int startCx, int startCy, int startCz, const Vector3& origin, Vector3 rotation, Color color, float scale, SimplexNoise* noise, Planetoid* planetoid, Scene& world) {
    struct QueueEntry {
        int cx, cy, cz;
    };
    std::queue<QueueEntry> q;
    q.push({startCx, startCy, startCz});
    constexpr int dirs[6][3] = {
        {1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}
    };
    while (!q.empty()) {
        time_t genStart = clock();
        auto [cx, cy, cz] = q.front();
        q.pop();
        Int3 ChunkKey{cx, cy, cz};
        // logger.logf("[iterativeChunk] Attempting chunk (%d, %d, %d) at planetoid origin (%.2f, %.2f, %.2f)\n", cx, cy, cz, origin.x, origin.y, origin.z);
        if (planetoid->generatedChunks.count(ChunkKey) > 0) {
            // logger.logf("[iterativeChunk] Skipping already generated chunk (%d, %d, %d)\n", cx, cy, cz);
            continue;
        }
        planetoid->generatedChunks[ChunkKey] = true;
        Vector3 chunkWorldPos = {
            static_cast<float>(cx * CHUNK_SIZE),
            static_cast<float>(cy * CHUNK_SIZE),
            static_cast<float>(cz * CHUNK_SIZE)
        };
        // chunkLocalPos is the chunk's offset within the planetoid (chunk-local)
        Vector3 chunkLocalPos = {
            static_cast<float>(cx * CHUNK_SIZE),
            static_cast<float>(cy * CHUNK_SIZE),
            static_cast<float>(cz * CHUNK_SIZE)
        };
        logger.logf("[iterativeChunk] Creating ChunkObject at chunkLocalPos (%.2f, %.2f, %.2f) for chunk (%d, %d, %d)",
            chunkLocalPos.x, chunkLocalPos.y, chunkLocalPos.z, cx, cy, cz);
        std::string chunkName;
        chunkName.reserve(32);
        chunkName = "chunk_" + std::to_string(cx) + "_" + std::to_string(cy) + "_" + std::to_string(cz);
        time_t genChunk = clock();
        Chunk chunk = generateChunk(chunkLocalPos, origin, noise); // Pass chunk-local position only
        // Log the chunk noise values min and max
        // logger.logf("[generateChunk] Chunk (%d, %d, %d) noise values min: %.4f, max: %.4f\n",
        //     cx, cy, cz,
        //     *std::min_element(chunk.noiseValues.begin(), chunk.noiseValues.end()),
        //     *std::max_element(chunk.noiseValues.begin(), chunk.noiseValues.end()));
        // writeNoiseValuesToFile(chunk.noiseValues, CHUNK_SIZE + 1, "assets/noise/" + chunkName + ".txt");
        weightNoise(chunk.noiseValues, CHUNK_SIZE + 1, chunkLocalPos, origin, 0.5f, planetoid->size);
        // Log the weighted noise values min and max
        // logger.logf("[weightNoise] Chunk (%d, %d, %d) weighted noise values min: %.4f, max: %.4f\n",
        //     cx, cy, cz,
        //     *std::min_element(chunk.noiseValues.begin(), chunk.noiseValues.end()),
        //     *std::max_element(chunk.noiseValues.begin(), chunk.noiseValues.end()));
        // writeNoiseValuesToFile(chunk.noiseValues, CHUNK_SIZE + 1, "assets/noise/" + chunkName + "_weighted.txt");
        // Log the chunk generation time
        time_t genChunkEnd = clock();
        double genChunkTime = static_cast<double>(genChunkEnd - genChunk) / CLOCKS_PER_SEC;
        logger.logf("Chunk (%d, %d, %d) generated in %.2f seconds.\n", cx, cy, cz, genChunkTime);
        if (isAllBelowThreshold(chunk.noiseValues, 0.01f)) { // Lowered threshold
            planetoid->chunkChildren.emplace(std::make_pair(Int3{cx, cy, cz}, std::make_unique<ChunkObject>("chunk", chunkName, chunkLocalPos, rotation, color, scale, chunk)));
            continue;
        }
        // --- Shared edge cache pointers for this chunk ---
        std::vector<EdgeCacheEntry> localEdgeCache(CHUNK_SIZE * CHUNK_SIZE * CHUNK_SIZE * 12);
        std::vector<EdgeCacheEntry>* edgeCacheXPos = nullptr;
        std::vector<EdgeCacheEntry>* edgeCacheYPos = nullptr;
        std::vector<EdgeCacheEntry>* edgeCacheZPos = nullptr;
        std::vector<EdgeCacheEntry>* edgeCacheXNeg = nullptr;
        std::vector<EdgeCacheEntry>* edgeCacheYNeg = nullptr;
        std::vector<EdgeCacheEntry>* edgeCacheZNeg = nullptr;
        // Canonical edge cache keying for all 6 faces
        // +X face: between (cx,cy,cz) and (cx+1,cy,cz), canonical is min(cx,cx+1)
        auto keyXPos = std::make_tuple(std::min(cx, cx+1), cy, cz, 0); // 0 = X face
        if (planetoid->sharedEdgeCaches.count(keyXPos) == 0) {
            logger.logf("[EdgeCache] Creating edgeCacheXPos for chunk (%d,%d,%d) key (%d,%d,%d,0)\n", cx, cy, cz, std::min(cx, cx+1), cy, cz);
            planetoid->sharedEdgeCaches[keyXPos] = std::vector<EdgeCacheEntry>(CHUNK_SIZE * CHUNK_SIZE * 12);
        } else {
            logger.logf("[EdgeCache] Reusing edgeCacheXPos for chunk (%d,%d,%d) key (%d,%d,%d,0)\n", cx, cy, cz, std::min(cx, cx+1), cy, cz);
        }
        edgeCacheXPos = &planetoid->sharedEdgeCaches[keyXPos];
        // -X face: between (cx,cy,cz) and (cx-1,cy,cz)
        auto keyXNeg = std::make_tuple(std::min(cx, cx-1), cy, cz, 0);
        if (planetoid->sharedEdgeCaches.count(keyXNeg) == 0) {
            logger.logf("[EdgeCache] Creating edgeCacheXNeg for chunk (%d,%d,%d) key (%d,%d,%d,0)\n", cx, cy, cz, std::min(cx, cx-1), cy, cz);
            planetoid->sharedEdgeCaches[keyXNeg] = std::vector<EdgeCacheEntry>(CHUNK_SIZE * CHUNK_SIZE * 12);
        } else {
            logger.logf("[EdgeCache] Reusing edgeCacheXNeg for chunk (%d,%d,%d) key (%d,%d,%d,0)\n", cx, cy, cz, std::min(cx, cx-1), cy, cz);
        }
        edgeCacheXNeg = &planetoid->sharedEdgeCaches[keyXNeg];
        // +Y face
        auto keyYPos = std::make_tuple(cx, std::min(cy, cy+1), cz, 1);
        if (planetoid->sharedEdgeCaches.count(keyYPos) == 0) {
            logger.logf("[EdgeCache] Creating edgeCacheYPos for chunk (%d,%d,%d) key (%d,%d,%d,1)\n", cx, cy, cz, cx, std::min(cy, cy+1), cz);
            planetoid->sharedEdgeCaches[keyYPos] = std::vector<EdgeCacheEntry>(CHUNK_SIZE * CHUNK_SIZE * 12);
        } else {
            logger.logf("[EdgeCache] Reusing edgeCacheYPos for chunk (%d,%d,%d) key (%d,%d,%d,1)\n", cx, cy, cz, cx, std::min(cy, cy+1), cz);
        }
        edgeCacheYPos = &planetoid->sharedEdgeCaches[keyYPos];
        // -Y face
        auto keyYNeg = std::make_tuple(cx, std::min(cy, cy-1), cz, 1);
        if (planetoid->sharedEdgeCaches.count(keyYNeg) == 0) {
            logger.logf("[EdgeCache] Creating edgeCacheYNeg for chunk (%d,%d,%d) key (%d,%d,%d,1)\n", cx, cy, cz, cx, std::min(cy, cy-1), cz);
            planetoid->sharedEdgeCaches[keyYNeg] = std::vector<EdgeCacheEntry>(CHUNK_SIZE * CHUNK_SIZE * 12);
        } else {
            logger.logf("[EdgeCache] Reusing edgeCacheYNeg for chunk (%d,%d,%d) key (%d,%d,%d,1)\n", cx, cy, cz, cx, std::min(cy, cy-1), cz);
        }
        edgeCacheYNeg = &planetoid->sharedEdgeCaches[keyYNeg];
        // +Z face
        auto keyZPos = std::make_tuple(cx, cy, std::min(cz, cz+1), 2);
        if (planetoid->sharedEdgeCaches.count(keyZPos) == 0) {
            logger.logf("[EdgeCache] Creating edgeCacheZPos for chunk (%d,%d,%d) key (%d,%d,%d,2)\n", cx, cy, cz, cx, cy, std::min(cz, cz+1));
            planetoid->sharedEdgeCaches[keyZPos] = std::vector<EdgeCacheEntry>(CHUNK_SIZE * CHUNK_SIZE * 12);
        } else {
            logger.logf("[EdgeCache] Reusing edgeCacheZPos for chunk (%d,%d,%d) key (%d,%d,%d,2)\n", cx, cy, cz, cx, cy, std::min(cz, cz+1));
        }
        edgeCacheZPos = &planetoid->sharedEdgeCaches[keyZPos];
        // -Z face
        auto keyZNeg = std::make_tuple(cx, cy, std::min(cz, cz-1), 2);
        if (planetoid->sharedEdgeCaches.count(keyZNeg) == 0) {
            logger.logf("[EdgeCache] Creating edgeCacheZNeg for chunk (%d,%d,%d) key (%d,%d,%d,2)\n", cx, cy, cz, cx, cy, std::min(cz, cz-1));
            planetoid->sharedEdgeCaches[keyZNeg] = std::vector<EdgeCacheEntry>(CHUNK_SIZE * CHUNK_SIZE * 12);
        } else {
            logger.logf("[EdgeCache] Reusing edgeCacheZNeg for chunk (%d,%d,%d) key (%d,%d,%d,2)\n", cx, cy, cz, cx, cy, std::min(cz, cz-1));
        }
        edgeCacheZNeg = &planetoid->sharedEdgeCaches[keyZNeg];
        // For now, pass +X, +Y, +Z caches to marchCube (legacy interface)
        chunk.vertices.clear();
        chunk.indices.clear();
        std::vector<Vector3> normals; // NEW: store per-vertex normals
        chunk.vertices.reserve(CHUNK_SIZE * CHUNK_SIZE * CHUNK_SIZE * 3 / 2);
        chunk.indices.reserve(CHUNK_SIZE * CHUNK_SIZE * CHUNK_SIZE * 6);
        normals.reserve(CHUNK_SIZE * CHUNK_SIZE * CHUNK_SIZE * 3 / 2);
        // time_t cubemarchStart = clock();
        
        // --- Border noise debug logging ---
        // Define neighbor chunk keys for border noise debug
        Int3 neighborX{cx+1, cy, cz};
        Int3 neighborY{cx, cy+1, cz};
        Int3 neighborZ{cx, cy, cz+1};
        if (planetoid->generatedChunks.count(neighborX) > 0) {
            ChunkObject* tempchunk = planetoid->chunkChildren[neighborX].get();
            checkChunkBorderNoise(chunk.noiseValues, tempchunk->chunk.noiseValues, CHUNK_SIZE + 1, 'x');
        }
        if (planetoid->generatedChunks.count(neighborY) > 0) {
            ChunkObject* tempchunk = planetoid->chunkChildren[neighborY].get();
            checkChunkBorderNoise(chunk.noiseValues, tempchunk->chunk.noiseValues, CHUNK_SIZE + 1, 'y');
        }
        if (planetoid->generatedChunks.count(neighborZ) > 0) {
            ChunkObject* tempchunk = planetoid->chunkChildren[neighborZ].get();
            checkChunkBorderNoise(chunk.noiseValues, tempchunk->chunk.noiseValues, CHUNK_SIZE + 1, 'z');
        }

        for(int k = 0; k < CHUNK_SIZE; ++k) {
            for(int j = 0; j < CHUNK_SIZE; ++j) {
                for(int i = 0; i < CHUNK_SIZE; ++i) {
                    // Select correct edge cache for each face
                    std::vector<EdgeCacheEntry>* edgeCacheX = nullptr;
                    std::vector<EdgeCacheEntry>* edgeCacheY = nullptr;
                    std::vector<EdgeCacheEntry>* edgeCacheZ = nullptr;
                    if (i == 0) edgeCacheX = edgeCacheXNeg;
                    else if (i == CHUNK_SIZE-1) edgeCacheX = edgeCacheXPos;
                    if (j == 0) edgeCacheY = edgeCacheYNeg;
                    else if (j == CHUNK_SIZE-1) edgeCacheY = edgeCacheYPos;
                    if (k == 0) edgeCacheZ = edgeCacheZNeg;
                    else if (k == CHUNK_SIZE-1) edgeCacheZ = edgeCacheZPos;
                    // Pass chunkLocalPos (not chunkWorldPos) as chunkOrigin to marchCube
                    marchCube(
                        i, j, k,
                        chunk.noiseValues, CHUNK_SIZE + 1,
                        chunk.vertices, chunk.indices, normals, // pass normals
                        &localEdgeCache, edgeCacheX, edgeCacheY, edgeCacheZ,
                        0.5f,
                        cx, cy, cz,
                        chunkLocalPos, // <-- fix: use chunkLocalPos, not chunkWorldPos
                        noise,
                        planetoid
                    );
                }
            }
        }
        // time_t cubemarchEnd = clock();
        // double cubemarchTime = static_cast<double>(cubemarchEnd - cubemarchStart) / CLOCKS_PER_SEC;
        // logger.logf("Chunk (%d, %d, %d) cubemarched in %.2f\n",cx, cy, cz, cubemarchTime);
        // normalizeNoise(chunk.noiseValues, CHUNK_SIZE + 1, 16.0f); // Normalize noise values
        for (size_t j = 0; j < std::min(chunk.vertices.size(), size_t(5)); ++j) {
            logger.logf("Chunk (%d,%d,%d) vertex %zu: (%.2f, %.2f, %.2f)", cx, cy, cz, j, chunk.vertices[j].x, chunk.vertices[j].y, chunk.vertices[j].z);
        }
        Mesh mesh = { 0 };
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
        // --- Upload normals ---
        mesh.normals = new float[mesh.vertexCount * 3];
        for (size_t j = 0; j < normals.size(); ++j) {
            mesh.normals[j * 3 + 0] = normals[j].x;
            mesh.normals[j * 3 + 1] = normals[j].y;
            mesh.normals[j * 3 + 2] = normals[j].z;
        }
        // --- Check for mismatch before uploading mesh ---
        if (normals.size() != chunk.vertices.size()) {
            logger.logf("ERROR: normals.size() (%zu) != vertices.size() (%zu)\n", normals.size(), chunk.vertices.size());
            // Optionally: abort or return here to avoid hang
            return;
        }
        UploadMesh(&mesh, false); // Possibly make this dynamic if i want to update the mesh later
        chunk.mesh = mesh;
        chunk.model = LoadModelFromMesh(chunk.mesh);
        chunk.model.materials[0].shader = world.lightingShader; // Use the lighting shader for the chunk
        auto inserted = planetoid->chunkChildren.emplace(std::make_pair(Int3{cx, cy, cz}, std::make_unique<ChunkObject>("chunk", chunkName, chunkLocalPos, rotation, color, scale, chunk)));
        if (inserted.second) {
            inserted.first->second->isActive = true;
            inserted.first->second->parent = planetoid; // Set parent pointer for correct transform
        }
        for (int d = 0; d < 6; ++d) {
            int nx = cx + dirs[d][0];
            int ny = cy + dirs[d][1];
            int nz = cz + dirs[d][2];
            Int3 nkey{nx, ny, nz};
            if (planetoid->generatedChunks.count(nkey) == 0) {
                q.push({nx, ny, nz});
            }
        }
        time_t genEnd = clock();
        double genTime = static_cast<double>(genEnd - genStart) / CLOCKS_PER_SEC;
        logger.logf("Chunk (%d, %d, %d) generated in %.2f seconds.\n", cx, cy, cz, genTime);
    }
}


void generatePlanetoid(float randScale,std::string name, Scene& world, Vector3 position, Vector3 rotation,Color color, size_t size, float scale) {
    Planetoid* planetoid = new Planetoid(name, position, rotation, color, scale, size);
    planetoid->parent = &world.rootObject; // Set the planetoid's parent to the root object of the scene
    // Initialize shared edge caches with correct type
    planetoid->sharedEdgeCaches = std::unordered_map<std::tuple<int, int, int, int>, std::vector<EdgeCacheEntry>, Tuple4Hash>(); // Initialize shared edge caches
    world.objects.push_back(std::unique_ptr<Object>(planetoid)); // Add the planetoid to the scene's object list

    clock_t planetoidGenStart = clock();
    
    float frequencyNoise = static_cast<float>(distrib(gen)/randScale)-0.5f; // Random frequency noise to add some variation from -0.5 to 0.5
    float seedX = static_cast<float>(distrib(gen)/10.f); // Random seed for X coordinate
    float seedY = static_cast<float>(distrib(gen)/10.f); // Random seed for Y coordinate
    float seedZ = static_cast<float>(distrib(gen)/10.f); // Random seed for Z coordinate
    planetoid->seed = {seedX, seedY, seedZ}; // Set the planetoid's seed for noise generation
    logger.logf("Planetoid generation parameters: frequencyNoise=%.4f, seedX=%.4f, seedY=%.4f, seedZ=%.4f\n", 
        frequencyNoise, seedX, seedY, seedZ);
    SimplexNoise* noise = new SimplexNoise((1.5f+frequencyNoise)/static_cast<float>(size), 2.0f, 2.0f, 0.5f,seedX,seedY,seedZ); // Create a new instance of SimplexNoise
    if (!noise) {
        std::cerr << "Failed to create SimplexNoise instance." << std::endl;
        return; // Skip this planetoid if noise generation fails;
    }
    logger.logf("Generating planetoid at (%f, %f, %f) with size %zu and scale %.2f with noise frequency %.4f...\n", position.x, position.y, position.z, size, scale,(1.5f+frequencyNoise)/static_cast<float>(size));

    // Use iterative chunk generation
    iterativeChunk(0, 0, 0, position, rotation, color, scale, noise, planetoid, world); 
    

    planetoid->isActive = true; // Set the planetoid as active
    delete noise; // Clean up the noise instance
    time_t planetoidGenEnd = clock();
    double planetoidGenTime = static_cast<double>(planetoidGenEnd - planetoidGenStart) / CLOCKS_PER_SEC;
    logger.logf("Planetoid generated in %.2f seconds at position (%f, %f, %f) with size %zu.\n", 
        planetoidGenTime, position.x, position.y, position.z, size);
}

int generatePlanetoids(float randScale, Scene& world, float genRange, float minSize, float maxSize) {
    int planetoidsGenerated = 0;

    //Generate planetoids
    int numPlanetoids = ((distrib(gen)/randScale) * 5) + 1; // Number of planetoids to generate
    numPlanetoids = 1; // For testing purposes, generate only one planetoid
    logger.logf("Generating %d planetoids...\n", numPlanetoids);

    for(int i = 0; i < numPlanetoids; ++i) {
        clock_t planetoidGenStart = clock();
        // Generate random position, rotation, color, size, and scale for the planetoid

        // Position is centered around (0,0,0) with a range of genRange
        Vector3 position = {    static_cast<float>(((distrib(gen)/randScale)*genRange)-(genRange/2)), 
                                static_cast<float>(((distrib(gen)/randScale)*genRange)-(genRange/2)), 
                                static_cast<float>(((distrib(gen)/randScale)*genRange)-(genRange/2))};

        position = {64.0f, 0.0f, 0.0f}; // Fixed position for testing

        // Rotation is random in degrees, scaled by randScale
        Vector3 rotation = {    static_cast<float>(distrib(gen)/randScale * 360), 
                                static_cast<float>(distrib(gen)/randScale * 360), 
                                static_cast<float>(distrib(gen)/randScale * 360)};
        rotation = {0.0f, 0.0f, 0.0f}; // Fixed rotation for testing

        // Using a random color generator with a range of 0-255 for RGB values
        Color color = { static_cast<unsigned char>(distrib(gen)/randScale * 255), 
                        static_cast<unsigned char>(distrib(gen)/randScale * 255), 
                        static_cast<unsigned char>(distrib(gen)/randScale * 255), 255};

        color = {255,255,255,255}; // Fixed color for testing

        // Size is a random float between minSize and maxSize
        size_t size = static_cast<size_t>((distrib(gen)/randScale)*(maxSize-minSize)) + static_cast<size_t>(minSize);
        size = 100; // Fixed size for testing

        // Scale is a fixed value for now, can be adjusted later
        float scale = 1.0f;

        logger.logf("Generating planetoid %d...\n", i);
        std::string name = "planetoid_" + std::to_string(i);
        generatePlanetoid(randScale, name, world, position, rotation, color, size, scale);
        
        clock_t planetoidGenEnd = clock();
        double planetoidGenTime = static_cast<double>(planetoidGenEnd - planetoidGenStart) / CLOCKS_PER_SEC;
        logger.logf("Planetoid %d generated in %.2f seconds.\n\n", i, planetoidGenTime);
        printf("Planetoid %d generated in %.2f seconds.\n\n", i, planetoidGenTime);
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
                model.materials[0].shader = world.lightingShader; // Ensure correct shader for loaded models
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

