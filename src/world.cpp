#include "world.h"
#include "utils/logger.h"
#include "worldtypes.h"
#include <unordered_map>
#include <vector>
#include <cmath>
#include <tuple>
#include <algorithm>

//Setup RNG
float randScale = 100000.0f; // Scale for random number generation
std::random_device rd; // Obtain a seed from the system
std::mt19937 gen(rd()); // Initialize the Mersenne Twister engine with the seed
std::uniform_int_distribution<> distrib(1, randScale);

// Helpers for 3D integer coordinates as chunk keys


const int CHUNK_SIZE = 32; // You can adjust this for performance/memory

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
            return chunkObject->noiseValues[static_cast<int>(localX) + static_cast<int>(localY) * size + static_cast<int>(localZ) * size * size];
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



// Helper to get noise value from this chunk or a neighbor
static float getNoiseAt(const Chunk* chunk, int x, int y, int z, int size) {
    if (x < 0 || x >= size || y < 0 || y >= size || z < 0 || z >= size) {
        logger.logf("[getNoiseAt WARNING] Out of bounds access: (%d, %d, %d) size=%d\n", x, y, z, size);
        return 0.0f;
    }
    return chunk->noiseValues[x + y * size + z * size * size];
}

template<typename T>
T clamp(T v, T lo, T hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

static Vector3 trilerpVec3(
    const Vector3& c000, const Vector3& c100, const Vector3& c010, const Vector3& c110,
    const Vector3& c001, const Vector3& c101, const Vector3& c011, const Vector3& c111,
    float tx, float ty, float tz
) {
    float x =
        c000.x * (1 - tx) * (1 - ty) * (1 - tz) +
        c100.x * tx * (1 - ty) * (1 - tz) +
        c010.x * (1 - tx) * ty * (1 - tz) +
        c110.x * tx * ty * (1 - tz) +
        c001.x * (1 - tx) * (1 - ty) * tz +
        c101.x * tx * (1 - ty) * tz +
        c011.x * (1 - tx) * ty * tz +
        c111.x * tx * ty * tz;

    float y =
        c000.y * (1 - tx) * (1 - ty) * (1 - tz) +
        c100.y * tx * (1 - ty) * (1 - tz) +
        c010.y * (1 - tx) * ty * (1 - tz) +
        c110.y * tx * ty * (1 - tz) +
        c001.y * (1 - tx) * (1 - ty) * tz +
        c101.y * tx * (1 - ty) * tz +
        c011.y * (1 - tx) * ty * tz +
        c111.y * tx * ty * tz;

    float z =
        c000.z * (1 - tx) * (1 - ty) * (1 - tz) +
        c100.z * tx * (1 - ty) * (1 - tz) +
        c010.z * (1 - tx) * ty * (1 - tz) +
        c110.z * tx * ty * (1 - tz) +
        c001.z * (1 - tx) * (1 - ty) * tz +
        c101.z * tx * (1 - ty) * tz +
        c011.z * (1 - tx) * ty * tz +
        c111.z * tx * ty * tz;

    return {x, y, z};
}


void Chunk::calculateNormals() {
    if (vertices.empty() || mesh.vertexCount == 0) {
        logger.logf("[deferredNormals]: Skipping normal calculation for chunk at (%d, %d, %d) because it has no vertices.\n", position.x, position.y, position.z);
        return;
    }
    logger.logf("[deferredNormals]: Calculating normals for chunk at (%d, %d, %d): vertices.size() = %zu, mesh.vertexCount = %d\n", position.x, position.y, position.z, vertices.size(), mesh.vertexCount);
    const int size = CHUNK_SIZE + 1;
    logger.logf("[deferredNormals]: newNormals size: %d, expected size: %d\n", size*size*size, mesh.vertexCount*3);
    std::vector<Vector3> newNormals(size * size * size);

    for (int z = 0; z < size; ++z) {
        for (int y = 0; y < size; ++y) {
            for (int x = 0; x < size; ++x) {
                float dx, dy, dz;
                // X
                if (x > 0 && x < size - 1) {
                    dx = getNoiseAt(this, x + 1, y, z, size) - getNoiseAt(this, x - 1, y, z, size);
                } else {
                    int nx = x == 0 ? -1 : 1;
                    int neighborIdx = -1;
                    for (int i = 0; i < 26; ++i) {
                        if (neighborOffsets[i][0] == nx && neighborOffsets[i][1] == 0 && neighborOffsets[i][2] == 0) {
                            neighborIdx = i;
                            break;
                        }
                    }
                    if (neighborIdx != -1 && neighbors[neighborIdx]) {
                        if (x == 0)
                            dx = getNoiseAt(this, x + 1, y, z, size) - getNoiseAt(neighbors[neighborIdx], size - 2, y, z, size);
                        else
                            dx = getNoiseAt(neighbors[neighborIdx], 1, y, z, size) - getNoiseAt(this, x - 1, y, z, size);
                    } else {
                        dx = 0.0f;
                    }
                }
                // Y
                if (y > 0 && y < size - 1) {
                    dy = getNoiseAt(this, x, y + 1, z, size) - getNoiseAt(this, x, y - 1, z, size);
                } else {
                    int ny = y == 0 ? -1 : 1;
                    int neighborIdx = -1;
                    for (int i = 0; i < 26; ++i) {
                        if (neighborOffsets[i][1] == ny && neighborOffsets[i][0] == 0 && neighborOffsets[i][2] == 0) {
                            neighborIdx = i;
                            break;
                        }
                    }
                    if (neighborIdx != -1 && neighbors[neighborIdx]) {
                        if (y == 0)
                            dy = getNoiseAt(this, x, y + 1, z, size) - getNoiseAt(neighbors[neighborIdx], x, size - 2, z, size);
                        else
                            dy = getNoiseAt(neighbors[neighborIdx], x, 1, z, size) - getNoiseAt(this, x, y - 1, z, size);
                    } else {
                        dy = 0.0f;
                    }
                }
                // Z
                if (z > 0 && z < size - 1) {
                    dz = getNoiseAt(this, x, y, z + 1, size) - getNoiseAt(this, x, y, z - 1, size);
                } else {
                    int nz = z == 0 ? -1 : 1;
                    int neighborIdx = -1;
                    for (int i = 0; i < 26; ++i) {
                        if (neighborOffsets[i][2] == nz && neighborOffsets[i][0] == 0 && neighborOffsets[i][1] == 0) {
                            neighborIdx = i;
                            break;
                        }
                    }
                    if (neighborIdx != -1 && neighbors[neighborIdx]) {
                        if (z == 0)
                            dz = getNoiseAt(this, x, y, z + 1, size) - getNoiseAt(neighbors[neighborIdx], x, y, size - 2, size);
                        else
                            dz = getNoiseAt(neighbors[neighborIdx], x, y, 1, size) - getNoiseAt(this, x, y, z - 1, size);
                    } else {
                        dz = 0.0f;
                    }
                }
                Vector3 n = { -dx, -dy, -dz };
                float len = sqrtf(n.x * n.x + n.y * n.y + n.z * n.z);
                if (len > 1e-6f) {
                    n.x /= len;
                    n.y /= len;
                    n.z /= len;
                }
                newNormals[x + y * size + z * size * size] = n;
                // Debug log for the produced normal and the noise values used
                // logger.logf("[deferredNormals] Normal at (%d,%d,%d): (%.3f, %.3f, %.3f) | dx=%.3f dy=%.3f dz=%.3f\n", x, y, z, n.x, n.y, n.z, dx, dy, dz);
            }
        }
    }
    // Copy new normals to mesh normals
    if (mesh.vertexCount > 0 && mesh.normals != nullptr) {
        if (mesh.vertexCount != vertices.size()) {
            logger.logf("[deferredNormals][ERROR] mesh.vertexCount (%d) != vertices.size() (%zu) in chunk (%d, %d, %d)\n", mesh.vertexCount, vertices.size(), position.x, position.y, position.z);
        }
        float minLx = 1e9f, maxLx = -1e9f, minLy = 1e9f, maxLy = -1e9f, minLz = 1e9f, maxLz = -1e9f;
        for (int i = 0; i < mesh.vertexCount && i < vertices.size(); ++i) {
            const Vector3& v = vertices[i];
            float fx = clamp(v.x, 0.0f, (float)(size - 1));
            float fy = clamp(v.y, 0.0f, (float)(size - 1));
            float fz = clamp(v.z, 0.0f, (float)(size - 1));
            int x0 = (int)floorf(fx), x1 = std::min(x0 + 1, size - 1);
            int y0 = (int)floorf(fy), y1 = std::min(y0 + 1, size - 1);
            int z0 = (int)floorf(fz), z1 = std::min(z0 + 1, size - 1);
            float tx = fx - x0, ty = fy - y0, tz = fz - z0;
            const Vector3& c000 = newNormals[x0 + y0 * size + z0 * size * size];
            const Vector3& c100 = newNormals[x1 + y0 * size + z0 * size * size];
            const Vector3& c010 = newNormals[x0 + y1 * size + z0 * size * size];
            const Vector3& c110 = newNormals[x1 + y1 * size + z0 * size * size];
            const Vector3& c001 = newNormals[x0 + y0 * size + z1 * size * size];
            const Vector3& c101 = newNormals[x1 + y0 * size + z1 * size * size];
            const Vector3& c011 = newNormals[x0 + y1 * size + z1 * size * size];
            const Vector3& c111 = newNormals[x1 + y1 * size + z1 * size * size];
            Vector3 n = trilerpVec3(c000, c100, c010, c110, c001, c101, c011, c111, tx, ty, tz);
            float len = sqrtf(n.x * n.x + n.y * n.y + n.z * n.z);
            if (len > 1e-6f) {
                n.x /= len; n.y /= len; n.z /= len;
            }
            mesh.normals[i * 3 + 0] = n.x;
            mesh.normals[i * 3 + 1] = n.y;
            mesh.normals[i * 3 + 2] = n.z;
            // Debug log for the interpolation process
            if (i < 10) {
                logger.logf("[deferredNormals] Vertex[%d] pos=(%.3f,%.3f,%.3f) local=(%.3f,%.3f,%.3f) fx=%.3f fy=%.3f fz=%.3f x0=%d y0=%d z0=%d tx=%.3f ty=%.3f tz=%.3f\n", i, v.x, v.y, v.z, lx, ly, lz, fx, fy, fz, x0, y0, z0, tx, ty, tz);
                logger.logf("[deferredNormals]  c000=(%.3f,%.3f,%.3f) c100=(%.3f,%.3f,%.3f) c010=(%.3f,%.3f,%.3f) c110=(%.3f,%.3f,%.3f)\n", c000.x, c000.y, c000.z, c100.x, c100.y, c100.z, c010.x, c010.y, c010.z, c110.x, c110.y, c110.z);
                logger.logf("[deferredNormals]  c001=(%.3f,%.3f,%.3f) c101=(%.3f,%.3f,%.3f) c011=(%.3f,%.3f,%.3f) c111=(%.3f,%.3f,%.3f)\n", c001.x, c001.y, c001.z, c101.x, c101.y, c101.z, c011.x, c011.y, c011.z, c111.x, c111.y, c111.z);
                logger.logf("[deferredNormals]  Interpolated normal: (%.3f, %.3f, %.3f)\n", n.x, n.y, n.z);
            }
        }
        // --- Additional Debugging: Check for invalid normals ---
        int invalidNormalCount = 0;
        for (int i = 0; i < mesh.vertexCount; ++i) {
            float nx = mesh.normals[i * 3 + 0];
            float ny = mesh.normals[i * 3 + 1];
            float nz = mesh.normals[i * 3 + 2];
            bool isZero = (fabs(nx) < 1e-6f && fabs(ny) < 1e-6f && fabs(nz) < 1e-6f);
            bool isNaN = (isnan(nx) || isnan(ny) || isnan(nz));
            if (isZero || isNaN) {
                logger.logf("[deferredNormals][INVALID] Normal[%d]: (%.3f, %.3f, %.3f)\n", i, nx, ny, nz);
                ++invalidNormalCount;
            }
        }
        if (invalidNormalCount > 0) {
            logger.logf("[deferredNormals][SUMMARY] Found %d invalid normals in chunk (%d, %d, %d)\n", invalidNormalCount, position.x, position.y, position.z);
        }
        for (int i = 0; i < std::min(5, mesh.vertexCount); ++i) {
            logger.logf("[deferredNormals]: Normal[%d]: (%.3f, %.3f, %.3f)\n", i, mesh.normals[i*3], mesh.normals[i*3+1], mesh.normals[i*3+2]);
        }
        
        if (mesh.normals == nullptr) {
            logger.logf("[deferredNormals][ERROR] mesh.normals is nullptr after assignment!\n");
        }
        UpdateMeshBuffer(mesh, 2, mesh.normals, mesh.vertexCount * 3, 0); // 2 = RLGL_ATTRIBUTE_NORMAL
    }
}

// --- Chunk neighborOffsets definition ---
const int Chunk::neighborOffsets[26][3] = {
    {-1,-1,-1},{ 0,-1,-1},{ 1,-1,-1},
    {-1, 0,-1},{ 0, 0,-1},{ 1, 0,-1},
    {-1, 1,-1},{ 0, 1,-1},{ 1, 1,-1},
    {-1,-1, 0},{ 0,-1, 0},{ 1,-1, 0},
    {-1, 0, 0},           { 1, 0, 0},
    {-1, 1, 0},{ 0, 1, 0},{ 1, 1, 0},
    {-1,-1, 1},{ 0,-1, 1},{ 1,-1, 1},
    {-1, 0, 1},{ 0, 0, 1},{ 1, 0, 1},
    {-1, 1, 1},{ 0, 1, 1},{ 1, 1, 1}
};

void Chunk::tryCalculateNormals() {
    logger.logf("[Chunk] Trying to calculate normals for chunk at position (%d, %d, %d)\n", position.x, position.y, position.z);
    if (allNeighborsPresent()) {
        logger.logf("[Chunk] All neighbors present for chunk at position (%d, %d, %d), calculating normals\n", position.x, position.y, position.z);
        calculateNormals();
        normalsPending = false;
    } else {
        logger.logf("[Chunk] Not all neighbors present for chunk at position (%d, %d, %d), deferring normal calculation\n", position.x, position.y, position.z);
        normalsPending = true;
    }
}

// --- Chunk::assignNeighborsAndNotify implementation ---
void Chunk::assignNeighborsAndNotify(std::unordered_map<Int3, std::unique_ptr<Chunk>>& chunkChildren) {
    logger.logf("[Chunk] Assigning neighbors for chunk at position (%d, %d, %d)\n", position.x, position.y, position.z);
    
    int cx = position.x / CHUNK_SIZE;
    int cy = position.y / CHUNK_SIZE;
    int cz = position.z / CHUNK_SIZE;
    // Get planetoid center and size (radius)
    // Assumes parent pointer is set to Planetoid, otherwise you may need to pass planetoid as argument
    Planetoid* planetoid = dynamic_cast<Planetoid*>(parent);
    Vector3 planetoidCenter = planetoid ? planetoid->position : Vector3{0,0,0};
    float planetoidRadius = planetoid ? planetoid->size : 0.0f;

    for (int i = 0; i < 26; ++i) {
        Int3 npos = { cx + neighborOffsets[i][0], cy + neighborOffsets[i][1], cz + neighborOffsets[i][2] };
        auto it = chunkChildren.find(npos);
        if (it != chunkChildren.end() && it->second) {
            Chunk* neighbor = it->second.get();
            setNeighbor(i, neighbor);
            // Find the reverse index for this neighbor
            for (int j = 0; j < 26; j++) {
                if (neighborOffsets[j][0] == -neighborOffsets[i][0] &&
                    neighborOffsets[j][1] == -neighborOffsets[i][1] &&
                    neighborOffsets[j][2] == -neighborOffsets[i][2]) {
                    neighbor->onNeighborAdded(j, this);
                    break;
                }
            }
        } else {
            // Calculate world position of neighbor chunk center
            Vector3 neighborCenter = {
                (float)(npos.x * CHUNK_SIZE) + CHUNK_SIZE / 2.0f + planetoidCenter.x,
                (float)(npos.y * CHUNK_SIZE) + CHUNK_SIZE / 2.0f + planetoidCenter.y,
                (float)(npos.z * CHUNK_SIZE) + CHUNK_SIZE / 2.0f + planetoidCenter.z
            };
            float dist = sqrtf(
                (neighborCenter.x - planetoidCenter.x) * (neighborCenter.x - planetoidCenter.x) +
                (neighborCenter.y - planetoidCenter.y) * (neighborCenter.y - planetoidCenter.y) +
                (neighborCenter.z - planetoidCenter.z) * (neighborCenter.z - planetoidCenter.z)
            );
            if (dist > planetoidRadius) {
                // Out of bounds, treat as present for normal calculation
                neighborMask |= (1u << i);
            }
        }
    }
    tryCalculateNormals();
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

bool isAllBelowThreshold(const std::vector<float>& vec, float threshold) {
    return std::all_of(vec.begin(), vec.end(), [threshold](float i){return i < threshold; });
}

// Utility to compute dot product of two Vector3
static float dot(const Vector3& a, const Vector3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Utility to compute length of a Vector3
static float length(const Vector3& v) {
    return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

// Utility to compute angle (in degrees) between two normals
static float angleBetween(const Vector3& a, const Vector3& b) {
    float d = dot(a, b) / (length(a) * length(b) + 1e-8f);
    d = clamp(d, -1.0f, 1.0f);
    return acosf(d) * (180.0f / 3.14159265f);
}

// Check for normal discontinuities at chunk borders
void logNormalDiscontinuitiesAtBorders(Chunk& chunkA, Chunk& chunkB, char axis, float angleThresholdDeg = 5.0f) {
    // axis: 'x', 'y', or 'z'. chunkA is at lower coord, chunkB is at higher coord.
    int size = CHUNK_SIZE;

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            int idxA, idxB;
            if (axis == 'x') {
                idxA = (i * size + j) * size + (size - 1); // x = size-1 face of A
                idxB = (i * size + j) * size + 0;           // x = 0 face of B
            } else if (axis == 'y') {
                idxA = ((size - 1) * size + i) * size + j; // y = size-1 face of A
                idxB = (0 * size + i) * size + j;          // y = 0 face of B
            } else { // 'z'
                idxA = (i * size + (size - 1)) * size + j; // z = size-1 face of A
                idxB = (i * size + 0) * size + j;          // z = 0 face of B
            }
            if (idxA < 0 || idxA >= chunkA.mesh.vertexCount || idxB < 0 || idxB >= chunkB.mesh.vertexCount) continue;
            Vector3 nA = {
                chunkA.mesh.normals[idxA * 3 + 0],
                chunkA.mesh.normals[idxA * 3 + 1],
                chunkA.mesh.normals[idxA * 3 + 2]
            };
            Vector3 nB = {
                chunkB.mesh.normals[idxB * 3 + 0],
                chunkB.mesh.normals[idxB * 3 + 1],
                chunkB.mesh.normals[idxB * 3 + 2]
            };
            float angle = angleBetween(nA, nB);
            if (angle > angleThresholdDeg) {
                logger.logf("Normal discontinuity at border (%c): idxA=%d, idxB=%d, angle=%.2f deg\n", axis, idxA, idxB, angle);
            }
        }
    }
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
    // logger.logf("[weightNoise] minScale=%.4f maxScale=%.4f\n", minScale, maxScale);
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

std::unique_ptr<Chunk> generateChunk(const Vector3& chunkWorldPos, const Vector3& origin, SimplexNoise* noise, Object* parent) {
    auto chunk = std::make_unique<Chunk>(); // Default-construct chunk
    chunk->position = {0, 0, 0}; // Initialize chunk position with zero
    chunk->parent = parent; // Set the parent object for the chunk
    // logger.logf("[generateChunk] Generating chunk at chunkWorldPos (%.2f, %.2f, %.2f) and origin (%.2f, %.2f, %.2f)\n", chunkWorldPos.x, chunkWorldPos.y, chunkWorldPos.z, origin.x, origin.y, origin.z);
    // Allocate noise for this chunk with a 1-voxel border
    const int chunkNoiseSize = (CHUNK_SIZE + 1);
    const int totalNoise = chunkNoiseSize * chunkNoiseSize * chunkNoiseSize;
    chunk->noiseValues.resize(totalNoise);
    chunk->position = {(int)chunkWorldPos.x, (int)chunkWorldPos.y, (int)chunkWorldPos.z}; // Set to relative position only
    // Use local chunk coordinates for mesh, but sample noise at world position
    size_t idx = 0;
    for (int z = 0; z <= CHUNK_SIZE; ++z) {
        for (int y = 0; y <= CHUNK_SIZE; ++y) {
            for (int x = 0; x <= CHUNK_SIZE; ++x, ++idx) {
                float wx = origin.x + chunkWorldPos.x + x;
                float wy = origin.y + chunkWorldPos.y + y;
                float wz = origin.z + chunkWorldPos.z + z;
                chunk->noiseValues[idx] = noise->fractal(4, wx, wy, wz) + 1.0f;
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
        // chunkLocalPos is the chunk's offset within the planetoid (chunk-local)
        Vector3 chunkLocalPos = {
            static_cast<float>(cx * CHUNK_SIZE),
            static_cast<float>(cy * CHUNK_SIZE),
            static_cast<float>(cz * CHUNK_SIZE)
        };
        // logger.logf("[iterativeChunk] Creating ChunkObject at chunkLocalPos (%.2f, %.2f, %.2f) for chunk (%d, %d, %d)",
        //    chunkLocalPos.x, chunkLocalPos.y, chunkLocalPos.z, cx, cy, cz);
        std::string chunkName;
        chunkName.reserve(32);
        chunkName = "chunk_" + std::to_string(cx) + "_" + std::to_string(cy) + "_" + std::to_string(cz);
        time_t genChunk = clock();
        auto chunk = generateChunk(chunkLocalPos, origin, noise, planetoid); // Pass chunk-local position only
        // Log the chunk noise values min and max
        // logger.logf("[generateChunk] Chunk (%d, %d, %d) noise values min: %.4f, max: %.4f\n",
        //     cx, cy, cz,
        //     *std::min_element(chunk->noiseValues.begin(), chunk->noiseValues.end()),
        //     *std::max_element(chunk->noiseValues.begin(), chunk->noiseValues.end()));
        // writeNoiseValuesToFile(chunk->noiseValues, CHUNK_SIZE + 1, "assets/noise/" + chunkName + ".txt");
        weightNoise(chunk->noiseValues, CHUNK_SIZE + 1, chunkLocalPos, origin, 0.5f, planetoid->size);
        // Log the weighted noise values min and max
        // logger.logf("[weightNoise] Chunk (%d, %d, %d) weighted noise values min: %.4f, max: %.4f\n",
        //     cx, cy, cz,
        //     *std::min_element(chunk->noiseValues.begin(), chunk->noiseValues.end()),
        //     *std::max_element(chunk->noiseValues.begin(), chunk->noiseValues.end()));
        // Log the chunk generation time
        time_t genChunkEnd = clock();
        double genChunkTime = static_cast<double>(genChunkEnd - genChunk) / CLOCKS_PER_SEC;
        logger.logf("Chunk noise (%d, %d, %d) generated in %.2f seconds.\n", cx, cy, cz, genChunkTime);
        if (isAllBelowThreshold(chunk->noiseValues, 0.01f)) { // Lowered threshold
            planetoid->chunkChildren.emplace(Int3{cx, cy, cz}, std::move(chunk));
            planetoid->chunkChildren[Int3{cx, cy, cz}]->assignNeighborsAndNotify(planetoid->chunkChildren);
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
        // For each face, always use the lower chunk coordinate for the cache key
        // +X face: shared between (cx,cy,cz) and (cx+1,cy,cz), key is (cx,cy,cz,0)
        auto keyXPos = std::make_tuple(cx, cy, cz, 0);
        if (planetoid->sharedEdgeCaches.count(keyXPos) == 0) {
            // logger.logf("[EdgeCache] Creating edgeCacheXPos for chunk (%d,%d,%d) key (%d,%d,%d,0)\n", cx, cy, cz, cx, cy, cz);
            planetoid->sharedEdgeCaches[keyXPos] = std::vector<EdgeCacheEntry>(CHUNK_SIZE * CHUNK_SIZE * 12);
        } else {
            // logger.logf("[EdgeCache] Reusing edgeCacheXPos for chunk (%d,%d,%d) key (%d,%d,%d,0)\n", cx, cy, cz, cx, cy, cz);
        }
        edgeCacheXPos = &planetoid->sharedEdgeCaches[keyXPos];
        // -X face: shared between (cx-1,cy,cz) and (cx,cy,cz), key is (cx-1,cy,cz,0)
        auto keyXNeg = std::make_tuple(cx-1, cy, cz, 0);
        if (planetoid->sharedEdgeCaches.count(keyXNeg) == 0) {
            // logger.logf("[EdgeCache] Creating edgeCacheXNeg for chunk (%d,%d,%d) key (%d,%d,%d,0)\n", cx, cy, cz, cx-1, cy, cz);
            planetoid->sharedEdgeCaches[keyXNeg] = std::vector<EdgeCacheEntry>(CHUNK_SIZE * CHUNK_SIZE * 12);
        } else {
            // logger.logf("[EdgeCache] Reusing edgeCacheXNeg for chunk (%d,%d,%d) key (%d,%d,%d,0)\n", cx, cy, cz, cx-1, cy, cz);
        }
        edgeCacheXNeg = &planetoid->sharedEdgeCaches[keyXNeg];
        // +Y face: shared between (cx,cy,cz) and (cx,cy+1,cz), key is (cx,cy,cz,1)
        auto keyYPos = std::make_tuple(cx, cy, cz, 1);
        if (planetoid->sharedEdgeCaches.count(keyYPos) == 0) {
            // logger.logf("[EdgeCache] Creating edgeCacheYPos for chunk (%d,%d,%d) key (%d,%d,%d,1)\n", cx, cy, cz, cx, cy, cz);
            planetoid->sharedEdgeCaches[keyYPos] = std::vector<EdgeCacheEntry>(CHUNK_SIZE * CHUNK_SIZE * 12);
        } else {
            // logger.logf("[EdgeCache] Reusing edgeCacheYPos for chunk (%d,%d,%d) key (%d,%d,%d,1)\n", cx, cy, cz, cx, cy, cz);
        }
        edgeCacheYPos = &planetoid->sharedEdgeCaches[keyYPos];
        // -Y face: shared between (cx,cy-1,cz) and (cx,cy,cz), key is (cx,cy-1,cz,1)
        auto keyYNeg = std::make_tuple(cx, cy-1, cz, 1);
        if (planetoid->sharedEdgeCaches.count(keyYNeg) == 0) {
            // logger.logf("[EdgeCache] Creating edgeCacheYNeg for chunk (%d,%d,%d) key (%d,%d,%d,1)\n", cx, cy, cz, cx, cy-1, cz);
            planetoid->sharedEdgeCaches[keyYNeg] = std::vector<EdgeCacheEntry>(CHUNK_SIZE * CHUNK_SIZE * 12);
        } else {
            // logger.logf("[EdgeCache] Reusing edgeCacheYNeg for chunk (%d,%d,%d) key (%d,%d,%d,1)\n", cx, cy, cz, cx, cy-1, cz);
        }
        edgeCacheYNeg = &planetoid->sharedEdgeCaches[keyYNeg];
        // +Z face: shared between (cx,cy,cz) and (cx,cy,cz+1), key is (cx,cy,cz,2)
        auto keyZPos = std::make_tuple(cx, cy, cz, 2);
        if (planetoid->sharedEdgeCaches.count(keyZPos) == 0) {
            // logger.logf("[EdgeCache] Creating edgeCacheZPos for chunk (%d,%d,%d) key (%d,%d,%d,2)\n", cx, cy, cz, cx, cy, cz);
            planetoid->sharedEdgeCaches[keyZPos] = std::vector<EdgeCacheEntry>(CHUNK_SIZE * CHUNK_SIZE * 12);
        } else {
            // logger.logf("[EdgeCache] Reusing edgeCacheZPos for chunk (%d,%d,%d) key (%d,%d,%d,2)\n", cx, cy, cz, cx, cy, cz);
        }
        edgeCacheZPos = &planetoid->sharedEdgeCaches[keyZPos];
        // -Z face: shared between (cx,cy,cz-1) and (cx,cy,cz), key is (cx,cy,cz-1,2)
        auto keyZNeg = std::make_tuple(cx, cy, cz-1, 2);
        if (planetoid->sharedEdgeCaches.count(keyZNeg) == 0) {
            // logger.logf("[EdgeCache] Creating edgeCacheZNeg for chunk (%d,%d,%d) key (%d,%d,%d,2)\n", cx, cy, cz, cx, cy, cz-1);
            planetoid->sharedEdgeCaches[keyZNeg] = std::vector<EdgeCacheEntry>(CHUNK_SIZE * CHUNK_SIZE * 12);
        } else {
            // logger.logf("[EdgeCache] Reusing edgeCacheZNeg for chunk (%d,%d,%d) key (%d,%d,%d,2)\n", cx, cy, cz, cx, cy, cz-1);
        }
        edgeCacheZNeg = &planetoid->sharedEdgeCaches[keyZNeg];
        // For now, pass +X, +Y, +Z caches to marchCube (legacy interface)
        chunk->vertices.clear();
        chunk->indices.clear();
        chunk->vertices.reserve(CHUNK_SIZE * CHUNK_SIZE * CHUNK_SIZE * 3 / 2);
        chunk->indices.reserve(CHUNK_SIZE * CHUNK_SIZE * CHUNK_SIZE * 6);
        // time_t cubemarchStart = clock();
        

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
                        chunk->noiseValues, CHUNK_SIZE + 1,
                        chunk->vertices, chunk->indices,
                        &localEdgeCache, edgeCacheX, edgeCacheY, edgeCacheZ,
                        0.5f,
                        cx, cy, cz,
                        chunkLocalPos,
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
        Mesh mesh = { 0 };
        mesh.vertexCount = static_cast<int>(chunk->vertices.size());
        if (mesh.vertexCount != chunk->vertices.size()) {
            logger.logf("[ERROR] mesh.vertexCount (%d) != vertices.size() (%zu)\n", mesh.vertexCount, chunk->vertices.size());
        }
        mesh.vertices = new float[mesh.vertexCount * 3];
        for (size_t j = 0; j < chunk->vertices.size(); ++j) {
            const Vector3& v = chunk->vertices[j];
            mesh.vertices[j * 3 + 0] = v.x;
            mesh.vertices[j * 3 + 1] = v.y;
            mesh.vertices[j * 3 + 2] = v.z;
        }
        mesh.triangleCount = static_cast<int>(chunk->indices.size() / 3);
        mesh.indices = new unsigned short[chunk->indices.size()];
        for (size_t j = 0; j < chunk->indices.size(); ++j) {
            mesh.indices[j] = static_cast<unsigned short>(chunk->indices[j]);
        }
        // Allocate and zero normals
        mesh.normals = new float[mesh.vertexCount * 3]();

        if (mesh.vertexCount != chunk->vertices.size()) {
            logger.logf("[ERROR] mesh.vertexCount (%d) != vertices.size() (%zu)\n", mesh.vertexCount, chunk->vertices.size());
        }

        UploadMesh(&mesh, true); // Make mesh dynamic so normals can be updated later
        chunk->mesh = mesh;
        chunk->model = LoadModelFromMesh(chunk->mesh);
        chunk->model.materials[0].shader = world.lightingShader; // Use the lighting shader for the chunk
        logger.logf("[iterativeChunk] Chunk (%d, %d, %d) model position (%.2f, %.2f, %.2f)\n",
            cx, cy, cz, chunk->model.transform.m0, chunk->model.transform.m5, chunk->model.transform.m10);
        auto inserted = planetoid->chunkChildren.emplace(Int3{cx, cy, cz}, std::move(chunk));
        planetoid->chunkChildren[Int3{cx, cy, cz}]->assignNeighborsAndNotify(planetoid->chunkChildren);
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

    clock_t genStart = clock();
    // Use iterative chunk generation
    iterativeChunk(0, 0, 0, position, rotation, color, scale, noise, planetoid, world); 
    clock_t genEnd = clock();
    double genTime = static_cast<double>(genEnd - genStart) / CLOCKS_PER_SEC;
    logger.logf("Chunks generated in %.2f seconds at position (%f, %f, %f) with size %zu and scale %.2f.\n", 
        genTime, position.x, position.y, position.z, size, scale);

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
        size = 50; // Fixed size for testing

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