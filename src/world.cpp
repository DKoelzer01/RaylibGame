#include "world.h"
#include <unordered_map>
#include <tuple>

// Helper for 3D integer coordinates as chunk keys
namespace std {
    template<>
    struct hash<Int3> {
        size_t operator()(const Int3& k) const {
            return ((hash<int>()(k.x) ^ (hash<int>()(k.y) << 1)) >> 1) ^ (hash<int>()(k.z) << 2);
        }
    };
}

constexpr size_t CHUNK_SIZE = 64;

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

void weightNoise(std::vector<float>& noiseValues, int size, float weight) {
    if (size <= 0) return;

    int dim = static_cast<int>(std::cbrt(size));
    if (dim * dim * dim != size) return; // Ensure it's a perfect cube

    float center = (dim - 1) / 2.0f;
    for (int z = 0; z < dim; ++z) {
        for (int y = 0; y < dim; ++y) {
            for (int x = 0; x < dim; ++x) {
                int idx = x + y * dim + z * dim * dim;
                float dx = x - center;
                float dy = y - center;
                float dz = z - center;
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                float maxDist = std::sqrt(3 * center * center);
                float scale = 2.0f - (dist / maxDist) * weight;
                if (scale < 0.0f) scale = 0.0f;
                noiseValues[idx] *= scale;
            }
        }
    }
}

void generateWorld(Scene& world) {
    printf("Generating world scene...\n");
    clock_t worldGenStart = clock();

    //Delete old planetoids
    //TODO: Implement a proper cleanup of old planetoids

    float randScale = 100000.0f; // Scale for random number generation
    float genRange = 500.0f; // Range for planetoid generation
    float minSize = 50.0f; // Minimum size of planetoids
    float maxSize = 250.0f; // Maximum size of planetoids
    // Initialize the global edgeCacheFlat with -1
    edgeCacheFlat.assign(static_cast<size_t>((maxSize-1)*(maxSize-1)*(maxSize-1)*12), -1); // Corrected edge cache size
    printf("Generated %d planetoids.\n", generatePlanetoids(randScale, world, genRange, minSize, maxSize)); // Generate planetoids in the world

    clock_t worldGenEnd = clock();
    double worldGenTime = static_cast<double>(worldGenEnd - worldGenStart) / CLOCKS_PER_SEC;
    printf("World generation completed in %.2f seconds.\n", worldGenTime);
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
    printf("Generating %d planetoids...\n", numPlanetoids);

    for(int i = 0; i < numPlanetoids; ++i) {
        // Generate random position, rotation, color, size, and scale for the planetoid

        // Position is centered around (0,0,0) with a range of genRange
        Vector3 position = {    static_cast<float>(((distrib(gen)/randScale)*genRange)-(genRange/2)), 
                                static_cast<float>(((distrib(gen)/randScale)*genRange)-(genRange/2)), 
                                static_cast<float>(((distrib(gen)/randScale)*genRange)-(genRange/2))};

        //DEBUG: Fixed position for testing
        position = { 0.0f, 0.0f, 0.0f }; // Centered position for testing

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
        size = 300; // Fixed size for testing

        // Initialize the global edgeCacheFlat for this planetoid's size
        edgeCacheFlat.assign((size-1)*(size-1)*(size-1)*12, -1); // Corrected edge cache size

        // Scale is a fixed value for now, can be adjusted later
        float scale = 1.0f;

        printf("\nGenerating planetoid %d...\n", i);
        clock_t planetoidGenStart = clock();

        float frequencyNoise = static_cast<float>(distrib(gen)/randScale)-0.5f; // Random frequency noise to add some variation from -0.5 to 0.5

        SimplexNoise* noise = new SimplexNoise((1.5f+frequencyNoise)/static_cast<float>(size), 2.0f, 2.0f, 0.5f); // Create a new instance of SimplexNoise
        if (!noise) {
            std::cerr << "Failed to create SimplexNoise instance." << std::endl;
            continue; // Skip this planetoid if noise generation fails;
        }

        printf("Generating noise for planetoid %d at position (%f, %f, %f) with size %zu and frequency %f\n", i, position.x, position.y, position.z, size, 1.5f+frequencyNoise);
        clock_t noiseGenStart = clock();

        // Generate noise value for the position
        std::vector<float> noiseValues(size * size * size);
        for (size_t z = 0; z < size; ++z) {
            for (size_t y = 0; y < size; ++y) {
                for (size_t x = 0; x < size; ++x) {
                    size_t idx = x + y * size + z * size * size;
                    noiseValues[idx] = (noise->fractal(6, position.x + static_cast<float>(x), position.y + static_cast<float>(y), position.z + static_cast<float>(z))+1.0f); // Ensure noise values are positive
                }
            }
        }
        clock_t noiseGenEnd = clock();
        double noiseGenTime = static_cast<double>(noiseGenEnd - noiseGenStart) / CLOCKS_PER_SEC;
        printf("Noise generation completed in %.2f seconds for planetoid %d.\n", noiseGenTime, i);

        printf("Weighing noise for planetoid %d at position (%f, %f, %f) with size %f\n", i, position.x, position.y, position.z, size);
        clock_t weightNoiseStart = clock();

        // Weigh the noise values based on their distance from the center of the planetoid
        // This will create a spherical shape for the planetoid
        weightNoise(noiseValues, pow(size,3), 4.0f); // Weight factor to control the shape of the planetoid

        clock_t weightNoiseEnd = clock();
        double weightNoiseTime = static_cast<double>(weightNoiseEnd - weightNoiseStart) / CLOCKS_PER_SEC;
        printf("Noise weighing completed in %.2f seconds for planetoid %d.\n", weightNoiseTime, i);

        // Marching cubes algorithm to generate the mesh
        printf("Marching cubes for planetoid %d at position (%f, %f, %f) with size %f\n", i, position.x, position.y, position.z, size);
        clock_t marchCubeStart = clock();
        std::vector<Vector3> vertices;
        std::vector<size_t> indices;
        std::unordered_map<EdgeKey, size_t> edgeCache;

        for(size_t z = 0; z < size - 1; ++z) {
            for (size_t y = 0; y < size - 1; ++y) {
                for (size_t x = 0; x < size - 1; ++x) {
                    //if( noiseValues[x + y * size + z * size * size] > 0.5f) world.objects.push_back(std::make_unique<GameObject>("cube", "testCube", Vector3{static_cast<float>(x+100), static_cast<float>(y), static_cast<float>(z)}, rotation, color, scale));
                    marchCube(x, y, z, noiseValues, static_cast<size_t>(size), vertices, indices, edgeCache, 0.5f);
                }
            }
        }

        printf("Vertices generated: %zu, Indices generated: %zu\n", vertices.size(), indices.size());
        clock_t marchCubeEnd = clock();
        double marchCubeTime = static_cast<double>(marchCubeEnd - marchCubeStart) / CLOCKS_PER_SEC; 
        printf("Marching cubes completed in %.2f seconds for planetoid %d.\n", marchCubeTime, i);


        printf("Creating model for planetoid %d at position (%f, %f, %f) with size %f\n", i, position.x, position.y, position.z, size);
        clock_t modelGenStart = clock();
        // Create a mesh from the generated vertices and indices
        if (vertices.empty() || indices.empty()) {
            std::cerr << "No vertices or indices generated for planetoid " << i << ". Skipping model generation." << std::endl;
            continue;   // Skip this planetoid if no vertices or indices were generated
        }
        Mesh modelMesh = { 0 };
        // Each triangle is 3 vertices, so flatten the triangles into a vertex array
        modelMesh.triangleCount = static_cast<int>(indices.size() / 3);
        modelMesh.vertexCount = static_cast<int>(indices.size());
        modelMesh.vertices = new float[modelMesh.vertexCount * 3];
        for (size_t j = 0; j < indices.size(); ++j) {
            const Vector3& v = vertices[indices[j]];
            modelMesh.vertices[j * 3 + 0] = v.x;
            modelMesh.vertices[j * 3 + 1] = v.y;
            modelMesh.vertices[j * 3 + 2] = v.z;
        }
        // No indices, so set to nullptr
        modelMesh.indices = nullptr;
        UploadMesh(&modelMesh, false);
        Model model = LoadModelFromMesh(modelMesh);

        clock_t modelGenEnd = clock();
        double modelGenTime = static_cast<double>(modelGenEnd - modelGenStart) / CLOCKS_PER_SEC;
        printf("Model for planetoid %d created in %.2f seconds.\n", i, modelGenTime);


        clock_t modelFileStart = clock();
        printf("Saving model for planetoid %d at position (%f, %f, %f) with size %f\n", i, position.x, position.y, position.z, size);
        std::ofstream modelFile;
        std::string modelPath = "assets/models/planetoid/planetoid" + std::to_string(i) + ".obj";
        modelFile.open(modelPath);
        if (!modelFile.is_open()) {
            std::cerr << "Failed to open model file for writing." << std::endl;
        } else {
            modelFile << "#planetoid" << i << "\n";
            // Write vertices
            for (const auto& v : vertices) {
                modelFile << "v " << v.x << " " << v.y << " " << v.z << "\n";
            }
            // Write faces using indices (OBJ uses 1-based indexing)
            for (size_t j = 0; j + 2 < indices.size(); j += 3) {
                modelFile   << "f "
                            << (indices[j] + 1) << " "
                            << (indices[j + 1] + 1) << " "
                            << (indices[j + 2] + 1) << "\n";
            }
            modelFile.close();
        }
        clock_t modelFileEnd = clock();
        double modelFileTime = static_cast<double>(modelFileEnd - modelFileStart) / CLOCKS_PER_SEC;
        printf("Model for planetoid %d saved in %.2f seconds to %s\n", i, modelFileTime, modelPath.c_str());


        world.objects.push_back(std::make_unique<GameObject>("gameobject", "planetoid" + std::to_string(i), position, rotation, color, scale, model));
        clock_t planetoidGenEnd = clock();
        double planetoidGenTime = static_cast<double>(planetoidGenEnd - planetoidGenStart) / CLOCKS_PER_SEC;
        printf("Planetoid %d generated in %.2f seconds.\n\n", i, planetoidGenTime);
        delete noise; // Clean up the noise instance
        planetoidsGenerated++;
    }
    return planetoidsGenerated; // Return the number of planetoids generated
}


void loadWorld(Scene& world) {
    // Load the world scene from a file or initialize it
    printf("Loading world scene...\n");
    std::ifstream file("assets/scenes/world");
    if (!file.is_open()) {
        std::cerr << "Failed to open world scene file." << std::endl;
        return;
    }

    std::string line;
    while (getline(file, line)) {
        if (line.empty()) continue; // Skip empty lines

        printf("Processing line: %s\n", line.c_str());
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
                printf("Creating GameObject: %s\n", name.c_str());
                std::string modelPath = "assets/models/" + name + "/" + name + ".obj";
                std::string texturePath = "assets/models/" + name + "/diffuse.png";
                Model model;
                // Model model = LoadModelFromMesh(GenMeshCube(1.0f, 1.0f, 1.0f)); // Default model if loading fails
                if (!FileExists(modelPath.c_str())) {
                    std::cerr << "Model file does not exist: " << modelPath << std::endl;
                    continue;
                } else {
                    model = LoadModel(modelPath.c_str());
                    printf("meshCount: %d, materialCount: %d, boneCount: %d\n", model.meshCount, model.materialCount, model.boneCount);
                    printf("bones: %p, bindPose: %p\n", model.bones, model.bindPose);
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
                    printf("Texture file does not exist or materials/maps not allocated: %s\n", texturePath.c_str());
                }

                world.objects.push_back(std::make_unique<GameObject>(type, name, position, rotation, color, scale, model));
                printf("GameObject created: %s\n", name.c_str());
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