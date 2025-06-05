#include "world.h"

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
    // Generate a simple world with a ground plane and some objects
    printf("Generating world scene...\n");

    //Generate planetoids
    std::random_device rd; // Obtain a seed from the system
    std::mt19937 gen(rd()); // Initialize the Mersenne Twister engine with the seed
    float randScale = 1000.0f;
    std::uniform_int_distribution<> distrib(1, randScale); // Define a uniform distribution between 1 and 10

    // for(int i = 0; i < distrib(gen); ++i) {
        Vector3 position = { static_cast<float>(distrib(gen)/10), static_cast<float>(distrib(gen)/10), static_cast<float>(distrib(gen)/10) };
        Vector3 rotation = { static_cast<float>(distrib(gen)/randScale * 360), static_cast<float>(distrib(gen)/randScale * 360), static_cast<float>(distrib(gen)/randScale * 360) };
        Color color = { static_cast<unsigned char>(distrib(gen)/randScale * 255), static_cast<unsigned char>(distrib(gen)/randScale * 255), 
                        static_cast<unsigned char>(distrib(gen)/randScale * 255), 255 };
        float size = static_cast<int>((distrib(gen)/randScale)*50 + 1); // Ensure size is at least 1.0
        size = 50.0f; // Fixed size for testing
        float scale = 1.0f;

        SimplexNoise* noise = new SimplexNoise(0.01f, 4.0f, 2.0f, 0.5f); // Create a new instance of SimplexNoise
        if (!noise) {
            std::cerr << "Failed to create SimplexNoise instance." << std::endl;
            return;
        }

        int i = 0;

        printf("Generating noise for planetoid %d at position (%f, %f, %f) with size %f\n", i, position.x, position.y, position.z, size);

        // Generate noise value for the position
        std::vector<float> noiseValues(size * size * size);
        for (int z = 0; z < size; ++z) {
            for (int y = 0; y < size; ++y) {
                for (int x = 0; x < size; ++x) {
                    int idx = x + y * size + z * size * size;
                    noiseValues[idx] = noise->fractal(4, position.x + static_cast<float>(x), position.y + static_cast<float>(y), position.z + static_cast<float>(z));
                }
            }
        }

        printf("Weighing noise for planetoid %d at position (%f, %f, %f) with size %f\n", i, position.x, position.y, position.z, size);

        weightNoise(noiseValues, pow(size,3), 3.5f);
        
        for(int z = 0; z < size; ++z) {
            for (int y = 0; y < size; ++y) {
                for (int x = 0; x < size; ++x) {
                    int idx = x + y * size + z * size * size;
                    float noiseVal = noiseValues[idx];
                    if (noiseVal > 0.02f) { // Threshold to determine if a voxel should be created
                        Vector3 voxelPosition = { position.x + static_cast<float>(x), position.y + static_cast<float>(y), position.z + static_cast<float>(z) };
                        //printf("Creating voxel at (%f, %f, %f) with noise value %f\n", voxelPosition.x, voxelPosition.y, voxelPosition.z, noiseVal);
                        world.objects.push_back(std::make_unique<GameObject>("cube", "voxel" + std::to_string(i) + "_" + std::to_string(x) + "_" + std::to_string(y) + "_" + std::to_string(z),voxelPosition, rotation, color, scale));
                    }
                }
            }
        }

        // world.objects.push_back(std::make_unique<GameObject>("sphere", "planetoid" + std::to_string(i), position, rotation, color, scale));
    // }


    printf("World scene generated with %zu objects.\n", world.objects.size());
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