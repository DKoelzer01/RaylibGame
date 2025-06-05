#include "world.h"

void worldHandler(Scene& world) {
    // Initialize the camera for the world scene
    world.camera.position = { 0.0f, 10.0f, 10.0f };
    world.camera.target = { 0.0f, 0.0f, 0.0f };
    world.camera.up = { 0.0f, 1.0f, 0.0f };
    world.camera.fovy = 45.0f;
    world.camera.projection = CAMERA_PERSPECTIVE;

    loadWorld(world);
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
            // Process each cell in the line
            std::cout << "Cell: " << cell << std::endl;
            result.push_back(cell);
        }
        if (result.size() == 12) {
            std::string type = result[0];
            Vector3 position = { std::stof(result[1]), std::stof(result[2]), std::stof(result[3]) };
            Vector3 rotation = { std::stof(result[4]), std::stof(result[5]), std::stof(result[6]) };
            Color color = { static_cast<unsigned char>(std::stoi(result[7])),
                            static_cast<unsigned char>(std::stoi(result[8])),
                            static_cast<unsigned char>(std::stoi(result[9])),
                            static_cast<unsigned char>(std::stoi(result[10])) };
            float scale = std::stof(result[11]);

            world.objects.push_back(std::make_unique<GameObject>(type, position, rotation, color, scale));
        } else {
            std::cerr << "Invalid line format: " << line << std::endl;
            continue; // Skip invalid lines
        }
    }
    file.close();
}