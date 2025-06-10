#include "scene.h"

#if defined(PLATFORM_DESKTOP)
    #define GLSL_VERSION            330
#endif


Scene::Scene(std::string name, bool isActive)
    : name(name), isActive(isActive),
      rootObject("root", "root", {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, WHITE, 1.0f)
{
    Mesh cube = GenMeshCube(1.0f, 1.0f, 1.0f);
    skybox = LoadModelFromMesh(cube);

    skybox.materials[0].shader = LoadShader(TextFormat("resources/skybox.vs", GLSL_VERSION),
                                            TextFormat("resources/skybox.fs", GLSL_VERSION));

    SetShaderValue(skybox.materials[0].shader, GetShaderLocation(skybox.materials[0].shader, "environmentMap"), (int[1]){ MATERIAL_MAP_CUBEMAP }, SHADER_UNIFORM_INT);
    SetShaderValue(skybox.materials[0].shader, GetShaderLocation(skybox.materials[0].shader, "doGamma"), (int[1]){0}, SHADER_UNIFORM_INT);
    SetShaderValue(skybox.materials[0].shader, GetShaderLocation(skybox.materials[0].shader, "vflipped"), (int[1]){0}, SHADER_UNIFORM_INT);

    // Load cubemap shader and setup required shader locations
    Shader shdrCubemap = LoadShader(TextFormat("resources/cubemap.vs", GLSL_VERSION),
                                    TextFormat("resources/cubemap.fs", GLSL_VERSION));

    SetShaderValue(shdrCubemap, GetShaderLocation(shdrCubemap, "equirectangularMap"), (int[1]){ 0 }, SHADER_UNIFORM_INT);

    Image img = LoadImage("resources/skybox.png");
    skybox.materials[0].maps[MATERIAL_MAP_CUBEMAP].texture = LoadTextureCubemap(img, CUBEMAP_LAYOUT_AUTO_DETECT);    // CUBEMAP_LAYOUT_PANORAMA
    UnloadImage(img);
    std::cout << "Skybox created: " << name << std::endl;

    // Initialize lighting
    lightingShader = LoadShader(TextFormat("resources/lighting.vs", GLSL_VERSION),
                     TextFormat("resources/lighting.fs", GLSL_VERSION));
    int ambientLoc = GetShaderLocation(lightingShader, "ambient");
    SetShaderValue(lightingShader, ambientLoc, (float[4]){ 0.1f, 0.1f, 0.1f, 1.0f }, SHADER_UNIFORM_VEC4);

    lights.push_back(CreateLight(LIGHT_DIRECTIONAL, (Vector3){ 0.0f, 0.5f, 0.0f }, (Vector3){ -14.0f, -0.4f, -0.5f }, WHITE, lightingShader));

    std::cout << "Scene created: " << name << std::endl;
}

void Scene::drawScene(int gamestate) {
    if (!isActive) return; // Skip drawing if the scene is not active
    BeginMode3D(camera);
    if(gamestate != 0 && gamestate != 2) { // If not in main menu or pause menu
        customUpdateCamera(&camera);
    }

    rlDisableBackfaceCulling();
    rlDisableDepthMask();
    DrawModel(skybox, (Vector3){0, 0, 0}, 1.0f, WHITE);
    rlEnableBackfaceCulling();
    rlEnableDepthMask();

    float cameraPos[3] = { camera.position.x, camera.position.y, camera.position.z };
    SetShaderValue(lightingShader, lightingShader.locs[SHADER_LOC_VECTOR_VIEW], cameraPos, SHADER_UNIFORM_VEC3);

    for (const auto& light : lights) { 
        UpdateLightValues(lightingShader, light);
    }

    BeginShaderMode(lightingShader);
    for (const auto& objPtr : objects) { objPtr->draw(); }
    for (const auto& objPtr : rootObject.children) { objPtr->draw(); }
    EndShaderMode();
    EndMode3D();
}

void Scene::drawUI(int gamestate) {
    if (!isActive) return; // Skip drawing if the scene is not active
    DrawFPS(10,10);
    std::string camPosStr = "Camera: X=" + std::to_string(camera.position.x) +
                            " Y=" + std::to_string(camera.position.y) +
                            " Z=" + std::to_string(camera.position.z);
    DrawText(camPosStr.c_str(), 10, 30, 20, GREEN);

    std::string camFacingVector = "Facing: X=" + std::to_string(camera.target.x) +
                                  " Y=" + std::to_string(camera.target.y) +
                                  " Z=" + std::to_string(camera.target.z);
    DrawText(camFacingVector.c_str(), 10, 50, 20, GREEN);
    for (const auto& objPtr : uiObjects) { objPtr->draw(); }
}

Scene::~Scene() {
    // Unload skybox model and its shader
    if (skybox.meshCount > 0) {
        UnloadModel(skybox);
        skybox.meshCount = 0;
    }
    // Unload lighting shader
    if (lightingShader.id > 0) {
        UnloadShader(lightingShader);
        lightingShader.id = 0;
    }
    // If you have any other dynamically loaded resources, unload them here
}


