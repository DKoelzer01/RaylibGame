#include "scene.h"

#if defined(PLATFORM_DESKTOP)
    #define GLSL_VERSION            330
#endif


Scene::Scene(std::string name, bool isActive) : name(name), isActive(isActive) { 
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

    std::cout << "Scene created: " << name << std::endl;
}

void Scene::drawScene(int gamestate) {
    if (!isActive) return; // Skip drawing if the scene is not active
    BeginMode3D(camera);
    if(gamestate != 0 && gamestate != 2) { // If not in main menu or pause menu
        UpdateCamera(&camera, CAMERA_FREE);
    }
    
    rlDisableBackfaceCulling();
    rlDisableDepthMask();
    DrawModel(skybox, (Vector3){0, 0, 0}, 1.0f, WHITE);
    rlEnableBackfaceCulling();
    rlEnableDepthMask();

    for (const auto& objPtr : objects) { objPtr->draw(); }
    EndMode3D();
}

void Scene::drawUI(int gamestate) {
    if (!isActive) return; // Skip drawing if the scene is not active
    for (const auto& objPtr : uiObjects) { objPtr->draw(); }
}


