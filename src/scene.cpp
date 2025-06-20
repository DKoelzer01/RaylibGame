#include "scene.h"
#include <raylib.h>
#include <rlgl.h>

#if defined(PLATFORM_DESKTOP)
    #define GLSL_VERSION            330
#endif

// --- Shadow mapping resources ---
#define SHADOW_MAP_SIZE 2048
Matrix lightSpaceMatrix;
Matrix cameraView;
Matrix cameraProj;
Texture2D shadowMapTexture; // Texture for shadow map

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
    printf("Lighting shader loaded: %d\n", lightingShader.id);
    int ambientLoc = GetShaderLocation(lightingShader, "ambient");
    SetShaderValue(lightingShader, ambientLoc, (float[4]){ 0.1f, 0.1f, 0.1f, 1.0f }, SHADER_UNIFORM_VEC4);

    lights.push_back(CreateLight(LIGHT_DIRECTIONAL, (Vector3){ 0.0f, 0.5f, 0.0f }, (Vector3){ -14.0f, -0.4f, -0.5f }, WHITE, lightingShader));

    // --- Shadow map FBO/texture ---
    shadowMap = LoadRenderTexture(SHADOW_MAP_SIZE, SHADOW_MAP_SIZE);
    SetTextureFilter(shadowMap.texture, TEXTURE_FILTER_POINT); // Ensure shadow map uses point sampling
    // Depth buffer is automatically handled by LoadRenderTexture in raylib
    // --- Load depth-only shader for shadow mapping ---
    depthShader = LoadShader("resources/depth.vs", "resources/depth.fs");
    printf("Depth shader loaded: %s\n", depthShader.id > 0 ? "Success" : "Failed");

    std::cout << "Scene created: " << name << std::endl;
}

void Scene::updateAllChunkShaders() {
    for (const auto& objPtr : objects) {
        auto planetoid = dynamic_cast<Planetoid*>(objPtr.get());
        if (planetoid) {
            for (auto& chunkPair : planetoid->chunkChildren) {
                chunkPair.second->model.materials[0].shader = lightingShader;
            }
        }
    }
}

void Scene::drawScene(int gamestate) {
    if (!isActive) return; // Skip drawing if the scene is not active

    // This is for debugging purposes only; in production, shaders should be loaded once and reused
    lightingShader = LoadShader(TextFormat("resources/lighting.vs", GLSL_VERSION),
                     TextFormat("resources/lighting.fs", GLSL_VERSION));
    depthShader = LoadShader("resources/depth.vs", "resources/depth.fs");
    // --- Compute light view/projection matrix for shadow mapping ---
    // Make the directional light follow the camera/player
    Vector3 cameraPosVec = camera.position;
    // cameraPosVec = {0.0f, 0.0f, 0.0f}; // DEBUG: Reset camera position to origin for testing
    // Choose a sun direction (normalized)
    Vector3 sunDir = Vector3Normalize((Vector3){ 26.0f, 0.0f, 0.0f }); // Example: from above and behind
    float sunDistance = 500.0f; // Larger offset for debugging
    Vector3 lightPos = Vector3Add(cameraPosVec, Vector3Scale(sunDir, sunDistance));
    Vector3 lightTarget = cameraPosVec;

    // Debug: Set light position to fixed position for testing
    lightPos = (Vector3){ sunDistance, 0.0f, 0.0f }; // DEBUG: Set light position to fixed position for testing
    // Debug set lightPos to orbit around the origin over time
    // lightPos = (Vector3){ sin(GetTime()) * sunDistance, 0.0f, -cos(GetTime()) * sunDistance };
    lightTarget = (Vector3){ 0.0f, 0.0f, 0.0f }; // DEBUG: Set light target to origin for testing
    // Update the first light's position/target
    if (!lights.empty()) {
        lights[0].position = lightPos;
        lights[0].target = lightTarget;
    }
    Camera lightCamera = {0};
    lightCamera.position = lightPos;
    lightCamera.target = lightTarget;
    lightCamera.up = (Vector3){0,1,0};
    lightCamera.fovy = 90.0f; // or appropriate value
    lightCamera.projection = CAMERA_ORTHOGRAPHIC;

    Matrix lightView = MatrixLookAt(lightPos, lightTarget, lightCamera.up);
    float orthoSize = sunDistance/2.0f; // Size of the orthographic projection box
    float nearPlane = 100.0f; // Near plane for shadow mapping
    float farPlane = orthoSize + 100.0f; // Ensure far plane is beyond the light position
    Matrix lightProj = MatrixOrtho(-orthoSize, orthoSize, -orthoSize, orthoSize, nearPlane, farPlane); 
    lightSpaceMatrix = MatrixMultiply(lightProj, lightView);   

    // --- Shadow pass ---
    BeginTextureMode(shadowMap);
    ClearBackground(BLACK);
    BeginMode3D(lightCamera);
    // Set global camera matrices for shadow pass
    extern Matrix cameraProj, cameraView;
    cameraProj = lightProj;
    cameraView = lightView;
    for (const auto& objPtr : objects) { objPtr->drawDepthOnly(lightSpaceMatrix, &depthShader); }
    for (const auto& objPtr : rootObject.children) { objPtr->drawDepthOnly(lightSpaceMatrix, &depthShader); }
    EndMode3D();
    EndTextureMode();
    shadowMapTexture = shadowMap.texture; // Update the global shadow map texture

    // --- Main pass ---
    Matrix proj = GetCameraProjectionMatrix(&camera, CAMERA_PERSPECTIVE);
    Matrix view = GetCameraMatrix(camera);
    
    extern Matrix cameraProj, cameraView;
    cameraProj = proj;
    cameraView = view;

    int lightSpaceLoc = GetShaderLocation(lightingShader, "lightSpaceMatrix");
    SetShaderValueMatrix(lightingShader, lightSpaceLoc, lightSpaceMatrix);

    // --- Main scene render ---
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
    int viewPosLoc = GetShaderLocation(lightingShader, "viewPos");
    if (viewPosLoc != -1) {
        SetShaderValue(lightingShader, viewPosLoc, cameraPos, SHADER_UNIFORM_VEC3);
    }

    int lightCount = 0;
    for (const auto& light : lights) { 
    if (!light.enabled) continue; // Skip disabled lights
        // Draw light source as a sphere at the light position
        // logger.logf("Drawing light: %s at position (%f, %f, %f) color (%d, %d, %d)\n",
        //     light.type == LIGHT_DIRECTIONAL ? "Directional" : "Point",
        //     light.position.x, light.position.y, light.position.z,
        //     light.color.r, light.color.g, light.color.b);
        DrawSphere(light.position, 20.0f, ColorAlpha(RED, 0.5f));
        UpdateLightValues(lightingShader, light, lightCount); // Update light values in shader
        lightCount++;
    }
    

    BeginShaderMode(lightingShader);
    // logger.logf("Drawing %zu objects\n", objects.size());
    for (const auto& objPtr : objects) {
        if (!objPtr) continue;
        // logger.logf("[Scene] Drawing object: %s at ptr %p\n", objPtr->name.c_str(), objPtr.get());
        objPtr->draw(&lightingShader);
    }
    for (const auto& objPtr : rootObject.children) {
        if(!objPtr) continue;
        // logger.logf("[Scene] Drawing root object: %s at ptr %p\n", objPtr->name.c_str(), objPtr.get());
        objPtr->draw(&lightingShader); 
    }
    EndShaderMode();
    // --- Debug: Draw orthographic projection box ---
    DrawCubeWires(lightTarget, orthoSize * 2, orthoSize * 2, farPlane - nearPlane, RED);
    EndMode3D();

    // --- Debug: Draw shadow map as quad ---
    DrawTexturePro(
        shadowMap.texture,
        (Rectangle){ 0, 0, (float)shadowMap.texture.width, -(float)shadowMap.texture.height },
        (Rectangle){ 0, 60, 400, 400 }, // Draw in a 400x400 box at top-left
        (Vector2){ 0, 0 },
        0.0f,
        WHITE
    );
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
    for (const auto& objPtr : uiObjects) { objPtr->draw(&lightingShader); }
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
    // Unload shadow map resources
    UnloadRenderTexture(shadowMap);
    // Unload depth shader
    if (depthShader.id > 0) {
        UnloadShader(depthShader);
        depthShader.id = 0;
    }
    // If you have any other dynamically loaded resources, unload them here
}


