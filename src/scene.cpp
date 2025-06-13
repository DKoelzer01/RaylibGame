#include "scene.h"
#include <raylib.h>
#include <rlgl.h>

#if defined(PLATFORM_DESKTOP)
    #define GLSL_VERSION            330
#endif

// --- Shadow mapping resources ---
#define SHADOW_MAP_SIZE 2048
Matrix lightSpaceMatrix;

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

    // --- Shadow map FBO/texture ---
    shadowMap = LoadRenderTexture(SHADOW_MAP_SIZE, SHADOW_MAP_SIZE);
    // Depth buffer is automatically handled by LoadRenderTexture in raylib
    // --- Load depth-only shader for shadow mapping ---
    depthShader = LoadShader("resources/depth.vs", "resources/depth.fs");

    std::cout << "Scene created: " << name << std::endl;
}

void Scene::updateAllChunkShaders() {
    for (const auto& objPtr : objects) {
        auto planetoid = dynamic_cast<Planetoid*>(objPtr.get());
        if (planetoid) {
            for (auto& chunkPair : planetoid->chunkChildren) {
                chunkPair.second->chunk.model.materials[0].shader = lightingShader;
            }
        }
    }
}

void Scene::drawScene(int gamestate) {
    if (!isActive) return; // Skip drawing if the scene is not active

    // --- Compute light view/projection matrix for shadow mapping ---
    // Make the directional light follow the camera/player
    Vector3 cameraPosVec = camera.position;
    // Choose a sun direction (normalized)
    Vector3 sunDir = Vector3Normalize((Vector3){ -14.0f, 0.0f, 0.0f }); // Example: from above and behind
    float sunDistance = 500.0f; // Larger offset for debugging
    Vector3 lightPos = Vector3Add(cameraPosVec, Vector3Scale(sunDir, sunDistance));
    Vector3 lightTarget = cameraPosVec;
    // Update the first light's position/target
    if (!lights.empty()) {
        lights[0].position = lightPos;
        lights[0].target = lightTarget;
    }
    Matrix lightView = MatrixLookAt(lightPos, lightTarget, (Vector3){0,1,0});
    float orthoSize = 3000.0f; // Larger ortho size for debugging
    Matrix lightProj = MatrixOrtho(-orthoSize, orthoSize, -orthoSize, orthoSize, 1.0f, 1000.0f);
    lightSpaceMatrix = MatrixMultiply(lightProj, lightView);    

    // logger.logf("lightSpaceMatrix: \n"
    //             "  m0: %f, m1: %f, m2: %f, m3: %f\n"
    //             "  m4: %f, m5: %f, m6: %f, m7: %f\n"
    //             "  m8: %f, m9: %f, m10: %f, m11: %f\n"
    //             "  m12: %f, m13: %f, m14: %f, m15: %f\n",
    //             lightSpaceMatrix.m0, lightSpaceMatrix.m1, lightSpaceMatrix.m2, lightSpaceMatrix.m3,
    //             lightSpaceMatrix.m4, lightSpaceMatrix.m5, lightSpaceMatrix.m6, lightSpaceMatrix.m7,
    //             lightSpaceMatrix.m8, lightSpaceMatrix.m9, lightSpaceMatrix.m10, lightSpaceMatrix.m11,
    //             lightSpaceMatrix.m12, lightSpaceMatrix.m13, lightSpaceMatrix.m14, lightSpaceMatrix.m15);
    // --- Shadow map render pass ---
    BeginTextureMode(shadowMap);
    ClearBackground(BLACK);
    // Set camera to lightView/lightProj, render depth only
    for (const auto& objPtr : objects) { objPtr->drawDepthOnly(lightSpaceMatrix, &depthShader); }
    for (const auto& objPtr : rootObject.children) { objPtr->drawDepthOnly(lightSpaceMatrix, &depthShader); }
    EndTextureMode();

    int shadowMapLoc = GetShaderLocation(lightingShader, "shadowMap");
    SetShaderValueTexture(lightingShader, shadowMapLoc, shadowMap.texture);
    int lightSpaceLoc = GetShaderLocation(lightingShader, "lightSpaceMatrix");
    SetShaderValueMatrix(lightingShader, lightSpaceLoc, lightSpaceMatrix);

    // Calculate projection and view matrices
    Matrix proj = GetCameraProjectionMatrix(&camera, CAMERA_PERSPECTIVE); // or MatrixPerspective(...)
    Matrix view = GetCameraMatrix(camera);

    // MVP = Projection * View
    Matrix mvp = MatrixMultiply(proj, view);

    // Set the uniform on your shader
    int mvpLoc = GetShaderLocation(lightingShader, "mvp");
    SetShaderValueMatrix(lightingShader, mvpLoc, mvp);
    
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

    for (const auto& light : lights) { 
    if (!light.enabled) continue; // Skip disabled lights
        // Draw light source as a sphere at the light position
        // logger.logf("Drawing light: %s at position (%f, %f, %f) color (%d, %d, %d)\n",
        //     light.type == LIGHT_DIRECTIONAL ? "Directional" : "Point",
        //     light.position.x, light.position.y, light.position.z,
        //     light.color.r, light.color.g, light.color.b);
        DrawSphere(light.position, 20.0f, ColorAlpha(RED, 0.5f));
        UpdateLightValues(lightingShader, light);
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
    EndMode3D();
    // DrawTextureRec(shadowMap.texture, (Rectangle){0, 0, shadowMap.texture.width, -shadowMap.texture.height}, (Vector2){10, 10}, WHITE);
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


