#include "scene.h"
#include "utils/logger.h"
#include "profiler.h"
#include "worldtypes.h"
#include <raylib.h>
#include <rlgl.h>

#define SHADOW_MAP_SIZE 2048

// ---------------------------------------------------------------------------------------------
// HotShader
// ---------------------------------------------------------------------------------------------

bool HotShader::reloadIfChanged(bool force) {
    const long vt = GetFileModTime(vsPath.c_str());
    const long ft = GetFileModTime(fsPath.c_str());
    if (!force && vt == vsTime && ft == fsTime) return false;
    vsTime = vt;
    fsTime = ft;

    Shader s = LoadShader(vsPath.c_str(), fsPath.c_str());
    // On a compile/link error raylib hands back its default shader instead.
    if (!IsShaderValid(s) || s.id == rlGetShaderIdDefault()) {
        logger.logf("[ERROR] Shader %s / %s failed to compile; keeping previous version\n", vsPath.c_str(), fsPath.c_str());
        TraceLog(LOG_WARNING, "Shader %s failed to compile; keeping previous version", fsPath.c_str());
        return false;
    }
    unload();
    shader = s;
    TraceLog(LOG_WARNING, "Shader loaded: %s (id %u)", fsPath.c_str(), shader.id);
    return true;
}

void HotShader::unload() {
    if (shader.id > 0 && shader.id != rlGetShaderIdDefault()) UnloadShader(shader);
    shader = Shader{};
}

// ---------------------------------------------------------------------------------------------
// Shadow map: depth-texture-only framebuffer (same approach as raylib's shaders_shadowmap example).
// A real depth texture gives 24-bit precision instead of the 8 bits a colour channel would.
// ---------------------------------------------------------------------------------------------

static RenderTexture2D LoadShadowmapRenderTexture(int width, int height) {
    RenderTexture2D target = { 0 };
    target.id = rlLoadFramebuffer();
    target.texture.width = width;
    target.texture.height = height;
    if (target.id > 0) {
        rlEnableFramebuffer(target.id);
        target.depth.id = rlLoadTextureDepth(width, height, false);
        target.depth.width = width;
        target.depth.height = height;
        target.depth.format = 19; // DEPTH_COMPONENT_24BIT
        target.depth.mipmaps = 1;
        rlFramebufferAttach(target.id, target.depth.id, RL_ATTACHMENT_DEPTH, RL_ATTACHMENT_TEXTURE2D, 0);
        if (!rlFramebufferComplete(target.id)) TraceLog(LOG_WARNING, "Shadow map FBO incomplete");
        rlDisableFramebuffer();
    }
    return target;
}

static void UnloadShadowmapRenderTexture(RenderTexture2D target) {
    if (target.id > 0) rlUnloadFramebuffer(target.id); // also deletes attached depth texture
}

// ---------------------------------------------------------------------------------------------
// Scene
// ---------------------------------------------------------------------------------------------

Scene::Scene(std::string name, bool isActive)
    : name(name), isActive(isActive),
      rootObject("root", "root", {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, WHITE, 1.0f)
{
    // Skybox
    Mesh cube = GenMeshCube(1.0f, 1.0f, 1.0f);
    skybox = LoadModelFromMesh(cube);
    skybox.materials[0].shader = LoadShader("resources/skybox.vs", "resources/skybox.fs");
    Shader& sky = skybox.materials[0].shader;
    SetShaderValue(sky, GetShaderLocation(sky, "environmentMap"), (int[1]){ MATERIAL_MAP_CUBEMAP }, SHADER_UNIFORM_INT);
    SetShaderValue(sky, GetShaderLocation(sky, "doGamma"), (int[1]){ 0 }, SHADER_UNIFORM_INT);
    SetShaderValue(sky, GetShaderLocation(sky, "vflipped"), (int[1]){ 0 }, SHADER_UNIFORM_INT);
    Image img = LoadImage("resources/skybox.png");
    skybox.materials[0].maps[MATERIAL_MAP_CUBEMAP].texture = LoadTextureCubemap(img, CUBEMAP_LAYOUT_AUTO_DETECT);
    UnloadImage(img);

    // Materials. LoadMaterialDefault allocates the maps array once; we own it for the scene's lifetime.
    terrainMaterial = LoadMaterialDefault();
    depthMaterial = LoadMaterialDefault();

    shadowMap = LoadShadowmapRenderTexture(SHADOW_MAP_SIZE, SHADOW_MAP_SIZE);
    terrainMaterial.maps[MATERIAL_MAP_BRDF].texture = shadowMap.depth;

    depth.reloadIfChanged(true);
    depthMaterial.shader = depth.shader;
    lighting.reloadIfChanged(true);

    lights.push_back(CreateLight(LIGHT_DIRECTIONAL, Vector3Scale(sunDirection, sunDistance), Vector3Zero(), WHITE, lighting.shader));
    applyLightingShader();

    std::cout << "Scene created: " << name << std::endl;
}

void Scene::applyLightingShader() {
    Shader& s = lighting.shader;
    // Map raylib's material slots to our sampler names. DrawMesh binds maps[i] to texture unit i
    // and sets locs[SHADER_LOC_MAP_DIFFUSE + i] = i, so textures can never end up on the wrong unit.
    s.locs[SHADER_LOC_MAP_ALBEDO]    = GetShaderLocation(s, "colorTexture");
    s.locs[SHADER_LOC_MAP_NORMAL]    = GetShaderLocation(s, "normalTexture");
    s.locs[SHADER_LOC_MAP_ROUGHNESS] = GetShaderLocation(s, "roughAoTexture");
    s.locs[SHADER_LOC_MAP_BRDF]      = GetShaderLocation(s, "shadowMap");
    locViewPos    = GetShaderLocation(s, "viewPos");
    locLightSpace = GetShaderLocation(s, "lightSpaceMatrix");
    locAmbient    = GetShaderLocation(s, "ambient");
    locShadowRes  = GetShaderLocation(s, "shadowMapResolution");
    locTexScale   = GetShaderLocation(s, "texScale");

    const int res = SHADOW_MAP_SIZE;
    SetShaderValue(s, locShadowRes, &res, SHADER_UNIFORM_INT);

    terrainMaterial.shader = s;
}

void Scene::reloadShaders() {
    if (lighting.reloadIfChanged(true)) applyLightingShader();
    if (depth.reloadIfChanged(true)) depthMaterial.shader = depth.shader;
}

// Mipmapped, trilinear + anisotropic. Without mipmaps a 2K texture minified onto distant terrain
// samples scattered texels: it shimmers and thrashes the texture cache.
static Texture2D UploadFiltered(Image& img) {
    Texture2D t = LoadTextureFromImage(img);
    GenTextureMipmaps(&t);
    SetTextureFilter(t, TEXTURE_FILTER_TRILINEAR);
    SetTextureFilter(t, TEXTURE_FILTER_ANISOTROPIC_8X);
    SetTextureWrap(t, TEXTURE_WRAP_REPEAT);
    return t;
}

void Scene::loadTerrainTextures() {
    if (terrainTexturesLoaded) return;

    Image color = LoadImage("assets/textures/ground_0010_color_2k.jpg");
    ImageFormat(&color, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    terrainMaterial.maps[MATERIAL_MAP_ALBEDO].texture = UploadFiltered(color);
    UnloadImage(color);

    // The source normal map is 16-bit; 8 bits per channel is plenty and halves VRAM.
    Image normal = LoadImage("assets/textures/ground_0010_normal_opengl_2k.png");
    ImageFormat(&normal, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    terrainMaterial.maps[MATERIAL_MAP_NORMAL].texture = UploadFiltered(normal);
    UnloadImage(normal);

    // Pack roughness (R) and ambient occlusion (G) into one texture: one fetch instead of two.
    Image rough = LoadImage("assets/textures/ground_0010_roughness_2k.jpg");
    Image ao = LoadImage("assets/textures/ground_0010_ao_2k.jpg");
    ImageFormat(&rough, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    ImageFormat(&ao, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    if (ao.width != rough.width || ao.height != rough.height) ImageResize(&ao, rough.width, rough.height);
    unsigned char* r = (unsigned char*)rough.data;
    const unsigned char* a = (const unsigned char*)ao.data;
    for (int i = 0; i < rough.width * rough.height; ++i) {
        r[i * 4 + 1] = a[i * 4 + 0];
        r[i * 4 + 2] = 0;
        r[i * 4 + 3] = 255;
    }
    terrainMaterial.maps[MATERIAL_MAP_ROUGHNESS].texture = UploadFiltered(rough);
    UnloadImage(rough);
    UnloadImage(ao);

    terrainTexturesLoaded = true;
}

void Scene::drawScene(int gamestate) {
    if (!isActive) return;
    PROFILE_FUNCTION();

    // --- Hot reload: only touches the GPU when a shader file actually changed ---
    if (lighting.reloadIfChanged()) applyLightingShader();
    if (depth.reloadIfChanged()) depthMaterial.shader = depth.shader;

    if (IsKeyPressed(KEY_F3)) debugDraw = !debugDraw;

    // --- Update camera first so both passes use this frame's position ---
    extern bool benchMode;
    if (gamestate == 1 && !benchMode && cameraInputEnabled) customUpdateCamera(&camera);

    // --- Directional light: an orthographic box that follows the camera ---
    const Vector3 centre = camera.position;
    const Vector3 lightPos = Vector3Add(centre, Vector3Scale(sunDirection, sunDistance));
    if (!lights.empty()) {
        lights[0].position = lightPos;
        lights[0].target = centre;
    }

    Camera3D lightCam = { 0 };
    lightCam.position = lightPos;
    lightCam.target = centre;
    lightCam.up = (fabsf(sunDirection.y) > 0.99f) ? Vector3{ 0, 0, 1 } : Vector3{ 0, 1, 0 };
    lightCam.projection = CAMERA_ORTHOGRAPHIC;
    lightCam.fovy = shadowHalfExtent * 2.0f; // raylib ortho: half-height = fovy/2

    // --- Shadow pass ---
    // Read the light's matrices back from rlgl inside BeginMode3D so the shader's lightSpaceMatrix is
    // exactly what the depth pass rendered with. (DrawMesh overwrites any hand-set "mvp" uniform.)
    const double defaultNear = rlGetCullDistanceNear(), defaultFar = rlGetCullDistanceFar();
    {
    PROFILE_SCOPE("Shadow pass");
    rlSetClipPlanes(1.0, sunDistance * 2.0);
    BeginTextureMode(shadowMap);
        ClearBackground(WHITE);
        BeginMode3D(lightCam);
            const Matrix lightView = rlGetMatrixModelview();
            const Matrix lightProj = rlGetMatrixProjection();
            lightSpaceMatrix = MatrixMultiply(lightView, lightProj); // raymath order: view, then projection
            const Frustum lightFrustum = Frustum::fromViewProj(lightSpaceMatrix);
            RenderContext shadowCtx;
            shadowCtx.depthMaterial = &depthMaterial;
            shadowCtx.frustum = &lightFrustum;
            for (const auto& obj : objects) obj->drawDepthOnly(shadowCtx);
            for (const auto& obj : rootObject.children) obj->drawDepthOnly(shadowCtx);
        EndMode3D();
    EndTextureMode();
    rlSetClipPlanes(defaultNear, defaultFar);
    }

    // --- Per-frame lighting uniforms ---
    PROFILE_SCOPE("Main pass");
    Shader& s = lighting.shader;
    SetShaderValueMatrix(s, locLightSpace, lightSpaceMatrix);
    SetShaderValue(s, locTexScale, &texScale, SHADER_UNIFORM_FLOAT);
    const float ambient[4] = { ambientColor.x, ambientColor.y, ambientColor.z, 1.0f };
    SetShaderValue(s, locAmbient, ambient, SHADER_UNIFORM_VEC4);
    SetShaderValue(s, locViewPos, &camera.position, SHADER_UNIFORM_VEC3);
    for (size_t i = 0; i < lights.size(); ++i) {
        if (lights[i].enabled) UpdateLightValues(s, lights[i], (int)i);
    }

    // --- Main pass ---
    int drawCalls = 0;
    BeginMode3D(camera);
        rlDisableBackfaceCulling();
        rlDisableDepthMask();
        DrawModel(skybox, (Vector3){0, 0, 0}, 1.0f, WHITE); // skybox.vs strips translation and scales the cube
        rlEnableBackfaceCulling();
        rlEnableDepthMask();

        const Frustum viewFrustum = Frustum::fromViewProj(MatrixMultiply(rlGetMatrixModelview(), rlGetMatrixProjection()));
        RenderContext ctx;
        ctx.terrainMaterial = &terrainMaterial;
        ctx.depthMaterial = &depthMaterial;
        ctx.frustum = &viewFrustum;
        ctx.drawCalls = &drawCalls;
        for (const auto& obj : objects) if (obj) obj->draw(ctx);
        for (const auto& obj : rootObject.children) if (obj) obj->draw(ctx);

        if (debugDraw) {
            DrawSphere(lightPos, 10.0f, ColorAlpha(RED, 0.5f));
            DrawCubeWires(centre, shadowHalfExtent * 2, shadowHalfExtent * 2, shadowHalfExtent * 2, RED);
            // LOD node bounds, coloured by level (green = finest).
            static const Color levelColors[] = { GREEN, LIME, YELLOW, ORANGE, RED, MAGENTA, PURPLE, VIOLET };
            for (const auto& obj : objects) {
                const Planetoid* p = dynamic_cast<const Planetoid*>(obj.get());
                if (!p) continue;
                for (const Chunk* c : p->drawList) {
                    if (!c->hasMesh()) continue;
                    const float h = c->key.size() * 0.5f * p->scale;
                    const Vector3 ctr = Vector3Transform(Vector3AddValue(c->key.origin(), c->key.size() * 0.5f), p->localTransform());
                    DrawCubeWires(ctr, 2 * h, 2 * h, 2 * h, levelColors[std::min(c->key.level, 7)]);
                }
            }
        }
        if (onDraw3D) onDraw3D();
    EndMode3D();
    lastDrawCalls = drawCalls;

    if (debugDraw) {
        // Depth texture preview (red channel = depth).
        DrawTexturePro(shadowMap.depth,
            (Rectangle){ 0, 0, (float)shadowMap.depth.width, -(float)shadowMap.depth.height },
            (Rectangle){ 10, 100, 300, 300 }, (Vector2){ 0, 0 }, 0.0f, WHITE);
    }
}

void Scene::drawUI(int gamestate) {
    if (!isActive) return;
    RenderContext none;
    if (gamestate != 0) {
        DrawFPS(10, 10);
        const DVec3 cw = cameraWorld();
        DrawText(TextFormat("Camera: %.1f %.1f %.1f   meshes drawn: %d   [F3] debug",
                 cw.x, cw.y, cw.z, lastDrawCalls), 10, 30, 20, GREEN);
    }
    for (const auto& obj : uiObjects) obj->draw(none);
}

Scene::~Scene() {
    if (skybox.meshCount > 0) {
        UnloadShader(skybox.materials[0].shader);
        UnloadTexture(skybox.materials[0].maps[MATERIAL_MAP_CUBEMAP].texture);
        UnloadModel(skybox);
    }
    if (terrainTexturesLoaded) {
        UnloadTexture(terrainMaterial.maps[MATERIAL_MAP_ALBEDO].texture);
        UnloadTexture(terrainMaterial.maps[MATERIAL_MAP_NORMAL].texture);
        UnloadTexture(terrainMaterial.maps[MATERIAL_MAP_ROUGHNESS].texture);
    }
    // Free the maps arrays only; shaders and textures are released explicitly here.
    MemFree(terrainMaterial.maps);
    MemFree(depthMaterial.maps);
    lighting.unload();
    depth.unload();
    UnloadShadowmapRenderTexture(shadowMap);
}

// ---------------------------------------------------------------------------------------------
// Floating origin
// ---------------------------------------------------------------------------------------------

void Scene::syncRenderPositions() {
    auto sync = [this](Object& o) {
        if (!o.hasWorldPos) { o.worldPos = renderOrigin + o.position; o.hasWorldPos = true; }
        o.position = (o.worldPos - renderOrigin).toFloat();
    };
    for (auto& obj : objects) sync(*obj);
    for (auto& obj : rootObject.children) sync(*obj);
}

void Scene::updateFloatingOrigin() {
    const Vector3 c = camera.position;
    if (Vector3Length(c) < originRebaseDistance) { syncRenderPositions(); return; }
    // Shift the origin to the camera; every render-space position moves by -c.
    renderOrigin += DVec3(c);
    camera.position = Vector3Subtract(camera.position, c);
    camera.target = Vector3Subtract(camera.target, c);
    syncRenderPositions();
}

void Scene::setCameraWorld(DVec3 position, DVec3 target) {
    renderOrigin = position;
    camera.position = { 0, 0, 0 };
    camera.target = (target - position).toFloat();
    syncRenderPositions();
}
