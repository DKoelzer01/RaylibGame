#include "scene.h"
#include "gputimer.h"
#include "utils/logger.h"
#include "profiler.h"
#include "worldtypes.h"
#include "world.h"
#include <raylib.h>
#include <rlgl.h>
#include <cmath>
#include <algorithm>

#define SHADOW_MAP_SIZE 2048

bool benchMode = false; // set by the game executables (--bench)

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

// Depth-only clear between the far and near passes (rlgl only clears colour and depth together).
extern "C" { typedef void (*PfnClear)(unsigned int); extern PfnClear glad_glClear; }
static void ClearDepthBuffer() {
    rlDrawRenderBatchActive();
    if (glad_glClear) glad_glClear(0x00000100); // GL_DEPTH_BUFFER_BIT
}

// Translucent water: overlapping LOD copies of the sea must blend only once per pixel, so the
// stencil lets the first (nearest, see Planetoid::drawTransparent) water fragment claim it.
extern "C" {
    typedef void (*PfnEnable)(unsigned int);
    typedef void (*PfnStencilFunc)(unsigned int, int, unsigned int);
    typedef void (*PfnStencilOp)(unsigned int, unsigned int, unsigned int);
    extern PfnEnable glad_glEnable;
    extern PfnEnable glad_glDisable;
    extern PfnStencilFunc glad_glStencilFunc;
    extern PfnStencilOp glad_glStencilOp;
}
static void BeginWaterStencil() {
    rlDrawRenderBatchActive();
    if (!glad_glStencilFunc || !glad_glStencilOp) return;
    glad_glClear(0x00000400);                   // GL_STENCIL_BUFFER_BIT
    glad_glEnable(0x0B90);                      // GL_STENCIL_TEST
    glad_glStencilFunc(0x0202, 0, 0xFF);        // GL_EQUAL 0: unclaimed pixels only
    glad_glStencilOp(0x1E00, 0x1E00, 0x1E02);   // keep, keep, GL_INCR where drawn
}
static void EndWaterStencil() {
    rlDrawRenderBatchActive();
    if (glad_glDisable) glad_glDisable(0x0B90);
}

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
    waterMaterial = LoadMaterialDefault();
    water.reloadIfChanged(true);
    grass.reloadIfChanged(true);
    trees.reloadIfChanged(true);
    farSphere.reloadIfChanged(true);
    parts.reloadIfChanged(true);
    atmosphereShader.reloadIfChanged(true);
    atmosphereDome = LoadModelFromMesh(GenMeshCube(1.0f, 1.0f, 1.0f));

    shadowMap = LoadShadowmapRenderTexture(SHADOW_MAP_SIZE, SHADOW_MAP_SIZE);
    terrainMaterial.maps[MATERIAL_MAP_BRDF].texture = shadowMap.depth;

    depth.reloadIfChanged(true);
    depthMaterial.shader = depth.shader;
    lighting.reloadIfChanged(true);

    // Slot 0 is the sun. (Not rlights' CreateLight: it numbers lights with a global counter across
    // all scenes, so the second scene's sun was written into lights[1] and left there as a frozen
    // phantom light while slot 0 was the one being updated.)
    Light sun{};
    sun.type = LIGHT_DIRECTIONAL;
    sun.enabled = true;
    sun.position = Vector3Scale(sunDirection, sunDistance);
    sun.target = Vector3Zero();
    sun.color = WHITE;
    lights.push_back(sun);
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
    locPalette    = GetShaderLocation(s, "palette");
    locFogColor   = GetShaderLocation(s, "fogColor");
    locFogDensity = GetShaderLocation(s, "fogDensity");
    locFogShell = GetShaderLocation(s, "fogShell");
    locChunkOrigin = GetShaderLocation(s, "chunkTexOrigin");
    locChunkClip = GetShaderLocation(s, "chunkClipMax");
    locTerrainLayers = GetShaderLocation(s, "terrainLayers");
    locSunOccluder = GetShaderLocation(s, "sunOccluder");
    // Always on its own unit: two samplers of different types (sampler2D / sampler2DArray) on one unit
    // make draws fail, and an unset sampler would default to unit 0, the albedo map's.
    { const int unit = ChunkBatch::kLayersUnit; if (locTerrainLayers >= 0) SetShaderValue(s, locTerrainLayers, &unit, SHADER_UNIFORM_INT); }
    const Vector3 noClip = { 1e30f, 1e30f, 1e30f }; // default for anything that is not a chunk
    SetShaderValue(s, locChunkClip, &noClip, SHADER_UNIFORM_VEC3);

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

    terrainLayers = LoadTerrainLayers(512);

    terrainTexturesLoaded = true;
}

void Scene::drawScene(int gamestate) {
    if (!isActive) return;
    PROFILE_FUNCTION();

    // --- Hot reload: only touches the GPU when a shader file actually changed ---
    if (lighting.reloadIfChanged()) applyLightingShader();
    if (depth.reloadIfChanged()) depthMaterial.shader = depth.shader;
    water.reloadIfChanged();
    grass.reloadIfChanged();
    trees.reloadIfChanged();
    farSphere.reloadIfChanged();
    parts.reloadIfChanged();
    atmosphereShader.reloadIfChanged();

    if (IsKeyPressed(KEY_F3)) debugDraw = !debugDraw;

    // --- Update camera first so both passes use this frame's position ---

    if (gamestate == 1 && !benchMode && cameraInputEnabled) customUpdateCamera(&camera);

    // --- Directional light: an orthographic box that follows the camera ---
    // The box centre is snapped to whole shadow-map texels in light space, so moving the camera
    // never slides the texel grid across the terrain (which made shadow edges shimmer/crawl).
    // On a spinning body the ground moves through world space (~105 m/s on Terra's equator): a grid
    // snapped in world space let the terrain slide across the texels every frame and shadow edges
    // crawl. So the grid is snapped in the frame of the body under the camera, and the light's "up"
    // turns with it: the texels then stay put on the ground (bar the sun's own slow motion).
    Planetoid* ground = NearestPlanetoid(*this, camera.position);
    if (ground && ground->flat) ground = nullptr;
    Vector3 upHint = (fabsf(sunDirection.y) > 0.99f) ? Vector3{ 0, 0, 1 } : Vector3{ 0, 1, 0 };
    if (ground) {
        upHint = ground->dirToWorld({ 0, 1, 0 });
        if (fabsf(Vector3DotProduct(upHint, Vector3Normalize(sunDirection))) > 0.99f) upHint = ground->dirToWorld({ 0, 0, 1 });
    }
    // Snapping happens in double precision so floating-origin rebases don't jump it.
    Vector3 centre = camera.position;
    {
        const Vector3 f = Vector3Normalize(Vector3Negate(sunDirection));
        const Vector3 r = Vector3Normalize(Vector3CrossProduct(f, upHint));
        const Vector3 u = Vector3CrossProduct(r, f);
        const double texel = 2.0 * shadowHalfExtent / (double)SHADOW_MAP_SIZE;
        // Relative to an anchor fixed on the body near the camera (a 256 m grid): about the body's
        // centre, the sun's slow turn (0.1 deg/s on Terra) still swept the grid 60 km out at 105 m/s.
        DVec3 cw = ground ? ground->toLocalD(cameraWorld()) : cameraWorld();
        DVec3 anchor{ 0, 0, 0 };
        if (ground) {
            anchor = { std::floor(cw.x / 256.0) * 256.0, std::floor(cw.y / 256.0) * 256.0, std::floor(cw.z / 256.0) * 256.0 };
            cw = cw - anchor;
        }
        auto snapAxis = [&](Vector3 worldAxis) {
            const Vector3 a = ground ? ground->dirToLocal(worldAxis) : worldAxis;
            const double d = cw.x * a.x + cw.y * a.y + cw.z * a.z;
            const double s = std::floor(d / texel + 0.5) * texel;
            return DVec3{ a.x * s, a.y * s, a.z * s };
        };
        const DVec3 snapped = snapAxis(r) + snapAxis(u) + snapAxis(f) + anchor;
        centre = ((ground ? ground->toWorldD(snapped) : snapped) - renderOrigin).toFloat();
    }
    const Vector3 lightPos = Vector3Add(centre, Vector3Scale(Vector3Normalize(sunDirection), sunDistance));
    if (!lights.empty()) {
        lights[0].position = lightPos;
        lights[0].target = centre;
        const unsigned char k = (unsigned char)(255.0f * (1.0f - 0.7f * std::clamp(weatherCloud, 0.0f, 1.0f))); // cloud cover dims the sun
        lights[0].color = { k, k, (unsigned char)std::min(255, k + 6), 255 };
    }

    Camera3D lightCam = { 0 };
    lightCam.position = lightPos;
    lightCam.target = centre;
    lightCam.up = upHint;
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
            if (shadowProbe) {
                const Vector3 wp = (shadowProbe() - renderOrigin).toFloat();
                const Matrix& m = lightSpaceMatrix;
                const float x = m.m0 * wp.x + m.m4 * wp.y + m.m8 * wp.z + m.m12, y = m.m1 * wp.x + m.m5 * wp.y + m.m9 * wp.z + m.m13;
                const Vector2 tx = { (x * 0.5f + 0.5f) * SHADOW_MAP_SIZE, (y * 0.5f + 0.5f) * SHADOW_MAP_SIZE };
                if (shadowProbeHave) {
                    // Whole-texel jumps are the box re-snapping; the fraction is the ground sliding.
                    const float dx = tx.x - shadowProbeLast.x, dy = tx.y - shadowProbeLast.y;
                    shadowProbeSlide += sqrtf((dx - roundf(dx)) * (dx - roundf(dx)) + (dy - roundf(dy)) * (dy - roundf(dy)));
                }
                shadowProbeLast = tx; shadowProbeHave = true;
            }
            const Frustum lightFrustum = Frustum::fromViewProj(lightSpaceMatrix);
            RenderContext shadowCtx;
            shadowCtx.depthMaterial = &depthMaterial;
            shadowCtx.frustum = &lightFrustum;
            shadowCtx.renderOrigin = renderOrigin;
            for (const auto& obj : objects) obj->drawDepthOnly(shadowCtx);
            for (const auto& obj : rootObject.children) obj->drawDepthOnly(shadowCtx);
        EndMode3D();
    EndTextureMode();
    GpuTimer_Mark("shadow");
    rlSetClipPlanes(defaultNear, defaultFar);
    }

    // --- Per-frame lighting uniforms ---
    PROFILE_SCOPE("Main pass");
    Shader& s = lighting.shader;
    SetShaderValueMatrix(s, locLightSpace, lightSpaceMatrix);
    SetShaderValue(s, locTexScale, &texScale, SHADER_UNIFORM_FLOAT);
    SetShaderValue(s, locPalette, &materialPalette, SHADER_UNIFORM_INT);
    // Atmosphere: by day near a planet the air hazes distant terrain toward the horizon sky colour and
    // lights shadows a little (sky ambient); both fade with altitude and at night. Under the sea the
    // camera sees a dense dark-blue fog instead.
    const Vector3 sun = Vector3Normalize(sunDirection);
    Vector3 effFog = fogColor, effAmbient = ambientColor, skyHorizon = { 0, 0, 0 };
    float effFogDensity = fogDensity, atmDensity = 0.0f;
    Vector3 atmUp = { 0, 1, 0 };
    if (atmosphere.enabled) {
        const Vector3 rel = Vector3Subtract(camera.position, atmosphere.centre);
        const float d = Vector3Length(rel);
        atmUp = d > 1e-3f ? Vector3Scale(rel, 1.0f / d) : Vector3{ 0, 1, 0 };
        const float alt = d - atmosphere.radius;
        atmDensity = powf(std::clamp(1.0f - alt / atmosphere.height, 0.0f, 1.0f), 1.5f);
        const float s = Vector3DotProduct(atmUp, sun);
        const float t = std::clamp((s + 0.18f) / 0.30f, 0.0f, 1.0f);
        const float day = t * t * (3.0f - 2.0f * t);
        skyHorizon = Vector3Scale({ 0.35f, 0.53f, 0.83f }, day); // linear horizon sky
        effFog = Vector3Add(effFog, Vector3Scale(skyHorizon, atmDensity));
        effFogDensity += atmDensity / 35000.0f;                    // ~35 km to fade a third of the way
        effAmbient = Vector3Add(effAmbient, Vector3Scale(skyHorizon, 0.06f * atmDensity));
        // Weather: cloud cover greys the haze and lifts the (diffuse) ambient a little; ground fog
        // thickens the haze near the ground.
        const Vector3 grey = Vector3Scale({ 0.30f, 0.32f, 0.35f }, day);
        effFog = Vector3Lerp(effFog, Vector3Scale(grey, atmDensity), weatherCloud * 0.6f);
        effFogDensity += atmDensity * (weatherCloud * weatherCloud / 9000.0f + weatherFog / 60.0f);
        effAmbient = Vector3Add(effAmbient, Vector3Scale(grey, 0.05f * weatherCloud * atmDensity));
        if (alt < 0.0f && IsUnderwater(*this, camera.position)) { effFog = { 0.004f, 0.02f, 0.035f }; effFogDensity = 0.06f; } // in the sea (not a dry cave below sea level)
    }
    SetShaderValue(s, locFogColor, &effFog, SHADER_UNIFORM_VEC3);
    SetShaderValue(s, locFogDensity, &effFogDensity, SHADER_UNIFORM_FLOAT);
    {
        const bool underwater = atmosphere.enabled && effFogDensity >= 0.06f;
        const float shell[4] = { atmosphere.centre.x, atmosphere.centre.y, atmosphere.centre.z,
                                 atmosphere.enabled && !underwater ? atmosphere.radius + atmosphere.hazeHeight : 0.0f };
        SetShaderValue(s, locFogShell, shell, SHADER_UNIFORM_VEC4);
    }
    const float ambient[4] = { effAmbient.x, effAmbient.y, effAmbient.z, 1.0f };
    SetShaderValue(s, locAmbient, ambient, SHADER_UNIFORM_VEC4);
    SetShaderValue(s, locViewPos, &camera.position, SHADER_UNIFORM_VEC3);
    // Upload every slot the shader has, so unused slots are explicitly disabled.
    for (int i = 0; i < MAX_LIGHTS; ++i) {
        if (i < (int)lights.size()) UpdateLightValues(s, lights[i], i);
        else { Light off{}; UpdateLightValues(s, off, i); }
    }

    // --- Main pass ---
    int drawCalls = 0;
    BeginMode3D(camera);
        rlDisableBackfaceCulling();
        rlDisableDepthMask();
        if (drawSkybox) DrawModel(skybox, (Vector3){0, 0, 0}, 1.0f, WHITE); // skybox.vs strips translation and scales the cube
        if (drawSkybox) {
            // The sky (and the sun disc, also in space) blended over the stars.
            Shader& as = atmosphereShader.shader;
            SetShaderValue(as, GetShaderLocation(as, "sunDir"), &sun, SHADER_UNIFORM_VEC3);
            SetShaderValue(as, GetShaderLocation(as, "upDir"), &atmUp, SHADER_UNIFORM_VEC3);
            SetShaderValue(as, GetShaderLocation(as, "density"), &atmDensity, SHADER_UNIFORM_FLOAT);
            const float overcast = std::clamp(weatherCloud, 0.0f, 1.0f) * std::min(1.0f, atmDensity * 3.0f);
            SetShaderValue(as, GetShaderLocation(as, "overcast"), &overcast, SHADER_UNIFORM_FLOAT);
            {   // Clouds: the camera in the planet's frame, and the layer's radius (none without air).
                const Vector3 rel = Vector3Subtract(camera.position, atmosphere.centre);
                const Vector3 camLocal = Vector3Add(Vector3Add(Vector3Scale(atmosphere.toLocalX, rel.x), Vector3Scale(atmosphere.toLocalY, rel.y)), Vector3Scale(atmosphere.toLocalZ, rel.z));
                const float shell = atmosphere.enabled ? atmosphere.radius + atmosphere.cloudHeight : 0.0f;
                const float cover = std::clamp(weatherCloud, 0.0f, 1.0f);
                SetShaderValue(as, GetShaderLocation(as, "camLocal"), &camLocal, SHADER_UNIFORM_VEC3);
                SetShaderValue(as, GetShaderLocation(as, "toLocalX"), &atmosphere.toLocalX, SHADER_UNIFORM_VEC3);
                SetShaderValue(as, GetShaderLocation(as, "toLocalY"), &atmosphere.toLocalY, SHADER_UNIFORM_VEC3);
                SetShaderValue(as, GetShaderLocation(as, "toLocalZ"), &atmosphere.toLocalZ, SHADER_UNIFORM_VEC3);
                SetShaderValue(as, GetShaderLocation(as, "cloudShell"), &shell, SHADER_UNIFORM_FLOAT);
                SetShaderValue(as, GetShaderLocation(as, "cloudDrift"), &cloudDrift, SHADER_UNIFORM_VEC3);
                SetShaderValue(as, GetShaderLocation(as, "cloudCover"), &cover, SHADER_UNIFORM_FLOAT);
            }
            atmosphereDome.materials[0].shader = as;
            DrawModel(atmosphereDome, (Vector3){ 0, 0, 0 }, 1.0f, WHITE);
        }
        rlEnableBackfaceCulling();
        rlEnableDepthMask();
        GpuTimer_Mark("sky");

        // Far pass: bodies whose nearest possible surface is beyond kFarSplit (the moon, the ring,
        // the planet seen from orbit) are drawn first with clip planes fitted to them, then the depth
        // buffer is cleared for the near pass. One depth range from 0.5 m to 1000 km would leave
        // the near terrain no precision.
        const double kFarSplit = 20000.0;
        std::vector<Object*> farBodies;
        double farMin = 1e30, farMax = 0.0;
        for (const auto& obj : objects) {
            const Planetoid* p = dynamic_cast<const Planetoid*>(obj.get());
            if (!p || !p->field) continue;
            const double d = Vector3Distance(camera.position, p->position);
            const double rMax = p->field->maxRadius() * p->scale, rOcc = p->field->occluderRadius() * p->scale;
            const double dMin = d - rMax;
            if (dMin < kFarSplit) continue;
            // Farthest visible point: over the body's horizon (its occluder) to the tallest peaks.
            const double horizon = sqrt(std::max(0.0, d * d - rOcc * rOcc)) + sqrt(std::max(0.0, rMax * rMax - rOcc * rOcc));
            farBodies.push_back(obj.get());
            farMin = std::min(farMin, dMin);
            farMax = std::max(farMax, std::min(d + rMax, horizon));
        }
        bool farExtra = false;
        if (farExtraRange) { double a, b; if (farExtraRange(a, b)) { farExtra = true; farMin = std::min(farMin, a); farMax = std::max(farMax, b); } }
        RenderContext ctx;
        ctx.terrainMaterial = &terrainMaterial;
        ctx.depthMaterial = &depthMaterial;
        ctx.drawCalls = &drawCalls;
        ctx.chunkOriginLoc = locChunkOrigin;
        ctx.chunkClipLoc = locChunkClip;
        ctx.renderOrigin = renderOrigin;
        ctx.texScale = texScale;
        ctx.terrainLayers = &terrainLayers;
        ctx.terrainLayersLoc = locTerrainLayers;
        ctx.paletteLoc = locPalette;
        ctx.sunOccluderLoc = locSunOccluder;
        ctx.scenePalette = materialPalette;
        ctx.cameraPos = camera.position;
        ctx.sunDir = sun;
        ctx.ambient = effAmbient;
        if (farSphere.shader.id > 0) {
            Shader& fs = farSphere.shader;
            SetShaderValue(fs, GetShaderLocation(fs, "sunDir"), &sun, SHADER_UNIFORM_VEC3);
            SetShaderValue(fs, GetShaderLocation(fs, "viewPos"), &camera.position, SHADER_UNIFORM_VEC3);
            SetShaderValue(fs, GetShaderLocation(fs, "ambient"), &effAmbient, SHADER_UNIFORM_VEC3);
            const float fshell[4] = { atmosphere.centre.x, atmosphere.centre.y, atmosphere.centre.z,
                                      atmosphere.enabled && effFogDensity < 0.06f ? atmosphere.radius + atmosphere.hazeHeight : 0.0f };
            SetShaderValue(fs, GetShaderLocation(fs, "fogColor"), &effFog, SHADER_UNIFORM_VEC3);
            SetShaderValue(fs, GetShaderLocation(fs, "fogDensity"), &effFogDensity, SHADER_UNIFORM_FLOAT);
            SetShaderValue(fs, GetShaderLocation(fs, "fogShell"), fshell, SHADER_UNIFORM_VEC4);
            ctx.farShader = &fs;
        }
        if (parts.shader.id > 0) {
            Shader& ps = parts.shader;
            const float pshell[4] = { atmosphere.centre.x, atmosphere.centre.y, atmosphere.centre.z,
                                      atmosphere.enabled && effFogDensity < 0.06f ? atmosphere.radius + atmosphere.hazeHeight : 0.0f };
            SetShaderValue(ps, GetShaderLocation(ps, "sunDir"), &sun, SHADER_UNIFORM_VEC3);
            SetShaderValue(ps, GetShaderLocation(ps, "viewPos"), &camera.position, SHADER_UNIFORM_VEC3);
            SetShaderValue(ps, GetShaderLocation(ps, "ambient"), &effAmbient, SHADER_UNIFORM_VEC3);
            SetShaderValue(ps, GetShaderLocation(ps, "fogColor"), &effFog, SHADER_UNIFORM_VEC3);
            SetShaderValue(ps, GetShaderLocation(ps, "fogDensity"), &effFogDensity, SHADER_UNIFORM_FLOAT);
            SetShaderValue(ps, GetShaderLocation(ps, "fogShell"), pshell, SHADER_UNIFORM_VEC4);
            ctx.partShader = &ps;
            ctx.partOccluderLoc = GetShaderLocation(ps, "sunOccluder");
        }
        auto isFar = [&](const Object* o) { return std::find(farBodies.begin(), farBodies.end(), o) != farBodies.end(); };
        const double nearDefault = rlGetCullDistanceNear(), farDefault = rlGetCullDistanceFar();
        // Sea surface material, set up once for the far, mid and near passes.
        auto prepareWater = [&] {
            Shader& ws = water.shader;
            const float now = (float)fmod(GetTime(), 3600.0);
            SetShaderValue(ws, GetShaderLocation(ws, "viewPos"), &camera.position, SHADER_UNIFORM_VEC3);
            SetShaderValue(ws, GetShaderLocation(ws, "sunDir"), &sun, SHADER_UNIFORM_VEC3);
            SetShaderValue(ws, GetShaderLocation(ws, "skyColor"), &skyHorizon, SHADER_UNIFORM_VEC3);
            SetShaderValue(ws, GetShaderLocation(ws, "fogColor"), &effFog, SHADER_UNIFORM_VEC3);
            SetShaderValue(ws, GetShaderLocation(ws, "fogDensity"), &effFogDensity, SHADER_UNIFORM_FLOAT);
            SetShaderValue(ws, GetShaderLocation(ws, "time"), &now, SHADER_UNIFORM_FLOAT);
            waterMaterial.shader = ws;
            ctx.waterMaterial = &waterMaterial;
            ctx.waterTexOriginLoc = GetShaderLocation(ws, "chunkTexOrigin");
            ctx.waterClipLoc = GetShaderLocation(ws, "chunkClipMax");
            ctx.waterOccluderLoc = GetShaderLocation(ws, "sunOccluder");
        };
        prepareWater();
        if (!farBodies.empty() || farExtra) {
            EndMode3D();
            rlSetClipPlanes(std::max(1.0, farMin * 0.9), farMax * 1.05 + 100.0);
            BeginMode3D(camera);
            const Frustum farFrustum = Frustum::fromViewProj(MatrixMultiply(rlGetMatrixModelview(), rlGetMatrixProjection()));
            ctx.frustum = &farFrustum;
            for (Object* o : farBodies) o->draw(ctx);
            BeginWaterStencil(); // the sea too (a planet seen from orbit is a far body)
            rlDisableBackfaceCulling();
            for (Object* o : farBodies) o->drawTransparent(ctx);
            rlEnableBackfaceCulling();
            EndWaterStencil();
            if (drawFarExtra) drawFarExtra(ctx); // translucent ring band, depth-tested against the far bodies
            EndMode3D();
            ClearDepthBuffer();
            rlSetClipPlanes(nearDefault, farDefault);
            BeginMode3D(camera);
            GpuTimer_Mark("far");
        }


        // Trees and bushes: one instanced draw per variant per body. The mid pass draws the simple
        // far meshes (trees reach 1 km, past the near/mid split); the near pass draws them all.
        // Foliage lighting (grass, trees): the sun's shadow map and the headlamp, as the terrain gets.
        auto setFoliageLighting = [&](Shader& sh) {
            const int slot = 6;
            rlActiveTextureSlot(slot);
            rlEnableTexture(shadowMap.depth.id);
            rlSetUniform(GetShaderLocation(sh, "shadowMap"), &slot, RL_SHADER_UNIFORM_INT, 1);
            rlSetUniformMatrix(GetShaderLocation(sh, "lightSpaceMatrix"), lightSpaceMatrix);
            const int res = shadowMap.depth.width;
            rlSetUniform(GetShaderLocation(sh, "shadowMapResolution"), &res, RL_SHADER_UNIFORM_INT, 1);
            float lamp[4] = { 0, 0, 0, 0 }, target[4] = { 0, 0, 0, 0 };
            if (lights.size() > 1 && lights[1].enabled && lights[1].type == LIGHT_SPOT) {
                lamp[0] = lights[1].position.x; lamp[1] = lights[1].position.y; lamp[2] = lights[1].position.z; lamp[3] = lights[1].color.a / 255.0f;
                target[0] = lights[1].target.x; target[1] = lights[1].target.y; target[2] = lights[1].target.z;
            }
            rlSetUniform(GetShaderLocation(sh, "lampPos"), lamp, RL_SHADER_UNIFORM_VEC4, 1);
            rlSetUniform(GetShaderLocation(sh, "lampTarget"), target, RL_SHADER_UNIFORM_VEC4, 1);
            rlActiveTextureSlot(0);
        };
        auto drawTrees = [&](int onlyLod) {
            if (trees.shader.id <= 0 || !foliageEnabled) return;
            rlDrawRenderBatchActive();
            Shader& ts = trees.shader;
            rlEnableShader(ts.id);
            const Matrix vp = MatrixMultiply(rlGetMatrixModelview(), rlGetMatrixProjection());
            rlSetUniformMatrix(GetShaderLocation(ts, "matViewProj"), vp);
            const float now = (float)fmod(GetTime(), 3600.0);
            const float shell[4] = { atmosphere.centre.x, atmosphere.centre.y, atmosphere.centre.z,
                                     atmosphere.enabled && effFogDensity < 0.06f ? atmosphere.radius + atmosphere.hazeHeight : 0.0f };
            rlSetUniform(GetShaderLocation(ts, "sunDir"), &sun, RL_SHADER_UNIFORM_VEC3, 1);
            rlSetUniform(GetShaderLocation(ts, "viewPos"), &camera.position, RL_SHADER_UNIFORM_VEC3, 1);
            rlSetUniform(GetShaderLocation(ts, "ambient"), &effAmbient, RL_SHADER_UNIFORM_VEC3, 1);
            rlSetUniform(GetShaderLocation(ts, "fogColor"), &effFog, RL_SHADER_UNIFORM_VEC3, 1);
            rlSetUniform(GetShaderLocation(ts, "fogDensity"), &effFogDensity, RL_SHADER_UNIFORM_FLOAT, 1);
            rlSetUniform(GetShaderLocation(ts, "fogShell"), shell, RL_SHADER_UNIFORM_VEC4, 1);
            rlSetUniform(GetShaderLocation(ts, "time"), &now, RL_SHADER_UNIFORM_FLOAT, 1);
            setFoliageLighting(ts);
            RenderContext tctx = ctx;
            tctx.treeShader = &ts;
            tctx.treeModelLoc = GetShaderLocation(ts, "matModel");
            tctx.treeOccluderLoc = GetShaderLocation(ts, "sunOccluder");
            tctx.treeFadeLoc = GetShaderLocation(ts, "fadeEnd");
            for (const auto& obj : objects)
                if (auto* p = dynamic_cast<Planetoid*>(obj.get()); p && !isFar(p) && p->trees) p->trees->draw(*p, tctx, onlyLod);
            rlDisableShader();
        };

        // Mid pass: near bodies' terrain and sea beyond kNearSplit, in their own depth range, then a
        // depth clear for the close range. One 0.5 m .. 60 km range resolved only ~12 m of depth at
        // 10 km, so far LOD seams and shorelines z-fought. Chunks straddling the split are drawn in
        // both passes; each pass clips to its range (and frustum-culls chunks outside it).
        const double kNearSplit = 2000.0;
        {
            EndMode3D();
            rlSetClipPlanes(kNearSplit * 0.9, farDefault);
            BeginMode3D(camera);
            const Frustum midFrustum = Frustum::fromViewProj(MatrixMultiply(rlGetMatrixModelview(), rlGetMatrixProjection()));
            ctx.frustum = &midFrustum;
            for (const auto& obj : objects) if (obj && !isFar(obj.get())) obj->draw(ctx);
            if (drawPassExtra) drawPassExtra(ctx, true);  // ring rocks (instances)
            rlDisableBackfaceCulling();
            BeginWaterStencil();
            for (const auto& obj : objects) if (obj && !isFar(obj.get())) obj->drawTransparent(ctx);
            EndWaterStencil();
            rlEnableBackfaceCulling();
            EndMode3D();
            ClearDepthBuffer();
            rlSetClipPlanes(nearDefault, kNearSplit);
            BeginMode3D(camera);
            GpuTimer_Mark("mid");
        }

        const Frustum viewFrustum = Frustum::fromViewProj(MatrixMultiply(rlGetMatrixModelview(), rlGetMatrixProjection()));
        ctx.terrainMaterial = &terrainMaterial;
        ctx.depthMaterial = &depthMaterial;
        ctx.frustum = &viewFrustum;
        for (const auto& obj : objects) if (obj && !isFar(obj.get())) obj->draw(ctx);
        for (const auto& obj : rootObject.children) if (obj) obj->draw(ctx);
        if (drawPassExtra) drawPassExtra(ctx, false);
        GpuTimer_Mark("near");

        // Foliage: grass tufts on near chunks (instanced; one draw per chunk with grass in range).
        if (grass.shader.id > 0 && foliageEnabled) {
            rlDrawRenderBatchActive();
            Shader& gs = grass.shader;
            rlEnableShader(gs.id);
            const Matrix vp = MatrixMultiply(rlGetMatrixModelview(), rlGetMatrixProjection());
            rlSetUniformMatrix(GetShaderLocation(gs, "matViewProj"), vp);
            const float now = (float)fmod(GetTime(), 3600.0);
            const float fadeEnd = kGrassDrawDistance;
            const float shell[4] = { atmosphere.centre.x, atmosphere.centre.y, atmosphere.centre.z,
                                     atmosphere.enabled && effFogDensity < 0.06f ? atmosphere.radius + atmosphere.hazeHeight : 0.0f };
            rlSetUniform(GetShaderLocation(gs, "sunDir"), &sun, RL_SHADER_UNIFORM_VEC3, 1);
            rlSetUniform(GetShaderLocation(gs, "viewPos"), &camera.position, RL_SHADER_UNIFORM_VEC3, 1);
            rlSetUniform(GetShaderLocation(gs, "ambient"), &effAmbient, RL_SHADER_UNIFORM_VEC3, 1);
            rlSetUniform(GetShaderLocation(gs, "fogColor"), &effFog, RL_SHADER_UNIFORM_VEC3, 1);
            rlSetUniform(GetShaderLocation(gs, "fogDensity"), &effFogDensity, RL_SHADER_UNIFORM_FLOAT, 1);
            rlSetUniform(GetShaderLocation(gs, "fogShell"), shell, RL_SHADER_UNIFORM_VEC4, 1);
            rlSetUniform(GetShaderLocation(gs, "time"), &now, RL_SHADER_UNIFORM_FLOAT, 1);
            rlSetUniform(GetShaderLocation(gs, "fadeEnd"), &fadeEnd, RL_SHADER_UNIFORM_FLOAT, 1);
            setFoliageLighting(gs);
            RenderContext gctx = ctx;
            gctx.grassShader = &gs;
            gctx.grassModelLoc = GetShaderLocation(gs, "matModel");
            gctx.grassOccluderLoc = GetShaderLocation(gs, "sunOccluder");
            gctx.cameraPos = camera.position;
            rlDisableBackfaceCulling();
            for (const auto& obj : objects)
                if (auto* p = dynamic_cast<Planetoid*>(obj.get()); p && !isFar(p)) p->drawFoliage(gctx);
            rlEnableBackfaceCulling();
            rlDisableShader();
        }
        GpuTimer_Mark("grass");
        drawTrees(-1);
        GpuTimer_Mark("trees");

        // Transparent pass (sea surfaces) after everything opaque. Depth writes stay on: where LOD
        // nodes overlap, the second copy of the same surface then fails the depth test instead of
        // blending twice. Both faces are drawn (seen from under the water too).
        {
            rlDisableBackfaceCulling();
            BeginWaterStencil();
            for (const auto& obj : objects) if (obj && !isFar(obj.get())) obj->drawTransparent(ctx);
            EndWaterStencil();
            rlEnableBackfaceCulling();
        }

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
        GpuTimer_Mark("water");
        if (onDraw3D) onDraw3D();
    EndMode3D();
    rlSetClipPlanes(nearDefault, farDefault); // the passes above narrowed them
    lastDrawCalls = drawCalls;
    GpuTimer_Mark("draw3d");

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
        UnloadTextureArray(terrainLayers);
    }
    // Free the maps arrays only; shaders and textures are released explicitly here.
    MemFree(terrainMaterial.maps);
    MemFree(depthMaterial.maps);
    lighting.unload();
    depth.unload();
    water.unload();
    grass.unload();
    trees.unload();
    farSphere.unload();
    parts.unload();
    atmosphereShader.unload();
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
