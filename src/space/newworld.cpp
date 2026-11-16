#include "newworld.h"
#include "voxel.h"
#include "world.h"

#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include <vector>

// The New World screen (#504). See newworld.h.

namespace {
    constexpr int kTexW = 1024, kTexH = 512; // the preview's equirectangular map (a build takes ~0.05 s)
    constexpr int kSphereRT = 640;           // the turning planet is drawn into this square, then onto the screen
    const Color kPanel = { 18, 24, 34, 230 }, kButton = { 44, 58, 80, 255 }, kHover = { 70, 96, 132, 255 }, kText = { 230, 236, 245, 255 };

    // ---- Screen state (main thread) ----
    unsigned int seed = 0;
    std::string seedText;
    bool seedFocus = false;
    WorldGenParams params;
    double changedAt = 0.0;       // the last input change (the preview rebuilds 0.35 s after it)
    bool dirty = true;
    float spin = 0.0f;

    // ---- Preview (a worker thread builds the map; the main thread uploads it) ----
    std::thread worker;
    std::atomic<bool> building{ false }, cancel{ false }, ready{ false };
    std::atomic<double> lastBuildSeconds{ -1.0 };
    std::vector<Color> pixels(kTexW * kTexH);
    std::mutex pixelMx;
    Texture2D mapTex{};
    RenderTexture2D sphereRT{};
    Model sphere{};
    Shader sphereShader{};
    int lightLoc = -1;
    bool gpuReady = false;
    bool haveMap = false;

    float Scale() { return std::max(1.0f, GetScreenHeight() / 1080.0f); }

    bool Button(Rectangle r, const char* label) {
        const bool hover = CheckCollisionPointRec(GetMousePosition(), r);
        DrawRectangleRounded(r, 0.18f, 6, hover ? kHover : kButton);
        if (hover) DrawRectangleRoundedLinesEx(r, 0.18f, 6, 2.0f, Fade(SKYBLUE, 0.8f));
        const int fs = (int)(r.height * 0.45f);
        DrawText(label, (int)(r.x + (r.width - MeasureText(label, fs)) * 0.5f), (int)(r.y + (r.height - fs) * 0.5f), fs, kText);
        return hover && IsMouseButtonReleased(MOUSE_BUTTON_LEFT);
    }
    // A slider row (label, bar, value), as the settings menu's. True while it changes.
    bool Slider(Rectangle r, const char* label, float& v, float lo, float hi, const char* fmt) {
        const int fs = (int)(r.height * 0.5f);
        DrawText(label, (int)r.x, (int)(r.y + (r.height - fs) * 0.5f), fs, kText);
        const Rectangle bar = { r.x + r.width * 0.42f, r.y + r.height * 0.4f, r.width * 0.4f, r.height * 0.2f };
        DrawRectangleRec(bar, kButton);
        const float f = (v - lo) / (hi - lo);
        DrawRectangle((int)bar.x, (int)bar.y, (int)(bar.width * f), (int)bar.height, SKYBLUE);
        DrawCircle((int)(bar.x + bar.width * f), (int)(bar.y + bar.height * 0.5f), r.height * 0.22f, kText);
        DrawText(TextFormat(fmt, v), (int)(bar.x + bar.width + 16), (int)(r.y + (r.height - fs) * 0.5f), fs, kText);
        const Rectangle hit = { bar.x - 10, r.y, bar.width + 20, r.height };
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) && CheckCollisionPointRec(GetMousePosition(), hit)) {
            const float nf = std::clamp((GetMousePosition().x - bar.x) / bar.width, 0.0f, 1.0f);
            const float nv = lo + nf * (hi - lo);
            if (nv != v) { v = nv; return true; }
        }
        return false;
    }

    // The colour of a point of the map: the sea by depth; land by biome, lighter with height, rock and snow up high.
    Color MapColour(const PlanetField& pf, Vector3 dir, float e) {
        if (e < 0.0f) {
            const float t = std::clamp(-e / 600.0f, 0.0f, 1.0f);
            return ColorLerp({ 46, 120, 160, 255 }, { 10, 30, 70, 255 }, t);
        }
        static const Color kBiome[(int)Biome::Count] = {
            { 46, 120, 160, 255 },  // Ocean (a low shore: the shallows' colour)
            { 30, 92, 34, 255 },    // Jungle
            { 168, 160, 78, 255 },  // Savanna
            { 222, 196, 132, 255 }, // Desert
            { 178, 96, 58, 255 },   // Badlands
            { 52, 110, 46, 255 },   // Forest
            { 106, 160, 72, 255 },  // Meadow
            { 76, 96, 54, 255 },    // Swamp
            { 46, 84, 62, 255 },    // Taiga
            { 140, 136, 112, 255 }, // Tundra
            { 132, 130, 126, 255 }, // Alpine
            { 236, 242, 248, 255 }, // Polar
        };
        Color c = kBiome[std::clamp((int)pf.biomeAt(dir, e), 0, (int)Biome::Count - 1)];
        c = ColorLerp(c, { 120, 116, 110, 255 }, std::clamp((e - 500.0f) / 700.0f, 0.0f, 0.8f)); // rock up high
        if (e > 1150.0f) c = ColorLerp(c, { 245, 248, 252, 255 }, std::clamp((e - 1150.0f) / 250.0f, 0.0f, 1.0f)); // snow caps
        const float lift = std::clamp(e / 1500.0f, 0.0f, 1.0f) * 0.25f;
        return { (unsigned char)std::min(255.0f, c.r * (0.85f + lift)), (unsigned char)std::min(255.0f, c.g * (0.85f + lift)), (unsigned char)std::min(255.0f, c.b * (0.85f + lift)), 255 };
    }

}

// --seed-survey (#532): the preview's map at any size, single threaded (the caller runs seeds in parallel).
void NewWorld_RenderMap(unsigned int terrainSeed, float radius, const WorldGenParams& p, int w, int h, std::vector<Color>& out, std::vector<float>* heights) {
    PlanetField pf(terrainSeed, radius, false);
    pf.setGenParams(p);
    out.assign((size_t)w * h, BLACK);
    if (heights) heights->assign((size_t)w * h, 0.0f);
    const float fp = 2.0f * PI * radius / kTexW; // the preview's footprint, whatever the size
    for (int y = 0; y < h; ++y) {
        const float lat = (0.5f - (y + 0.5f) / h) * PI;
        for (int x = 0; x < w; ++x) {
            const float lon = (x + 0.5f) / w * 2.0f * PI;
            const Vector3 dir = { cosf(lat) * cosf(lon), sinf(lat), cosf(lat) * sinf(lon) };
            const float e = pf.elevationBase(dir, fp);
            out[(size_t)y * w + x] = MapColour(pf, dir, e);
            if (heights) (*heights)[(size_t)y * w + x] = e;
        }
    }
}

namespace {
    // The worker: the start planet's base terrain (no rivers, no bake) sampled on the map, in a few threads.
    void BuildMap(unsigned int terrainSeed, float radius, WorldGenParams p) {
        const auto t0 = std::chrono::steady_clock::now();
        PlanetField pf(terrainSeed, radius, false);
        pf.setGenParams(p);
        std::vector<Color> out(kTexW * kTexH);
        const float fp = 2.0f * PI * radius / kTexW; // metres a texel: coarse noise only
        const int nThreads = std::clamp((int)std::thread::hardware_concurrency() / 2, 1, 6);
        std::vector<std::thread> rows;
        for (int t = 0; t < nThreads; ++t)
            rows.emplace_back([&, t] {
                for (int y = t; y < kTexH && !cancel.load(); y += nThreads) {
                    const float lat = (0.5f - (y + 0.5f) / kTexH) * PI;
                    for (int x = 0; x < kTexW; ++x) {
                        const float lon = (x + 0.5f) / kTexW * 2.0f * PI;
                        const Vector3 dir = { cosf(lat) * cosf(lon), sinf(lat), cosf(lat) * sinf(lon) };
                        out[y * kTexW + x] = MapColour(pf, dir, pf.elevationBase(dir, fp));
                    }
                }
            });
        for (auto& r : rows) r.join();
        if (cancel.load()) { building.store(false); return; }
        { std::lock_guard<std::mutex> lk(pixelMx); pixels.swap(out); }
        lastBuildSeconds.store(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        ready.store(true);
        building.store(false);
    }

    void StartBuild() {
        if (worker.joinable()) worker.join();
        unsigned int tseed = 0; float radius = 0.0f;
        if (!WorldGen_StartPlanet(seed, tseed, radius)) return; // no planet in the system: nothing to show
        cancel.store(false);
        building.store(true);
        worker = std::thread(BuildMap, tseed, radius, params);
        dirty = false;
    }

    void InitGpu() {
        if (gpuReady) return;
        Image img = GenImageColor(kTexW, kTexH, { 10, 30, 70, 255 });
        mapTex = LoadTextureFromImage(img);
        UnloadImage(img);
        SetTextureFilter(mapTex, TEXTURE_FILTER_BILINEAR);
        sphereRT = LoadRenderTexture(kSphereRT, kSphereRT);
        SetTextureFilter(sphereRT.texture, TEXTURE_FILTER_BILINEAR);
        // par_shapes' sphere has its poles on z, texcoord.x = the polar angle / pi and texcoord.y = the longitude / 2 pi:
        // the map is the other way round (x longitude, y latitude from the north), so its cold rows ran along a meridian
        // (#513). Drawn turned z -> y (below), a vertex shows the direction (x, z, -y): latitude pi/2 - phi, longitude
        // -theta, so the map's uv is (1 - t.y, t.x). The seam vertices par_shapes duplicates stay seamless.
        Mesh globe = GenMeshSphere(1.0f, 64, 128);
        for (int i = 0; i < globe.vertexCount; ++i) {
            float* t = &globe.texcoords[i * 2];
            const float polar = t[0], around = t[1];
            t[0] = 1.0f - around;
            t[1] = polar;
        }
        UpdateMeshBuffer(globe, 1, globe.texcoords, globe.vertexCount * 2 * (int)sizeof(float), 0);
        sphere = LoadModelFromMesh(globe);
        static const char* vs =
            "#version 330\n"
            "in vec3 vertexPosition; in vec2 vertexTexCoord; in vec3 vertexNormal;\n"
            "uniform mat4 mvp; uniform mat4 matModel;\n"
            "out vec2 uv; out vec3 n;\n"
            "void main() { uv = vertexTexCoord; n = normalize(mat3(matModel) * vertexNormal); gl_Position = mvp * vec4(vertexPosition, 1.0); }\n";
        static const char* fs =
            "#version 330\n"
            "in vec2 uv; in vec3 n; uniform sampler2D texture0; uniform vec3 lightDir; out vec4 finalColor;\n"
            "void main() {\n"
            "  vec3 c = texture(texture0, uv).rgb;\n"
            "  float l = max(dot(normalize(n), normalize(lightDir)), 0.0);\n"
            "  float rim = pow(1.0 - abs(normalize(n).z), 3.0);\n"         // a thin bright limb (air)
            "  finalColor = vec4(c * (0.12 + 0.88 * l) + vec3(0.25, 0.45, 0.8) * rim * 0.5, 1.0);\n"
            "}\n";
        sphereShader = LoadShaderFromMemory(vs, fs);
        lightLoc = GetShaderLocation(sphereShader, "lightDir");
        sphere.materials[0].shader = sphereShader;
        sphere.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = mapTex;
        gpuReady = true;
    }

    void DrawPreview(Rectangle area) {
        if (ready.exchange(false)) {
            std::lock_guard<std::mutex> lk(pixelMx);
            UpdateTexture(mapTex, pixels.data());
            haveMap = true;
        }
        spin += GetFrameTime() * 0.12f; // a turn in ~50 s
        BeginTextureMode(sphereRT);
        ClearBackground(BLANK);
        Camera3D cam{};
        cam.position = { 0.0f, 0.0f, 3.2f }; cam.target = { 0, 0, 0 }; cam.up = { 0, 1, 0 };
        cam.fovy = 40.0f; cam.projection = CAMERA_PERSPECTIVE;
        BeginMode3D(cam);
        const Vector3 light = Vector3Normalize({ -0.6f, 0.35f, 0.75f });
        SetShaderValue(sphereShader, lightLoc, &light, SHADER_UNIFORM_VEC3);
        rlDisableDepthTest(); // a convex sphere: back faces culled, no depth needed (the menu's GL state stays as it was)
        // GenMeshSphere puts its poles on the z axis: turned onto y first (the map's poles at the top and bottom, its
        // pinch out of sight), then spun about y and tilted toward the viewer a little.
        sphere.transform = MatrixMultiply(MatrixMultiply(MatrixRotateX(-PI * 0.5f), MatrixRotateY(spin)), MatrixRotateX(0.35f));
        DrawModel(sphere, { 0, 0, 0 }, 1.0f, WHITE);
        rlEnableDepthTest();
        EndMode3D();
        EndTextureMode();
        // The render texture is upside down in screen space.
        const float side = std::min(area.width, area.height);
        const Rectangle dst = { area.x + (area.width - side) * 0.5f, area.y + (area.height - side) * 0.5f, side, side };
        DrawTexturePro(sphereRT.texture, { 0, 0, (float)kSphereRT, -(float)kSphereRT }, dst, { 0, 0 }, 0.0f, WHITE);
        const float s = Scale();
        if (!haveMap || building.load()) {
            const char* msg = haveMap ? "updating..." : "shaping the preview...";
            DrawText(msg, (int)(dst.x + (dst.width - MeasureText(msg, (int)(20 * s))) * 0.5f), (int)(dst.y + dst.height - 30 * s), (int)(20 * s), Fade(kText, 0.7f));
        }
    }

    void Changed() { changedAt = GetTime(); dirty = true; if (building.load()) cancel.store(true); }

    void RollSeed() {
        std::random_device rd;
        seed = rd();
        if (seed == 0) seed = 1; // 0 means "random" to the world generator
        seedText = std::to_string(seed);
        Changed();
    }
}

void NewWorld_Open() {
    InitGpu();
    params = WorldGenParams{};
    seedFocus = false;
    haveMap = false;
    RollSeed();
    changedAt = 0.0; // build at once
}

NewWorldChoice NewWorld_Draw() {
    NewWorldChoice r;
    const float s = Scale();
    const float W = (float)GetScreenWidth(), H = (float)GetScreenHeight();
    DrawRectangleGradientV(0, 0, (int)W, (int)H, { 4, 6, 14, 255 }, { 14, 22, 40, 255 });
    // A few fixed stars behind the planet.
    for (int i = 0; i < 160; ++i) {
        const uint32_t h = (uint32_t)i * 2654435761u;
        DrawPixel((int)((h % 10007) / 10007.0f * W), (int)(((h >> 8) % 9973) / 9973.0f * H), Fade(WHITE, 0.25f + 0.5f * ((h >> 3) % 7) / 7.0f));
    }
    const int titleFs = (int)(52 * s);
    DrawText("NEW WORLD", (int)(W * 0.06f), (int)(H * 0.07f), titleFs, kText);

    // The panel: seed, parameters, buttons.
    const float px = W * 0.06f, pw = std::min(W * 0.42f, 720.0f * s), rowH = 46.0f * s, gap = 12.0f * s;
    float y = H * 0.18f;
    DrawRectangleRounded({ px - 20 * s, y - 20 * s, pw + 40 * s, H * 0.74f }, 0.04f, 6, kPanel);
    {   // Seed: digits, editable; Random rolls a new one.
        const int fs = (int)(rowH * 0.5f);
        DrawText("Seed", (int)px, (int)(y + (rowH - fs) * 0.5f), fs, kText);
        const Rectangle field = { px + pw * 0.42f, y, pw * 0.36f, rowH };
        const bool hover = CheckCollisionPointRec(GetMousePosition(), field);
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) seedFocus = hover;
        DrawRectangleRec(field, seedFocus ? kHover : kButton);
        if (seedFocus) DrawRectangleLinesEx(field, 2.0f, SKYBLUE);
        const bool caret = seedFocus && fmodf((float)GetTime(), 1.0f) < 0.5f;
        DrawText(TextFormat("%s%s", seedText.c_str(), caret ? "_" : ""), (int)(field.x + 10 * s), (int)(field.y + (rowH - fs) * 0.5f), fs, kText);
        if (seedFocus) {
            for (int ch = GetCharPressed(); ch > 0; ch = GetCharPressed())
                if (ch >= '0' && ch <= '9' && seedText.size() < 10) seedText.push_back((char)ch);
            if (IsKeyPressed(KEY_BACKSPACE) && !seedText.empty()) seedText.pop_back();
            const unsigned long long v = seedText.empty() ? 0ull : std::stoull(seedText);
            const unsigned int ns = (unsigned int)std::min<unsigned long long>(v, 4294967295ull);
            if (ns != seed && ns != 0) { seed = ns; Changed(); }
            if (IsKeyPressed(KEY_ENTER)) seedFocus = false;
        }
        if (Button({ px + pw * 0.8f, y, pw * 0.2f, rowH }, "Random")) RollSeed();
        y += rowH + gap * 2;
    }
    const Rectangle row = { px, y, pw, rowH };
    auto next = [&]() { Rectangle r = row; r.y = y; y += rowH + gap; return r; };
    if (Slider(next(), "Ocean", params.ocean, -1.0f, 1.0f, "%+.2f")) Changed();
    if (Slider(next(), "Mountain height", params.mountains, 0.3f, 2.0f, "%.2fx")) Changed();
    if (Slider(next(), "Roughness", params.roughness, 0.3f, 2.0f, "%.2fx")) Changed();
    if (Slider(next(), "Biome size", params.biomeSize, 0.5f, 2.0f, "%.2fx")) Changed();
    if (Slider(next(), "Resources", params.resources, 0.5f, 2.0f, "%.2fx")) Changed();
    {   // Erosion on/off and its strength (not in the preview: the bake runs when the world is made).
        const Rectangle r = next();
        const int fs = (int)(rowH * 0.5f);
        DrawText("Erosion", (int)r.x, (int)(r.y + (rowH - fs) * 0.5f), fs, kText);
        if (Button({ r.x + pw * 0.42f, r.y + rowH * 0.1f, rowH * 2.2f, rowH * 0.8f }, params.erosion ? "On" : "Off")) { params.erosion = !params.erosion; Changed(); }
    }
    if (params.erosion && Slider(next(), "Erosion strength", params.erosionStrength, 0.5f, 2.0f, "%.2fx")) Changed();
    if (params.erosion && Slider(next(), "Cliffs and canyons", params.cliffs, 0.0f, 2.0f, "%.2fx")) Changed(); // #589: resistant rock layers (Plan 9)
    else if (!params.erosion) y += rowH + gap;
    if (params.erosion) {   // Plan 9 landforms (#592, #593), one row
        const Rectangle r = next();
        const int fs = (int)(rowH * 0.5f);
        const float ty = r.y + (rowH - fs) * 0.5f;
        DrawText("Karst", (int)r.x, (int)ty, fs, kText);
        if (Button({ r.x + pw * 0.42f, r.y + rowH * 0.1f, rowH * 2.2f, rowH * 0.8f }, params.karst ? "On" : "Off")) { params.karst = !params.karst; Changed(); }
        DrawText("Volcanoes", (int)(r.x + pw * 0.42f + rowH * 2.6f), (int)ty, fs, kText);
        if (Button({ r.x + pw * 0.42f + rowH * 2.6f + MeasureText("Volcanoes", fs) + rowH * 0.4f, r.y + rowH * 0.1f, rowH * 2.2f, rowH * 0.8f },
                   params.volcanoes ? "On" : "Off")) { params.volcanoes = !params.volcanoes; Changed(); }
    } else y += rowH + gap;
    const int noteFs = (int)(17 * s);
    DrawText("The preview shows the base terrain; erosion shapes it when the world is created (about a minute).", (int)px, (int)y, noteFs, Fade(kText, 0.6f));
    y += noteFs + gap * 2;
    if (Button({ px, y, pw * 0.3f, rowH }, "Defaults")) { params = WorldGenParams{}; Changed(); }
    const float by = H * 0.18f + H * 0.74f - rowH - 40 * s;
    if (Button({ px, by, pw * 0.46f, rowH * 1.15f }, "Create World")) { r.kind = NewWorldChoice::Create; r.seed = seed; r.params = params; }
    if (Button({ px + pw * 0.54f, by, pw * 0.46f, rowH * 1.15f }, "Cancel") || (IsKeyPressed(KEY_ESCAPE) && !seedFocus)) r.kind = NewWorldChoice::Cancel;

    // The preview, rebuilt 0.35 s after the last change.
    if (dirty && !building.load() && GetTime() - changedAt > 0.35) StartBuild();
    DrawPreview({ W * 0.5f, H * 0.1f, W * 0.46f, H * 0.82f });
    const int capFs = (int)(18 * s);
    const char* cap = TextFormat("seed %u", seed);
    DrawText(cap, (int)(W * 0.73f - MeasureText(cap, capFs) * 0.5f), (int)(H * 0.93f), capFs, Fade(kText, 0.6f));

    if (r.kind != NewWorldChoice::None && building.load()) cancel.store(true); // leaving: stop the preview
    return r;
}

void NewWorld_Apply(const NewWorldChoice& c) {
    worldGen.seed = c.seed;
    WorldGenParams_Save(WorldSavePath(), c.params); // read back by generateWorld (and by every later load of this save)
}

double NewWorld_PreviewBuildSeconds() { return lastBuildSeconds.load(); }
void NewWorld_TestSetParams(const WorldGenParams& p) { params = p; Changed(); }

void NewWorld_Shutdown() {
    cancel.store(true);
    if (worker.joinable()) worker.join();
}
