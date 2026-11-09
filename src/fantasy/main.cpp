// Fantasy open-world survival prototype: heightmap + cave terrain on the shared voxel engine,
// first-person exploration, and drawn-seal magic (see spells.h).
#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>

#include "audio.h"
#include "debugui.h"
#include "player.h"
#include "profiler.h"
#include "rlImGui.h"
#include "scene.h"
#include "terraintool.h"
#include "utils/logger.h"
#include "world.h"
#include "spells.h"
#include "sealrecognizer.h"

Logger logger("log_fantasy.txt");

namespace {
    // Spell name -> (core, modifiers) for tests.
    const std::map<std::string, std::pair<std::string, std::vector<std::string>>> kSpellRecipes = {
        { "fireball", { "flame", { "arrow" } } }, { "flare", { "flame", { "orb" } } },     { "ember", { "flame", {} } },
        { "wall", { "earth", { "line" } } },      { "wallbreaker", { "earth", { "arrow" } } }, { "mound", { "earth", {} } },
        { "mending", { "water", {} } },           { "gale", { "wind", { "arrow" } } },     { "featherfall", { "wind", { "orb" } } },
        { "updraft", { "wind", {} } },            { "lantern", { "light", { "orb" } } },   { "sunlance", { "light", { "arrow" } } },
        { "glow", { "light", {} } },
    };

    // Recogniser accuracy on synthetic seals at increasing sloppiness (no input needed).
    void RunRecognizerTest() {
        std::mt19937 rng(12345);
        const float levels[] = { 0.0f, 0.4f, 0.8f, 1.2f };
        for (float sl : levels) {
            int total = 0, ringOk = 0, coreOk = 0, modsOk = 0, valid = 0;
            float qualitySum = 0.0f;
            for (const auto& kv : kSpellRecipes) {
                for (int rep = 0; rep < 8; ++rep) {
                    const int dots = rep % 3;
                    auto strokes = SynthesizeSeal(kv.second.first, kv.second.second, dots, sl, { 960, 540 }, 300.0f, rng);
                    SealAnalysis a = AnalyzeSeal(strokes);
                    total++;
                    if (a.ring.found) ringOk++;
                    if (a.core >= 0 && a.glyphs[a.core].sign == kv.second.first) coreOk++;
                    std::vector<std::string> mods;
                    for (int m : a.modifiers) mods.push_back(a.glyphs[m].sign);
                    if (mods == kv.second.second) modsOk++;
                    if (a.problem.empty() && a.quality >= 0.45f) valid++;
                    qualitySum += a.quality;
                }
            }
            printf("bench: recognizer sloppiness %.1f: ring %d/%d, core %d/%d, modifiers %d/%d, castable %d/%d, avg quality %.2f\n",
                   sl, ringOk, total, coreOk, total, modsOk, total, valid, total, qualitySum / total);
        }
    }
}

int main(int argc, char** argv) {
    double benchSeconds = 0.0;
    unsigned int seed = 0;
    std::string testSpell;
    float sloppiness = 0.2f;
    bool testRecognizer = false;
    std::string grimoireDemo;   // --grimoire-demo SPELL: screenshot the canvas with that seal drawn
    Vector3 camOffset = { 0, 0, 0 };
    bool camSet = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--bench") benchSeconds = (i + 1 < argc) ? atof(argv[i + 1]) : 10.0;
        if (a == "--seed" && i + 1 < argc) seed = (unsigned int)strtoul(argv[i + 1], nullptr, 10);
        if (a == "--test-spell" && i + 1 < argc) testSpell = argv[i + 1];
        if (a == "--sloppy" && i + 1 < argc) sloppiness = (float)atof(argv[i + 1]);
        if (a == "--test-recognizer") testRecognizer = true;
        if (a == "--grimoire-demo" && i + 1 < argc) grimoireDemo = argv[i + 1];
        if (a == "--ui") DebugUI_Toggle();
        if (a == "--cam" && i + 3 < argc) { camOffset = { (float)atof(argv[i + 1]), (float)atof(argv[i + 2]), (float)atof(argv[i + 3]) }; camSet = true; }
    }
    const bool bench = benchSeconds > 0.0;
    benchMode = bench;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    if (bench) { InitWindow(1920, 1080, "Inkwild bench"); SetTargetFPS(0); }
    else {
        SetConfigFlags(FLAG_VSYNC_HINT | FLAG_WINDOW_RESIZABLE);
        InitWindow(1600, 900, "Inkwild");
        SetTargetFPS(144);
    }
    SetExitKey(KEY_NULL);
    SetTraceLogLevel(LOG_WARNING);
    rlSetClipPlanes(0.3, 6000.0);
    DebugUI_Init();
    if (!bench) Audio_Init();
    if (testRecognizer) RunRecognizerTest();

    {
        Scene world("fantasy", true);
        world.materialPalette = 1;
        world.drawSkybox = false;
        world.ambientColor = { 0.30f, 0.33f, 0.40f };
        world.sunDirection = Vector3Normalize({ 0.45f, 0.8f, 0.3f });
        world.shadowHalfExtent = 90.0f;
        world.texScale = 0.08f;
        world.fogColor = { 0.55f, 0.68f, 0.85f };
        world.fogDensity = 0.0022f;
        world.camera.fovy = 65.0f;
        world.camera.up = { 0, 1, 0 };
        world.camera.projection = CAMERA_PERSPECTIVE;
        world.loadTerrainTextures();

        if (seed == 0) seed = (unsigned int)GetRandomValue(1, 1 << 30);
        auto terrain = std::make_shared<TerrainField>(seed);
        auto body = std::make_unique<Planetoid>("terrain", DVec3(0, 0, 0), terrain);
        Planetoid* ground = body.get();
        body->parent = &world.rootObject;
        world.objects.push_back(std::move(body));
        world.syncRenderPositions();
        logger.logf("[Fantasy] terrain seed %u\n", seed);

        player.spawnOn(world, *ground);
        player.applyToCamera(world);
        terrainTool.creative = true; // building is free here; magic is the focus
        world.onDraw3D = [&world] { TerrainTool_Draw3D(); Spells_Draw3D(world); };
        if (bench && camSet) { player.flyMode = true; world.camera.target = world.camera.position; world.camera.position = Vector3Add(world.camera.position, camOffset); }

        if (!grimoireDemo.empty()) {
            auto it = kSpellRecipes.find(grimoireDemo);
            std::mt19937 rng(seed);
            if (it != kSpellRecipes.end()) Spells_ShowInGrimoire(SynthesizeSeal(it->second.first, it->second.second, 2, sloppiness, { 960, 540 }, 300.0f, rng));
        }
        double benchStart = -1.0, castTime = -1.0;
        long benchFrames = 0;
        double benchMaxFrame = 0.0;
        bool spellCast = false;
        int cratersBefore = 0, editsBefore = 0;

        while (!WindowShouldClose()) {
            const float dt = GetFrameTime();
            if (bench) {
                if (benchStart < 0.0 && WorldIdle()) {
                    if (!testSpell.empty() && !spellCast) {
                        // Scripted cast: synthesise a seal for the spell, face some terrain, cast it.
                        auto it = kSpellRecipes.find(testSpell);
                        if (it == kSpellRecipes.end()) { printf("bench: unknown spell %s\n", testSpell.c_str()); testSpell.clear(); }
                        else {
                            std::mt19937 rng(seed);
                            const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                            player.forward = Vector3Normalize(Vector3Add(flat, Vector3Scale(player.up, -0.25f)));
                            player.applyToCamera(world);
                            cratersBefore = spells.craters;
                            editsBefore = (int)ground->edits->size();
                            const Vector3 posBefore = (player.position - world.renderOrigin).toFloat();
                            auto strokes = SynthesizeSeal(it->second.first, it->second.second, 1, sloppiness, { 960, 540 }, 300.0f, rng);
                            SpellResult r = Spells_Cast(world, strokes);
                            printf("bench: cast '%s' (sloppiness %.1f): %s | %s\n", testSpell.c_str(), sloppiness, r.success ? "ok" : "failed", r.message.c_str());
                            (void)posBefore;
                            spellCast = true;
                            castTime = GetTime();
                        }
                    } else if (spellCast && castTime > 0.0 && (Spells_ActiveProjectiles() == 0 || GetTime() - castTime > 6.0) && GetTime() - castTime > 0.5) {
                        printf("bench: after cast: %d new terrain edits, %d craters, mana %.0f, HP %.0f, lights on %d\n",
                               (int)ground->edits->size() - editsBefore, spells.craters - cratersBefore, spells.mana, spells.health,
                               [&] { int n = 0; for (size_t i = 1; i < world.lights.size(); ++i) n += world.lights[i].enabled; return n; }());
                        castTime = -1.0;
                    } else if (castTime < 0.0) {
                        benchStart = GetTime();
                        printf("bench: world ready (%d meshes, %d nodes, avg build %.2f ms)\n", GetWorldStats().chunks, GetWorldStats().nodes, GetWorldStats().avgBuildMs);
                    }
                }
                const double t = benchStart < 0.0 ? 0.0 : GetTime() - benchStart;
                if (benchStart >= 0.0 && t > 1.0) { benchFrames++; benchMaxFrame = std::max(benchMaxFrame, (double)dt); }
                if (benchStart >= 0.0 && t > 1.0 + benchSeconds) {
                    printf("bench: %ld frames in %.1fs | avg %.1f FPS | worst %.2f ms\n", benchFrames, benchSeconds, benchFrames / benchSeconds, benchMaxFrame * 1000.0);
                    TakeScreenshot("bench.png");
                    break;
                }
                if (GetTime() > 120.0) { printf("bench: timeout\n"); break; }
            }

            if (!bench && IsKeyPressed(KEY_F1)) DebugUI_Toggle();
            if (!bench && IsKeyPressed(KEY_ESCAPE) && !spells.grimoireOpen) break;
            const bool uiOpen = DebugUI_IsOpen() || spells.grimoireOpen;
            const bool gameInput = !bench && !uiOpen && IsWindowFocused();
            if (!bench) { if (uiOpen) ShowCursor(); else HideCursor(); }
            world.cameraInputEnabled = player.flyMode && !uiOpen;

            player.update(world, gameInput, dt);
            if (player.impactSpeed > 12.0f) spells.health = std::max(0.0f, spells.health - (player.impactSpeed - 12.0f) * 6.0f); // fall damage
            if (spells.health <= 0.0f) { spells.health = 100.0f; player.spawnOn(world, *ground); spells.last.message = "You blacked out and woke at the start"; spells.lastTime = GetTime(); }
            player.applyToCamera(world);
            Spells_Update(world, dt, !bench && !DebugUI_IsOpen());
            world.updateFloatingOrigin();
            worldUpdate(world);
            TerrainTool_Update(world, gameInput);

            BeginDrawing();
                ClearBackground({ 140, 175, 215, 255 });
                world.drawScene(1);
                world.drawUI(1);
                const int cx = GetScreenWidth() / 2, cy = GetScreenHeight() / 2;
                DrawLine(cx - 8, cy, cx + 8, cy, RAYWHITE);
                DrawLine(cx, cy - 8, cx, cy + 8, RAYWHITE);
                Spells_DrawOverlay(world);
                rlImGuiBegin();
                DebugUI_Draw(world);
                rlImGuiEnd();
                if (!bench && !uiOpen) SetMousePosition(GetScreenWidth() / 2, GetScreenHeight() / 2);
            EndDrawing();
            PROFILE_FRAME();
        }
    }
    DebugUI_Shutdown();
    Audio_Shutdown();
    CloseWindow();
    return 0;
}
