#include "debugui.h"
#include "scene.h"
#include "world.h"
#include "profiler.h"
#include "terraintool.h"
#include "audio.h"

#include "imgui.h"
#include "rlImGui.h"

#include <cmath>

namespace {
    bool open = false;
    float frameTimes[240] = {};
    int frameIndex = 0;

    // Sun direction edited as azimuth/elevation (degrees); initialised from the scene on first draw.
    bool sunInit = false;
    float sunAzimuth = 0.0f, sunElevation = 0.0f;
}

void DebugUI_Init() {
    rlImGuiSetup(true);
    ImGui::GetIO().IniFilename = "imgui.ini";
}

void DebugUI_Shutdown() { rlImGuiShutdown(); }
bool DebugUI_IsOpen() { return open; }
void DebugUI_Toggle() { open = !open; }

void DebugUI_Draw(Scene& world) {
    // Always record frame times so the graph has history when the panel opens.
    frameTimes[frameIndex] = GetFrameTime() * 1000.0f;
    frameIndex = (frameIndex + 1) % IM_ARRAYSIZE(frameTimes);
    if (!open) return;
    PROFILE_FUNCTION();

    if (!sunInit) {
        const Vector3 d = world.sunDirection;
        sunAzimuth = atan2f(d.z, d.x) * RAD2DEG;
        sunElevation = asinf(d.y) * RAD2DEG;
        sunInit = true;
    }

    ImGui::SetNextWindowPos(ImVec2(10, 60), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(380, 0), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Debug (F1)")) {
        // --- Frame ---
        float worst = 0.0f, sum = 0.0f;
        for (float t : frameTimes) { worst = fmaxf(worst, t); sum += t; }
        ImGui::Text("%d FPS   avg %.2f ms   worst %.2f ms", GetFPS(), sum / IM_ARRAYSIZE(frameTimes), worst);
        ImGui::PlotLines("##frametimes", frameTimes, IM_ARRAYSIZE(frameTimes), frameIndex, "frame ms", 0.0f, 33.0f, ImVec2(-1, 60));
        ImGui::Text("Chunks drawn: %d", world.lastDrawCalls);

        // --- Rendering ---
        if (ImGui::CollapsingHeader("Rendering", ImGuiTreeNodeFlags_DefaultOpen)) {
            bool sunChanged = ImGui::SliderFloat("Sun azimuth", &sunAzimuth, -180.0f, 180.0f, "%.0f deg");
            sunChanged |= ImGui::SliderFloat("Sun elevation", &sunElevation, -89.0f, 89.0f, "%.0f deg");
            if (sunChanged) {
                const float az = sunAzimuth * DEG2RAD, el = sunElevation * DEG2RAD;
                world.sunDirection = { cosf(el) * cosf(az), sinf(el), cosf(el) * sinf(az) };
            }
            ImGui::SliderFloat("Shadow box half-size", &world.shadowHalfExtent, 25.0f, 500.0f, "%.0f");
            ImGui::SliderFloat("Texture scale", &world.texScale, 0.005f, 0.5f, "%.3f", ImGuiSliderFlags_Logarithmic);
            ImGui::ColorEdit3("Ambient", &world.ambientColor.x);
            ImGui::Checkbox("Debug overlay (F3)", &world.debugDraw);
            ImGui::SameLine();
            ImGui::Checkbox("Frustum culling", &frustumCullingEnabled);
            if (ImGui::Button("Reload shaders")) world.reloadShaders();
            ImGui::SameLine();
            ImGui::TextDisabled("(also reload automatically on save)");
            if (world.shadowMap.depth.id > 0 && ImGui::TreeNode("Shadow map")) {
                rlImGuiImageSize(&world.shadowMap.depth, 256, 256);
                ImGui::TreePop();
            }
        }

        // --- Camera ---
        if (ImGui::CollapsingHeader("Audio")) {
            float vol = Audio_GetMasterVolume();
            if (ImGui::SliderFloat("Master volume", &vol, 0.0f, 1.0f, "%.2f")) Audio_SetMasterVolume(vol);
        }

        if (ImGui::CollapsingHeader("Camera", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::DragFloat3("Position", &world.camera.position.x, 0.5f);
            ImGui::SliderFloat("Move speed", &cameraMoveSpeed, 1.0f, 10000.0f, "%.0f u/s", ImGuiSliderFlags_Logarithmic);
            const DVec3 cw = world.cameraWorld();
            ImGui::Text("World position: %.1f %.1f %.1f", cw.x, cw.y, cw.z);
            ImGui::Text("Render origin:  %.1f %.1f %.1f", world.renderOrigin.x, world.renderOrigin.y, world.renderOrigin.z);
            ImGui::SliderFloat("FOV", &world.camera.fovy, 20.0f, 120.0f, "%.0f deg");
        }

        // --- Chunk streaming ---
        if (ImGui::CollapsingHeader("Chunk streaming", ImGuiTreeNodeFlags_DefaultOpen)) {
            const WorldStats& s = GetWorldStats();
            ImGui::Text("Meshes drawn: %d   LOD nodes: %d   edits: %d", s.chunks, s.nodes, s.edits);
            ImGui::SliderFloat("LOD factor", &lodSettings.lodFactor, 0.5f, 6.0f, "%.2f");
            ImGui::Text("Jobs queued: %d   running: %d / %d threads", s.queuedJobs, s.runningJobs, s.workerThreads);
            ImGui::Text("Uploads: %d pending, %d last frame", s.pendingUploads, s.uploadsLastFrame);
            ImGui::Text("Avg chunk build: %.2f ms   last world load: %.2f s", s.avgBuildMs, s.lastGenerationSeconds);
        }

        // --- Terrain tool ---
        if (ImGui::CollapsingHeader("Terrain tool")) {
            ImGui::SliderFloat("Brush radius", &terrainTool.brushRadius, 1.0f, 24.0f, "%.1f");
            ImGui::SliderFloat("Strength", &terrainTool.digStrength, 0.02f, 1.0f, "%.2f");
            ImGui::SliderFloat("Strokes / s", &terrainTool.strokesPerSecond, 1.0f, 30.0f, "%.0f");
            if (ImGui::Button("Save (F5)")) SaveWorld(world);
            ImGui::SameLine();
            if (ImGui::Button("Load (F9)")) LoadWorldSave(world);
        }

        // --- World generation ---
        if (ImGui::CollapsingHeader("World generation", ImGuiTreeNodeFlags_DefaultOpen)) {
            int seed = (int)worldGen.seed;
            if (ImGui::InputInt("Seed (0 = random)", &seed)) worldGen.seed = (unsigned int)std::max(0, seed);
            ImGui::Text("Last seed used: %u", worldGen.lastSeed);
            ImGui::SameLine();
            if (ImGui::SmallButton("Reuse")) worldGen.seed = worldGen.lastSeed;
            ImGui::SliderInt("Planetoids", &worldGen.planetoidCount, 1, 32);
            ImGui::SliderFloat("System radius", &worldGen.systemRadius, 500.0f, 50000.0f, "%.0f", ImGuiSliderFlags_Logarithmic);
            ImGui::SliderFloat("Max spin (deg/s)", &worldGen.maxSpinDegPerSec, 0.0f, 20.0f, "%.1f");
            ImGui::SliderInt("Planetoid size", &worldGen.planetoidSize, 16, 250);
            ImGui::TextDisabled("Chunks build on worker threads, nearest first.");
            ImGui::DragFloat3("First position", &worldGen.firstPosition.x, 1.0f);
            if (ImGui::Button("Regenerate world")) regenerateWorld(world);
        }
    }
    ImGui::End();
}
