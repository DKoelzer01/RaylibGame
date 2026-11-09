#include "terraintool.h"
#include "world.h"
#include "inventory.h"
#include "audio.h"
#include "player.h"

#include <raymath.h>
#include <algorithm>
#include <string>

TerrainToolSettings terrainTool;

namespace {
    TerrainHit hit;
    bool hasHit = false;
    bool digging = false, adding = false;
    float strokeTimer = 0.0f;
    std::string toast;
    double toastUntil = 0.0;

    void Toast(const std::string& msg) {
        toast = msg;
        toastUntil = GetTime() + 2.5;
    }
}

void TerrainTool_Update(Scene& world, bool inputEnabled) {
    PROFILE_FUNCTION();
    const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
    hasHit = RaycastTerrain(world, Ray{ world.camera.position, dir }, player.flyMode ? 150.0f : terrainTool.reach, hit);

    if (!inputEnabled) { digging = adding = false; return; }

    if (IsKeyPressed(KEY_F5)) Toast(SaveWorld(world) ? "Saved to saves/world.sav" : "Save failed");
    if (IsKeyPressed(KEY_F9)) {
        if (LoadWorldSave(world)) { Toast("Loaded saves/world.sav"); hasHit = false; return; }
        Toast("No save found");
    }

    const float wheel = GetMouseWheelMove();
    if (wheel != 0.0f) terrainTool.brushRadius = std::clamp(terrainTool.brushRadius * (wheel > 0 ? 1.15f : 1.0f / 1.15f), 1.0f, 24.0f);

    digging = IsMouseButtonDown(MOUSE_BUTTON_LEFT);
    adding = !digging && IsMouseButtonDown(MOUSE_BUTTON_RIGHT);
    if (!(digging || adding) || !hasHit) { strokeTimer = 0.0f; return; }

    // Fixed stroke rate while held, so editing speed does not depend on frame rate.
    strokeTimer -= GetFrameTime();
    if (strokeTimer > 0.0f) return;
    strokeTimer = 1.0f / terrainTool.strokesPerSecond;
    TerrainTool_Stroke(digging);
}

// One brush stroke at the current crosshair hit (also used by scripted tests).
bool TerrainTool_Stroke(bool dig) {
    if (!hasHit) return false;
    const bool digging = dig;
    Planetoid& p = *hit.planetoid;
    const float r3 = terrainTool.brushRadius * terrainTool.brushRadius * terrainTool.brushRadius;
    if (digging) {
        // Yield scales with the volume the stroke removes (~strength * r^3).
        inventory.addFractional(ItemForMaterial(hit.material), terrainTool.digStrength * r3 * 0.15f);
        Audio_Play(Sfx::Dig, 0.5f, 0.85f + 0.3f * (float)GetRandomValue(0, 100) / 100.0f);
    } else if (!terrainTool.creative) {
        // Building consumes bulk material, best first.
        static const char* fillers[] = { "concrete", "stone", "regolith" };
        bool paid = false;
        for (const char* f : fillers) if (inventory.remove(f, 1)) { paid = true; break; }
        if (!paid) {
            Toast("Need stone, regolith or lunacrete to build");
            Audio_Play(Sfx::Error, 0.5f);
            strokeTimer = 0.5f;
            return false;
        }
        Audio_Play(Sfx::Build, 0.6f);
    }
    TerrainEdit e;
    e.radius = terrainTool.brushRadius / p.scale;
    e.strength = digging ? -terrainTool.digStrength : terrainTool.digStrength;
    // Centre the brush slightly into (dig) or out of (add) the surface so strokes bite immediately.
    const Vector3 localNormal = Vector3Normalize(Vector3Subtract(p.toLocal(Vector3Add(hit.point, hit.normal)), hit.local));
    e.center = Vector3Add(hit.local, Vector3Scale(localNormal, (digging ? -0.35f : 0.35f) * e.radius));
    ApplyTerrainEdit(p, e);
    return true;
}

void TerrainTool_Draw3D() {
    if (!hasHit || hit.distance < terrainTool.brushRadius * 1.1f) return; // don't engulf the camera
    const Color c = digging ? RED : adding ? GREEN : Fade(WHITE, 0.6f);
    DrawSphereWires(hit.point, terrainTool.brushRadius, 8, 12, c);
}

void TerrainTool_DrawHUD() {
    const int cx = GetScreenWidth() / 2, cy = GetScreenHeight() / 2;
    DrawLine(cx - 8, cy, cx + 8, cy, RAYWHITE);
    DrawLine(cx, cy - 8, cx, cy + 8, RAYWHITE);

    const int y = GetScreenHeight() - 60;
    if (hasHit) {
        DrawText(TextFormat("%s  %.1f m   brush %.1f", VoxelMaterialName(hit.material), hit.distance, terrainTool.brushRadius),
                 cx + 14, cy + 10, 20, RAYWHITE);
    }
    DrawText("LMB dig  RMB build  wheel brush  Tab inventory  L lamp  F fly  F5 save  F9 load  F1 debug", 10, y + 30, 20, Fade(RAYWHITE, 0.7f));
    if (GetTime() < toastUntil) DrawText(toast.c_str(), 10, y, 20, YELLOW);
}
