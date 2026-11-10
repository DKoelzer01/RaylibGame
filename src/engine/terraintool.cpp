#include "terraintool.h"
#include "world.h"
#include "inventory.h"
#include "audio.h"
#include "player.h"
#include "worldtypes.h"

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
    // The plant under the crosshair (closer than the terrain hit), and progress chopping it.
    bool hasPlant = false;
    PlantHit plantHit{};
    Planetoid* plantBody = nullptr;
    uint64_t chopId = 0;
    int chopHits = 0;

    void Toast(const std::string& msg);
    int HitsToFell(const Plant& p) { return p.tree() ? 6 + (int)(p.scale * 4.0f) : 1; }

    // One swing at the aimed plant; true when it came down (items into the inventory).
    bool Chop() {
        if (!hasPlant || !plantBody || !plantBody->trees) return false;
        const Plant& p = plantHit.plant;
        if (p.id != chopId) { chopId = p.id; chopHits = 0; }
        Audio_Play(Sfx::Dig, 0.45f, p.tree() ? 0.6f : 1.4f);
        if (++chopHits < HitsToFell(p)) return false;
        plantBody->trees->remove(p.id);
        if (p.tree()) {
            const int wood = 3 + (int)(p.scale * 4.0f);
            inventory.add("wood", wood);
            Toast(TextFormat("+%d wood", wood));
            Audio_Play(Sfx::Build, 0.6f, 0.5f);
        } else {
            inventory.add("fibre", 2);
            const bool berries = (p.id & 1) == 0;
            if (berries) inventory.add("berries", 1);
            Toast(berries ? "+2 fibre  +1 berries" : "+2 fibre");
            Audio_Play(Sfx::Pickup, 0.5f);
        }
        chopId = 0; chopHits = 0;
        hasPlant = false;
        return true;
    }

    void Toast(const std::string& msg) {
        toast = msg;
        toastUntil = GetTime() + 2.5;
    }
}

void TerrainTool_Update(Scene& world, bool inputEnabled, bool editEnabled) {
    PROFILE_FUNCTION();
    const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
    hasHit = RaycastTerrain(world, Ray{ world.camera.position, dir }, player.flyMode ? 150.0f : terrainTool.reach, hit);
    hasPlant = false;
    plantBody = player.nearest;
    if (plantBody && plantBody->trees && !player.flyMode) {
        const Vector3 lo = plantBody->toLocal(world.camera.position), ld = plantBody->dirToLocal(dir);
        const float reach = std::min(terrainTool.reach, 6.0f); // arm's length for an axe
        hasPlant = plantBody->trees->raycast(lo, ld, hasHit ? std::min(reach, hit.distance) : reach, plantHit);
    }
    if (!hasPlant) chopHits = 0;

    if (!inputEnabled) { digging = adding = false; return; }

    if (IsKeyPressed(KEY_F5)) Toast(SaveWorld(world) ? TextFormat("Saved to %s", WorldSavePath()) : "Save failed");
    if (IsKeyPressed(KEY_F9)) {
        if (LoadWorldSave(world)) { Toast("Loaded saves/world.sav"); hasHit = false; return; }
        Toast("No save found");
    }
    if (!editEnabled) { digging = adding = false; return; } // e.g. piloting: save/load only

    const float wheel = GetMouseWheelMove();
    if (wheel != 0.0f) terrainTool.brushRadius = std::clamp(terrainTool.brushRadius * (wheel > 0 ? 1.15f : 1.0f / 1.15f), 1.0f, 24.0f);

    digging = IsMouseButtonDown(MOUSE_BUTTON_LEFT);
    adding = !digging && IsMouseButtonDown(MOUSE_BUTTON_RIGHT);
    if (!(digging || adding) || !(hasHit || (digging && hasPlant))) { strokeTimer = 0.0f; return; }

    // Fixed stroke rate while held, so editing speed does not depend on frame rate: keep the
    // timer's remainder and run every interval that elapsed this frame (re-aiming between catch-up
    // strokes so they don't hit an already-dug spot). Catch-up is capped after a long hitch.
    strokeTimer -= GetFrameTime();
    const float interval = 1.0f / terrainTool.strokesPerSecond;
    for (int n = 0; strokeTimer <= 0.0f && n < 8; ++n) {
        if (n > 0) hasHit = RaycastTerrain(world, Ray{ world.camera.position, dir }, player.flyMode ? 150.0f : terrainTool.reach, hit);
        if (!hasHit && !hasPlant) break;
        if (digging && hasPlant) { // ~8 swings a second; a short pause after it comes down
            if (Chop()) { strokeTimer = 0.4f; break; }
            strokeTimer += 0.12f;
            continue;
        }
        TerrainTool_Stroke(digging);
        strokeTimer += interval;
    }
    strokeTimer = std::max(strokeTimer, 0.0f);
}

// Sand, gravel or dirt on a river or pond bed (planets).
static bool WetBed(Planetoid& p, Vector3 local, VoxelMaterial m) {
    if (m != VoxelMaterial::Sand && m != VoxelMaterial::Gravel && m != VoxelMaterial::Regolith) return false;
    const auto* pf = dynamic_cast<const PlanetField*>(p.field.get());
    return pf && pf->wetBed(local);
}

// One brush stroke at the current crosshair hit (also used by scripted tests).
bool TerrainTool_Stroke(bool dig) {
    if (!hasHit) return false;
    const bool digging = dig;
    Planetoid& p = *hit.planetoid;
    if (!p.field->editable(hit.local)) {
        // A planet's mantle: edits there would have no effect, so refuse before charging anything.
        Toast("The mantle is too hard to dig into");
        Audio_Play(Sfx::Error, 0.5f);
        strokeTimer = 0.5f;
        return false;
    }
    const float r3 = terrainTool.brushRadius * terrainTool.brushRadius * terrainTool.brushRadius;
    if (digging) {
        // Yield scales with the volume the stroke removes (~strength * r^3).
        const char* item = hit.material == VoxelMaterial::Ore ? OreKindAt(hit.local, hit.planetoid->seed) : ItemForMaterial(hit.material, HasSoil(*hit.planetoid));
        inventory.addFractional(item, terrainTool.digStrength * r3 * 0.15f);
        if (WetBed(p, hit.local, hit.material)) inventory.addFractional("clay", terrainTool.digStrength * r3 * 0.06f); // river and pond beds hold clay
        Audio_Play(Sfx::Dig, 0.5f, 0.85f + 0.3f * (float)GetRandomValue(0, 100) / 100.0f);
    } else if (!terrainTool.creative) {
        // Building consumes bulk material, best first.
        static const char* fillers[] = { "concrete", "stone", "dirt", "regolith", "gravel", "sand" };
        bool paid = false;
        for (const char* f : fillers) if (inventory.remove(f, 1)) { paid = true; break; }
        if (!paid) {
            Toast("Need stone, dirt, regolith, gravel, sand or concrete to build");
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
    if (!hasHit || hasPlant || hit.distance < terrainTool.brushRadius * 1.1f) return; // don't engulf the camera
    const Color c = digging ? RED : adding ? GREEN : Fade(WHITE, 0.6f);
    DrawSphereWires(hit.point, terrainTool.brushRadius, 8, 12, c);
}

void TerrainTool_DrawHUD() {
    const int cx = GetScreenWidth() / 2, cy = GetScreenHeight() / 2;
    DrawLine(cx - 8, cy, cx + 8, cy, RAYWHITE);
    DrawLine(cx, cy - 8, cx, cy + 8, RAYWHITE);

    const int y = GetScreenHeight() - 60;
    if (hasPlant) {
        const Plant& p = plantHit.plant;
        DrawText(TextFormat("%s  %.1f m   LMB %s", p.tree() ? "tree" : "bush", plantHit.distance, p.tree() ? "chop" : "harvest"), cx + 14, cy + 10, 20, RAYWHITE);
        if (chopId == p.id && chopHits > 0) {
            const float f = (float)chopHits / HitsToFell(p);
            DrawRectangle(cx + 14, cy + 34, 120, 6, Fade(BLACK, 0.5f));
            DrawRectangle(cx + 14, cy + 34, (int)(120 * f), 6, ORANGE);
        }
    } else if (hasHit) {
        DrawText(TextFormat("%s  %.1f m   brush %.1f", hit.material == VoxelMaterial::Regolith && HasSoil(*hit.planetoid) ? "dirt" : VoxelMaterialName(hit.material), hit.distance, terrainTool.brushRadius),
                 cx + 14, cy + 10, 20, RAYWHITE);
    }
    DrawText("LMB dig  RMB fill  wheel brush  G build  1-0 hotbar  Tab inventory  T telescope  L lamp  F fly  F5 save  F1 debug", 10, y + 30, 20, Fade(RAYWHITE, 0.7f));
    if (GetTime() < toastUntil) DrawText(toast.c_str(), 10, y, 20, YELLOW);
}

bool TerrainTool_Chop() { return Chop(); }
bool TerrainTool_AimingPlant(uint64_t* id) {
    if (hasPlant && id) *id = plantHit.plant.id;
    return hasPlant;
}

// Scripted gathering (tests): dig at a body-local point with the tool's yield rule (one item of
// whatever is there per call), and chop / harvest the nearest plant.
bool TerrainTool_DigAt(Planetoid& p, Vector3 local) {
    if (!p.field->editable(local)) return false;
    const float base = p.field->baseDensity(local);
    if (base < ISO_LEVEL) return false; // air
    const VoxelMaterial m = p.field->material(local, base);
    const char* item = m == VoxelMaterial::Ore ? OreKindAt(local, p.seed) : ItemForMaterial(m, HasSoil(p));
    inventory.add(item, 1);
    if (WetBed(p, local, m)) inventory.addFractional("clay", 0.4f); // beds hold clay
    ApplyTerrainEdit(p, TerrainEdit{ local, 1.2f, -0.6f });
    return true;
}

bool TerrainTool_HarvestNearest(Planetoid& p, Vector3 local, bool tree) {
    if (!p.trees) return false;
    Plant plant;
    if (!p.trees->nearestPlant(local, 400.0f, tree, plant)) return false;
    hasPlant = true;
    plantBody = &p;
    plantHit.plant = plant;
    plantHit.distance = 1.0f;
    for (int i = 0; i < 64; ++i) if (Chop()) return true;
    return false;
}
