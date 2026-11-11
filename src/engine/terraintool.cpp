#include "terraintool.h"
#include "world.h"
#include "inventory.h"
#include "audio.h"
#include "player.h"
#include "worldtypes.h"
#include "tools.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cstring>
#include <map>
#include <string>

// The depth test is flipped to draw the mining sphere only where it lies behind the terrain (inside
// the material): opengl32 exports glDepthFunc itself (GL 1.1), so no loader is needed.
extern "C" __declspec(dllimport) void __stdcall glDepthFunc(unsigned int func);
static constexpr unsigned int kGlGreater = 0x0204, kGlLequal = 0x0203;

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
    // The tool in hand and what it does to the aimed material.
    ToolSpec tool;
    DigRule rule;
    std::string hitItem;
    double whyShownAt = -10.0;

    void Toast(const std::string& msg);
    int HitsToFell(const Plant& p) {
        if (!p.tree()) return 1;
        return std::max(1, (int)lroundf((6 + (int)(p.scale * 4.0f)) * ChopFactor(tool)));
    }

    ToolSpec HeldTool() {
        if (terrainTool.heldTool.empty() || inventory.count(terrainTool.heldTool) <= 0) return ToolSpec{};
        return ToolOf(terrainTool.heldTool);
    }

    const char* ItemAt(Planetoid& p, Vector3 local, VoxelMaterial m) {
        return m == VoxelMaterial::Ore ? OreKindOf(p, local) : ItemForMaterial(m, HasSoil(p));
    }

    // One swing at the aimed plant; true when it came down (items into the inventory).
    bool Chop() {
        if (!hasPlant || !plantBody || !plantBody->trees) return false;
        const Plant& p = plantHit.plant;
        if (p.id != chopId) { chopId = p.id; chopHits = 0; }
        Audio_Play(Sfx::Dig, 0.45f, p.tree() ? 0.6f : 1.4f);
        ++terrainTool.swings;
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
    tool = HeldTool();
    if (hasHit) {
        hitItem = ItemAt(*hit.planetoid, hit.local, hit.material);
        rule = DigRuleFor(tool, hit.material, hitItem.c_str());
    }

    if (!inputEnabled) { digging = adding = false; return; }

    if (IsKeyPressed(KEY_F5)) Toast(SaveWorld(world) ? TextFormat("Saved to %s", WorldSavePath()) : "Save failed");
    if (IsKeyPressed(KEY_F9)) {
        if (LoadWorldSave(world)) { Toast("Loaded saves/world.sav"); hasHit = false; return; }
        Toast("No save found");
    }
    if (!editEnabled) { digging = adding = false; return; } // e.g. piloting: save/load only

    if (terrainTool.terraform) {
        const float wheel = GetMouseWheelMove();
        if (wheel != 0.0f) terrainTool.brushRadius = std::clamp(terrainTool.brushRadius * (wheel > 0 ? 1.15f : 1.0f / 1.15f), 1.0f, 24.0f);
    }

    digging = IsMouseButtonDown(MOUSE_BUTTON_LEFT);
    adding = !digging && IsMouseButtonDown(MOUSE_BUTTON_RIGHT);
    if (!(digging || adding) || !(hasHit || (digging && hasPlant))) { strokeTimer = 0.0f; return; }

    // Fixed stroke rate while held, so editing speed does not depend on frame rate: keep the
    // timer's remainder and run every interval that elapsed this frame (re-aiming between catch-up
    // strokes so they don't hit an already-dug spot). Catch-up is capped after a long hitch.
    strokeTimer -= GetFrameTime();
    const float rate = terrainTool.terraform ? terrainTool.strokesPerSecond : std::max(rule.strokesPerSecond, 2.0f);
    const float interval = 1.0f / rate;
    for (int n = 0; strokeTimer <= 0.0f && n < 8; ++n) {
        if (n > 0) {
            hasHit = RaycastTerrain(world, Ray{ world.camera.position, dir }, player.flyMode ? 150.0f : terrainTool.reach, hit);
            if (hasHit) { hitItem = ItemAt(*hit.planetoid, hit.local, hit.material); rule = DigRuleFor(tool, hit.material, hitItem.c_str()); }
        }
        if (!hasHit && !hasPlant) break;
        if (digging && hasPlant) { // ~8 swings a second; a short pause after it comes down
            if (Chop()) { strokeTimer = 0.4f; break; }
            strokeTimer += 0.12f;
            continue;
        }
        if (!TerrainTool_Stroke(digging)) { strokeTimer = std::max(strokeTimer, 0.4f); break; }
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

// Building consumes bulk material, best first.
static bool PayFiller() {
    if (terrainTool.creative) return true;
    static const char* fillers[] = { "concrete", "stone", "dirt", "regolith", "gravel", "sand" };
    for (const char* f : fillers) if (inventory.remove(f, 1)) return true;
    Toast("Need stone, dirt, regolith, gravel, sand or concrete to build");
    Audio_Play(Sfx::Error, 0.5f);
    return false;
}

// The old brush (debug): any material, soft-edged, yield by the brush volume.
static bool BrushStroke(bool digging) {
    Planetoid& p = *hit.planetoid;
    const float r3 = terrainTool.brushRadius * terrainTool.brushRadius * terrainTool.brushRadius;
    if (digging) {
        inventory.addFractional(hitItem, terrainTool.digStrength * r3 * 0.15f);
        if (WetBed(p, hit.local, hit.material)) inventory.addFractional("clay", terrainTool.digStrength * r3 * 0.06f); // river and pond beds hold clay
        Audio_Play(Sfx::Dig, 0.5f, 0.85f + 0.3f * (float)GetRandomValue(0, 100) / 100.0f);
        ++terrainTool.swings;
    } else {
        if (!PayFiller()) return false;
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

// The mining sphere for the current hit: centred on the surface point, the tool's radius.
static float MiningRadius() { return rule.radius > 0.0f ? rule.radius : DigRuleFor(tool, VoxelMaterial::Regolith, nullptr).radius; }

// One stroke at the current crosshair hit (also used by scripted tests). Tools carve exactly the
// sphere shown (what lies inside it, of what the tool can dig, is collected); RMB places soil.
bool TerrainTool_Stroke(bool dig) {
    if (!hasHit) return false;
    Planetoid& p = *hit.planetoid;
    if (!p.field->editable(hit.local)) {
        // A planet's mantle: edits there would have no effect, so refuse before charging anything.
        Toast("The mantle is too hard to dig into");
        Audio_Play(Sfx::Error, 0.5f);
        return false;
    }
    if (terrainTool.terraform) return BrushStroke(dig);

    const float radius = MiningRadius();
    if (!dig) {
        if (!PayFiller()) return false;
        TerrainEdit e{};
        e.radius = radius * 1.3f / p.scale;
        e.strength = 0.6f;
        const Vector3 localNormal = Vector3Normalize(Vector3Subtract(p.toLocal(Vector3Add(hit.point, hit.normal)), hit.local));
        e.center = Vector3Add(hit.local, Vector3Scale(localNormal, 0.3f * e.radius));
        ApplyTerrainEdit(p, e);
        Audio_Play(Sfx::Build, 0.6f);
        return true;
    }
    if (!rule.can) {
        if (GetTime() - whyShownAt > 1.5) { Toast(rule.why); whyShownAt = GetTime(); }
        Audio_Play(Sfx::Error, 0.3f);
        return false;
    }
    // Collect what the sphere holds: sample it on a grid (solid samples by material, two items per
    // cubic metre), counting only what the tool can dig - hands and shovels leave rock in place.
    const bool softOnly = tool.type != ToolType::Pickaxe;
    const float r = radius / p.scale;
    const float step = r / 4.0f;
    const float cell = step * step * step * p.scale * p.scale * p.scale;
    std::map<std::string, float> got;
    float wet = 0.0f;
    for (int z = -4; z <= 4; ++z)
        for (int y = -4; y <= 4; ++y)
            for (int x = -4; x <= 4; ++x) {
                if (x * x + y * y + z * z > 16) continue;
                const Vector3 q = Vector3Add(hit.local, { x * step, y * step, z * step });
                const float d = p.densityAt(q);
                if (d < ISO_LEVEL || !p.field->editable(q)) continue;
                const VoxelMaterial m = p.field->material(q, d);
                if (softOnly && !IsSoft(m)) continue;
                const char* item = ItemAt(p, q, m);
                got[item] += cell * 2.0f;
                if (strcmp(item, "dirt") == 0) got["gravel"] += cell * 0.4f; // stones in the soil (enough for the first stone tools)
                if (WetBed(p, q, m)) wet += cell;
            }
    for (auto& kv : got) inventory.addFractional(kv.first, kv.second);
    if (wet > 0.0f) inventory.addFractional("clay", wet * 0.8f); // river and pond beds hold clay
    TerrainEdit e{};
    e.center = hit.local;
    e.radius = r;
    e.carve = true;
    e.softOnly = softOnly;
    ApplyTerrainEdit(p, e);
    Audio_Play(Sfx::Dig, 0.5f, 0.85f + 0.3f * (float)GetRandomValue(0, 100) / 100.0f);
    ++terrainTool.swings;
    return true;
}

void TerrainTool_Draw3D() {
    if (!hasHit || hasPlant) return;
    if (terrainTool.terraform) {
        if (hit.distance < terrainTool.brushRadius * 1.1f) return; // don't engulf the camera
        const Color c = digging ? RED : adding ? GREEN : Fade(WHITE, 0.6f);
        DrawSphereWires(hit.point, terrainTool.brushRadius, 8, 12, c);
        return;
    }
    // The mining sphere, only where it lies inside the material: drawn where it is behind the
    // terrain surface (depth test flipped), so it shows through the ground and never in the air.
    const float radius = MiningRadius();
    if (hit.distance < radius * 1.05f) return;
    const Color c = !rule.can ? Color{ 255, 80, 60, 255 } : digging ? Color{ 255, 220, 120, 255 } : Color{ 120, 220, 255, 255 };
    rlDrawRenderBatchActive();
    glDepthFunc(kGlGreater);
    rlDisableDepthMask();
    rlDisableBackfaceCulling(); // the bowl inside the ground is the sphere's far side
    DrawSphereEx(hit.point, radius, 12, 16, Fade(c, 0.35f));
    DrawSphereWires(hit.point, radius, 8, 12, Fade(c, 0.8f));
    rlDrawRenderBatchActive();
    glDepthFunc(kGlLequal);
    rlEnableDepthMask();
    rlEnableBackfaceCulling();
}

void TerrainTool_DrawHUD() {
    const int cx = GetScreenWidth() / 2, cy = GetScreenHeight() / 2;
    DrawLine(cx - 8, cy, cx + 8, cy, RAYWHITE);
    DrawLine(cx, cy - 8, cx, cy + 8, RAYWHITE);

    const int y = GetScreenHeight() - 60;
    const std::string held = tool.type == ToolType::Hands ? std::string("hands") : std::string(ToolTierName(tool.tier)) + " " + ToolTypeName(tool.type);
    if (hasPlant) {
        const Plant& p = plantHit.plant;
        DrawText(TextFormat("%s  %.1f m   LMB %s (%s)", p.tree() ? "tree" : "bush", plantHit.distance, p.tree() ? "chop" : "harvest", held.c_str()), cx + 14, cy + 10, 20, RAYWHITE);
        if (chopId == p.id && chopHits > 0) {
            const float f = (float)chopHits / HitsToFell(p);
            DrawRectangle(cx + 14, cy + 34, 120, 6, Fade(BLACK, 0.5f));
            DrawRectangle(cx + 14, cy + 34, (int)(120 * f), 6, ORANGE);
        }
    } else if (hasHit) {
        const char* name = hit.material == VoxelMaterial::Regolith && HasSoil(*hit.planetoid) ? "dirt" : VoxelMaterialName(hit.material);
        if (terrainTool.terraform)
            DrawText(TextFormat("%s  %.1f m   brush %.1f (debug terraform)", name, hit.distance, terrainTool.brushRadius), cx + 14, cy + 10, 20, RAYWHITE);
        else {
            DrawText(TextFormat("%s  %.1f m   %s", hit.material == VoxelMaterial::Ore ? hitItem.c_str() : name, hit.distance, held.c_str()), cx + 14, cy + 10, 20, RAYWHITE);
            if (!rule.can) DrawText(rule.why.c_str(), cx + 14, cy + 34, 18, Color{ 255, 120, 90, 255 });
        }
    }
    DrawText(terrainTool.terraform ? "LMB dig  RMB fill  wheel brush  G build  1-0 hotbar  Tab inventory  T telescope  L lamp  F fly  F5 save  F1 debug"
                                   : "LMB mine/chop  RMB place soil  G build  1-0 hotbar  Tab inventory  T telescope  L lamp  F fly  F5 save  F1 debug",
             10, y + 30, 20, Fade(RAYWHITE, 0.7f));
    if (GetTime() < toastUntil) DrawText(toast.c_str(), 10, y, 20, YELLOW);
}

bool TerrainTool_Chop() { tool = HeldTool(); return Chop(); }
bool TerrainTool_AimingPlant(uint64_t* id) {
    if (hasPlant && id) *id = plantHit.plant.id;
    return hasPlant;
}

// Scripted gathering (tests): dig at a body-local point with the held tool's rules (one item of
// whatever is there per call, if the tool can dig it), and chop / harvest the nearest plant.
bool TerrainTool_DigAt(Planetoid& p, Vector3 local) {
    if (!p.field->editable(local)) return false;
    const float base = p.field->baseDensity(local);
    if (base < ISO_LEVEL) return false; // air
    const VoxelMaterial m = p.field->material(local, base);
    const char* item = ItemAt(p, local, m);
    if (!terrainTool.terraform && !DigRuleFor(HeldTool(), m, item).can) return false;
    inventory.add(item, 1);
    if (strcmp(item, "dirt") == 0) inventory.addFractional("gravel", 0.2f); // stones in the soil
    if (WetBed(p, local, m)) inventory.addFractional("clay", 0.4f); // beds hold clay
    ApplyTerrainEdit(p, TerrainEdit{ local, 1.2f, -0.6f });
    return true;
}

bool TerrainTool_HarvestNearest(Planetoid& p, Vector3 local, bool tree) {
    if (!p.trees) return false;
    Plant plant;
    if (!p.trees->nearestPlant(local, 400.0f, tree, plant)) return false;
    tool = HeldTool();
    hasPlant = true;
    plantBody = &p;
    plantHit.plant = plant;
    plantHit.distance = 1.0f;
    for (int i = 0; i < 64; ++i) if (Chop()) return true;
    return false;
}
