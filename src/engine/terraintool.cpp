#include "terraintool.h"
#include "world.h"
#include "inventory.h"
#include "audio.h"
#include "particles.h"
#include "player.h"
#include "worldtypes.h"
#include "tools.h"
#include "plots.h"
#include "fluid.h"

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
    double nextStrokeAt = 0.0;   // the next swing may start then (GetTime)
    Vector3 aimLocal{};          // the view direction, body-local, at the last hit (digs go into the ground along it)
    float strokeLength = 0.5f;   // the current cooldown's length (the HUD bar)
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
        if (m == VoxelMaterial::Ore) return OreKindOf(p, local);
        // The moon's rock is anorthite (M24: smelted into aluminium on Terra).
        if (m == VoxelMaterial::Rock) if (const auto* pf = dynamic_cast<const PlanetField*>(p.field.get()); pf && pf->moon) return "anorthite";
        return ItemForMaterial(m, HasSoil(p));
    }

    // One swing at the aimed plant; true when it came down (items into the inventory).
    bool Chop() {
        if (!hasPlant || !plantBody || !plantBody->trees) return false;
        const Plant& p = plantHit.plant;
        if (p.id != chopId) { chopId = p.id; chopHits = 0; }
        Audio_Play(Sfx::Dig, 0.45f, p.tree() ? 0.6f : 1.4f);
        ++terrainTool.swings;
        {   // Chips off the trunk (leaves off a bush, grit off a stone).
            ParticleBurst chips;
            chips.kind = ParticleKind::Debris;
            chips.color = p.tree() ? Color{ 196, 160, 110, 255 } : p.stone() ? Color{ 130, 128, 124, 255 } : Color{ 70, 110, 45, 255 };
            chips.count = 6; chips.size = p.tree() ? 0.06f : 0.05f; chips.life = 1.4f; chips.speed = 2.5f; chips.spread = 0.8f;
            Particles_Emit(*plantBody, Vector3Add(p.base, Vector3Scale(p.up, p.tree() ? std::min(1.2f, p.height * 0.3f) : p.height * 0.5f)), chips);
        }
        if (++chopHits < HitsToFell(p)) return false;
        if (!p.stone()) {
            // It comes down: a burst of leaves from the crown.
            ParticleBurst leaves;
            leaves.kind = ParticleKind::Debris; leaves.color = { 60, 105, 40, 255 };
            leaves.count = p.tree() ? 40 : 14; leaves.size = 0.09f; leaves.life = 2.5f; leaves.speed = 2.0f; leaves.spread = 1.0f;
            Particles_Emit(*plantBody, Vector3Add(p.base, Vector3Scale(p.up, p.height * (p.tree() ? 0.8f : 0.5f))), leaves);
        }
        plantBody->trees->remove(p.id);
        // The biome it grew in (Plan 5 B resources): desert scrub is cactus, jungle trees give rubber.
        Biome biome = Biome::Forest;
        if (const auto* pf = dynamic_cast<const PlanetField*>(plantBody->field.get())) {
            const float r = Vector3Length(p.base);
            if (r > 1.0f) biome = pf->biomeAt(Vector3Scale(p.base, 1.0f / r), r - pf->radius);
        }
        if (p.tree() && biome == Biome::Jungle && itemDB.find("rubber")) { inventory.add("rubber", 1); }
        if (!p.tree() && !p.stone() && (biome == Biome::Desert || biome == Biome::Badlands) && itemDB.find("cactus")) {
            inventory.add("cactus", 2);
            Toast("+2 cactus");
            Audio_Play(Sfx::Pickup, 0.5f);
            chopId = 0; chopHits = 0; hasPlant = false;
            return true;
        }
        if (p.tree()) {
            // Logs of its species (saw them into planks: recipes take any planks).
            const WoodSpeciesDef& sp = WoodSpeciesInfo(SpeciesOfVariant(p.variant));
            const int logs = 1 + (int)(p.scale * 2.0f);
            inventory.add(sp.log, logs);
            Toast(TextFormat("+%d %s log%s%s", logs, sp.name, logs > 1 ? "s" : "", biome == Biome::Jungle ? "  +1 rubber" : ""));
            Audio_Play(Sfx::Build, 0.6f, 0.5f);
        } else if (p.stone()) {
            const int n = 1 + (int)((p.id >> 3) & 1); // a stone or two, picked up off the ground
            inventory.add("stone", n);
            Toast(TextFormat("+%d stone", n));
            Audio_Play(Sfx::Pickup, 0.5f, 0.8f);
        } else {
            inventory.add("fibre", 2);
            const bool berries = (p.id & 1) == 0;
            if (berries) inventory.add("berries", 1);
            // Seeds for farming: one of grain, roots or flax from half the bushes.
            static const char* seeds[3] = { "grain_seeds", "root_seeds", "flax_seeds" };
            const char* seed = ((p.id >> 1) & 1) == 0 ? seeds[(p.id >> 2) % 3] : nullptr;
            if (seed && itemDB.find(seed)) inventory.add(seed, 1);
            Toast(TextFormat("+2 fibre%s%s", berries ? "  +1 berries" : "", seed && itemDB.find(seed) ? "  +1 seeds" : ""));
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

// A built part nearer along the ray than the terrain hit blocks it: no digging, tilling or placing
// ground through walls and floors (breaking parts is the build tool's).
static bool PartInFront(Scene& world, Vector3 dir, const TerrainHit& h) {
    const Planetoid* body = h.planetoid;
    if (!body || !body->plots) return false;
    const Vector3 lo = body->toLocal(world.camera.position), ld = body->dirToLocal(dir);
    for (auto& pl : body->plots->plots) {
        Int3 c, f; float d;
        if (pl->raycast(lo, ld, h.distance / std::max(body->scale, 1e-6f), c, f, d)) return true;
    }
    return false;
}

// The bucket in hand (Plan 4, F5): empty, it scoops 1 m^3 from the first water along the view within
// reach (static or flowing); full, it pours where you look. Uses the aim of the last TerrainTool_Update.
bool TerrainTool_Bucket(Scene& world) {
    Planetoid* bodyP = player.nearest ? player.nearest : (hasHit ? hit.planetoid : nullptr); // nearest is set by the next update after a teleport
    if (!bodyP || (terrainTool.heldTool != "bucket" && terrainTool.heldTool != "water_bucket")) return false;
    const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
    {
    Planetoid& body = *bodyP;
    const std::string from = terrainTool.heldTool;
    const int kindHere = from == "bucket" && hasHit ? Fluid::KindAt(body, hit.local, 1.5f) : -1;
    if (kindHere > 0 && Fluid::FluidAt(body, hit.local, 1.5f, kindHere) >= 0.5f) {
        // Another fluid (fluids.txt): a bucketful yields the kind's bucket item; the bucket comes back empty.
        const Fluid::Kind& fk = Fluid::KindInfo(kindHere);
        if (!fk.bucketItem.empty() && itemDB.find(fk.bucketItem) && Fluid::Take(body, hit.local, 1.0f, 1.5f, kindHere) >= 0.5f) {
            inventory.add(fk.bucketItem, fk.bucketCount);
            Audio_Play(Sfx::Pickup, 0.6f, 0.5f);
            Toast(TextFormat("+%d %s", fk.bucketCount, fk.name.c_str()));
        } else Toast("Not enough to scoop here");
        return true;
    }
    if (from == "bucket") {
        const float maxD = hasHit ? hit.distance + 1.0f : terrainTool.reach; // the water may start just past a shallow bed
        bool found = false;
        for (float s = 0.4f; s <= maxD && !found; s += 0.25f) {
            const Vector3 q = Vector3Add(world.camera.position, Vector3Scale(dir, s));
            if (!IsUnderwater(world, q)) continue;
            found = true;
            if (Fluid::Take(body, body.toLocal(q), 1.0f, 1.5f) >= 0.5f && inventory.remove("bucket", 1)) {
                inventory.add("water_bucket", 1);
                terrainTool.heldTool = "water_bucket";
                Audio_Play(Sfx::Pickup, 0.6f, 0.6f);
                Toast("Bucket filled (1 m^3): RMB pours it");
                if (terrainTool.swapHeld) terrainTool.swapHeld(from, "water_bucket", "scoop");
            } else Toast("Not enough water here");
        }
        if (!found) Toast("No water within reach");
    } else if (!hasHit) {
        Toast("Look at the ground to pour");
    } else {
        const Vector3 up = Vector3Normalize(hit.local);
        Fluid::Pour(*hit.planetoid, Vector3Add(hit.local, Vector3Scale(up, 0.6f)), 1.0f);
        ParticleBurst splash;
        splash.kind = ParticleKind::Splash; splash.color = { 150, 190, 235, 220 };
        splash.count = 18; splash.size = 0.08f; splash.life = 0.9f; splash.speed = 2.5f; splash.spread = 0.7f;
        Particles_Emit(*hit.planetoid, Vector3Add(hit.local, Vector3Scale(up, 0.3f)), splash);
        if (inventory.remove("water_bucket", 1)) inventory.add("bucket", 1);
        terrainTool.heldTool = "bucket";
        Audio_Play(Sfx::Dig, 0.5f, 1.6f);
        if (terrainTool.swapHeld) terrainTool.swapHeld(from, "bucket", "pour");
    }
    }
    return true;
}

void TerrainTool_Update(Scene& world, bool inputEnabled, bool editEnabled) {
    PROFILE_FUNCTION();
    const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
    hasHit = RaycastTerrain(world, Ray{ world.camera.position, dir }, player.flyMode ? 150.0f : terrainTool.reach, hit) && !PartInFront(world, dir, hit);
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
        aimLocal = hit.planetoid->dirToLocal(dir);
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

    // Buckets (Plan 4, F5): RMB with an empty bucket scoops 1 m^3 of water (static or flowing) from
    // the first water along the view within reach; RMB with a full one pours it where you look.
    const bool bucket = terrainTool.heldTool == "bucket" || terrainTool.heldTool == "water_bucket";
    if (bucket && IsMouseButtonPressed(MOUSE_BUTTON_RIGHT)) TerrainTool_Bucket(world);

    digging = IsMouseButtonDown(MOUSE_BUTTON_LEFT);
    adding = !bucket && !digging && IsMouseButtonDown(MOUSE_BUTTON_RIGHT);
    // One swing per cooldown, and the cooldown runs whether the button is held or not, so clicking
    // fast is no faster than holding (it was: a release reset the timer). The swing animation lasts
    // exactly the cooldown (terrainTool.swingTime), and the HUD bar under the crosshair shows it.
    const double now = GetTime();
    if (!(digging || adding) || !(hasHit || (digging && hasPlant))) return;
    if (now < nextStrokeAt) return;
    const double base = now - nextStrokeAt < 0.1 ? nextStrokeAt : now; // held on: keep the beat (frame-rate independent)
    if (digging && hasPlant) {
        const float chop = 0.3f; // an axe (or hand) swing at a trunk or bush
        const bool felled = Chop();
        terrainTool.swingTime = chop;
        nextStrokeAt = base + (felled ? 0.45 : chop);
        strokeLength = (float)(nextStrokeAt - now);
        return;
    }
    const float rate = terrainTool.terraform ? terrainTool.strokesPerSecond : std::max(rule.strokesPerSecond, 2.0f);
    const float interval = 1.0f / rate;
    if (!TerrainTool_Stroke(digging)) { nextStrokeAt = now + 0.4; strokeLength = 0.4f; return; }
    terrainTool.swingTime = interval;
    nextStrokeAt = base + interval;
    strokeLength = (float)(nextStrokeAt - now);
}

// Sand, gravel or dirt on a river or pond bed (planets).
static bool WetBed(Planetoid& p, Vector3 local, VoxelMaterial m) {
    if (m != VoxelMaterial::Sand && m != VoxelMaterial::Gravel && m != VoxelMaterial::Regolith) return false;
    const auto* pf = dynamic_cast<const PlanetField*>(p.field.get());
    return pf && pf->wetBed(local);
}

// The ground a filler item makes when placed (-1: not a filler).
static int FillerMaterial(const std::string& item) {
    static const struct { const char* item; VoxelMaterial m; } kFillers[] = {
        { "dirt", VoxelMaterial::Regolith }, { "regolith", VoxelMaterial::Regolith }, { "gravel", VoxelMaterial::Gravel },
        { "sand", VoxelMaterial::Sand }, { "stone", VoxelMaterial::Rock }, { "concrete", VoxelMaterial::Rock },
        { "granite", VoxelMaterial::Rock }, { "basalt", VoxelMaterial::Rock }, { "limestone", VoxelMaterial::Rock },
        { "marble", VoxelMaterial::Rock }, { "sandstone", VoxelMaterial::Rock },
    };
    for (const auto& f : kFillers) if (item == f.item) return (int)f.m;
    return -1;
}

// Building consumes bulk material: the filler in hand when there is one (it is what gets placed),
// else the best one carried. `paint` = the material the placed ground is made of.
static bool PayFiller(int& paint) {
    const std::string& held = terrainTool.heldTool;
    if (FillerMaterial(held) >= 0 && (terrainTool.creative || inventory.remove(held, 1))) { paint = FillerMaterial(held); return true; }
    if (terrainTool.creative) { paint = -1; return true; }
    static const char* fillers[] = { "concrete", "stone", "dirt", "regolith", "gravel", "sand" };
    for (const char* f : fillers) if (inventory.remove(f, 1)) { paint = FillerMaterial(f); return true; }
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
        int paint = -1;
        if (!PayFiller(paint)) return false;
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

// Where a dig's sphere sits: half its radius into the ground along the view (or straight down), so
// the hole is where the crosshair is and the circle drawn on the surface is exactly its mouth
// (radius sqrt(1 - 0.25) = 0.87 of the sphere's). Centred on the surface point, half of it was air
// and how deep it bit depended on the angle.
static constexpr float kDigSink = 0.5f;
static Vector3 DigCentre(float r) {
    const Vector3 d = Vector3LengthSqr(aimLocal) > 0.5f ? Vector3Normalize(aimLocal) : Vector3Negate(Vector3Normalize(hit.local));
    return Vector3Add(hit.local, Vector3Scale(d, kDigSink * r));
}

// Particles for a stroke (particles.h): a puff of dust in the ground's colour and a few clods thrown
// out of the hole, toward the camera side.
static Color DustColor(VoxelMaterial m) {
    switch (m) {
        case VoxelMaterial::Regolith: return { 122, 92, 62, 255 };   // dirt on planets
        case VoxelMaterial::Ice: return { 200, 225, 245, 255 };
        case VoxelMaterial::Ore: return { 120, 96, 80, 255 };
        case VoxelMaterial::Sand: return { 214, 192, 146, 255 };
        case VoxelMaterial::Gravel: return { 128, 124, 118, 255 };
        case VoxelMaterial::Grass: return { 96, 84, 52, 255 };
        case VoxelMaterial::Snow: return { 240, 244, 250, 255 };
        default: return { 128, 126, 124, 255 };
    }
}
static void StrokeParticles(const Planetoid& p, Vector3 local, VoxelMaterial m, bool dig) {
    const Color c = DustColor(m);
    ParticleBurst dust;
    dust.kind = ParticleKind::Dust; dust.color = { c.r, c.g, c.b, 200 }; dust.count = dig ? 10 : 7;
    dust.size = 0.35f; dust.life = 1.1f; dust.speed = 1.2f; dust.spread = 0.9f;
    Particles_Emit(p, local, dust);
    if (!dig) return;
    ParticleBurst clods;
    clods.kind = ParticleKind::Debris; clods.color = c; clods.count = 7;
    clods.size = 0.07f; clods.life = 1.6f; clods.speed = 3.2f; clods.spread = 0.5f;
    Particles_Emit(p, local, clods);
}

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
    if (!dig && (tool.type == ToolType::Shovel || tool.type == ToolType::Hoe) && IsSoft(hit.material) && terrainTool.till) {
        // A shovel or a hoe tills soft ground into farmland (no soil placed); an iron or better hoe
        // tills the 3x3 metres around the aimed square at once.
        bool any = terrainTool.till(p, hit.local);
        if (tool.type == ToolType::Hoe && tool.tier >= 2) {
            const Vector3 up = Vector3Normalize(hit.local);
            Vector3 east = Vector3CrossProduct({ 0, 1, 0 }, up);
            east = Vector3Length(east) < 1e-3f ? Vector3{ 1, 0, 0 } : Vector3Normalize(east);
            const Vector3 north = Vector3CrossProduct(up, east);
            for (int a = -1; a <= 1; ++a)
                for (int b = -1; b <= 1; ++b) {
                    if (a == 0 && b == 0) continue;
                    const Vector3 d = Vector3Normalize(Vector3Add(hit.local, Vector3Add(Vector3Scale(east, (float)a), Vector3Scale(north, (float)b))));
                    const auto* pf = dynamic_cast<const PlanetField*>(p.field.get());
                    const Vector3 at = pf ? Vector3Scale(d, pf->radius + pf->elevation(d)) : Vector3Add(hit.local, Vector3Add(Vector3Scale(east, (float)a), Vector3Scale(north, (float)b)));
                    any = terrainTool.till(p, at) || any;
                }
        }
        if (any) { Audio_Play(Sfx::Dig, 0.5f, 0.7f); ++terrainTool.swings; return true; }
        return false;
    }
    if (!dig) {
        int paint = -1;
        if (!PayFiller(paint)) return false;
        TerrainEdit e{};
        e.radius = radius * 1.3f / p.scale;
        e.strength = 0.6f;
        e.paint = paint;
        const Vector3 localNormal = Vector3Normalize(Vector3Subtract(p.toLocal(Vector3Add(hit.point, hit.normal)), hit.local));
        e.center = Vector3Add(hit.local, Vector3Scale(localNormal, 0.3f * e.radius));
        StrokeParticles(p, hit.local, hit.material, false);
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
    const Vector3 centre = DigCentre(r);
    const float cell = step * step * step * p.scale * p.scale * p.scale;
    std::map<std::string, float> got;
    float wet = 0.0f;
    bool grass = false;
    const char* stoneKindAt = nullptr;
    bool swampHere = false; // Plan 5 B: digging swamp ground turns up peat
    if (const auto* pfb = dynamic_cast<const PlanetField*>(p.field.get()); pfb && itemDB.find("peat")) {
        const float r0 = Vector3Length(hit.local);
        if (r0 > 1.0f) swampHere = pfb->biomeAt(Vector3Scale(hit.local, 1.0f / r0), pfb->elevation(Vector3Scale(hit.local, 1.0f / r0))) == Biome::Swamp;
    }
    for (int z = -4; z <= 4; ++z)
        for (int y = -4; y <= 4; ++y)
            for (int x = -4; x <= 4; ++x) {
                if (x * x + y * y + z * z > 16) continue;
                const Vector3 q = Vector3Add(centre, { x * step, y * step, z * step });
                const float d = p.densityAt(q);
                if (d < ISO_LEVEL || !p.field->editable(q)) continue;
                const VoxelMaterial m = p.materialAt(q);
                if (softOnly && !IsSoft(m)) continue;
                const char* item = ItemAt(p, q, m);
                got[item] += cell * 2.0f;
                if (strcmp(item, "stone") == 0) { // the kind of stone the rock is here (building variety), one look per stroke
                    if (!stoneKindAt) { const auto* pfk = dynamic_cast<const PlanetField*>(p.field.get()); stoneKindAt = pfk ? pfk->stoneKind(hit.local) : "stone"; }
                    if (strcmp(stoneKindAt, "stone") != 0) got[stoneKindAt] += cell * 0.8f;
                }
                if (strcmp(item, "dirt") == 0) got["gravel"] += cell * 0.4f; // stones in the soil (enough for the first stone tools)
                if (WetBed(p, q, m)) wet += cell;
                if (m == VoxelMaterial::Grass) grass = true;
                if (swampHere && (m == VoxelMaterial::Grass || m == VoxelMaterial::Regolith)) got["peat"] += cell * 0.8f; // Plan 5 B: swamp ground is peat
            }
    for (auto& kv : got) inventory.addFractional(kv.first, kv.second);
    if (wet > 0.0f) inventory.addFractional("clay", wet * 0.8f); // river and pond beds hold clay
    if (stoneKindAt) { // rock now and then holds a gem: quartz anywhere, and the gem of its stone kind (M22)
        const char* gem = GetRandomValue(0, 29999) < 10 * terrainTool.gemBoost ? GemForStone(stoneKindAt) : GetRandomValue(0, 29999) < 16 * terrainTool.gemBoost ? "quartz" : nullptr; // rare: 1 in 3000 strokes a gem, 1 in 1900 quartz (veins later)
        if (gem && itemDB.find(gem)) { inventory.add(gem, 1); Toast(TextFormat("+1 %s!", gem)); }
    }
    if (grass && GetRandomValue(0, 99) < 6) { // grass now and then turns up a seed
        static const char* seeds[3] = { "grain_seeds", "root_seeds", "flax_seeds" };
        const char* s = seeds[GetRandomValue(0, 2)];
        if (itemDB.find(s)) { inventory.add(s, 1); Toast(TextFormat("+1 %s from the grass", s)); }
    }
    TerrainEdit e{};
    e.center = centre;
    e.radius = r;
    e.carve = true;
    e.softOnly = softOnly;
    StrokeParticles(p, hit.local, hit.material, true);
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
    // The dig's mouth: a circle on the surface where the sphere cuts it (DigCentre), depth-tested like
    // anything else, so it never shows through structures, foliage or hills.
    const float radius = MiningRadius() * sqrtf(1.0f - kDigSink * kDigSink);
    const Color c = !rule.can ? Color{ 255, 80, 60, 255 } : digging ? Color{ 255, 220, 120, 255 } : Color{ 120, 220, 255, 255 };
    const Vector3 n = Vector3Normalize(hit.normal);
    Vector3 u = Vector3CrossProduct(n, fabsf(n.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 });
    u = Vector3Normalize(u);
    const Vector3 v = Vector3CrossProduct(n, u);
    const Vector3 o = Vector3Add(hit.point, Vector3Scale(n, 0.05f)); // just off the ground (no z-fighting)
    constexpr int kSeg = 40;
    auto at = [&](int k) { const float a = 2.0f * PI * k / kSeg; return Vector3Add(o, Vector3Add(Vector3Scale(u, cosf(a) * radius), Vector3Scale(v, sinf(a) * radius))); };
    rlDrawRenderBatchActive();
    rlDisableDepthMask();
    rlDisableBackfaceCulling();
    for (int k = 0; k < kSeg; ++k) DrawTriangle3D(o, at(k), at(k + 1), Fade(c, 0.16f));
    for (int k = 0; k < kSeg; ++k) DrawLine3D(at(k), at(k + 1), Fade(c, 0.9f));
    rlDrawRenderBatchActive();
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
        DrawText(TextFormat("%s  %.1f m   LMB %s (%s)", p.tree() ? "tree" : p.stone() ? "stones" : "bush", plantHit.distance, p.tree() ? "chop" : p.stone() ? "pick up" : "harvest", held.c_str()), cx + 14, cy + 10, 20, RAYWHITE);
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
    if (strcmp(item, "stone") == 0) if (const auto* pfk = dynamic_cast<const PlanetField*>(p.field.get())) { // its kind of stone, and now and then a gem
        const char* k = pfk->stoneKind(local);
        if (strcmp(k, "stone") != 0) inventory.addFractional(k, 0.4f);
        const char* gem = GetRandomValue(0, 29999) < 10 * terrainTool.gemBoost ? GemForStone(k) : GetRandomValue(0, 29999) < 16 * terrainTool.gemBoost ? "quartz" : nullptr; // rare: 1 in 3000 strokes a gem, 1 in 1900 quartz (veins later)
        if (gem && itemDB.find(gem)) inventory.add(gem, 1);
    }
    if (WetBed(p, local, m)) inventory.addFractional("clay", 0.4f); // beds hold clay
    ApplyTerrainEdit(p, TerrainEdit{ local, 1.2f, -0.6f });
    return true;
}

bool TerrainTool_HarvestNearest(Planetoid& p, Vector3 local, bool tree, bool stone) {
    if (!p.trees) return false;
    Plant plant;
    if (!p.trees->nearestPlant(local, 400.0f, tree, plant, stone)) return false;
    tool = HeldTool();
    hasPlant = true;
    plantBody = &p;
    plantHit.plant = plant;
    plantHit.distance = 1.0f;
    for (int i = 0; i < 64; ++i) if (Chop()) return true;
    return false;
}

bool TerrainTool_IsFiller(const std::string& item) { return FillerMaterial(item) >= 0; }

float TerrainTool_Cooldown() {
    const double left = nextStrokeAt - GetTime();
    return left <= 0.0 ? 0.0f : std::min(1.0f, (float)(left / std::max(strokeLength, 0.01f)));
}
