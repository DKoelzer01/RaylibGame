#include "terraintool.h"
#include "gamerng.h"
#include "world.h"
#include "playerinput.h"
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
#include <unordered_map>
#include <string>


// (terrainTool: the local PlayerSession's, src/space/playersession.cpp, #313)

namespace {
    TerrainHit hit;
    bool hasHit = false;
    bool digging = false, adding = false;
    double nextStrokeAt = 0.0;   // the next swing may start then (SimTime, #315)
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
        return std::max(1, (int)lroundf((4.0f + 9.0f * p.trunkVolume()) * ChopFactor(tool))); // the trunk's volume (#274, user 2026-11-15): a giant oak ~190 stone-axe hits (~1 min), a medium tree 10-25, a sapling ~5; better axes proportionally faster
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
            const int logs = std::clamp((int)lroundf(3.0f * p.chopSize()), 1, 48); // the trunk (#274): 3-5 from a medium tree, ~26 from a giant
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

static int FillerMaterial(const std::string& item); // below
void TerrainTool_Update(Scene& world, const PlayerInput& in, bool editEnabled) {
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

    if (!in.enabled) { digging = adding = false; return; }

    if (in.hit(Act::Save)) Toast(SaveWorld(world) ? TextFormat("Saved to %s", WorldSavePath()) : "Save failed");
    if (in.hit(Act::Load)) {
        if (LoadWorldSave(world)) { Toast("Loaded saves/world.sav"); hasHit = false; return; }
        Toast("No save found");
    }
    if (!editEnabled) { digging = adding = false; return; } // e.g. piloting: save/load only

    if (terrainTool.terraform) {
        const float wheel = in.wheel;
        const bool modWheel = in.down(Act::Down) || in.down(Act::FreeLook); // Ctrl or Alt (#227: the plain wheel scrolls the hotbar)
        if (modWheel && wheel != 0.0f) terrainTool.brushRadius = std::clamp(terrainTool.brushRadius * (wheel > 0 ? 1.15f : 1.0f / 1.15f), 1.0f, 24.0f);
    }

    // Fill shapes (#228, user 2026-11-15): with a filler in hand, Ctrl+wheel steps through the shapes RMB places, as it
    // does through build shapes: the soft blob, then exact spheres of 0.5-3 m (paid by the volume they fill).
    if (!terrainTool.terraform && FillerMaterial(terrainTool.heldTool) >= 0) {
        static const float kFillShapes[] = { 0.0f, 0.5f, 1.0f, 1.5f, 2.0f, 3.0f }; // 0: the blob
        constexpr int kN = (int)(sizeof(kFillShapes) / sizeof(kFillShapes[0]));
        const bool modWheel = in.down(Act::Down) || in.down(Act::FreeLook);
        if (modWheel && in.wheel != 0.0f) {
            int k = 0;
            if (terrainTool.sphereFill) for (int i = 1; i < kN; ++i) if (fabsf(kFillShapes[i] - terrainTool.fillRadius) < 0.01f) k = i;
            k = (k + (in.wheel > 0 ? 1 : kN - 1)) % kN;
            terrainTool.sphereFill = k > 0;
            if (k > 0) terrainTool.fillRadius = kFillShapes[k];
            Toast(k > 0 ? TextFormat("Fill shape: sphere %.1f m", kFillShapes[k]) : "Fill shape: soft blob");
        }
    }
    // Buckets (Plan 4, F5): RMB with an empty bucket scoops 1 m^3 of water (static or flowing) from
    // the first water along the view within reach; RMB with a full one pours it where you look.
    const bool bucket = terrainTool.heldTool == "bucket" || terrainTool.heldTool == "water_bucket";
    if (bucket && in.hit(Act::Secondary)) TerrainTool_Bucket(world);

    digging = in.down(Act::Primary);
    adding = !bucket && !terrainTool.holdingItem && !digging && in.down(Act::Secondary); // food / feed in hand never builds (#201)
    // One swing per cooldown, and the cooldown runs whether the button is held or not, so clicking
    // fast is no faster than holding (it was: a release reset the timer). The swing animation lasts
    // exactly the cooldown (terrainTool.swingTime), and the HUD bar under the crosshair shows it.
    const double now = SimTime();
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

// The ground a filler item makes when placed (-1: not a filler).
static int FillerMaterial(const std::string& item) {
    static const struct { const char* item; VoxelMaterial m; } kFillers[] = {
        { "dirt", VoxelMaterial::Regolith }, { "regolith", VoxelMaterial::Regolith }, { "gravel", VoxelMaterial::Gravel },
        { "sand", VoxelMaterial::Sand }, { "clay", VoxelMaterial::Clay }, { "stone", VoxelMaterial::Rock }, { "concrete", VoxelMaterial::Rock },
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

// Ore veins (#207): once a deposit is ~92% mined, the rest comes out by itself (into the inventory, no tool wear). A deposit
// is sampled once (1.5 m grid over its blob: the points that were ore, and their original density); later checks only
// add the edits, so a check is ~1 ms. Noise pockets (small) are left as they are.
namespace {
    struct VeinSamples { std::vector<Vector3> pts; std::vector<float> base; bool done = false; };
    std::unordered_map<uint64_t, VeinSamples> veinCache;
}
static void CheckOreVein(Planetoid& p, Vector3 local) {
    const auto* pf = dynamic_cast<const PlanetField*>(p.field.get());
    Vector3 c; float rad; int kind;
    if (!pf || !pf->depositAt(local, c, rad, kind)) return;
    const uint64_t key = p.id * 0x9e3779b97f4a7c15ull ^ ((uint64_t)(int64_t)floorf(c.x) * 73856093u) ^ ((uint64_t)(int64_t)floorf(c.y) * 19349663u) ^ ((uint64_t)(int64_t)floorf(c.z) * 83492791u);
    VeinSamples& v = veinCache[key];
    if (v.done) return;
    constexpr float kStep = 1.5f;
    if (v.pts.empty()) {
        const float reach = rad * 1.4f;
        for (float z = -reach; z <= reach; z += kStep)
            for (float y = -reach; y <= reach; y += kStep)
                for (float x = -reach; x <= reach; x += kStep) {
                    if (x * x + y * y + z * z > reach * reach) continue;
                    const Vector3 q = Vector3Add(c, { x, y, z });
                    const float d = Vector3Length(q);
                    if (d < 1.0f) continue;
                    const float depth = pf->radius + pf->elevation(Vector3Scale(q, 1.0f / d)) - d;
                    if (pf->oreKind(q, depth) != kind) continue;
                    const float base = p.field->density(q, nullptr, 0.0f);
                    if (base < ISO_LEVEL) continue; // a cave through it: never ore to mine
                    v.pts.push_back(q);
                    v.base.push_back(base);
                }
        if (v.pts.empty()) { v.done = true; return; }
    }
    std::vector<Vector3> left;
    for (size_t i = 0; i < v.pts.size(); ++i)
        if (v.base[i] + (p.edits ? p.edits->sample(v.pts[i]) : 0.0f) >= ISO_LEVEL) left.push_back(v.pts[i]);
    if (v.pts.size() < 12 || (float)left.size() > 0.08f * (float)v.pts.size()) return;
    v.done = true;
    if (left.empty()) return;
    std::vector<TerrainEdit> edits;
    for (const Vector3& q : left) { TerrainEdit e{}; e.center = q; e.radius = 1.3f; e.carve = true; edits.push_back(e); } // 1.3 m: the whole 1.5 m cell (its diagonal), no ore left between
    ApplyTerrainEdits(p, edits.data(), edits.size());
    const char* item = PlanetField::OreKindName(kind);
    const float amount = (float)left.size() * kStep * kStep * kStep * 2.0f * p.scale * p.scale * p.scale; // 2 items a m^3, as digging
    if (item) {
        inventory.addFractional(item, amount);
        Toast(TextFormat("The rest of the vein comes loose: %.0f %s", amount, item));
    }
    Audio_Play(Sfx::Dig, 0.7f, 0.7f);
}

bool TerrainTool_TestVein(Planetoid& p, Vector3 home) {
    // #207: the first deposit within ~400 m of `home` (a body-local direction) with enough ore: mine 93% of its ore samples
    // directly, run the check: the rest must come loose and land in the inventory.
    const auto* pf = dynamic_cast<const PlanetField*>(p.field.get());
    if (!pf) return false;
    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
    for (int a = -4; a <= 4; ++a)
        for (int b = -4; b <= 4; ++b) {
            const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, a * 96.0f / pf->radius), Vector3Scale(t2, b * 96.0f / pf->radius))));
            for (float depth = 6.0f; depth < 90.0f; depth += 12.0f) {
                const Vector3 q = Vector3Scale(d, pf->radius + pf->elevation(d) - depth);
                Vector3 c; float rad; int kind;
                if (!pf->depositAt(q, c, rad, kind)) continue;
                const float cd = Vector3Length(c);
                const float cdepth = pf->radius + pf->elevation(Vector3Scale(c, 1.0f / cd)) - cd;
                if (pf->oreKind(c, cdepth) != kind) continue; // its middle is ore (deep enough for the kind)
                veinCache.clear();
                CheckOreVein(p, c); // samples it (nothing mined: nothing happens)
                const uint64_t key = p.id * 0x9e3779b97f4a7c15ull ^ ((uint64_t)(int64_t)floorf(c.x) * 73856093u) ^ ((uint64_t)(int64_t)floorf(c.y) * 19349663u) ^ ((uint64_t)(int64_t)floorf(c.z) * 83492791u);
                VeinSamples& v = veinCache[key];
                if (v.pts.size() < 40) continue;
                std::vector<TerrainEdit> mine;
                for (size_t i = 0; i < v.pts.size() * 93 / 100; ++i) { TerrainEdit e{}; e.center = v.pts[i]; e.radius = 0.8f; e.carve = true; mine.push_back(e); }
                ApplyTerrainEdits(p, mine.data(), mine.size());
                const char* item = PlanetField::OreKindName(kind);
                const double before = inventory.amount(item);
                CheckOreVein(p, c);
                const double got = inventory.amount(item) - before;
                size_t left = 0;
                for (size_t i = 0; i < v.pts.size(); ++i) left += v.base[i] + (p.edits ? p.edits->sample(v.pts[i]) : 0.0f) >= ISO_LEVEL;
                const bool ok = v.done && got > 0.0 && left == 0;
                printf("bench: vein (#207): a %s deposit (r %.1f m, %zu ore samples), 93%% mined by hand: the rest came loose %s (+%.1f %s, %zu left) -> %s\n",
                       item, rad, v.pts.size(), v.done ? "yes" : "NO", got, item, left, ok ? "PASS" : "FAIL");
                return ok;
            }
        }
    printf("bench: vein (#207): no deposit near the start -> FAIL\n");
    return false;
}

// The old brush (debug): any material, soft-edged, yield by the brush volume.
static bool BrushStroke(bool digging) {
    Planetoid& p = *hit.planetoid;
    const float r3 = terrainTool.brushRadius * terrainTool.brushRadius * terrainTool.brushRadius;
    if (digging) {
        inventory.addFractional(hitItem, terrainTool.digStrength * r3 * 0.15f);
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
        case VoxelMaterial::Clay: return { 170, 120, 95, 255 };
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

// Sphere fill (#228): the sphere sits on the aimed surface (its centre 0.4 radius out along the normal, so it joins
// the ground), and costs what it adds: the lattice points inside it that are air now, 2 items a cubic metre (as digging
// yields).
static constexpr float kFillOut = 0.4f;
static Vector3 FillCentre(const Planetoid& p, float r) {
    const Vector3 localNormal = Vector3Normalize(Vector3Subtract(p.toLocal(Vector3Add(hit.point, hit.normal)), hit.local));
    return Vector3Add(hit.local, Vector3Scale(localNormal, kFillOut * r));
}
static float FillCost(const Planetoid& p, Vector3 c, float r, int stride = 1) { // stride 2: every other lattice point (the HUD estimate, 8x cheaper)
    int air = 0;
    const int x0 = (int)floorf((c.x - r) / kVoxelSize), x1 = (int)ceilf((c.x + r) / kVoxelSize);
    const int y0 = (int)floorf((c.y - r) / kVoxelSize), y1 = (int)ceilf((c.y + r) / kVoxelSize);
    const int z0 = (int)floorf((c.z - r) / kVoxelSize), z1 = (int)ceilf((c.z + r) / kVoxelSize);
    for (int z = z0; z <= z1; z += stride)
        for (int y = y0; y <= y1; y += stride)
            for (int x = x0; x <= x1; x += stride) {
                const Vector3 q = { x * kVoxelSize, y * kVoxelSize, z * kVoxelSize };
                if (Vector3DistanceSqr(q, c) >= r * r) continue;
                if (p.densityAt(q) < ISO_LEVEL) ++air;
            }
    const float s = p.scale;
    return air * (float)(stride * stride * stride) * kVoxelSize * kVoxelSize * kVoxelSize * s * s * s * 2.0f;
}
static bool FillStroke(Planetoid& p) {
    const float r = terrainTool.fillRadius / p.scale;
    const Vector3 c = FillCentre(p, r);
    // Pay with the filler in hand (what gets placed), by the volume added.
    const std::string held = terrainTool.heldTool;
    const int paint = FillerMaterial(held);
    if (paint < 0) return false;
    const float cost = FillCost(p, c, r);
    if (cost < 0.05f) { Toast("Already solid there"); return false; }
    if (!terrainTool.creative && !inventory.removeAmount(held, cost)) {
        Toast(TextFormat("Need %.1f %s for this fill (have %.1f)", cost, held.c_str(), inventory.amount(held)));
        Audio_Play(Sfx::Error, 0.5f);
        return false;
    }
    TerrainEdit e{};
    e.center = c; e.radius = r; e.fill = true; e.paint = paint;
    StrokeParticles(p, hit.local, hit.material, false);
    ApplyTerrainEdit(p, e);
    Audio_Play(Sfx::Build, 0.6f, 0.8f);
    return true;
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
    if (!dig && terrainTool.sphereFill && FillerMaterial(terrainTool.heldTool) >= 0) return FillStroke(p); // #228
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
                if (m == VoxelMaterial::Grass) grass = true;
                if (swampHere && (m == VoxelMaterial::Grass || m == VoxelMaterial::Regolith)) got["peat"] += cell * 0.8f; // Plan 5 B: swamp ground is peat
            }
    for (auto& kv : got) inventory.addFractional(kv.first, kv.second);
    if (stoneKindAt) { // rock now and then holds a gem: quartz anywhere, and the gem of its stone kind (M22)
        const char* gem = GameRng_Int(0, 29999) < 10 * terrainTool.gemBoost ? GemForStone(stoneKindAt) : GameRng_Int(0, 29999) < 16 * terrainTool.gemBoost ? "quartz" : nullptr; // rare: 1 in 3000 strokes a gem, 1 in 1900 quartz (veins later)
        if (gem && itemDB.find(gem)) { inventory.add(gem, 1); Toast(TextFormat("+1 %s!", gem)); }
    }
    if (grass && GameRng_Int(0, 99) < 6) { // grass now and then turns up a seed
        static const char* seeds[3] = { "grain_seeds", "root_seeds", "flax_seeds" };
        const char* s = seeds[GameRng_Int(0, 2)];
        if (itemDB.find(s)) { inventory.add(s, 1); Toast(TextFormat("+1 %s from the grass", s)); }
    }
    TerrainEdit e{};
    e.center = centre;
    e.radius = r;
    e.carve = true;
    e.softOnly = softOnly;
    StrokeParticles(p, hit.local, hit.material, true);
    ApplyTerrainEdit(p, e);
    if (!softOnly) for (int k = 0; k < PlanetField::kOreKindCount; ++k) // a stroke that took ore: the vein may be nearly done (#207)
        if (const char* oreName = PlanetField::OreKindName(k); oreName && got.count(oreName)) { CheckOreVein(p, hit.local); break; }
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
    if (terrainTool.sphereFill && FillerMaterial(terrainTool.heldTool) >= 0 && !digging) { // #228: the sphere it will add
        Planetoid& p = *hit.planetoid;
        const float r = terrainTool.fillRadius / p.scale;
        const Vector3 c = p.toWorld(FillCentre(p, r));
        if (Vector3Distance(c, hit.point) + 0.3f > hit.distance && hit.distance < terrainTool.fillRadius * 1.2f) return; // round the camera
        DrawSphereWires(c, terrainTool.fillRadius, 8, 14, Fade(adding ? Color{ 120, 255, 160, 255 } : Color{ 120, 220, 255, 255 }, 0.8f));
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
            if (FillerMaterial(terrainTool.heldTool) >= 0) { // #228: the fill mode, and what a sphere here costs
                if (terrainTool.sphereFill) {
                    Planetoid& p = *hit.planetoid;
                    const float r = terrainTool.fillRadius / p.scale;
                    static double costAt = -1.0; static float cost = 0.0f;
                    if (GetTime() - costAt > 0.2) { costAt = GetTime(); cost = FillCost(p, FillCentre(p, r), r, 2); } // an estimate (the stroke counts exactly)
                    DrawText(TextFormat("sphere fill %.2f m: %.1f %s (have %.1f)   Ctrl+wheel: shape", terrainTool.fillRadius, cost, terrainTool.heldTool.c_str(), inventory.amount(terrainTool.heldTool)),
                             cx + 14, cy + 58, 18, Color{ 140, 230, 255, 255 });
                } else DrawText("fill shape: soft blob   Ctrl+wheel: shape", cx + 14, cy + 58, 18, Fade(RAYWHITE, 0.6f));
            }
        }
    }
    DrawText(terrainTool.terraform ? "LMB dig  RMB fill  Ctrl+wheel brush  G build  1-0 hotbar  Tab inventory  T telescope  L lamp  F fly  F5 save  F1 debug"
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
        const char* gem = GameRng_Int(0, 29999) < 10 * terrainTool.gemBoost ? GemForStone(k) : GameRng_Int(0, 29999) < 16 * terrainTool.gemBoost ? "quartz" : nullptr; // rare: 1 in 3000 strokes a gem, 1 in 1900 quartz (veins later)
        if (gem && itemDB.find(gem)) inventory.add(gem, 1);
    }
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
    const double left = nextStrokeAt - SimTime();
    return left <= 0.0 ? 0.0f : std::min(1.0f, (float)(left / std::max(strokeLength, 0.01f)));
}

bool TerrainTool_LastHit(Planetoid*& body, Vector3& local, Vector3& normalWorld) {
    if (!hasHit || !hit.planetoid) return false;
    body = hit.planetoid; local = hit.local; normalWorld = hit.normal;
    return true;
}

bool TerrainTool_TestFill() {
    // --test-fill (#228): a 1 m sphere of dirt on the aimed ground: it costs exactly the air it fills (2 a m^3), leaves no
    // air inside, and has a sharp edge (air 0.35 m outside it, solid 0.2 m inside).
    if (!hasHit || !hit.planetoid) { printf("bench: fill: nothing aimed at -> FAIL\n"); return false; }
    Planetoid& p = *hit.planetoid;
    terrainTool.heldTool = "dirt";
    inventory.add("dirt", 30);
    terrainTool.sphereFill = true;
    terrainTool.fillRadius = 1.0f;
    const float r = 1.0f / p.scale;
    const Vector3 c = FillCentre(p, r);
    const Vector3 n = Vector3Normalize(Vector3Subtract(c, hit.local));
    const float cost = FillCost(p, c, r);
    const double before = inventory.amount("dirt");
    const bool stroke = TerrainTool_Stroke(false);
    const double paid = before - inventory.amount("dirt");
    const float left = FillCost(p, c, r);
    const bool edge = p.densityAt(Vector3Add(c, Vector3Scale(n, r + 0.35f))) < ISO_LEVEL;
    const bool inside = p.densityAt(Vector3Add(c, Vector3Scale(n, r - 0.2f))) >= ISO_LEVEL;
    const bool ok = stroke && fabs(paid - cost) < 0.02 && cost > 0.5f && cost < 8.5f && left < 0.05f && edge && inside;
    printf("bench: fill: a 1 m sphere of dirt cost %.2f (paid %.2f; a whole sphere would be 8.38), air left inside %.2f, sharp edge %d, solid inside %d -> %s\n",
           cost, paid, left, (int)edge, (int)inside, ok ? "PASS" : "FAIL");
    terrainTool.sphereFill = false;
    return ok;
}
