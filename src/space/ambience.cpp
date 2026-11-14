#include "ambience.h"
#include "machinemodels.h"
#include "weather.h"
#include "particles.h"
#include "factory.h"
#include "player.h"
#include "plots.h"
#include "scene.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"
#include "fluid.h"
#include "trees.h"
#include <unordered_set>
#include <vector>
#include <raymath.h>
#include <algorithm>
#include <cmath>

// Environmental particles around the player (particles.h): splashes where rain lands, fireflies over
// grassland on clear nights, and smoke (with flame and embers for campfires) from fire-driven machines
// while they work.

namespace {
    float rainDue = 0.0f, flyDue = 0.0f, smokeClock = 0.0f;

    float Rand(float a, float b) { return a + (b - a) * (float)GetRandomValue(0, 10000) / 10000.0f; }

    // A random point on the ground within `range` m of the player (body-local); false if none.
    bool GroundPoint(const Planetoid& body, const PlanetField& pf, Vector3 up, float range, Vector3& out, bool& water) {
        Vector3 e = Vector3CrossProduct(up, { 0, 0, 1 });
        if (Vector3Length(e) < 1e-3f) e = Vector3CrossProduct(up, { 1, 0, 0 });
        e = Vector3Normalize(e);
        const Vector3 n = Vector3CrossProduct(up, e);
        const float a = Rand(0.0f, 6.2831853f), d = range * sqrtf(Rand(0.02f, 1.0f));
        const Vector3 p = Vector3Add(Vector3Scale(up, pf.radius), Vector3Add(Vector3Scale(e, cosf(a) * d), Vector3Scale(n, sinf(a) * d)));
        const Vector3 dir = Vector3Normalize(p);
        const float ground = pf.radius + pf.elevation(dir), sea = pf.seaLevelRadius();
        water = sea > 0.0f && ground < sea;
        out = Vector3Scale(dir, std::max(ground, water ? sea : ground));
        (void)body;
        return true;
    }

    bool SmokyMachine(PartType t) {
        switch (t) {
            case PartType::Campfire: case PartType::Kiln: case PartType::CokeOven: case PartType::BlastFurnace:
            case PartType::Boiler: case PartType::SteamMacerator: case PartType::ForgeHammer: case PartType::AlloySmelter:
                return true;
            default: return false;
        }
    }
}

// Natural fluid sources (Plan 5), all through Fluid::AddSource (a rate and a head):
// - oil seeps in warm biomes (desert, badlands, savanna, jungle, swamp): a shallow basin is cut (saved
//   with the terrain edits) and the seep's head lies just under its rim, so the pool fills and rests;
// - water springs on hillsides in wet temperate and cold biomes (forest, meadow, taiga, jungle, alpine):
//   no basin, the head above the ground, so the water overflows and runs downhill as a stream.
// Each type has a grid of cells on the sea-level sphere; a hashed share of the cells hold one, at a hashed
// point, where the biome and height suit it. Sources are registered when the player first comes within
// 160 m in a session (the fluid save keeps them and their fluid; re-registering replaces them).
namespace {
    struct SourceType { const char* fluid; float cell, share, rate, minElev, maxElev; uint32_t salt; bool basin; };
    const SourceType kSourceTypes[] = {
        { "crude_oil", 800.0f, 0.33f, 0.08f, 3.0f, 600.0f, 0x0115eedu, true },   // seeps: pools
        { "water", 900.0f, 0.22f, 0.06f, 30.0f, 800.0f, 0x5291e5u, false },     // springs: streams
    };
    constexpr int kSourceTypeCount = (int)(sizeof(kSourceTypes) / sizeof(kSourceTypes[0]));
    std::unordered_set<uint64_t> seededSeeps;
    double seepClock = 0.0;
    uint32_t SeepHash(int x, int y, int z, uint32_t seed) {
        uint32_t h = (uint32_t)x * 73856093u ^ (uint32_t)y * 19349663u ^ (uint32_t)z * 83492791u ^ seed * 2654435761u;
        h ^= h >> 15; h *= 2246822519u; h ^= h >> 13; h *= 3266489917u; h ^= h >> 16;
        return h;
    }
    bool Suits(int type, Biome b) {
        if (type == 0) return b == Biome::Desert || b == Biome::Badlands || b == Biome::Savanna || b == Biome::Jungle || b == Biome::Swamp;
        return b == Biome::Forest || b == Biome::Meadow || b == Biome::Taiga || b == Biome::Jungle || b == Biome::Alpine;
    }
    int SourcesNear(const Planetoid& body, Vector3 local, float range, int type, std::vector<Vector3>& out) {
        out.clear();
        const SourceType& st = kSourceTypes[type];
        const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
        if (!pf || pf->moon) return 0;
        const float r = Vector3Length(local);
        if (r < 1.0f) return 0;
        const Vector3 s = Vector3Scale(local, pf->radius / r); // on the sea-level sphere
        const int cx = (int)floorf(s.x / st.cell), cy = (int)floorf(s.y / st.cell), cz = (int)floorf(s.z / st.cell);
        const int R = (int)ceilf(range / st.cell) + 1;
        for (int dz = -R; dz <= R; ++dz)
            for (int dy = -R; dy <= R; ++dy)
                for (int dx = -R; dx <= R; ++dx) {
                    const uint32_t h = SeepHash(cx + dx, cy + dy, cz + dz, body.seed ^ st.salt);
                    if ((h & 1023) > (uint32_t)(st.share * 1023.0f)) continue;
                    const Vector3 c = { (cx + dx + ((h >> 10) & 255) / 255.0f) * st.cell, (cy + dy + ((h >> 18) & 255) / 255.0f) * st.cell, (cz + dz + ((h >> 26) & 63) / 63.0f) * st.cell };
                    const float cl = Vector3Length(c);
                    if (cl < 1.0f) continue;
                    const Vector3 dir = Vector3Scale(c, 1.0f / cl);
                    // The point must be in this cell's own patch of the sphere (each cell owns its source).
                    const Vector3 onSphere = Vector3Scale(dir, pf->radius);
                    if ((int)floorf(onSphere.x / st.cell) != cx + dx || (int)floorf(onSphere.y / st.cell) != cy + dy || (int)floorf(onSphere.z / st.cell) != cz + dz) continue;
                    const float e = pf->elevation(dir);
                    if (e < st.minElev || e > st.maxElev || !Suits(type, pf->biomeAt(dir, e))) continue;
                    const Vector3 at = Vector3Scale(dir, pf->radius + e);
                    if (Vector3Distance(at, local) <= range) out.push_back(at);
                }
        return (int)out.size();
    }
}

void Ambience_Reset() { seededSeeps.clear(); seepClock = 0.0; } // a new world or a load: sources re-check their basin, clearing and registration

int Ambience_SeepsNear(const Planetoid& body, Vector3 local, float range, std::vector<Vector3>& out) { return SourcesNear(body, local, range, 0, out); }
int Ambience_SpringsNear(const Planetoid& body, Vector3 local, float range, std::vector<Vector3>& out) { return SourcesNear(body, local, range, 1, out); }

static void UpdateSeeps(Planetoid& body, Vector3 local, float dt) {
    seepClock += dt;
    if (seepClock < 0.5) return;
    seepClock = 0.0;
    static std::vector<Vector3> found;
    for (int type = 0; type < kSourceTypeCount; ++type) {
        const SourceType& st = kSourceTypes[type];
        if (!SourcesNear(body, local, 160.0f, type, found)) continue;
        for (const Vector3& at : found) {
            const Vector3 up = Vector3Normalize(at);
            const uint64_t key = (((uint64_t)(int64_t)lroundf(at.x) * 73856093ull) ^ ((uint64_t)(int64_t)lroundf(at.y) * 19349663ull) ^ ((uint64_t)(int64_t)lroundf(at.z) * 83492791ull)) + (uint64_t)type; // by place (a load makes new body ids)
            if (!seededSeeps.count(key)) {
                seededSeeps.insert(key);
                const int kind = Fluid::KindByName(st.fluid);
                if (kind < 0) continue;
                if (st.basin) {
                    // A bowl ~5 m across and ~1 m deep: the sphere's bottom 1 m into the ground.
                    TerrainEdit e{};
                    e.center = Vector3Add(at, Vector3Scale(up, 2.6f)); e.radius = 3.6f; e.carve = true; e.softOnly = false;
                    ApplyTerrainEdit(body, e);
                }
                if (body.trees) body.trees->addClearing(at, st.basin ? 5.0f : 2.0f); // no plants floating over it
                // On the real ground (the mesh strays from the height map): the first solid going down.
                float floorH = -1.0f;
                for (float h = 8.0f; h > -5.0f; h -= 0.1f) if (body.densityAt(Vector3Add(at, Vector3Scale(up, h))) >= ISO_LEVEL) { floorH = h; break; }
                const Vector3 src = Vector3Add(at, Vector3Scale(up, floorH + 0.3f));
                // Basins: the head just under the rim (fills and rests). Springs: 1.2 m over the ground (overflows).
                const float head = st.basin ? std::max(0.3f, -0.15f - (floorH + 0.3f)) : 1.2f;
                Fluid::AddSource(body, src, kind, st.rate, head);
                continue;
            }
            if (GetRandomValue(0, 2) == 0) { // bubbles breaking the surface: gas in a seep, the welling of a spring
                ParticleBurst b;
                b.kind = ParticleKind::Splash; b.color = st.basin ? Color{ 40, 32, 26, 230 } : Color{ 190, 215, 240, 200 };
                b.count = 3; b.size = 0.07f; b.life = 0.5f; b.speed = 0.8f; b.spread = 0.4f;
                Particles_Emit(body, Vector3Add(at, Vector3Scale(up, 0.1f)), b);
            }
        }
    }
}

void Ambience_Update(Scene& world, float dt) {
    // The body underfoot (when the player isn't colliding with one yet, e.g. flying or just spawned: the nearest).
    Planetoid* body = player.nearest ? player.nearest : NearestPlanetoid(world, (player.position - world.renderOrigin).toFloat());
    if (!body || dt <= 0.0f) return;
    const auto* pf = dynamic_cast<const PlanetField*>(body->field.get());
    if (!pf) return;
    const Vector3 local = body->toLocalD(player.position).toFloat();
    const float r = Vector3Length(local);
    if (r < 1.0f) return;
    const Vector3 up = Vector3Scale(local, 1.0f / r);
    const float surface = pf->radius + pf->elevation(up);
    const bool outdoors = r > surface - 3.0f && r < surface + 60.0f;
    UpdateSeeps(*body, local, dt); // oil seeps in warm biomes
    const WeatherState& w = Weather_Here();

    // Rain splashes: ~150 a second in a downpour, on the ground (or the sea) within 14 m.
    if (outdoors && !pf->moon && (w.kind == WeatherKind::Rain || w.kind == WeatherKind::Storm) && w.precipitation > 0.1f) {
        rainDue = std::min(rainDue + dt * 150.0f * w.precipitation, 30.0f);
        while (rainDue >= 1.0f) {
            rainDue -= 1.0f;
            Vector3 at; bool water = false;
            if (!GroundPoint(*body, *pf, up, 14.0f, at, water)) continue;
            ParticleBurst s;
            s.kind = ParticleKind::Splash; s.color = water ? Color{ 200, 220, 245, 210 } : Color{ 175, 195, 225, 190 };
            s.count = 3; s.size = 0.045f; s.life = 0.35f; s.speed = 1.3f; s.spread = 0.45f;
            Particles_Emit(*body, at, s);
        }
    } else rainDue = 0.0f;

    // Fireflies: clear, mild nights over temperate grassland; ~40 around the player.
    const float sunUp = Vector3DotProduct(body->dirToWorld(up), Vector3Normalize(world.sunDirection));
    const float lat = fabsf(up.y);
    if (outdoors && !pf->moon && sunUp < -0.08f && w.precipitation < 0.05f && lat > 0.12f && lat < 0.6f && surface - pf->radius > 4.0f && surface - pf->radius < 700.0f) {
        flyDue = std::min(flyDue + dt * 7.0f, 10.0f);
        while (flyDue >= 1.0f) {
            flyDue -= 1.0f;
            Vector3 at; bool water = false;
            if (!GroundPoint(*body, *pf, up, 22.0f, at, water) || water) continue;
            at = Vector3Add(at, Vector3Scale(Vector3Normalize(at), Rand(0.4f, 2.4f)));
            ParticleBurst f;
            f.kind = ParticleKind::Firefly; f.color = { 200, 255, 110, 255 };
            f.count = 1; f.size = 0.1f; f.life = 6.0f; f.speed = 0.3f; f.spread = 1.0f;
            Particles_Emit(*body, at, f);
        }
    } else flyDue = 0.0f;

    // Smoke from working fire machines within 120 m (checked ten times a second).
    smokeClock += dt;
    if (smokeClock < 0.1f || !body->plots) return;
    const float tick = smokeClock;
    smokeClock = 0.0f;
    for (auto& plot : body->plots->plots) {
        if (Vector3Distance(plot->toBody({ 0, 0, 0 }), local) > 160.0f) continue;
        for (size_t i = 0; i < plot->parts.size(); ++i) {
            const Part& part = plot->parts[i];
            if (!SmokyMachine(part.type)) continue;
            // From the chimney tops of the model (#220: the base assumed the footprint ran +x/+z from the origin cell,
            // so rotated machines smoked from the opposite corner, in the air).
            Vector3 tops[4];
            const int nTops = Machine_SmokePoints(*plot, i, tops, 4);
            if (nTops == 0 || Vector3Distance(tops[0], local) > 120.0f) continue;
            if (!Factory::Working(body->id, plot->id, (int)i)) continue;
            const bool campfire = part.type == PartType::Campfire;
            ParticleBurst smoke;
            smoke.kind = ParticleKind::Smoke; smoke.color = campfire ? Color{ 190, 188, 185, 255 } : Color{ 150, 148, 145, 255 };
            smoke.count = std::max(1, (int)(tick * 8.0f)); smoke.size = campfire ? 0.35f : 0.5f; smoke.life = 4.5f; smoke.speed = 0.9f; smoke.spread = 0.25f;
            smoke.dir = plot->up;
            for (int k = 0; k < nTops; ++k) Particles_Emit(*body, tops[k], smoke);
            const Vector3 base = Vector3Subtract(tops[0], Vector3Scale(plot->up, 0.5f)); // the campfire's floor (its smoke point is 0.5 m up)
            if (campfire) {
                ParticleBurst fire;
                fire.kind = ParticleKind::Fire; fire.color = { 255, 120, 30, 255 };
                fire.count = std::max(2, (int)(tick * 30.0f)); fire.size = 0.35f; fire.life = 0.6f; fire.speed = 0.8f; fire.spread = 0.2f;
                fire.dir = plot->up;
                Particles_Emit(*body, Vector3Add(base, Vector3Scale(plot->up, 0.15f)), fire);
                if (GetRandomValue(0, 2) == 0) {
                    ParticleBurst ember;
                    ember.kind = ParticleKind::Spark; ember.color = { 255, 170, 60, 255 };
                    ember.count = 1; ember.size = 0.04f; ember.life = 1.6f; ember.speed = 2.0f; ember.spread = 0.3f;
                    ember.dir = plot->up;
                    Particles_Emit(*body, Vector3Add(base, Vector3Scale(plot->up, 0.3f)), ember);
                }
            }
        }
    }
}
