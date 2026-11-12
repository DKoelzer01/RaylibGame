#include "ambience.h"
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

// Oil seeps (Plan 5 B): in warm biomes (desert, badlands, savanna, jungle, swamp), about a third of the
// 800 m cells on the sea-level sphere hold one, at a hashed point. Near the player a seep wells up crude
// oil through the flowing-fluid sim (fluid.h: oil cells don't mix with water): the first time this session
// its basin is cut (a shallow bowl, saved with the terrain edits) and seeded with a few m^3, then oil
// seeps in slowly while the pool holds less than kSeepCap. Scoop it with a bucket for crude oil.
namespace {
    constexpr float kSeepCell = 800.0f, kSeepCap = 8.0f, kSeepRate = 0.08f; // m, m^3, m^3/s
    std::unordered_set<uint64_t> seededSeeps;
    double seepClock = 0.0;
    uint32_t SeepHash(int x, int y, int z, uint32_t seed) {
        uint32_t h = (uint32_t)x * 73856093u ^ (uint32_t)y * 19349663u ^ (uint32_t)z * 83492791u ^ seed * 2654435761u;
        h ^= h >> 15; h *= 2246822519u; h ^= h >> 13; h *= 3266489917u; h ^= h >> 16;
        return h;
    }
    bool WarmBiome(Biome b) { return b == Biome::Desert || b == Biome::Badlands || b == Biome::Savanna || b == Biome::Jungle || b == Biome::Swamp; }
}

int Ambience_SeepsNear(const Planetoid& body, Vector3 local, float range, std::vector<Vector3>& out) {
    out.clear();
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    if (!pf || pf->moon) return 0;
    const float r = Vector3Length(local);
    if (r < 1.0f) return 0;
    const Vector3 s = Vector3Scale(local, pf->radius / r); // on the sea-level sphere
    const int cx = (int)floorf(s.x / kSeepCell), cy = (int)floorf(s.y / kSeepCell), cz = (int)floorf(s.z / kSeepCell);
    const int R = (int)ceilf(range / kSeepCell) + 1;
    for (int dz = -R; dz <= R; ++dz)
        for (int dy = -R; dy <= R; ++dy)
            for (int dx = -R; dx <= R; ++dx) {
                uint32_t h = SeepHash(cx + dx, cy + dy, cz + dz, body.seed ^ 0x0115eedu);
                if ((h & 1023) > 340) continue; // ~a third of the cells
                const Vector3 c = { (cx + dx + ((h >> 10) & 255) / 255.0f) * kSeepCell, (cy + dy + ((h >> 18) & 255) / 255.0f) * kSeepCell, (cz + dz + ((h >> 26) & 63) / 63.0f) * kSeepCell };
                const float cl = Vector3Length(c);
                if (cl < 1.0f) continue;
                const Vector3 dir = Vector3Scale(c, 1.0f / cl);
                // The point must be in this cell's own patch of the sphere (each cell owns its seep).
                const Vector3 onSphere = Vector3Scale(dir, pf->radius);
                if ((int)floorf(onSphere.x / kSeepCell) != cx + dx || (int)floorf(onSphere.y / kSeepCell) != cy + dy || (int)floorf(onSphere.z / kSeepCell) != cz + dz) continue;
                const float e = pf->elevation(dir);
                if (e < 3.0f || e > 600.0f || !WarmBiome(pf->biomeAt(dir, e))) continue;
                const Vector3 at = Vector3Scale(dir, pf->radius + e);
                if (Vector3Distance(at, local) <= range) out.push_back(at);
            }
    return (int)out.size();
}

static void UpdateSeeps(Planetoid& body, Vector3 local, float dt) {
    seepClock += dt;
    if (seepClock < 0.5) return;
    const float step = (float)seepClock;
    seepClock = 0.0;
    static std::vector<Vector3> seeps;
    if (!Ambience_SeepsNear(body, local, 160.0f, seeps)) return;
    for (const Vector3& at : seeps) {
        const Vector3 up = Vector3Normalize(at);
        const uint64_t key = ((uint64_t)(int64_t)lroundf(at.x) * 73856093ull) ^ ((uint64_t)(int64_t)lroundf(at.y) * 19349663ull) ^ ((uint64_t)(int64_t)lroundf(at.z) * 83492791ull); // by place (a load makes new body ids)
        const Vector3 pool = Vector3Add(at, Vector3Scale(up, 0.3f));
        if (!seededSeeps.count(key)) {
            seededSeeps.insert(key);
            // The basin: a bowl ~5 m across and ~1 m deep (oil fills it from the next update).
            TerrainEdit e{};
            e.center = Vector3Add(at, Vector3Scale(up, 2.6f)); e.radius = 3.6f; e.carve = true; e.softOnly = false; // the sphere's bottom 1 m into the ground
            ApplyTerrainEdit(body, e);
            // A source with a rate and a head just under the rim: the bowl fills and the seep rests (Fluid::AddSource).
            const int oil = Fluid::KindByName("crude_oil");
            if (oil > 0) Fluid::AddSource(body, Vector3Subtract(at, Vector3Scale(up, 0.8f)), oil, kSeepRate, 0.65f);

            continue;
        }
        if (GetRandomValue(0, 2) == 0) { // a bubble of gas breaking the surface
            ParticleBurst b;
            b.kind = ParticleKind::Splash; b.color = { 40, 32, 26, 230 };
            b.count = 3; b.size = 0.07f; b.life = 0.5f; b.speed = 0.8f; b.spread = 0.4f;
            Particles_Emit(body, Vector3Add(at, Vector3Scale(up, 0.1f)), b);
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
            const Vector3 base = plot->toBody({ (part.origin.x + 1) * kCellSize, part.origin.y * kCellSize, (part.origin.z + 1) * kCellSize });
            if (Vector3Distance(base, local) > 120.0f) continue;
            if (Factory::View(body->id, plot->id, (int)i).status != MachineStatus::Working) continue;
            const bool campfire = part.type == PartType::Campfire;
            ParticleBurst smoke;
            smoke.kind = ParticleKind::Smoke; smoke.color = campfire ? Color{ 190, 188, 185, 255 } : Color{ 150, 148, 145, 255 };
            smoke.count = std::max(1, (int)(tick * 8.0f)); smoke.size = campfire ? 0.35f : 0.5f; smoke.life = 4.5f; smoke.speed = 0.9f; smoke.spread = 0.25f;
            smoke.dir = plot->up;
            const Vector3 top = Vector3Add(base, Vector3Scale(plot->up, campfire ? 0.5f : 1.05f));
            Particles_Emit(*body, top, smoke);
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
