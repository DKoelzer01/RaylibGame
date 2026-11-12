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

void Ambience_Update(Scene& world, float dt) {
    Planetoid* body = player.nearest;
    if (!body || dt <= 0.0f) return;
    const auto* pf = dynamic_cast<const PlanetField*>(body->field.get());
    if (!pf) return;
    const Vector3 local = body->toLocalD(player.position).toFloat();
    const float r = Vector3Length(local);
    if (r < 1.0f) return;
    const Vector3 up = Vector3Scale(local, 1.0f / r);
    const float surface = pf->radius + pf->elevation(up);
    const bool outdoors = r > surface - 3.0f && r < surface + 60.0f;
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
