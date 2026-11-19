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
#include "menus.h"
#include <unordered_map>
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
            if (part.type == PartType::Torch) { // a placed torch's flame (#558: particles, it was two glowing boxes)
                Vector3 at;
                if (Machine_SmokePoints(*plot, i, &at, 1) == 0 || Vector3Distance(at, local) > 60.0f) continue;
                ParticleBurst fire;
                fire.kind = ParticleKind::Fire; fire.color = { 255, 140, 40, 255 };
                fire.count = std::max(1, (int)(tick * 16.0f)); fire.size = 0.11f; fire.life = 0.45f; fire.speed = 0.45f; fire.spread = 0.12f;
                fire.dir = plot->up;
                Particles_Emit(*body, at, fire);
                if (GetRandomValue(0, 5) == 0) {
                    ParticleBurst wisp;
                    wisp.kind = ParticleKind::Smoke; wisp.color = { 120, 118, 115, 255 };
                    wisp.count = 1; wisp.size = 0.08f; wisp.life = 1.6f; wisp.speed = 0.5f; wisp.spread = 0.15f;
                    wisp.dir = plot->up;
                    Particles_Emit(*body, Vector3Add(at, Vector3Scale(plot->up, 0.2f)), wisp);
                }
                continue;
            }
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

// Local lights (#139) near the camera: lamp blocks (steady warm white), lit campfires (flickering orange) and working
// fire machines (a dim glow at the door; the arc furnace blue-white). Each plot's light spots are cached per version.
void Ambience_GatherLights(Scene& world) {
    world.localLights.clear();
    Planetoid* body = player.nearest;
    if (!body || !body->plots) return;
    const Vector3 cam = body->toLocal(world.camera.position);
    static int lampMat = -2;
    if (lampMat == -2) lampMat = BuildMaterialForItem("lamp");
    struct Spot { int part; Vector3 centre; uint8_t kind; }; // kind: 0 lamp, 1 campfire, 2 fire machine, 3 arc furnace, 4 placed torch
    struct Cached { uint32_t version = 0xFFFFFFFFu; std::vector<Spot> spots; };
    static std::unordered_map<uint64_t, Cached> cache;
    static std::vector<std::pair<Int3, CellShape>> fp;
    const float t = (float)fmod(GameClock(), 3600.0); // flicker stops behind the pause menu (#511)
    for (auto& plot : body->plots->plots) {
        if (Vector3Distance(plot->toBody({ 0, 0, 0 }), cam) > 300.0f) continue;
        Cached& c = cache[(body->id * 1315423911ull) ^ plot->id];
        if (c.version != plot->version) {
            c.version = plot->version;
            c.spots.clear();
            for (size_t i = 0; i < plot->parts.size(); ++i) {
                const Part& p = plot->parts[i];
                int kind = -1;
                if (p.type == PartType::Campfire) kind = 1;
                else if (p.type == PartType::Torch) kind = 4;
                else if (p.type == PartType::ArcFurnace) kind = 3;
                else if (p.type == PartType::Kiln || p.type == PartType::CokeOven || p.type == PartType::BlastFurnace || p.type == PartType::Furnace ||
                         p.type == PartType::AlloySmelter) kind = 2;
                else if (lampMat >= 0 && IsShape(p.type) && p.material == lampMat) kind = 0;
                if (kind < 0) continue;
                Plot::Footprint(p, fp);
                if (fp.empty()) continue;
                Vector3 sum{};
                for (auto& cell : fp) sum = Vector3Add(sum, { (cell.first.x + 0.5f) * kCellSize, (cell.first.y + 0.5f) * kCellSize, (cell.first.z + 0.5f) * kCellSize });
                c.spots.push_back({ (int)i, Vector3Scale(sum, 1.0f / fp.size()), (uint8_t)kind });
            }
        }
        for (const Spot& s : c.spots) {
            const Vector3 at = plot->toBody(s.centre);
            if (Vector3Distance(at, cam) > 110.0f) continue;
            LocalLight L;
            L.position = body->toWorld(Vector3Add(at, Vector3Scale(plot->up, s.kind == 1 ? 0.3f : 0.0f)));
            if (s.kind == 0) { L.range = 8.0f; L.color = { 1.0f, 0.92f, 0.76f }; L.intensity = 1.1f; }
            else if (s.kind == 4) { // a placed torch (#218): always lit, flickering, from its flame (top of the post)
                const float flicker = settings.fireFlicker ? 0.85f + 0.1f * sinf(t * 13.0f + s.part) + 0.05f * sinf(t * 31.0f + s.part) : 0.85f; // #558: steady when off
                L.position = body->toWorld(Vector3Add(at, Vector3Scale(plot->up, 0.4f)));
                L.range = 10.0f; L.color = { 1.0f, 0.62f, 0.3f }; L.intensity = 1.3f * flicker;
            }
            else {
                if (!Factory::Working(body->id, plot->id, s.part)) continue;
                const float flicker = settings.fireFlicker ? 0.82f + 0.12f * sinf(t * 11.0f + s.part) + 0.06f * sinf(t * 27.3f + 2.0f * s.part) : 0.82f;
                if (s.kind == 1) { L.range = 11.0f; L.color = { 1.0f, 0.55f, 0.22f }; L.intensity = 1.7f * flicker; }
                else if (s.kind == 2) { L.range = 5.0f; L.color = { 1.0f, 0.5f, 0.2f }; L.intensity = 0.9f * flicker; }
                else { L.range = 6.0f; L.color = { 0.75f, 0.85f, 1.0f }; L.intensity = 1.2f * flicker; }
            }
            world.localLights.push_back(L);
        }
    }
    // (#972, user 2026-11-19) Room bounce light: sunlit floor patches under a roof (the sun in through a door or a window)
    // light the room a little. Every half second, the plots within 40 m: each building-block cell with air over it, a roof
    // over that (a solid cell within 6 m straight up) and a clear line to the sun through the plot (12 m) is a sunlit
    // indoor floor cell; they are binned into 2 m squares, each a faint, wide point light 0.4 m over the patch in the sun's
    // colour times the floor's, strength by its sunlit area, the sun's height and the cloud cover. RAYLIBGAME_ROOM_BOUNCE=0:
    // off. The lights are cached between updates (they follow the sun slowly).
    static const bool bounceOn = !(getenv("RAYLIBGAME_ROOM_BOUNCE") && getenv("RAYLIBGAME_ROOM_BOUNCE")[0] == '0');
    static std::vector<LocalLight> bounce; // positions BODY-LOCAL: the body spins between updates (render space drifted off)
    static double nextBounce = 0.0;
    if (bounceOn && GetTime() >= nextBounce) {
        nextBounce = GetTime() + 0.5;
        const double bt0 = GetTime(); int bounceCells = 0;
        bounce.clear();
        const Vector3 sunL = body->dirToLocal(Vector3Normalize(world.sunDirection));
        const float cloud = std::clamp(world.weatherCloud, 0.0f, 1.0f);
        for (auto& plot : body->plots->plots) {
            if (plot->cells.empty() || Vector3Distance(plot->toBody({ 0, 0, 0 }), cam) > 40.0f + 0.5f * kCellSize * 64.0f) continue;
            const Vector3 s = { Vector3DotProduct(sunL, plot->east), Vector3DotProduct(sunL, plot->up), Vector3DotProduct(sunL, plot->north) };
            if (s.y < 0.05f) continue;                                      // the sun is down (for this plot)
            const float sunK = std::min(1.0f, s.y * 3.0f) * (1.0f - 0.8f * cloud);
            auto solid = [&](int x, int y, int z) {
                auto it = plot->cells.find(Int3{ x, y, z });
                return it != plot->cells.end() && PartBlocksSky(plot->parts[it->second.part]);
            };
            struct Bin { int n = 0; Vector3 sum{}, col{}; };
            std::unordered_map<long long, Bin> bins;
            for (const auto& [c, ref] : plot->cells) {
                if (!PartBlocksSky(plot->parts[ref.part]) || solid(c.x, c.y + 1, c.z)) continue; // a floor: solid, air over it
                const Vector3 top = { (c.x + 0.5f) * kCellSize, (c.y + 1.0f) * kCellSize, (c.z + 0.5f) * kCellSize };
                if (Vector3Distance(plot->toBody(top), cam) > 40.0f) continue;
                bool roofed = false;
                for (int k = 2; k <= 12 && !roofed; ++k) roofed = solid(c.x, c.y + k, c.z);
                if (!roofed) continue;                                       // outdoors: the sky's ambient covers it
                bool lit = true;
                for (float t = 0.6f; t <= 24.0f && lit; t += 0.5f) {         // cells: 12 m toward the sun
                    const Vector3 q = { c.x + 0.5f + s.x * t, c.y + 1.05f + s.y * t, c.z + 0.5f + s.z * t };
                    lit = !solid((int)floorf(q.x), (int)floorf(q.y), (int)floorf(q.z));
                }
                if (!lit) continue;
                Bin& b = bins[((long long)(c.x >> 2) * 73856093LL) ^ ((long long)(c.y) * 19349663LL) ^ ((long long)(c.z >> 2) * 83492791LL)];
                const Color mc = GetBuildMaterial(plot->parts[ref.part].material).color;
                ++bounceCells; b.n++; b.sum = Vector3Add(b.sum, top); b.col = Vector3Add(b.col, { mc.r / 255.0f, mc.g / 255.0f, mc.b / 255.0f });
            }
            for (const auto& [k, b] : bins) {
                if (bounce.size() >= 24) break;
                LocalLight L;
                const Vector3 at = Vector3Add(plot->toBody(Vector3Scale(b.sum, 1.0f / b.n)), Vector3Scale(plot->up, 0.4f));
                L.position = at; // body-local (to render space each frame, below)
                static const float brange = getenv("RAYLIBGAME_ROOM_BOUNCE_RANGE") ? (float)atof(getenv("RAYLIBGAME_ROOM_BOUNCE_RANGE")) : 4.0f; L.range = brange;
                const Vector3 alb = Vector3Scale(b.col, 1.0f / b.n);
                L.color = { 1.0f * alb.x, 0.96f * alb.y, 0.88f * alb.z };   // the sun's warm white off the floor
                static const float gain = getenv("RAYLIBGAME_ROOM_BOUNCE_GAIN") ? (float)atof(getenv("RAYLIBGAME_ROOM_BOUNCE_GAIN")) : 1.0f; // tuning
                L.intensity = 0.012f * std::min(b.n, 16) * sunK * gain;              // a full 2 m patch: ~0.19 at noon (small and faded: the user, #972)
                bounce.push_back(L);
            }
        }
        static const bool bounceLog = getenv("RAYLIBGAME_ROOM_BOUNCE_LOG") && getenv("RAYLIBGAME_ROOM_BOUNCE_LOG")[0] == '1';
        if (bounceLog) printf("room bounce: %zu lights from %d sunlit indoor floor cells in %.3f ms (cloud %.2f)\n", bounce.size(), bounceCells, (GetTime() - bt0) * 1000.0, cloud);
    }
    if (bounceOn) for (LocalLight L : bounce) { L.position = body->toWorld(L.position); world.localLights.push_back(L); }
}
