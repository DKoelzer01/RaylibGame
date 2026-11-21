#include "trapping.h"
#include "audio.h"
#include "gameui.h"
#include "inventory.h"
#include "player.h"
#include "plots.h"
#include "scene.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <cmath>
#include <cstring>
#include <map>
#include <string>
#include <tuple>

namespace {
    struct Snare { double setAt = 0.0; bool caught = false; };
    std::map<std::tuple<std::string, uint32_t, int>, Snare> snares; // (body name, plot id, part) -> state
    double lastCheck = -10.0;
    bool forceCatch = false;                                         // tests
    constexpr double kInterval = 1200.0;                             // world seconds between chances
    constexpr float kCollect = 2.5f;

    uint32_t Hash(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdull; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull; x ^= x >> 33; return (uint32_t)x; }

    float GameOdds(Biome b) {
        switch (b) {
            case Biome::Meadow: case Biome::Forest: case Biome::Savanna: case Biome::Taiga: return 0.7f;
            case Biome::Jungle: case Biome::Tundra: case Biome::Alpine: return 0.5f;
            case Biome::Swamp: case Biome::Badlands: return 0.35f;
            case Biome::Desert: return 0.25f;
            default: return 0.1f; // polar, ocean shelves
        }
    }
}

void Trapping_Reset() { snares.clear(); lastCheck = -10.0; }

bool Trapping_HasCatch(uint64_t bodyId, uint32_t plotId, int part) {
    Planetoid* b = Planetoid::find(bodyId);
    if (!b) return false;
    auto it = snares.find(std::make_tuple(b->name, plotId, part));
    return it != snares.end() && it->second.caught;
}

void Trapping_Update(Scene& world) {
    if (GetTime() - lastCheck < 1.0 && !forceCatch) return;
    lastCheck = GetTime();
    Planetoid* body = player.nearest;
    if (!body || !body->plots) return;
    const auto* pf = dynamic_cast<const PlanetField*>(body->field.get());
    const Vector3 me = body->toLocal((player.position - world.renderOrigin).toFloat());
    const double now = WorldTime();
    for (auto& plot : body->plots->plots) {
        if (Vector3Distance(plot->toBody({ 0, 0, 0 }), me) > 300.0f) continue;
        for (size_t i = 0; i < plot->parts.size(); ++i) {
            if (plot->parts[i].type != PartType::Snare) continue;
            const auto key = std::make_tuple(body->name, plot->id, (int)i);
            auto it = snares.find(key);
            if (it == snares.end()) { snares[key] = Snare{ now, false }; continue; } // just set
            Snare& s = it->second;
            const Int3 o = plot->parts[i].origin;
            const Vector3 at = plot->toBody({ (o.x + 0.5f) * kCellSize, (o.y + 0.25f) * kCellSize, (o.z + 0.5f) * kCellSize });
            if (!s.caught) {
                const uint32_t h = Hash(((uint64_t)plot->id << 32) ^ (uint64_t)i ^ (uint64_t)(int64_t)(s.setAt * 7.0));
                const double due = kInterval * (0.7 + 0.6 * ((h & 0xFFFF) / 65535.0));
                if (now - s.setAt >= due || forceCatch) {
                    const float d = Vector3Length(at);
                    const Biome b = pf && d > 1.0f ? pf->biomeAt(Vector3Scale(at, 1.0f / d), d - pf->radius) : Biome::Meadow;
                    if (forceCatch || ((h >> 16) & 0xFFFF) / 65535.0f < GameOdds(b)) s.caught = true;
                    else s.setAt = now; // nothing this time: another chance later
                }
            }
            if (s.caught && Vector3Distance(at, me) < kCollect + 1.0f) { // standing at it (feet ~1 m above the frame)
                inventory.add("raw_meat", 1);
                const bool hide = (Hash((uint64_t)(now * 3.0) ^ i) & 1) != 0;
                if (hide) inventory.add("hide", 1);
                GameUI_Toast(hide ? "The snare held small game: +1 raw meat, +1 hide" : "The snare held small game: +1 raw meat");
                Audio_Play(Sfx::Pickup, 0.6f, 0.9f);
                s.caught = false; s.setAt = now;
            }
        }
    }
    // Snares no longer standing are forgotten.
    for (auto it = snares.begin(); it != snares.end();) {
        bool alive = false;
        if (std::get<0>(it->first) != body->name) { ++it; continue; }
        for (auto& plot : body->plots->plots)
            if (plot->id == std::get<1>(it->first)) { const int p = std::get<2>(it->first); alive = p >= 0 && p < (int)plot->parts.size() && plot->parts[p].type == PartType::Snare; break; }
        if (!alive) it = snares.erase(it); else ++it;
    }
    forceCatch = false;
}

void Trapping_Save(FILE* f) {
    for (const auto& kv : snares)
        fprintf(f, "snare %s %u %d %.1f %d\n", std::get<0>(kv.first).c_str(), std::get<1>(kv.first), std::get<2>(kv.first), kv.second.setAt, kv.second.caught ? 1 : 0);
}

bool Trapping_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "snare") != 0) return false;
    char name[128]; unsigned plot = 0; int part = 0, caught = 0; double at = 0.0;
    if (fscanf(f, " %127s %u %d %lf %d", name, &plot, &part, &at, &caught) == 5) snares[std::make_tuple(std::string(name), (uint32_t)plot, part)] = Snare{ at, caught != 0 };
    return true;
}

int Trapping_TestAge(Scene& world, double seconds) {
    lastCheck = -10.0;
    Trapping_Update(world); // registers new snares
    int n = 0;
    for (auto& kv : snares) { kv.second.setAt -= seconds; ++n; }
    forceCatch = true;
    return n;
}
