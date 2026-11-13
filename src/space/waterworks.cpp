#include "waterworks.h"
#include "factory.h"
#include "fluid.h"
#include "plots.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <deque>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {
    double accum = 0.0;
    WaterworksStats stats;
    inline int64_t CellKey(Int3 c) { return ((int64_t)(c.x + 100000) << 40) | ((int64_t)(c.y + 100000) << 20) | (int64_t)(c.z + 100000); }
}

void Waterworks_Update(Scene& world, float dt) {
    accum += dt;
    if (accum < 0.5) return;
    accum = 0.0;
    stats.pumps = stats.networks = stats.outlets = stats.dry = 0;
    std::vector<std::pair<Int3, CellShape>> fp;
    for (auto& obj : world.objects) {
        Planetoid* body = dynamic_cast<Planetoid*>(obj.get());
        if (!body || !body->plots) continue;
        for (auto& pl : body->plots->plots) {
            // The plot's pumps, pipe cells and outlets.
            std::vector<int> pumps;
            std::unordered_map<int64_t, int> pipes, outlets, tanks; // cell -> part (tanks: machines that take water, boilers #229)
            for (int i = 0; i < (int)pl->parts.size(); ++i) {
                const Part& part = pl->parts[i];
                if (part.type == PartType::WaterPump) pumps.push_back(i);
                else if (part.type == PartType::Pipe || part.type == PartType::Outlet) {
                    Plot::Footprint(part, fp);
                    for (auto& c : fp) (part.type == PartType::Pipe ? pipes : outlets)[CellKey(c.first)] = i;
                } else if (Factory::TakesWater(part.type)) {
                    Plot::Footprint(part, fp);
                    for (auto& c : fp) tanks[CellKey(c.first)] = i;
                }
            }
            if (pumps.empty()) continue;
            for (int pump : pumps) {
                ++stats.pumps;
                // Its network: pipe cells joined face to face from the pump's own faces; outlets on it.
                Plot::Footprint(pl->parts[pump], fp);
                std::vector<Int3> pumpCells;
                Vector3 bottom = { 0, 0, 0 };
                int minY = INT32_MAX;
                for (auto& c : fp) { pumpCells.push_back(c.first); minY = std::min(minY, c.first.y); }
                int nb = 0;
                for (const Int3& c : pumpCells) if (c.y == minY) { bottom = Vector3Add(bottom, { (float)c.x, (float)c.y, (float)c.z }); ++nb; }
                bottom = Vector3Scale(bottom, 1.0f / std::max(1, nb));
                std::unordered_set<int64_t> seen, outs;
                std::unordered_set<int> tankParts; // boilers on the network: they fill from it
                std::deque<Int3> q;
                for (const Int3& c : pumpCells) { seen.insert(CellKey(c)); q.push_back(c); }
                const int d[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
                while (!q.empty() && seen.size() < 20000) {
                    const Int3 c = q.front(); q.pop_front();
                    for (auto& o : d) {
                        const Int3 n = { c.x + o[0], c.y + o[1], c.z + o[2] };
                        const int64_t k = CellKey(n);
                        if (seen.count(k)) continue;
                        if (outlets.count(k)) { outs.insert(k); seen.insert(k); continue; }
                        if (tanks.count(k)) { tankParts.insert(tanks[k]); seen.insert(k); continue; }
                        if (pipes.count(k)) { seen.insert(k); q.push_back(n); }
                    }
                }
                if (outs.empty() && tankParts.empty()) continue; // nowhere to send it: the pump fills up and stops
                ++stats.networks;
                // Each outlet part once, however many of its cells touch the network.
                std::unordered_set<int> outParts;
                for (int64_t k : outs) outParts.insert(outlets[k]);
                stats.outlets += (int)outParts.size();
                const int units = Factory::PumpWater(body->id, pl->id, pump); // less what earlier drains still have queued
                if (units <= 0) continue;
                // From under the pump; only the units the world actually gave leave the pump (a dry or short
                // intake keeps the rest).
                const Vector3 intake = Vector3Subtract(pl->toBody({ (bottom.x + 0.5f) * kCellSize, bottom.y * kCellSize, (bottom.z + 0.5f) * kCellSize }), Vector3Scale(pl->up, 0.7f));
                const float got = Fluid::Take(*body, intake, units * kWaterPerUnit, 1.6f);
                if (got <= 0.0f) { ++stats.dry; continue; }
                Factory::Drain(body->id, pl->id, pump, "water", std::min(units, (int)ceilf(got / kWaterPerUnit - 1e-4f)));
                stats.pumped += got;
                const float each = got / (float)(outParts.size() + tankParts.size());
                for (int part : tankParts) Factory::AddWater(body->id, pl->id, part, (int)(each * 1000.0f)); // m^3 to litres
                for (int part : outParts) {
                    Plot::Footprint(pl->parts[part], fp);
                    const Int3 c = fp.empty() ? pl->parts[part].origin : fp[0].first;
                    Fluid::Pour(*body, pl->toBody({ (c.x + 0.5f) * kCellSize, (c.y + 0.5f) * kCellSize, (c.z + 0.5f) * kCellSize }), each);
                }
            }
        }
    }
}

WaterworksStats Waterworks_Stats() { return stats; }
