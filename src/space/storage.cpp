#include "storage.h"
#include "factory.h"
#include "hotbar.h"
#include "inventory.h"
#include "player.h"
#include "plots.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <algorithm>

std::vector<NearChest> Storage_NearbyChests(Scene& world, float radius) {
    std::vector<NearChest> out;
    Planetoid* body = player.nearest;
    if (!body || !body->plots) return out;
    const Vector3 me = body->toLocal((player.position - world.renderOrigin).toFloat());
    for (const auto& pl : body->plots->plots)
        for (int i = 0; i < (int)pl->parts.size(); ++i) {
            const Part& p = pl->parts[i];
            if (p.type != PartType::Chest) continue;
            // A chest is 2x2x2 cells: its middle is one cell in from its origin.
            const Vector3 c = pl->toBody({ (p.origin.x + 1.0f) * kCellSize, (p.origin.y + 1.0f) * kCellSize, (p.origin.z + 1.0f) * kCellSize });
            const float d = Vector3Distance(c, me);
            if (d <= radius) out.push_back({ body->id, pl->id, i, d });
        }
    std::sort(out.begin(), out.end(), [](const NearChest& a, const NearChest& b) { return a.dist < b.dist; });
    return out;
}

int Storage_CountNearby(Scene& world, const std::string& item) {
    int n = 0;
    for (const NearChest& c : Storage_NearbyChests(world))
        for (const auto& s : Factory::ChestContents(c.body, c.plot, c.part)) if (s.first == item) n += s.second;
    return n;
}

bool Storage_Pull(Scene& world, const std::string& item, int count) {
    if (count <= 0) return true;
    int left = count;
    for (const NearChest& c : Storage_NearbyChests(world)) {
        for (const auto& s : Factory::ChestContents(c.body, c.plot, c.part)) {
            if (s.first != item || left <= 0) continue;
            const int k = std::min(left, s.second);
            Factory::Take(c.body, c.plot, c.part, item, k);
            left -= k;
        }
        if (left <= 0) break;
    }
    return left <= 0;
}

int Storage_Deposit(Scene& world) {
    const auto chests = Storage_NearbyChests(world);
    if (chests.empty()) return 0;
    int moved = 0;
    for (const auto& kv : inventory.all()) {
        const std::string& id = kv.first;
        if (kv.second <= 0 || KindOf(id) == ItemKind::Tool) continue;
        bool bound = false;
        for (int s = 0; s < 10; ++s) bound = bound || Hotbar_Slot(s) == id;
        if (bound) continue;
        for (const NearChest& c : chests) {
            bool holds = false;
            for (const auto& s : Factory::ChestContents(c.body, c.plot, c.part)) holds = holds || s.first == id;
            if (!holds) continue;
            const int n = kv.second;
            if (inventory.remove(id, n)) { Factory::Insert(c.body, c.plot, c.part, id, n); ++moved; } // what does not fit comes back
            break;
        }
    }
    return moved;
}
