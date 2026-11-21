#include "magicmachines.h"

#include "manafield.h"
#include "inventory.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <vector>

std::string MagicMachines_SiphonJar(const ManaField& f, Vector3 local) {
    const float r = Vector3Length(local);
    if (r < 1.0f) return "jar_none";
    const Vector3 dir = Vector3Scale(local, 1.0f / r);
    std::vector<ManaNode> nodes;
    f.nodesNear(dir, 8.0f, nodes);
    const ManaNode* best = nullptr; float bestD = 1e9f;
    for (const ManaNode& n : nodes) {
        const float d = acosf(std::clamp(Vector3DotProduct(dir, n.dir), -1.0f, 1.0f));
        if (d < bestD) { bestD = d; best = &n; }
    }
    if (!best) return "jar_none";
    int kind = best->kind;
    if (kind < 0 || !Mana_Kind(kind).slice) kind = Mana_Find(best->a == kManaEarth ? "earth" : best->a == kManaWater ? "water" : best->a == kManaFire ? "fire" : "air");
    return kind >= 0 ? "jar_" + Mana_Kind(kind).id : "jar_none";
}

// ---- The game ----
#include "factory.h"
#include "magic.h"
#include "net.h"
#include "player.h"
#include "plots.h"
#include "world.h"
#include "worldtypes.h"

namespace {
double nextAt = 0.0;
bool nearTable = false;
}

bool MagicMachines_NearResearchTable() { return nearTable; }

void MagicMachines_Update(Scene& world, double now) {
    if (!Magic_Enabled()) return;
    if (nextAt - now > 1.0) nextAt = now; // (the sim clock went back: a new world or a load)
    if (now < nextAt) return;
    nextAt = now + 0.5;
    bool near = false;
    for (auto& o : world.objects) {
        Planetoid* body = dynamic_cast<Planetoid*>(o.get());
        if (!body || !body->plots) continue;
        const ManaField* f = Magic_FieldFor(*body);
        const Vector3 me = body->toLocalOf(player.position, world.renderOrigin);
        for (auto& plot : body->plots->plots) {
            for (size_t i = 0; i < plot->parts.size(); ++i) {
                const Part& part = plot->parts[i];
                if (part.type != PartType::ManaSiphon && part.type != PartType::ResearchTable) continue;
                const Int3 c = part.origin;
                const Vector3 at = plot->toBody({ (c.x + 1) * kCellSize, (c.y + 1) * kCellSize, (c.z + 1) * kCellSize });
                if (part.type == PartType::ResearchTable) { near = near || Vector3Distance(at, me) <= 6.0f; continue; }
                if (gNetRole == NetRole::Client) continue; // (the server steers the machines)
                Factory::SetVariant(body->id, plot->id, (int)i, f ? MagicMachines_SiphonJar(*f, at) : std::string("jar_none"));
            }
        }
    }
    nearTable = near;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1252): a siphon fills the jar of the node it stands on, and idles elsewhere ----
#include <gtest/gtest.h>

TEST(MagicMachines, ASiphonFillsItsNodesJar) {
    ASSERT_TRUE(Mana_Load());
    ManaFieldParams p; p.seed = 42; p.radius = 47000.0f;
    ManaField f(p);
    std::vector<ManaNode> nodes;
    f.nodesNear(Vector3Normalize({ 1.0f, 0.2f, 0.3f }), 6000.0f, nodes);
    ASSERT_FALSE(nodes.empty());
    int checked = 0;
    for (const ManaNode& n : nodes) {
        const std::string jar = MagicMachines_SiphonJar(f, Vector3Scale(n.dir, p.radius + 30.0f));
        ASSERT_NE(jar, "jar_none") << "on a node";
        const std::string kind = jar.substr(4);
        const int k = Mana_Find(kind);
        ASSERT_GE(k, 0) << jar;
        EXPECT_TRUE(Mana_Kind(k).slice) << "only the first slice's jars: " << jar;
        if (n.kind >= 0 && Mana_Kind(n.kind).slice) { EXPECT_EQ(k, n.kind) << "the node's own kind"; }
        if (++checked >= 20) break;
    }
    // Far from any node: idle.
    int idle = 0;
    for (int i = 0; i < 50; ++i) {
        const Vector3 d = Vector3Normalize({ sinf(i * 1.7f), cosf(i * 0.9f), sinf(i * 0.4f + 2.0f) });
        std::vector<ManaNode> nn; f.nodesNear(d, 8.0f, nn);
        if (nn.empty()) { EXPECT_EQ(MagicMachines_SiphonJar(f, Vector3Scale(d, p.radius)), "jar_none"); ++idle; }
    }
    EXPECT_GT(idle, 30) << "most places are not on a node";
}
#endif
