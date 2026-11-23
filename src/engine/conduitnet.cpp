#include "conduitnet.h"

#include <algorithm>
#include <cstdlib>
#include <map>
#include <numeric>
#include <queue>
#include <string>
#include <utility>

int64_t Conduit_RunCapacity(const ConduitRun& run, const ConduitMaterial& mat) {
    const int64_t s2 = (int64_t)run.size * run.size;
    return (run.kind == ConduitKind::Cable ? (int64_t)mat.volts : (int64_t)mat.flow) * s2;
}

namespace {
// Closest centreline slot of the run to p, and its arc distance from the first point.
struct Proj { CInt3 at; int64_t arc = 0; int64_t d2 = 0; };

Proj ProjectOnRun(const ConduitRun& run, CInt3 p) {
    Proj best; bool have = false;
    int64_t base = 0;
    for (size_t i = 0; i + 1 < run.points.size(); i++) {
        const CInt3 a = run.points[i], b = run.points[i + 1];
        CInt3 q = a;
        const int64_t len = std::abs(b.x - a.x) + std::abs(b.y - a.y) + std::abs(b.z - a.z);
        int64_t off = 0;
        if (a.x != b.x)      { q.x = std::min(std::max(p.x, std::min(a.x, b.x)), std::max(a.x, b.x)); off = std::abs(q.x - a.x); }
        else if (a.y != b.y) { q.y = std::min(std::max(p.y, std::min(a.y, b.y)), std::max(a.y, b.y)); off = std::abs(q.y - a.y); }
        else if (a.z != b.z) { q.z = std::min(std::max(p.z, std::min(a.z, b.z)), std::max(a.z, b.z)); off = std::abs(q.z - a.z); }
        const int64_t dx = p.x - q.x, dy = p.y - q.y, dz = p.z - q.z;
        const int64_t d2 = dx * dx + dy * dy + dz * dz;
        if (!have || d2 < best.d2) { best.at = q; best.arc = base + off; best.d2 = d2; have = true; }
        base += len;
    }
    if (!have && !run.points.empty()) best.at = run.points[0];
    return best;
}

struct DSU {
    std::vector<int> p;
    explicit DSU(int n) : p(n) { std::iota(p.begin(), p.end(), 0); }
    int find(int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }
    void join(int a, int b) { a = find(a); b = find(b); if (a != b) p[std::max(a, b)] = std::min(a, b); }
};

typedef std::pair<std::string, CInt3> NodeKey; // (kind+material, position): different kinds/materials never share a node
} // namespace

ConduitGraph Conduit_BuildGraph(const std::vector<ConduitRun>& runsIn, const ConduitCatalog& cat, const std::vector<ConduitPortHit>& ports) {
    ConduitGraph g;
    g.portNode.assign(ports.size(), -1);

    std::vector<const ConduitRun*> runs;
    for (const ConduitRun& r : runsIn) if (r.points.size() >= 2) runs.push_back(&r);
    std::stable_sort(runs.begin(), runs.end(), [](const ConduitRun* a, const ConduitRun* b) { return a->id < b->id; });

    std::vector<std::vector<CInt3>> slots(runs.size());
    std::vector<std::string> key(runs.size());
    for (size_t i = 0; i < runs.size(); i++) {
        slots[i] = Conduit_Slots(*runs[i]);
        key[i] = std::string(runs[i]->kind == ConduitKind::Cable ? "c:" : "p:") + runs[i]->material;
    }

    std::map<NodeKey, int> nodeOf;
    auto nodeAt = [&](const std::string& k, CInt3 at) {
        auto it = nodeOf.find(NodeKey(k, at));
        if (it != nodeOf.end()) return it->second;
        const int idx = (int)g.nodes.size();
        ConduitNode n; n.at = at;
        g.nodes.push_back(n);
        nodeOf[NodeKey(k, at)] = idx;
        return idx;
    };

    // Each run's effective end positions: an end inside another run's slots moves onto that run's centreline (the joint).
    std::vector<std::pair<CInt3, CInt3>> endPos(runs.size());
    for (size_t i = 0; i < runs.size(); i++) {
        CInt3 e[2] = { runs[i]->points.front(), runs[i]->points.back() };
        for (int k = 0; k < 2; k++) {
            for (size_t j = 0; j < runs.size(); j++) {
                if (j == i || key[j] != key[i]) continue;
                if (!std::binary_search(slots[j].begin(), slots[j].end(), e[k])) continue;
                e[k] = ProjectOnRun(*runs[j], e[k]).at;
                break; // the first run (lowest id) that holds the end
            }
        }
        endPos[i] = { e[0], e[1] };
    }

    std::vector<std::pair<int, int>> endNodes(runs.size(), { -1, -1 }); // node of the run's first / last point
    for (size_t i = 0; i < runs.size(); i++) {
        const ConduitRun& R = *runs[i];
        struct Item { int64_t arc; CInt3 at; int end; }; // end: 0 first, 1 last, -1 other
        std::vector<Item> items;
        // corners and ends
        int64_t arc = 0;
        for (size_t p = 0; p < R.points.size(); p++) {
            if (p > 0) {
                const CInt3 a = R.points[p - 1], b = R.points[p];
                arc += std::abs(b.x - a.x) + std::abs(b.y - a.y) + std::abs(b.z - a.z);
            }
            Item it; it.arc = arc; it.at = R.points[p]; it.end = -1;
            if (p == 0) { it.end = 0; it.at = endPos[i].first; }
            else if (p + 1 == R.points.size()) { it.end = 1; it.at = endPos[i].second; }
            items.push_back(it);
        }
        // joints: other runs' ends inside this run's slots
        for (size_t j = 0; j < runs.size(); j++) {
            if (j == i || key[j] != key[i]) continue;
            const CInt3 oe[2] = { runs[j]->points.front(), runs[j]->points.back() };
            for (int k = 0; k < 2; k++) {
                if (!std::binary_search(slots[i].begin(), slots[i].end(), oe[k])) continue;
                const Proj pr = ProjectOnRun(R, oe[k]);
                Item it; it.arc = pr.arc; it.at = pr.at; it.end = -1;
                items.push_back(it);
            }
        }
        std::stable_sort(items.begin(), items.end(), [](const Item& a, const Item& b) { return a.arc < b.arc; });

        const ConduitMaterial* mat = cat.find(R.material);
        const int64_t cap = mat ? Conduit_RunCapacity(R, *mat) : 0;
        int prevNode = -1; int64_t prevArc = 0;
        for (const Item& it : items) {
            const int n = nodeAt(key[i], it.at);
            if (it.end == 0) endNodes[i].first = n;
            if (it.end == 1) endNodes[i].second = n;
            if (prevNode >= 0 && n != prevNode && it.arc > prevArc) {
                ConduitEdge e; e.a = prevNode; e.b = n; e.runId = R.id; e.capacity = cap; e.length = it.arc - prevArc;
                g.edges.push_back(e);
            }
            if (prevNode < 0 || n != prevNode) { prevNode = n; prevArc = it.arc; }
        }
    }

    // ports: the first run (by id) with an end exactly on the slot
    for (size_t q = 0; q < ports.size(); q++) {
        for (size_t i = 0; i < runs.size() && g.portNode[q] < 0; i++) {
            int n = -1;
            if (runs[i]->points.front() == ports[q].slot) n = endNodes[i].first;
            else if (runs[i]->points.back() == ports[q].slot) n = endNodes[i].second;
            if (n < 0) continue;
            g.portNode[q] = n;
            if (g.nodes[n].owner < 0) { g.nodes[n].owner = ports[q].owner; g.nodes[n].port = ports[q].port; }
        }
    }

    // components, numbered by the lowest node index
    DSU dsu((int)g.nodes.size());
    for (const ConduitEdge& e : g.edges) dsu.join(e.a, e.b);
    g.component.assign(g.nodes.size(), -1);
    std::map<int, int> idOfRoot;
    for (size_t n = 0; n < g.nodes.size(); n++) {
        const int r = dsu.find((int)n);
        auto it = idOfRoot.find(r);
        if (it == idOfRoot.end()) it = idOfRoot.emplace(r, g.networks++).first;
        g.component[n] = it->second;
    }
    return g;
}

namespace {
struct FlowArc { int to; int64_t cap; };
struct Dinic {
    std::vector<FlowArc> arcs;             // arc i and i^1 are a pair
    std::vector<std::vector<int>> adj;
    std::vector<int> level, it;
    explicit Dinic(int n) : adj(n), level(n), it(n) {}
    int add(int a, int b, int64_t cap, int64_t rev = 0) {
        const int id = (int)arcs.size();
        arcs.push_back({ b, cap }); adj[a].push_back(id);
        arcs.push_back({ a, rev }); adj[b].push_back(id + 1);
        return id;
    }
    bool bfs(int s, int t) {
        std::fill(level.begin(), level.end(), -1);
        std::queue<int> q; q.push(s); level[s] = 0;
        while (!q.empty()) {
            const int u = q.front(); q.pop();
            for (int id : adj[u]) if (arcs[id].cap > 0 && level[arcs[id].to] < 0) { level[arcs[id].to] = level[u] + 1; q.push(arcs[id].to); }
        }
        return level[t] >= 0;
    }
    int64_t dfs(int u, int t, int64_t f) {
        if (u == t) return f;
        for (int& i = it[u]; i < (int)adj[u].size(); i++) {
            const int id = adj[u][i];
            const int v = arcs[id].to;
            if (arcs[id].cap <= 0 || level[v] != level[u] + 1) continue;
            const int64_t d = dfs(v, t, std::min(f, arcs[id].cap));
            if (d > 0) { arcs[id].cap -= d; arcs[id ^ 1].cap += d; return d; }
        }
        return 0;
    }
    void run(int s, int t) {
        while (bfs(s, t)) {
            std::fill(it.begin(), it.end(), 0);
            while (dfs(s, t, INT64_MAX) > 0) {}
        }
    }
};
} // namespace

std::vector<int64_t> Conduit_MaxFlow(const ConduitGraph& g, const std::vector<ConduitFlowNode>& sources,
                                     const std::vector<ConduitFlowNode>& sinks, int ticks) {
    std::vector<int64_t> out(sinks.size(), 0);
    const int N = (int)g.nodes.size();
    if (N == 0 || ticks <= 0) return out;
    const int S = N, T = N + 1;
    Dinic d(N + 2);
    for (const ConduitEdge& e : g.edges) d.add(e.a, e.b, e.capacity * ticks, e.capacity * ticks); // undirected: both directions
    std::map<int, int64_t> supply, demand;
    for (const ConduitFlowNode& s : sources) if (s.node >= 0 && s.node < N && s.amount > 0) supply[s.node] += s.amount;
    for (const ConduitFlowNode& s : sinks) if (s.node >= 0 && s.node < N && s.amount > 0) demand[s.node] += s.amount;
    for (const auto& kv : supply) d.add(S, kv.first, kv.second);
    std::map<int, int> sinkArc;
    for (const auto& kv : demand) sinkArc[kv.first] = d.add(kv.first, T, kv.second);
    d.run(S, T);
    std::map<int, int64_t> delivered; // per node: what reached T through it
    for (const auto& kv : sinkArc) delivered[kv.first] = d.arcs[kv.second ^ 1].cap;
    for (size_t i = 0; i < sinks.size(); i++) {
        const ConduitFlowNode& s = sinks[i];
        if (s.node < 0 || s.node >= N || s.amount <= 0) continue;
        int64_t& left = delivered[s.node];
        const int64_t give = std::min(left, s.amount);
        out[i] = give; left -= give;
    }
    return out;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"

namespace {
ConduitRun NetRun(int id, ConduitKind kind, const std::string& mat, int size, std::vector<CInt3> pts) {
    ConduitRun r;
    r.id = id; r.kind = kind; r.material = mat; r.size = size; r.points = std::move(pts);
    return r;
}
ConduitRun NCable(int id, int size, std::vector<CInt3> pts) { return NetRun(id, ConduitKind::Cable, "copper", size, std::move(pts)); }
ConduitRun NPipe(int id, int size, std::vector<CInt3> pts) { return NetRun(id, ConduitKind::Pipe, "iron", size, std::move(pts)); }
const ConduitCatalog& NetCat() {
    static const ConduitCatalog cat = Conduit_ParseMaterials(
        "cable copper | item cable | units 8 | volts 32\n"
        "pipe iron | item pipe | units 8 | flow 40 | maxtemp 250\n", nullptr);
    return cat;
}
} // namespace

TEST(ConduitNet, StraightRunBetweenPorts) {
    for (int size : { 1, 2 }) {
        const std::vector<ConduitRun> runs = { NCable(1, size, { {0, 0, 0}, {40, 0, 0} }) };
        const std::vector<ConduitPortHit> ports = { { {0, 0, 0}, 7, 0 }, { {40, 0, 0}, 8, 1 } };
        const ConduitGraph g = Conduit_BuildGraph(runs, NetCat(), ports);
        ASSERT_EQ(g.nodes.size(), 2u);
        ASSERT_EQ(g.edges.size(), 1u);
        EXPECT_EQ(g.edges[0].capacity, size == 1 ? 32 : 128);
        EXPECT_EQ(g.edges[0].length, 40);
        EXPECT_EQ(g.edges[0].runId, 1);
        EXPECT_EQ(g.networks, 1);
        ASSERT_EQ(g.portNode.size(), 2u);
        EXPECT_GE(g.portNode[0], 0);
        EXPECT_GE(g.portNode[1], 0);
        EXPECT_NE(g.portNode[0], g.portNode[1]);
        EXPECT_EQ(g.nodes[g.portNode[0]].owner, 7);
        EXPECT_EQ(g.nodes[g.portNode[1]].port, 1);
    }
}

TEST(ConduitNet, LRunHasCornerNode) {
    const ConduitGraph g = Conduit_BuildGraph({ NCable(1, 1, { {0, 0, 0}, {30, 0, 0}, {30, 0, 20} }) }, NetCat(), {});
    EXPECT_EQ(g.nodes.size(), 3u);
    ASSERT_EQ(g.edges.size(), 2u);
    EXPECT_EQ(g.edges[0].length, 30);
    EXPECT_EQ(g.edges[1].length, 20);
    EXPECT_EQ(g.networks, 1);
}

TEST(ConduitNet, TJointSplitsTheRun) {
    const ConduitGraph g = Conduit_BuildGraph({ NCable(1, 1, { {0, 0, 0}, {40, 0, 0} }), NCable(2, 1, { {20, 0, 20}, {20, 0, 0} }) }, NetCat(), {});
    EXPECT_EQ(g.nodes.size(), 4u);
    EXPECT_EQ(g.edges.size(), 3u);
    EXPECT_EQ(g.networks, 1);
    int64_t onA = 0;
    for (const ConduitEdge& e : g.edges) if (e.runId == 1) onA += e.length;
    EXPECT_EQ(onA, 40);
}

TEST(ConduitNet, TJointProjectsOntoTheCentreline) {
    // Run B's end sits one slot off the centreline of the size-2 run A, inside A's slots (a size-2 cross-section spans
    // -1..0 around the centreline): it moves onto the centreline. (B comes from -z, so its end is at z = -1.)
    const ConduitGraph g = Conduit_BuildGraph({ NCable(1, 2, { {0, 0, 0}, {40, 0, 0} }), NCable(2, 1, { {20, 0, -20}, {20, 0, -1} }) }, NetCat(), {});
    EXPECT_EQ(g.nodes.size(), 4u);
    EXPECT_EQ(g.networks, 1);
    bool found = false;
    for (const ConduitNode& n : g.nodes) if (n.at == CInt3{ 20, 0, 0 }) found = true;
    EXPECT_TRUE(found);
}

TEST(ConduitNet, EndToEndSharesANode) {
    const ConduitGraph g = Conduit_BuildGraph({ NCable(1, 1, { {0, 0, 0}, {20, 0, 0} }), NCable(2, 1, { {20, 0, 0}, {20, 0, 30} }) }, NetCat(), {});
    EXPECT_EQ(g.nodes.size(), 3u);
    EXPECT_EQ(g.edges.size(), 2u);
    EXPECT_EQ(g.networks, 1);
}

TEST(ConduitNet, SeparateRunsAndKindsDontJoin) {
    ConduitGraph g = Conduit_BuildGraph({ NCable(1, 1, { {0, 0, 0}, {10, 0, 0} }), NCable(2, 1, { {0, 0, 100}, {10, 0, 100} }) }, NetCat(), {});
    EXPECT_EQ(g.networks, 2);
    ASSERT_EQ(g.component.size(), 4u);
    EXPECT_EQ(g.component[0], 0);
    EXPECT_EQ(g.component[1], 0);
    EXPECT_EQ(g.component[2], 1);
    EXPECT_EQ(g.component[3], 1);
    g = Conduit_BuildGraph({ NCable(1, 1, { {0, 0, 0}, {10, 0, 0} }), NPipe(2, 1, { {10, 0, 0}, {20, 0, 0} }) }, NetCat(), {});
    EXPECT_EQ(g.networks, 2);
    EXPECT_EQ(g.nodes.size(), 4u);
}

TEST(ConduitNet, OrderIndependent) {
    std::vector<ConduitRun> runs = { NCable(1, 1, { {0, 0, 0}, {40, 0, 0} }), NCable(2, 1, { {20, 0, 20}, {20, 0, 0} }),
                                     NCable(3, 2, { {40, 0, 0}, {40, 0, 30}, {60, 0, 30} }), NCable(4, 1, { {0, 0, 200}, {9, 0, 200} }) };
    const std::vector<ConduitPortHit> ports = { { {0, 0, 0}, 1, 0 }, { {60, 0, 30}, 2, 0 } };
    const ConduitGraph a = Conduit_BuildGraph(runs, NetCat(), ports);
    std::reverse(runs.begin(), runs.end());
    const ConduitGraph b = Conduit_BuildGraph(runs, NetCat(), ports);
    ASSERT_EQ(a.nodes.size(), b.nodes.size());
    ASSERT_EQ(a.edges.size(), b.edges.size());
    for (size_t i = 0; i < a.nodes.size(); i++) {
        EXPECT_TRUE(a.nodes[i].at == b.nodes[i].at);
        EXPECT_EQ(a.nodes[i].owner, b.nodes[i].owner);
        EXPECT_EQ(a.nodes[i].port, b.nodes[i].port);
    }
    for (size_t i = 0; i < a.edges.size(); i++) {
        EXPECT_EQ(a.edges[i].a, b.edges[i].a);
        EXPECT_EQ(a.edges[i].b, b.edges[i].b);
        EXPECT_EQ(a.edges[i].runId, b.edges[i].runId);
        EXPECT_EQ(a.edges[i].capacity, b.edges[i].capacity);
        EXPECT_EQ(a.edges[i].length, b.edges[i].length);
    }
    EXPECT_EQ(a.component, b.component);
    EXPECT_EQ(a.networks, b.networks);
    EXPECT_EQ(a.portNode, b.portNode);
}

TEST(ConduitNet, MaxFlowSingleCable) {
    const std::vector<ConduitPortHit> ports = { { {0, 0, 0}, 0, 0 }, { {40, 0, 0}, 1, 0 } };
    const ConduitGraph g = Conduit_BuildGraph({ NCable(1, 1, { {0, 0, 0}, {40, 0, 0} }) }, NetCat(), ports);
    const std::vector<ConduitFlowNode> src = { { g.portNode[0], 100 } }, snk = { { g.portNode[1], 100 } };
    EXPECT_EQ(Conduit_MaxFlow(g, src, snk, 1)[0], 32);
    EXPECT_EQ(Conduit_MaxFlow(g, src, snk, 4)[0], 100);
}

TEST(ConduitNet, MaxFlowBranchesAndSplits) {
    // Source at the origin; a size-2 branch along +x and a size-1 branch along +z.
    const std::vector<ConduitPortHit> ports = { { {0, 0, 0}, 0, 0 }, { {40, 0, 0}, 1, 0 }, { {0, 0, 40}, 2, 0 } };
    const ConduitGraph g = Conduit_BuildGraph({ NCable(1, 2, { {0, 0, 0}, {40, 0, 0} }), NCable(2, 1, { {0, 0, 0}, {0, 0, 40} }) }, NetCat(), ports);
    const std::vector<ConduitFlowNode> src = { { g.portNode[0], 1000 } };
    const std::vector<ConduitFlowNode> snk = { { g.portNode[2], 50 }, { g.portNode[1], 50 } };
    const std::vector<int64_t> got = Conduit_MaxFlow(g, src, snk, 1);
    ASSERT_EQ(got.size(), 2u);
    EXPECT_EQ(got[0], 32);
    EXPECT_EQ(got[1], 50);

    // Two sinks on one node share what arrives, in input order.
    const std::vector<ConduitPortHit> p2 = { { {0, 0, 0}, 0, 0 }, { {40, 0, 0}, 1, 0 } };
    const ConduitGraph g2 = Conduit_BuildGraph({ NCable(1, 1, { {0, 0, 0}, {40, 0, 0} }) }, NetCat(), p2);
    const std::vector<int64_t> two = Conduit_MaxFlow(g2, { { g2.portNode[0], 100 } }, { { g2.portNode[1], 20 }, { g2.portNode[1], 20 } }, 1);
    ASSERT_EQ(two.size(), 2u);
    EXPECT_EQ(two[0], 20);
    EXPECT_EQ(two[1], 12);
}

TEST(ConduitNet, MaxFlowAcrossNetworksIsZero) {
    const std::vector<ConduitPortHit> ports = { { {0, 0, 0}, 0, 0 }, { {40, 0, 100}, 1, 0 } };
    const ConduitGraph g = Conduit_BuildGraph({ NCable(1, 1, { {0, 0, 0}, {40, 0, 0} }), NCable(2, 1, { {0, 0, 100}, {40, 0, 100} }) }, NetCat(), ports);
    ASSERT_GE(g.portNode[0], 0);
    ASSERT_GE(g.portNode[1], 0);
    EXPECT_EQ(Conduit_MaxFlow(g, { { g.portNode[0], 100 } }, { { g.portNode[1], 100 } }, 4)[0], 0);
}

TEST(ConduitNet, UnreachedPortHasNoNode) {
    const std::vector<ConduitPortHit> ports = { { {0, 0, 0}, 0, 0 }, { {500, 0, 0}, 1, 0 }, { {20, 0, 0}, 2, 0 } };
    const ConduitGraph g = Conduit_BuildGraph({ NCable(1, 1, { {0, 0, 0}, {40, 0, 0} }) }, NetCat(), ports);
    EXPECT_GE(g.portNode[0], 0);
    EXPECT_EQ(g.portNode[1], -1);
    EXPECT_EQ(g.portNode[2], -1); // mid-run, not an end
}
#endif
