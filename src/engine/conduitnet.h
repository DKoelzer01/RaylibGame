#ifndef CONDUITNET_H
#define CONDUITNET_H

// (#1620, conduits infra PR C2a) Conduit graph networks and the capacity solve. Pure code: runs -> a graph of nodes
// (corners, ends, joints, ports) and edges (straight run stretches with a capacity), its networks, and a deterministic
// integer max-flow from sources to sinks. No game code calls it yet.

#include "conduit.h"

#include <cstdint>
#include <vector>

struct ConduitPortHit { CInt3 slot; int owner; int port; };   // a machine's port: the slot a conduit end must sit on; owner = machine index, port = its port index
struct ConduitNode { CInt3 at; int owner = -1, port = -1; };   // a centreline slot; owner/port set when a port sits there
struct ConduitEdge { int a, b; int runId; int64_t capacity; int64_t length; }; // a run's centreline stretch between two nodes

struct ConduitGraph {
    std::vector<ConduitNode> nodes;
    std::vector<ConduitEdge> edges;
    std::vector<int> component;                // per node: its network id (0..networks-1), numbered in order of the lowest node index
    int networks = 0;
    std::vector<int> portNode;                 // per ConduitPortHit given (same order): its node, -1 if no conduit reaches it
};

ConduitGraph Conduit_BuildGraph(const std::vector<ConduitRun>& runs, const ConduitCatalog& cat, const std::vector<ConduitPortHit>& ports);
int64_t Conduit_RunCapacity(const ConduitRun&, const ConduitMaterial&); // cable: volts * size^2 (EU/t); pipe: flow * size^2 (L/t)

struct ConduitFlowNode { int node; int64_t amount; };   // a source's supply or a sink's demand this step (>= 0)
// Max flow from sources to sinks over the graph's edges (undirected: each edge carries up to its capacity * ticks in
// either direction). Returns each sink's delivered amount (same order as sinks).
std::vector<int64_t> Conduit_MaxFlow(const ConduitGraph& g, const std::vector<ConduitFlowNode>& sources,
                                     const std::vector<ConduitFlowNode>& sinks, int ticks);

#endif
