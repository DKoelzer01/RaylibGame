#ifndef NETINTEREST_H
#define NETINTEREST_H

// Interest management (#340): what a server sends each client, and how often, within a per-client bandwidth budget.
// Each client keeps a priority accumulator per entity (players now; vehicles and animals later). Every snapshot each
// entity's accumulator grows by its relevancy rate (its share of snapshots): the client's own player 1 (always sent),
// players on the same body under 200 m 1, under 2 km 1/2, under 20 km 1/6, farther or on another body 1/60 (0.5 Hz at
// 30 Hz). Entities whose accumulator reached 1 are sent, highest first, until the snapshot's byte budget (the client's
// kbit/s over the snapshot rate, unused bytes carried one snapshot) runs out; the rest wait with their priority still
// rising, so nothing starves. Terrain edits use Interest_EditRelevant (same body, within the edit range).
// Standalone: the server's snapshot loop (hostserver.cpp SendSnapshots) asks Interest_Select; `--test-interest` checks
// the rates, the budget and starvation without a world or a network.

#include <cstdint>
#include <vector>
#include "dvec3.h"

struct InterestConfig {
    float budgetKbps = 128.0f;          // per client (RAYLIBGAME_NET_BUDGET_KBPS)
    int snapshotHz = 30;                // the server's snapshot rate
    float nearM = 200.0f, midM = 2000.0f, farM = 20000.0f;
    float rateNear = 1.0f, rateMid = 0.5f, rateFar = 1.0f / 6.0f, rateOther = 1.0f / 60.0f;
    float editRangeM = 3000.0f;         // terrain edits farther than this (or on another body) wait for the region's load
};

struct InterestItem {
    uint32_t id = 0;                    // e.g. the player's session
    uint64_t body = 0;                  // the body NetId its position is relative to
    DVec3 local{};                      // body-local position (body units; x scale = metres)
    uint32_t bytes = 0;                 // what sending it costs
    bool own = false;                   // the client's own player (always sent)
};

struct InterestStats { uint64_t snapshots = 0, bytes = 0, sent = 0, deferred = 0; };

void Interest_SetConfig(const InterestConfig& c);
const InterestConfig& Interest_Config();
// One snapshot for one client: the ids to send now, in priority order. `viewer*` is the client's own player; bodyScale
// converts body-local units to metres.
std::vector<uint32_t> Interest_Select(uint32_t client, uint64_t viewerBody, const DVec3& viewerLocal, double bodyScale, const std::vector<InterestItem>& items);
float Interest_Rate(uint64_t viewerBody, const DVec3& viewerLocal, double bodyScale, const InterestItem& item); // per snapshot
bool Interest_EditRelevant(uint64_t viewerBody, const DVec3& viewerLocal, double bodyScale, uint64_t editBody, const DVec3& editCentre);
void Interest_Forget(uint32_t client);              // a client left
InterestStats Interest_Stats(uint32_t client);
bool Interest_Test();                               // --test-interest

#endif
