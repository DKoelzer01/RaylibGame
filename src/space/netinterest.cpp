#include "netinterest.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>

// See netinterest.h (#340).
namespace {
InterestConfig config;
bool envRead = false;
struct ClientState {
    std::unordered_map<uint32_t, float> acc; // id -> priority
    double carry = 0.0;                       // unused bytes from the last snapshot
    InterestStats stats;
};
std::unordered_map<uint32_t, ClientState> clients;

void ReadEnv() {
    if (envRead) return;
    envRead = true;
    if (const char* s = getenv("RAYLIBGAME_NET_BUDGET_KBPS")) config.budgetKbps = std::max(1.0f, (float)atof(s));
}
}

void Interest_SetConfig(const InterestConfig& c) { config = c; envRead = true; }
const InterestConfig& Interest_Config() { ReadEnv(); return config; }

float Interest_Rate(uint64_t viewerBody, const DVec3& viewerLocal, double bodyScale, const InterestItem& it) {
    if (it.own) return 1.0f;
    if (it.body != viewerBody) return config.rateOther;
    const double d = (it.local - viewerLocal).length() * bodyScale;
    return d < config.nearM ? config.rateNear : d < config.midM ? config.rateMid : d < config.farM ? config.rateFar : config.rateOther;
}

bool Interest_EditRelevant(uint64_t viewerBody, const DVec3& viewerLocal, double bodyScale, uint64_t editBody, const DVec3& editCentre) {
    return editBody == viewerBody && (editCentre - viewerLocal).length() * bodyScale <= config.editRangeM;
}

std::vector<uint32_t> Interest_Select(uint32_t client, uint64_t viewerBody, const DVec3& viewerLocal, double bodyScale, const std::vector<InterestItem>& items) {
    ReadEnv();
    ClientState& cs = clients[client];
    ++cs.stats.snapshots;
    // Grow every accumulator; drop the ones for entities that are gone.
    std::unordered_map<uint32_t, float> next;
    next.reserve(items.size());
    struct Cand { const InterestItem* it; float pri; };
    std::vector<Cand> ready;
    for (const InterestItem& it : items) {
        auto old = cs.acc.find(it.id);
        const float a = (old == cs.acc.end() ? 1.0f : old->second) + Interest_Rate(viewerBody, viewerLocal, bodyScale, it); // new: sent soon
        next[it.id] = a;
        if (it.own || a >= 1.0f) ready.push_back({ &it, it.own ? 1e9f : a });
    }
    cs.acc.swap(next);
    std::sort(ready.begin(), ready.end(), [](const Cand& a, const Cand& b) { return a.pri > b.pri; });
    const double perSnapshot = config.budgetKbps * 1000.0 / 8.0 / std::max(1, config.snapshotHz);
    double budget = perSnapshot + cs.carry;
    std::vector<uint32_t> out;
    for (const Cand& c : ready) {
        if (!c.it->own && (double)c.it->bytes > budget) { ++cs.stats.deferred; continue; } // waits, still rising
        budget -= c.it->bytes;
        out.push_back(c.it->id);
        cs.stats.bytes += c.it->bytes; ++cs.stats.sent;
        cs.acc[c.it->id] = 0.0f; // sent: back to the start (the others keep rising, so they get their turn: no starvation)
    }
    cs.carry = std::clamp(budget, 0.0, perSnapshot); // at most one snapshot's worth carried
    return out;
}

void Interest_Forget(uint32_t client) { clients.erase(client); }
InterestStats Interest_Stats(uint32_t client) { auto it = clients.find(client); return it == clients.end() ? InterestStats{} : it->second.stats; }

// --test-interest: 32 players around one viewer (itself, 5 near, 8 at 1 km, 8 at 10 km, 10 on another body), 10 s of
// 30 Hz snapshots. With room in the budget each group is sent at its rate; on a tight budget the bytes stay within it,
// the own player is never skipped, nearer groups still get more, and no player waits longer than its bound.
bool Interest_Test() {
    const InterestConfig saved = config;
    bool ok = true;
    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: interest: %s FAILED\n", what); ok = false; } };
    const uint64_t home = 1, other = 2;
    const double scale = 1.0;
    std::vector<InterestItem> items;
    auto add = [&](uint32_t id, uint64_t body, double dist, bool own) { items.push_back({ id, body, DVec3(dist, 0.0, 0.0), 37 + 8, own }); };
    add(1, home, 0.0, true);
    for (int i = 0; i < 5; ++i) add(10 + i, home, 50.0 + 20.0 * i, false);
    for (int i = 0; i < 8; ++i) add(20 + i, home, 1000.0 + 50.0 * i, false);
    for (int i = 0; i < 8; ++i) add(30 + i, home, 10000.0 + 300.0 * i, false);
    for (int i = 0; i < 10; ++i) add(40 + i, other, 500.0, false);
    auto group = [](uint32_t id) { return id == 1 ? 0 : id < 20 ? 1 : id < 30 ? 2 : id < 40 ? 3 : 4; };
    const char* names[] = { "own", "near", "1 km", "10 km", "other body" };
    int kSnaps = 300; // 10 s at 30 Hz (60 s on the tight budget: the slowest group waits ~20 s there)
    for (int pass = 0; pass < 2; ++pass) {
        InterestConfig c; c.budgetKbps = pass == 0 ? 1000.0f : 24.0f; // roomy, then tight (~100 bytes a snapshot: 2 states)
        kSnaps = pass == 0 ? 300 : 1800;
        Interest_SetConfig(c);
        const uint32_t client = 100 + pass;
        Interest_Forget(client);
        double count[5] = {}; int lastSent[64]; int worstGap[5] = {};
        for (int& l : lastSent) l = -1;
        for (int s = 0; s < kSnaps; ++s) {
            const std::vector<uint32_t> sent = Interest_Select(client, home, DVec3(0.0, 0.0, 0.0), scale, items);
            for (uint32_t id : sent) {
                const int g = group(id); count[g] += 1.0;
                if (lastSent[id] >= 0) worstGap[g] = std::max(worstGap[g], s - lastSent[id]);
                lastSent[id] = s;
            }
        }
        const int per[5] = { 1, 5, 8, 8, 10 };
        double rate[5];
        for (int g = 0; g < 5; ++g) rate[g] = count[g] / (per[g] * (double)kSnaps);
        const InterestStats st = Interest_Stats(client);
        const double kbps = st.bytes * 8.0 / 1000.0 / (kSnaps / 30.0);
        printf("bench: interest: %s budget %.0f kbit/s: sent %.1f kbit/s; per-player share of snapshots: own %.2f, near %.2f, 1 km %.2f, 10 km %.3f, other body %.3f; "
               "longest wait (snapshots): near %d, 1 km %d, 10 km %d, other %d; %llu deferred\n",
               pass == 0 ? "roomy" : "tight", c.budgetKbps, kbps, rate[0], rate[1], rate[2], rate[3], rate[4], worstGap[1], worstGap[2], worstGap[3], worstGap[4],
               (unsigned long long)st.deferred);
        check(rate[0] > 0.999, "own player every snapshot");
        if (pass == 0) {
            check(fabs(rate[1] - 1.0) < 0.02 && fabs(rate[2] - 0.5) < 0.02 && fabs(rate[3] - 1.0 / 6.0) < 0.02 && fabs(rate[4] - 1.0 / 60.0) < 0.01, "rates with room in the budget");
            (void)names;
        } else {
            check(kbps <= c.budgetKbps * 1.05 + 1.0, "within the budget");
            check(rate[1] > rate[2] && rate[2] > rate[3] && rate[3] >= rate[4], "nearer players get more");
            // Each waits about (sum of rates / its rate) snapshots: near ~11, another body ~630 (21 s).
            check(worstGap[1] <= 15 && worstGap[4] > 0 && worstGap[4] <= 900, "no player starves (near within 0.5 s, another body within 30 s)");
        }
    }
    // Edits: same body within range only.
    check(Interest_EditRelevant(home, DVec3(0, 0, 0), 1.0, home, DVec3(2500, 0, 0)) && !Interest_EditRelevant(home, DVec3(0, 0, 0), 1.0, home, DVec3(3500, 0, 0)) &&
          !Interest_EditRelevant(home, DVec3(0, 0, 0), 1.0, other, DVec3(10, 0, 0)), "edit relevancy");
    config = saved;
    printf("bench: interest: rates, budget, starvation and edit relevancy -> %s\n", ok ? "PASS" : "FAIL");
    return ok;
}
