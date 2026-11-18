#include "mpbots.h"
#include "hostserver.h"
#include "net.h"
#include "netclient.h"

#include <raylib.h>
#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

// See mpbots.h (#347).
namespace {
struct Bot {
    NetClient c;
    std::string name;
    uint32_t session = 0;     // from the first join
    double joinMs = -1.0;     // connect -> Joined
    int lastPongs = 0;
    bool gone = false;        // left this churn churnRound
};
std::vector<Bot> bots;
uint16_t port = 0;
int phase = -1, churnRound = 0;
double t0 = 0.0, phaseT = 0.0, connectT = 0.0;
std::vector<double> joinMs, rtts, rejoinMs;
int sessionsKept = 0, rejoins = 0, onlineChecks = 0, onlineOk = 0;
const int kRounds = 3;

const char* Addr() { return TextFormat("127.0.0.1:%d", port); }
double Pct(std::vector<double> v, double p) {
    if (v.empty()) return -1.0;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, (size_t)(p * (v.size() - 1) + 0.5))];
}
bool AllJoined() { for (const Bot& b : bots) if (!b.gone && b.c.state != NetClient::State::Joined) return false; return true; }
void Finish(const char* why, bool pass) {
    const NetConditions nc = Net_GetConditions();
    const NetLinkStats ls = Net_GetLinkStats();
    printf("bench: mp-bots: conditions lag %.0f ms, jitter %.0f ms, loss %.1f%%, reorder %.1f%% (seed %u); link: %llu sent, %llu dropped, %llu delayed, %llu reordered, %llu retransmits\n",
           nc.lagMs, nc.jitterMs, nc.lossPct, nc.reorderPct, nc.seed, (unsigned long long)ls.sent, (unsigned long long)ls.dropped,
           (unsigned long long)ls.delayed, (unsigned long long)ls.reordered, (unsigned long long)ls.retransmits);
    printf("bench: mp-bots: timings: join p50 %.0f / max %.0f ms; rtt p50 %.0f / p95 %.0f / max %.0f ms (%zu pongs); rejoin p50 %.0f / max %.0f ms\n",
           Pct(joinMs, 0.5), Pct(joinMs, 1.0), Pct(rtts, 0.5), Pct(rtts, 0.95), Pct(rtts, 1.0), rtts.size(), Pct(rejoinMs, 0.5), Pct(rejoinMs, 1.0));
    printf("bench: mp-bots: %zu bots joined %zu; churn %d of %d rounds: %d rejoins kept their session %d, online count right %d of %d; "
           "position agreement SKIP (player replication is #327); %s after %.1f s -> %s\n",
           bots.size(), joinMs.size(), churnRound, kRounds, rejoins, sessionsKept, onlineOk, onlineChecks, why, GetTime() - t0, pass ? "PASS" : "FAIL");
    for (Bot& b : bots) if (!b.gone) NetClient_Disconnect(b.c);
    phase = -1;
}
}

bool MpBots_Active() { return phase >= 0; }

void MpBots_Start(int count, uint16_t p) {
    port = p;
    bots.assign((size_t)std::max(1, count), Bot{});
    for (size_t i = 0; i < bots.size(); ++i) bots[i].name = "bot" + std::to_string(i + 1);
    joinMs.clear(); rtts.clear(); rejoinMs.clear();
    phase = 0; churnRound = 0; sessionsKept = rejoins = onlineChecks = onlineOk = 0;
    t0 = GetTime();
}

bool MpBots_Update(Scene&) {
    if (phase < 0) return true;
    const double now = GetTime();
    for (Bot& b : bots) {
        if (b.gone) continue;
        NetClient_Update(b.c);
        if (b.c.pongs != b.lastPongs) { if (b.c.rttMs >= 0.0) rtts.push_back(b.c.rttMs); b.lastPongs = b.c.pongs; }
    }
    if (!Host_Listening()) { Finish(Net_Available() ? "the host is not listening" : "built without ENABLE_GNS", false); return true; }
    if (now - t0 > 120.0) { Finish("timed out", false); return true; }
    switch (phase) {
    case 0: // everyone connects at once
        for (Bot& b : bots) NetClient_Connect(b.c, Addr(), b.name);
        connectT = now; phase = 1;
        break;
    case 1: // joins
        for (Bot& b : bots)
            if (b.joinMs < 0.0 && b.c.state == NetClient::State::Joined) { b.joinMs = (now - connectT) * 1000.0; joinMs.push_back(b.joinMs); b.session = b.c.welcome.session; }
        for (const Bot& b : bots) if (b.c.state == NetClient::State::Rejected || b.c.state == NetClient::State::Closed) { Finish(TextFormat("%s was %s", b.name.c_str(), NetClient_StateName(b.c.state)), false); return true; }
        if (AllJoined()) { phase = 2; phaseT = now; }
        break;
    case 2: // pings under the conditions
        if (now - phaseT > 6.0) { phase = 3; phaseT = now; }
        break;
    case 3: // churn: half leave
        for (size_t i = 0; i < bots.size(); ++i) if ((int)(i % 2) == churnRound % 2) { NetClient_Disconnect(bots[i].c); bots[i].gone = true; }
        phase = 4; phaseT = now;
        break;
    case 4: { // the server sees them go
        int gone = 0; for (const Bot& b : bots) gone += b.gone ? 1 : 0;
        if (Host_OnlinePlayers() == (int)bots.size() - gone || now - phaseT > 20.0) {
            ++onlineChecks; onlineOk += Host_OnlinePlayers() == (int)bots.size() - gone ? 1 : 0;
            for (Bot& b : bots) if (b.gone) { b.c = NetClient{}; b.lastPongs = 0; NetClient_Connect(b.c, Addr(), b.name); b.gone = false; b.joinMs = -2.0; }
            phase = 5; phaseT = now;
        }
        break;
    }
    case 5: // and come back to their sessions
        for (Bot& b : bots)
            if (b.joinMs == -2.0 && b.c.state == NetClient::State::Joined) {
                b.joinMs = 0.0; ++rejoins; rejoinMs.push_back((now - phaseT) * 1000.0);
                sessionsKept += b.c.welcome.session == b.session && b.c.welcome.reconnected ? 1 : 0;
            }
        if (AllJoined()) {
            ++onlineChecks; onlineOk += Host_OnlinePlayers() == (int)bots.size() ? 1 : 0;
            if (++churnRound >= kRounds) {
                const bool pass = (int)joinMs.size() == (int)bots.size() && sessionsKept == rejoins && rejoins > 0 && onlineOk == onlineChecks;
                Finish("done", pass);
                return true;
            }
            phase = 3; phaseT = now;
        } else if (now - phaseT > 30.0) { Finish("a rejoin timed out", false); return true; }
        break;
    }
    return false;
}
