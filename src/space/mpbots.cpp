#include "mpbots.h"
#include "hostserver.h"
#include "net.h"
#include "netclient.h"
#include "player.h"
#include "playerinput.h"
#include "simclock.h"
#include "worldtypes.h"

#include <raylib.h>
#include <raymath.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

// See mpbots.h (#347).
namespace {
struct Bot {
    NetClient c;
    std::unique_ptr<Player> self;   // the predicted player (as the local player is in a game)
    std::string name;
    uint32_t session = 0;           // from the first join
    double joinMs = -1.0;           // connect -> Joined
    int lastPongs = 0, lastStates = 0, lastCorrections = 0;
    bool gone = false;              // left this churn round
    bool synced = false;
    double serverGap = -1.0;        // predicted vs the server's player at the end of the walk (m)
};
std::vector<Bot> bots;
uint16_t port = 0;
int phase = -1, churnRound = 0;
long long stepped = 0;
double t0 = 0.0, phaseT = 0.0, connectT = 0.0;
bool warm = false;
std::vector<double> joinMs, rtts, rejoinMs;
int sessionsKept = 0, rejoins = 0, onlineChecks = 0, onlineOk = 0;
float watchWorst = 0.0f; std::vector<DVec3> watchLast; std::vector<bool> watchSeen; // bot1 watching the others: speed / scripted speed
float ScriptSpeed(size_t k) { return k % 3 == 1 ? 12.0f : 6.0f; } // m/s: walkSpeed 6, sprint x2 (crouching is slower)
FILE* csv = nullptr;
const int kRounds = 3;
const double kWalkSeconds = 8.0, kWarmSeconds = 1.5;

const char* Addr() { return TextFormat("127.0.0.1:%d", port); }
double Pct(std::vector<double> v, double p) {
    if (v.empty()) return -1.0;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, (size_t)(p * (v.size() - 1) + 0.5))];
}
bool AllJoined() { for (const Bot& b : bots) if (!b.gone && b.c.state != NetClient::State::Joined) return false; return true; }

// The scripted inputs (one per sim step): bot k walks a circle (turning 30 deg/s), sprints in a wider one with a jump every
// second, or crouch-walks a zigzag, by k % 3.
PlayerInput Script(int k, long long step, Player& p) {
    PlayerInput in;
    in.hold(Act::Forward);
    const float turn = k % 3 == 1 ? 15.0f : k % 3 == 2 ? (step / 90 % 2 ? 40.0f : -40.0f) : 30.0f; // deg/s
    if (k % 3 == 1) { in.hold(Act::Sprint); if (step % 60 == 0) in.press(Act::Jump); }
    if (k % 3 == 2) in.hold(Act::Down);
    p.forward = Vector3RotateByAxisAngle(p.forward, p.up, turn * DEG2RAD / 60.0f);
    return in;
}

double GapLimit(size_t k);
void Finish(const char* why, bool movedOk, bool pass) {
    const NetConditions nc = Net_GetConditions();
    const NetLinkStats ls = Net_GetLinkStats();
    printf("bench: mp-bots: conditions lag %.0f ms, jitter %.0f ms, loss %.1f%%, reorder %.1f%% (seed %u); link: %llu sent, %llu dropped, %llu delayed, %llu reordered, %llu retransmits\n",
           nc.lagMs, nc.jitterMs, nc.lossPct, nc.reorderPct, nc.seed, (unsigned long long)ls.sent, (unsigned long long)ls.dropped,
           (unsigned long long)ls.delayed, (unsigned long long)ls.reordered, (unsigned long long)ls.retransmits);
    printf("bench: mp-bots: timings: join p50 %.0f / max %.0f ms; rtt p50 %.0f / p95 %.0f / max %.0f ms (%zu pongs); rejoin p50 %.0f / max %.0f ms\n",
           Pct(joinMs, 0.5), Pct(joinMs, 1.0), Pct(rtts, 0.5), Pct(rtts, 0.95), Pct(rtts, 1.0), rtts.size(), Pct(rejoinMs, 0.5), Pct(rejoinMs, 1.0));
    for (const Bot& b : bots)
        printf("bench: mp-bots: %s (session %u): %d own states, prediction off by %.3f m at worst (%.3f last), %d corrections; %.2f m from the server at the end\n",
               b.name.c_str(), b.session, b.c.states, b.c.worstError, b.c.lastError, b.c.corrections, b.serverGap, GapLimit(&b - bots.data()));
    printf("bench: mp-bots: %zu bots joined %zu; moved ok %d (bot1 saw the others move at most %.2f x their scripted speed); churn %d of %d rounds: %d rejoins kept their session %d, "
           "online count right %d of %d; %s after %.1f s -> %s\n",
           bots.size(), joinMs.size(), (int)movedOk, watchWorst, churnRound, kRounds, rejoins, sessionsKept, onlineOk, onlineChecks, why, GetTime() - t0, pass ? "PASS" : "FAIL");
    for (Bot& b : bots) if (!b.gone) NetClient_Disconnect(b.c);
    if (csv) { fclose(csv); csv = nullptr; }
    phase = -1;
}
// The predicted player runs ahead of the server's by about its speed x the round trip (its inputs are still on the way).
double GapLimit(size_t k) { const double rtt = std::max(0.0, Pct(rtts, 0.95)) / 1000.0; return ScriptSpeed(k) * rtt + 0.5; }
bool MovedOk() {
    for (size_t k = 0; k < bots.size(); ++k) {
        const Bot& b = bots[k];
        if (!b.synced || b.c.states < 80 || b.c.worstError >= 0.5 || b.c.corrections * 10 > b.c.states || b.serverGap < 0.0 || b.serverGap > GapLimit(k)) return false;
    }
    return watchWorst < 2.0f;
}
}

bool MpBots_Active() { return phase >= 0; }

void MpBots_Start(int count, uint16_t p) {
    port = p;
    bots.clear();
    bots.resize((size_t)std::max(1, count));
    for (size_t i = 0; i < bots.size(); ++i) bots[i].name = "bot" + std::to_string(i + 1);
    joinMs.clear(); rtts.clear(); rejoinMs.clear();
    watchWorst = 0.0f; watchLast.assign(bots.size(), DVec3{}); watchSeen.assign(bots.size(), false);
    phase = 0; churnRound = 0; warm = false; sessionsKept = rejoins = onlineChecks = onlineOk = 0;
    t0 = GetTime();
    if (const char* path = getenv("RAYLIBGAME_MP_CSV")) { // per-state trace of bot1 (#347: for #327's prediction work)
        csv = fopen(path, "w");
        if (csv) fprintf(csv, "client_tick,server_tick,ack_input,server_x,server_y,server_z,client_x,client_y,client_z,error_m,corrected\n");
    }
}

bool MpBots_Update(Scene& world) {
    if (phase < 0) return true;
    const double now = GetTime();
    for (Bot& b : bots) {
        if (b.gone) continue;
        NetClient_Update(b.c);
        if (b.c.pongs != b.lastPongs) { if (b.c.rttMs >= 0.0) rtts.push_back(b.c.rttMs); b.lastPongs = b.c.pongs; }
    }
    if (!Host_Listening()) { Finish(Net_Available() ? "the host is not listening" : "built without ENABLE_GNS: SKIP", false, false); return true; }
    if (now - t0 > 150.0) { Finish("timed out", false, false); return true; }
    switch (phase) {
    case 0: // everyone connects at once
        for (Bot& b : bots) NetClient_Connect(b.c, Addr(), b.name);
        connectT = now; phase = 1;
        break;
    case 1: { // joins; each bot's player starts where the server put it
        bool ready = true;
        for (Bot& b : bots) {
            if (b.joinMs < 0.0 && b.c.state == NetClient::State::Joined) { b.joinMs = (now - connectT) * 1000.0; joinMs.push_back(b.joinMs); b.session = b.c.welcome.session; }
            if (b.c.state == NetClient::State::Rejected || b.c.state == NetClient::State::Closed) { Finish(TextFormat("%s was %s", b.name.c_str(), NetClient_StateName(b.c.state)), false, false); return true; }
            if (!b.synced && b.c.state == NetClient::State::Joined && b.c.haveOwn)
                if (Planetoid* body = Planetoid::findNet(b.c.own.body)) {
                    b.self = std::make_unique<Player>();
                    b.self->spawnOn(world, *body, body->dirToWorld(Vector3Normalize(b.c.own.local.toFloat())));
                    b.self->applyNetState(world, NetClient_FromMsg(b.c.own));
                    b.synced = true;
                }
            ready = ready && b.synced;
        }
        if (ready) { phase = 2; phaseT = now; stepped = SimTick(); }
        else if (now - connectT > 30.0) { Finish("joins or first states timed out", false, false); return true; }
        break;
    }
    case 2: { // scripted movement, each bot predicting itself
        for (; stepped < SimTick(); ++stepped)
            for (size_t k = 0; k < bots.size(); ++k) {
                Bot& b = bots[k];
                NetClient_Step(b.c, world, *b.self, Script((int)k, stepped, *b.self));
                if (k == 0 && csv && b.c.states != b.lastStates) {
                    const DVec3 cl = b.self->nearest ? b.self->nearest->toLocalD(b.self->position) : DVec3{};
                    fprintf(csv, "%u,%u,%u,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d\n", b.c.tick, b.c.own.tick, b.c.own.ackInput, b.c.own.local.x, b.c.own.local.y, b.c.own.local.z,
                            cl.x, cl.y, cl.z, b.c.lastError, b.c.corrections != b.lastCorrections ? 1 : 0);
                }
                b.lastStates = b.c.states; b.lastCorrections = b.c.corrections;
            }
        if (!warm && now - phaseT > kWarmSeconds) { // the start settles first (as --test-net-move)
            warm = true; watchWorst = 0.0f;
            for (Bot& b : bots) { b.c.worstError = 0.0; b.c.corrections = 0; b.c.states = 0; b.lastStates = b.lastCorrections = 0; }
        }
        // bot1 draws the others: none may move faster than twice its scripted speed from frame to frame (as --test-net-move).
        for (size_t k = 1; k < bots.size(); ++k) {
            uint64_t vb = 0; DVec3 v{};
            if (NetClient_OtherView(bots[0].c, bots[k].c.welcome.session, vb, v)) {
                if (watchSeen[k] && warm && GetFrameTime() > 0.0f) watchWorst = std::max(watchWorst, (float)(v - watchLast[k]).length() / GetFrameTime() / ScriptSpeed(k));
                watchLast[k] = v; watchSeen[k] = true;
            }
        }
        if (now - phaseT > kWalkSeconds) {
            for (Bot& b : bots) {
                uint64_t sb = 0; DVec3 sl{};
                if (Host_PlayerLocal(world, b.c.welcome.session, sb, sl) && b.self->nearest) b.serverGap = (b.self->nearest->toLocalD(b.self->position) - sl).length();
            }
            phase = 3; phaseT = now;
        }
        break;
    }
    case 3: // churn: half leave
        for (size_t i = 0; i < bots.size(); ++i) if ((int)(i % 2) == churnRound % 2) { NetClient_Disconnect(bots[i].c); bots[i].gone = true; }
        phase = 4; phaseT = now;
        break;
    case 4: { // the server sees them go
        int gone = 0; for (const Bot& b : bots) gone += b.gone ? 1 : 0;
        if (Host_OnlinePlayers() == (int)bots.size() - gone || now - phaseT > 20.0) {
            ++onlineChecks; onlineOk += Host_OnlinePlayers() == (int)bots.size() - gone ? 1 : 0;
            for (Bot& b : bots) if (b.gone) {
                const int st = b.c.states, co = b.c.corrections; const double we = b.c.worstError, le = b.c.lastError; // keep the walk's numbers for the report
                b.c = NetClient{}; b.c.states = st; b.c.corrections = co; b.c.worstError = we; b.c.lastError = le;
                b.lastPongs = 0; NetClient_Connect(b.c, Addr(), b.name); b.gone = false; b.joinMs = -2.0;
            }
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
                const bool moved = MovedOk();
                const bool pass = moved && (int)joinMs.size() == (int)bots.size() && sessionsKept == rejoins && rejoins > 0 && onlineOk == onlineChecks;
                Finish("done", moved, pass);
                return true;
            }
            phase = 3; phaseT = now;
        } else if (now - phaseT > 30.0) { Finish("a rejoin timed out", MovedOk(), false); return true; }
        break;
    }
    return false;
}
