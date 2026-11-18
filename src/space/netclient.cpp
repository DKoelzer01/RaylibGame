#include "netclient.h"
#include "hostserver.h"
#include "net.h"
#include "player.h"
#include "scene.h"
#include "simclock.h"
#include "worldtypes.h"
#include <raymath.h>
#include <algorithm>

#include <raylib.h>
#include <cstdio>
#include <vector>

extern const char* kGameVersion; // "<branch> <commit>" (#501)

namespace {
    template <class M> void SendMsg(const NetClient& c, const M& m) {
        NetWriter w; NetWrite(w, m);
        Net_Send(c.conn, w.bytes.data(), (uint32_t)w.bytes.size(), true);
    }
    std::string Commit() { const std::string v = kGameVersion; const size_t sp = v.rfind(' '); return sp == std::string::npos ? v : v.substr(sp + 1); }
    void End(NetClient& c, NetClient::State s) {
        if (c.conn) Net_Close(c.conn);
        c.conn = 0;
        c.state = s;
    }
    void OnMessage(NetClient& c, const NetMessage& m) {
        c.lastHeard = GetTime();
        NetReader r(m.data.data() + (m.data.empty() ? 0 : 1), m.data.empty() ? 0 : m.data.size() - 1);
        switch (NetPeek(m.data)) {
            case NetMsg::Welcome: {
                MsgWelcome w;
                if (c.state != NetClient::State::AwaitingWelcome || !NetRead(r, w)) return;
                c.welcome = w; c.serverTick = w.serverTick; c.state = NetClient::State::Joined;
                return;
            }
            case NetMsg::Reject: {
                MsgReject j;
                if (!NetRead(r, j)) return;
                c.reject = j;
                End(c, NetClient::State::Rejected);
                return;
            }
            case NetMsg::Bye: { MsgBye b; if (NetRead(r, b)) { c.bye = b.reason; End(c, NetClient::State::Closed); } return; }
            case NetMsg::Pong: {
                MsgPong p;
                if (!NetRead(r, p)) return;
                c.rttMs = (GetTime() - p.clientTime) * 1000.0; c.serverTick = p.serverTick; ++c.pongs;
                return;
            }
            case NetMsg::Shower: { MsgShower s; if (NetRead(r, s)) c.shower = s; return; }
            case NetMsg::PlayerState: {
                MsgPlayerState s;
                if (c.state != NetClient::State::Joined || !NetRead(r, s)) return;
                if (s.session == c.welcome.session) {
                    if (!c.haveOwn || s.tick > c.own.tick) { c.own = s; c.haveOwn = true; } // (unreliable: may arrive out of order)
                    ++c.states;
                } else {
                    auto& q = c.others[s.session];
                    if (q.empty() || s.tick > q.back().tick) q.push_back(s);
                    while (q.size() > 32) q.pop_front();
                }
                return;
            }
            default: return;
        }
    }
}

const char* NetClient_StateName(NetClient::State s) {
    switch (s) {
        case NetClient::State::Idle: return "idle"; case NetClient::State::Connecting: return "connecting";
        case NetClient::State::AwaitingWelcome: return "awaiting welcome"; case NetClient::State::Joined: return "joined";
        case NetClient::State::Rejected: return "rejected"; case NetClient::State::Closed: return "closed";
    }
    return "?";
}

bool NetClient_Connect(NetClient& c, const char* address, const std::string& name, const std::string& password) {
    if (!Host_Listening()) { std::string err; if (!Net_Init(err)) { c.state = NetClient::State::Closed; return false; } }
    c.name = name; c.password = password;
    c.conn = Net_Connect(address);
    if (!c.conn) { c.state = NetClient::State::Closed; return false; }
    if (Host_Listening()) Host_SetLocalClient(c.conn); // a host in this process receives for us
    c.state = NetClient::State::Connecting;
    c.started = c.lastHeard = GetTime();
    return true;
}

void NetClient_Update(NetClient& c) {
    if (c.state == NetClient::State::Idle || c.state == NetClient::State::Rejected || c.state == NetClient::State::Closed) return;
    // The others' playback clock (server ticks): runs at real time and eases toward kInterpDelay behind the newest state
    // heard, so the others move smoothly between snapshots that arrive in bursts (keyed to the newest state, it jumped).
    {
        uint32_t newest = 0;
        for (const auto& [session, q] : c.others) if (!q.empty()) newest = std::max(newest, q.back().tick);
        if (newest) {
            const double target = newest - kInterpDelay * kSimHz;
            c.viewTick += GetFrameTime() * kSimHz;
            if (c.viewTick <= 0.0 || fabs(c.viewTick - target) > 30.0) c.viewTick = target; // (first, or 0.5 s off: jump)
            else c.viewTick += (target - c.viewTick) * (1.0 - exp(-GetFrameTime() / 1.0)); // a 1 s time constant (per frame it was frame-rate bound)
        }
    }
    std::vector<NetMessage> in;
    if (Host_Listening()) in = Host_TakeLocalClientMessages(c.conn);
    else { Net_Poll(); Net_Receive(in, 256); }
    for (const NetMessage& m : in) {
        if (m.conn != c.conn) continue;
        OnMessage(c, m);
        if (c.state == NetClient::State::Rejected || c.state == NetClient::State::Closed) return;
    }
    const double now = GetTime();
    if (c.state == NetClient::State::Connecting) {
        if (Net_IsConnected(c.conn)) {
            SendMsg(c, MsgHello{ c.protocol, c.build.empty() ? Commit() : c.build, c.name, c.password });
            c.state = NetClient::State::AwaitingWelcome;
            c.lastHeard = now;
        } else if (now - c.started > kClientConnectTimeout) { c.bye = NetBye::TimedOut; End(c, NetClient::State::Closed); }
        return;
    }
    if (c.state == NetClient::State::Joined && now - c.lastPing >= 1.0) { SendMsg(c, MsgPing{ ++c.pingSeq, now }); c.lastPing = now; }
    if (now - c.lastHeard > kClientSilentTimeout || !Net_IsConnected(c.conn)) { c.bye = NetBye::TimedOut; End(c, NetClient::State::Closed); }
}

void NetClient_Disconnect(NetClient& c) {
    if (!c.conn) return;
    SendMsg(c, MsgBye{ NetBye::Quit });
    c.bye = NetBye::Quit;
    End(c, NetClient::State::Closed);
}

// ---- Movement (#327) ----

MsgPlayerState NetClient_ToMsg(const Player::NetState& n) {
    MsgPlayerState s;
    s.body = n.body; s.local = n.local;
    s.vx = n.vel.x; s.vy = n.vel.y; s.vz = n.vel.z;
    s.jumpBuffer = n.jumpBuffer; s.coyote = n.coyote; s.jumpLock = n.jumpLock; s.jetpackFuel = n.jetpackFuel;
    s.flags = (n.grounded ? 1 : 0) | (n.hovering ? 16 : 0) | ((n.posture & 3) << 5);
    return s;
}
Player::NetState NetClient_FromMsg(const MsgPlayerState& s) {
    Player::NetState n;
    n.body = s.body; n.local = s.local;
    n.vel = { s.vx, s.vy, s.vz };
    n.grounded = (s.flags & 1) != 0; n.hovering = (s.flags & 16) != 0; n.posture = (uint8_t)((s.flags >> 5) & 3);
    n.jumpBuffer = s.jumpBuffer; n.coyote = s.coyote; n.jumpLock = s.jumpLock; n.jetpackFuel = s.jetpackFuel;
    return n;
}

void NetClient_Step(NetClient& c, Scene& world, Player& self, const PlayerInput& in) {
    if (c.state != NetClient::State::Joined) return;
    // Reconcile first: where were we, by our own reckoning, when the server applied that input?
    if (c.haveOwn) {
        c.haveOwn = false;
        Planetoid* body = Planetoid::findNet(c.own.body);
        while (!c.history.empty() && c.history.front().tick < c.own.ackInput) c.history.pop_front();
        if (body && !c.history.empty() && c.history.front().tick == c.own.ackInput) {
            const DVec3 d = c.history.front().local - c.own.local;
            c.lastError = d.length();
            c.worstError = std::max(c.worstError, c.lastError);
            c.history.pop_front(); // acknowledged
            if (c.lastError > kReconcileError) {
                ++c.corrections;
                // Back to the server's state, then replay what it hasn't applied yet, each step with the view it was
                // sent with (the view now is the newest; replaying with it turned every replayed step).
                const Vector3 viewLocal = body->dirToLocal(self.forward);
                self.applyNetState(world, NetClient_FromMsg(c.own));
                for (NetPredicted& h : c.history) {
                    self.forward = body->dirToWorld({ h.fx, h.fy, h.fz });
                    RenderOriginScope frame(world, self.position);
                    self.update(world, h.in, SimDt());
                    h.local = body->toLocalD(self.position);
                }
                self.forward = body->dirToWorld(viewLocal);
            }
        }
    }
    // Then this step: send it, and predict it.
    Planetoid* body = self.nearest;
    NetPredicted step;
    step.tick = ++c.tick; step.in = in;
    if (body) { const Vector3 f = body->dirToLocal(self.forward); step.fx = f.x; step.fy = f.y; step.fz = f.z; }
    {
        RenderOriginScope frame(world, self.position);
        self.update(world, in, SimDt());
    }
    if (body) step.local = body->toLocalD(self.position);
    c.history.push_back(step);
    while (c.history.size() > 600) c.history.pop_front(); // 10 s unacknowledged: something is wrong; don't grow forever
    // The newest kInputRedundancy steps, unreliable (a lost packet is covered by the next).
    MsgInputs m;
    for (auto it = c.history.rbegin(); it != c.history.rend() && (int)m.steps.size() < kInputRedundancy; ++it) {
        if (!m.steps.empty() && it->tick != m.steps.back().tick - 1) break; // (consecutive only)
        m.steps.push_back({ it->tick, it->in.held, it->in.pressed, it->fx, it->fy, it->fz });
    }
    NetWriter w; NetWrite(w, m);
    Net_Send(c.conn, w.bytes.data(), (uint32_t)w.bytes.size(), false);
}

bool NetClient_OtherView(const NetClient& c, uint32_t session, uint64_t& body, DVec3& local) {
    auto it = c.others.find(session);
    if (it == c.others.end() || it->second.empty()) return false;
    const auto& q = it->second;
    // kInterpDelay behind the newest state, in server ticks.
    const double t = c.viewTick > 0.0 ? c.viewTick : (double)q.back().tick - kInterpDelay * kSimHz;
    body = q.back().body;
    if (t <= q.front().tick) { local = q.front().local; return true; }
    for (size_t i = 1; i < q.size(); ++i) {
        if (q[i].tick < t) continue;
        const MsgPlayerState& a = q[i - 1];
        const MsgPlayerState& b = q[i];
        const double k = b.tick == a.tick ? 1.0 : (t - a.tick) / (double)(b.tick - a.tick);
        local = a.local + (b.local - a.local) * k;
        if (a.body != b.body) local = b.local; // (a new body: no blending across frames)
        return true;
    }
    local = q.back().local;
    return true;
}
