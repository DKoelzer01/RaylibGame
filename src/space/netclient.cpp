#include "netclient.h"
#include "hostserver.h"
#include "net.h"

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
