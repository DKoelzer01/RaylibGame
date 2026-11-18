#include "hostserver.h"
#include "events.h"
#include "joininfo.h"
#include "net.h"
#include "netproto.h"
#include "scene.h"
#include "simclock.h"
#include "survival.h"
#include "world.h"

#include <raylib.h>
#include <algorithm>
#include <cstdio>
#include <set>
#include <string>
#include <unordered_map>

extern const char* kGameVersion; // "<branch> <commit>" (cmake/game_version.cmake, #501)

namespace {
    bool listening = false;
    int joins = 0;
    std::set<uint32_t> localClients;                   // tests: client connections opened in this process
    std::unordered_map<uint32_t, std::vector<NetMessage>> localInbox;
    double autosaveSeconds = 0.0, lastSave = -1.0;
    int autosaves = 0;
    std::string password;
    int maxPlayers = 8;

    // The connection lifecycle (#321). A peer is a connection; a session is a player, kept by name for the server's
    // run, so a reconnect under the same name gets its session (and later its PlayerSession) back.
    struct Peer { uint32_t session = 0; bool joined = false; double lastHeard = 0.0; };
    std::unordered_map<uint32_t, Peer> peers;          // by connection
    struct Session { uint32_t id = 0; std::string name; uint32_t conn = 0; bool online = false; };
    std::vector<Session> sessions;
    uint32_t nextSession = 1;
    constexpr double kTimeoutSeconds = 15.0;           // nothing heard (pings come every second): gone

    void Send(uint32_t conn, const NetWriter& w) { Net_Send(conn, w.bytes.data(), (uint32_t)w.bytes.size(), true); }
    template <class M> void SendMsg(uint32_t conn, const M& m) { NetWriter w; NetWrite(w, m); Send(conn, w); }
    std::string BuildCommit() { const std::string v = kGameVersion; const size_t sp = v.rfind(' '); return sp == std::string::npos ? v : v.substr(sp + 1); }

    void Drop(uint32_t conn, NetBye why) {
        auto it = peers.find(conn);
        if (it == peers.end()) return;
        for (Session& s : sessions) if (s.id == it->second.session && s.conn == conn) { s.online = false; s.conn = 0; }
        if (why != NetBye::Quit) SendMsg(conn, MsgBye{ why });
        printf("host: connection %u (session %u) %s\n", conn, it->second.session, NetByeText(why));
        peers.erase(it);
        Net_Close(conn);
    }
    void Refuse(uint32_t conn, NetReject why, const std::string& detail) {
        SendMsg(conn, MsgReject{ why, detail });
        printf("host: refused connection %u: %s%s%s\n", conn, NetRejectText(why), detail.empty() ? "" : ": ", detail.c_str());
        peers.erase(conn);
        Net_Close(conn); // (lingers: the reject goes out first)
    }

    void OnHello(Scene& world, uint32_t conn, const MsgHello& h) {
        if (h.protocol != kNetProtocol) { Refuse(conn, NetReject::Protocol, TextFormat("server %u, client %u", kNetProtocol, h.protocol)); return; }
        if (!h.build.empty() && h.build != BuildCommit()) { Refuse(conn, NetReject::Build, "server " + BuildCommit() + ", client " + h.build); return; }
        if (!password.empty() && h.password != password) { Refuse(conn, NetReject::Password, ""); return; }
        Session* s = nullptr;
        for (Session& x : sessions) if (x.name == h.name) s = &x;
        if (s && s->online) { Refuse(conn, NetReject::NameTaken, h.name); return; }
        const int online = (int)std::count_if(sessions.begin(), sessions.end(), [](const Session& x) { return x.online; });
        if (online >= maxPlayers) { Refuse(conn, NetReject::Full, TextFormat("%d players", maxPlayers)); return; }
        const bool reconnected = s != nullptr;
        if (!s) { sessions.push_back({ nextSession++, h.name, 0, false }); s = &sessions.back(); }
        s->conn = conn; s->online = true;
        Peer& p = peers[conn];
        p.session = s->id; p.joined = true; p.lastHeard = GetTime();
        const JoinInfo info = JoinInfo_FromWorld(world);
        MsgWelcome w;
        w.session = s->id; w.reconnected = reconnected; w.seed = info.seed; w.bodies = info.bodies; w.fingerprint = info.fingerprint;
        w.serverTick = (uint32_t)SimTick();
        w.settings.difficulty = (uint8_t)std::clamp(worldDifficulty, 0, 2); w.settings.maxPlayers = (uint8_t)maxPlayers;
        SendMsg(conn, w);
        const EventsShower sh = Events_GetShower();
        SendMsg(conn, MsgShower{ sh.start, sh.end, sh.id });
        ++joins;
        printf("host: %s joined as session %u%s (connection %u, %d joins so far)\n", h.name.c_str(), s->id, reconnected ? " (reconnected)" : "", conn, joins);
    }

    void OnMessage(Scene& world, const NetMessage& m) {
        auto it = peers.find(m.conn);
        if (it == peers.end()) it = peers.emplace(m.conn, Peer{}).first; // first word on a new connection
        it->second.lastHeard = GetTime();
        NetReader r(m.data.data() + (m.data.empty() ? 0 : 1), m.data.empty() ? 0 : m.data.size() - 1);
        switch (NetPeek(m.data)) {
            case NetMsg::Hello: {
                MsgHello h;
                if (it->second.joined) return; // already in
                if (!NetRead(r, h)) { Refuse(m.conn, NetReject::Malformed, ""); return; }
                OnHello(world, m.conn, h);
                return;
            }
            case NetMsg::Ping: {
                MsgPing p;
                if (!it->second.joined || !NetRead(r, p)) return;
                SendMsg(m.conn, MsgPong{ p.seq, p.clientTime, (uint32_t)SimTick() });
                return;
            }
            case NetMsg::Bye: { MsgBye b; if (NetRead(r, b)) Drop(m.conn, NetBye::Quit); return; }
            default: return; // not joined, or not a message for the server yet: dropped
        }
    }
}

bool Host_Start(uint16_t port, std::string& status) {
    gNetRole = NetRole::Host; // the authority either way: the world runs and rolls here
    lastSave = GetTime();
    if (!Net_Available()) { status = "built without ENABLE_GNS: the world runs, nobody can connect"; return false; }
    std::string err;
    if (!Net_Init(err)) { status = "network init failed: " + err; return false; }
    if (!Net_Host(port)) { status = "could not listen on port " + std::to_string(port); return false; }
    listening = true;
    status = "listening on port " + std::to_string(port);
    return true;
}

void Host_SetAutosaveMinutes(double minutes) { autosaveSeconds = minutes > 0.0 ? minutes * 60.0 : 0.0; }
void Host_SetPassword(const std::string& p) { password = p; }
void Host_SetMaxPlayers(int n) { maxPlayers = std::clamp(n, 1, 64); }

void Host_Update(Scene& world) {
    // Autosave (user, 2026-11-15: every 30 minutes by default; --autosave-minutes N, 0 = off). Wall-clock time: a
    // server's saves protect against crashes, whatever the game clock does.
    if (autosaveSeconds > 0.0 && lastSave >= 0.0 && GetTime() - lastSave >= autosaveSeconds) {
        lastSave = GetTime();
        const bool ok = SaveWorld(world);
        if (ok) ++autosaves;
        printf("host: autosave %d to %s %s\n", autosaves, WorldSavePath(), ok ? "done" : "FAILED");
    }
    if (!listening) return;
    Net_Poll();
    std::vector<NetMessage> in;
    Net_Receive(in, 256);
    for (NetMessage& m : in) {
        if (localClients.count(m.conn)) { localInbox[m.conn].push_back(std::move(m)); continue; }
        OnMessage(world, m);
    }
    // Gone: closed under us (the transport noticed), or silent past the timeout.
    std::vector<uint32_t> gone, silent;
    for (const auto& [conn, p] : peers) {
        if (!Net_IsConnected(conn)) gone.push_back(conn);
        else if (GetTime() - p.lastHeard > kTimeoutSeconds) silent.push_back(conn);
    }
    for (uint32_t c : gone) Drop(c, NetBye::Quit);
    for (uint32_t c : silent) Drop(c, NetBye::TimedOut);
}

bool Host_SaveOnShutdown(Scene& world) {
    const bool ok = SaveWorld(world);
    printf("host: saved on shutdown to %s %s\n", WorldSavePath(), ok ? "done" : "FAILED");
    return ok;
}

void Host_Stop() {
    if (listening) {
        // Everyone hears why (#321): the server closed (each close lingers until the bye is out).
        std::vector<uint32_t> conns;
        for (const auto& [conn, p] : peers) conns.push_back(conn);
        for (uint32_t c : conns) Drop(c, NetBye::ServerClosed);
        Net_Shutdown();
    }
    listening = false;
    localClients.clear();
    localInbox.clear();
    peers.clear();
    sessions.clear();
    nextSession = 1;
    gNetRole = NetRole::Offline;
}

bool Host_Listening() { return listening; }
int Host_Joins() { return joins; }
int Host_Autosaves() { return autosaves; }
int Host_OnlinePlayers() { return (int)std::count_if(sessions.begin(), sessions.end(), [](const Session& s) { return s.online; }); }
void Host_SetLocalClient(uint32_t conn) { localClients.insert(conn); }
std::vector<NetMessage> Host_TakeLocalClientMessages(uint32_t conn) {
    std::vector<NetMessage> out;
    auto it = localInbox.find(conn);
    if (it != localInbox.end()) out.swap(it->second);
    return out;
}
