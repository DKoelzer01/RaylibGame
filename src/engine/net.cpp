#include "net.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <unordered_map>

// Simulated conditions (#347): shared by both builds (the stub build simply never sends).
namespace {
struct Held { double due = 0.0; uint64_t order = 0; uint32_t conn = 0; bool reliable = false; std::vector<uint8_t> data; };
NetConditions defaultCond;
std::unordered_map<uint32_t, NetConditions> linkCond;
std::unordered_map<uint32_t, double> lastReliableDue;   // keeps reliable messages in order per connection
std::deque<Held> held;                                   // sorted by (due, order)
NetLinkStats linkStats;
std::unordered_map<uint32_t, NetTraffic> traffic; // (#348)
uint64_t heldOrder = 0, rng = 0x9E3779B97F4A7C15ull;
bool envRead = false;

double NowMs() { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
double Rand01() { // splitmix64
    uint64_t z = (rng += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; z ^= z >> 31;
    return (double)(z >> 11) * (1.0 / 9007199254740992.0);
}
const NetConditions& CondFor(uint32_t conn) { auto it = linkCond.find(conn); return it == linkCond.end() ? defaultCond : it->second; }
void ReadEnv() {
    if (envRead) return;
    envRead = true;
    auto f = [](const char* n, float& v) { if (const char* s = getenv(n)) v = (float)atof(s); };
    NetConditions c = defaultCond;
    f("RAYLIBGAME_NET_LAG_MS", c.lagMs); f("RAYLIBGAME_NET_JITTER_MS", c.jitterMs);
    f("RAYLIBGAME_NET_LOSS", c.lossPct); f("RAYLIBGAME_NET_REORDER", c.reorderPct);
    if (const char* s = getenv("RAYLIBGAME_NET_SEED")) c.seed = (uint32_t)strtoul(s, nullptr, 10);
    if (c.active()) Net_SetConditions(c);
}
bool SendNow(uint32_t conn, const void* data, uint32_t size, bool reliable); // the real transport (below)
void FlushHeld() {
    const double now = NowMs();
    while (!held.empty() && held.front().due <= now) {
        Held h = std::move(held.front());
        held.pop_front();
        SendNow(h.conn, h.data.data(), (uint32_t)h.data.size(), h.reliable);
    }
}
}

void Net_SetConditions(const NetConditions& c, uint32_t conn) {
    envRead = true; // explicit settings win over the environment
    if (conn == 0) { defaultCond = c; rng = 0x9E3779B97F4A7C15ull ^ ((uint64_t)c.seed * 0xD1B54A32D192ED03ull); }
    else linkCond[conn] = c;
}
NetConditions Net_GetConditions(uint32_t conn) { return conn ? CondFor(conn) : defaultCond; }
NetLinkStats Net_GetLinkStats() { return linkStats; }
int Net_HeldMessages() { return (int)held.size(); }
std::vector<std::pair<uint32_t, NetTraffic>> Net_Traffic() { return { traffic.begin(), traffic.end() }; }
void Net_ResetTraffic() { traffic.clear(); }
// #327's test call: the default conditions, through the seeded layer (it was GNS's global fake lag and loss: not seedable).
void Net_SetFakeConditions(int lagMs, float lossPercent) { NetConditions c = defaultCond; c.lagMs = (float)lagMs; c.lossPct = lossPercent; Net_SetConditions(c); }

// Applies the conditions; true when the message was handled (held or dropped) instead of sent now.
[[maybe_unused]] static bool Net_Condition(uint32_t conn, const void* data, uint32_t size, bool reliable) {
    ReadEnv();
    const NetConditions& c = CondFor(conn);
    if (!c.active()) return false;
    ++linkStats.sent;
    const bool lost = Rand01() * 100.0 < c.lossPct;
    if (lost && !reliable) { ++linkStats.dropped; return true; }
    double delay = c.lagMs + (c.jitterMs > 0.0f ? (Rand01() * 2.0 - 1.0) * c.jitterMs : 0.0);
    if (lost) { delay += 2.0 * c.lagMs + 20.0; ++linkStats.retransmits; }                  // resent after a round trip
    if (Rand01() * 100.0 < c.reorderPct) { delay += c.lagMs * (1.0 + Rand01()) + 5.0; ++linkStats.reordered; }
    double due = NowMs() + std::max(0.0, delay);
    if (reliable) { double& last = lastReliableDue[conn]; due = std::max(due, last); last = due; }
    if (due > NowMs()) ++linkStats.delayed;
    Held h{ due, heldOrder++, conn, reliable, std::vector<uint8_t>((const uint8_t*)data, (const uint8_t*)data + size) };
    auto it = held.end();
    while (it != held.begin() && (it - 1)->due > due) --it;   // most land at the back
    held.insert(it, std::move(h));
    FlushHeld();
    return true;
}

#ifdef RAYLIBGAME_GNS

#include <steam/steamnetworkingsockets.h>
#include <steam/isteamnetworkingutils.h>
#include <algorithm>
#include <unordered_set>

namespace {
bool initialised = false;
HSteamListenSocket listenSocket = k_HSteamListenSocket_Invalid;
HSteamNetPollGroup pollGroup = k_HSteamNetPollGroup_Invalid;
std::unordered_set<uint32_t> connected;

void OnStatus(SteamNetConnectionStatusChangedCallback_t* info) {
    ISteamNetworkingSockets* s = SteamNetworkingSockets();
    switch (info->m_info.m_eState) {
    case k_ESteamNetworkingConnectionState_Connecting:
        // Incoming (it arrived on our listen socket): accept it into the poll group.
        if (info->m_info.m_hListenSocket != k_HSteamListenSocket_Invalid && s->AcceptConnection(info->m_hConn) == k_EResultOK)
            s->SetConnectionPollGroup(info->m_hConn, pollGroup);
        break;
    case k_ESteamNetworkingConnectionState_Connected:
        connected.insert(info->m_hConn);
        break;
    case k_ESteamNetworkingConnectionState_ClosedByPeer:
    case k_ESteamNetworkingConnectionState_ProblemDetectedLocally:
        connected.erase(info->m_hConn);
        s->CloseConnection(info->m_hConn, 0, nullptr, false);
        break;
    default: break;
    }
}

SteamNetworkingConfigValue_t StatusOption() {
    SteamNetworkingConfigValue_t opt;
    opt.SetPtr(k_ESteamNetworkingConfig_Callback_ConnectionStatusChanged, (void*)OnStatus);
    return opt;
}
}

bool Net_Available() { return true; }

bool Net_Init(std::string& error) {
    if (initialised) return true;
    SteamDatagramErrMsg msg;
    if (!GameNetworkingSockets_Init(nullptr, msg)) { error = msg; return false; }
    ReadEnv(); // simulated conditions from the environment (#347)
    // A bigger send buffer (#776): the default 512 KB refused reliable sends outright when it filled (a join sync streams
    // megabytes), and they were lost. Senders that stream also hold back on Net_PendingBytes.
    SteamNetworkingUtils()->SetGlobalConfigValueInt32(k_ESteamNetworkingConfig_SendBufferSize, 8 * 1024 * 1024);
    pollGroup = SteamNetworkingSockets()->CreatePollGroup();
    initialised = true;
    return true;
}

void Net_Shutdown() {
    if (!initialised) return;
    ISteamNetworkingSockets* s = SteamNetworkingSockets();
    for (uint32_t c : connected) s->CloseConnection(c, 0, "shutdown", false);
    connected.clear();
    if (listenSocket != k_HSteamListenSocket_Invalid) s->CloseListenSocket(listenSocket);
    listenSocket = k_HSteamListenSocket_Invalid;
    s->DestroyPollGroup(pollGroup);
    pollGroup = k_HSteamNetPollGroup_Invalid;
    GameNetworkingSockets_Kill();
    initialised = false;
}

bool Net_Host(uint16_t port) {
    if (!initialised || listenSocket != k_HSteamListenSocket_Invalid) return false;
    SteamNetworkingIPAddr addr;
    addr.Clear();
    addr.m_port = port;
    const SteamNetworkingConfigValue_t opt = StatusOption();
    listenSocket = SteamNetworkingSockets()->CreateListenSocketIP(addr, 1, &opt);
    return listenSocket != k_HSteamListenSocket_Invalid;
}

uint32_t Net_Connect(const char* address) {
    if (!initialised) return 0;
    SteamNetworkingIPAddr addr;
    addr.Clear();
    if (!addr.ParseString(address)) return 0;
    const SteamNetworkingConfigValue_t opt = StatusOption();
    const HSteamNetConnection c = SteamNetworkingSockets()->ConnectByIPAddress(addr, 1, &opt);
    if (c != k_HSteamNetConnection_Invalid) SteamNetworkingSockets()->SetConnectionPollGroup(c, pollGroup);
    return c;
}

void Net_Poll() { if (initialised) { FlushHeld(); SteamNetworkingSockets()->RunCallbacks(); } }
std::string Net_RemoteAddress(uint32_t conn) {
    SteamNetConnectionInfo_t info;
    if (!initialised || !SteamNetworkingSockets()->GetConnectionInfo(conn, &info)) return "";
    char buf[SteamNetworkingIPAddr::k_cchMaxString]; info.m_addrRemote.ToString(buf, sizeof buf, true);
    return buf;
}
bool Net_IsConnected(uint32_t conn) { return connected.count(conn) > 0; }
int Net_ConnectionCount() { return (int)connected.size(); }

namespace {
int reliableFailures = 0;
bool SendNow(uint32_t conn, const void* data, uint32_t size, bool reliable) {
    if (!initialised) return false;
    const bool ok = SteamNetworkingSockets()->SendMessageToConnection(conn, data, size, reliable ? k_nSteamNetworkingSend_Reliable : k_nSteamNetworkingSend_Unreliable, nullptr) == k_EResultOK;
    if (!ok && reliable && ++reliableFailures <= 5) printf("net: a reliable send of %u bytes on connection %u was refused (the send buffer is full?)\n", size, conn);
    return ok;
}
}

int Net_ReliableSendFailures() { return reliableFailures; }
uint32_t Net_PendingBytes(uint32_t conn) {
    uint32_t pending = 0; // (the simulated-conditions queue counts too)
    for (const Held& h : held) if (h.conn == conn) pending += (uint32_t)h.data.size();
    if (!initialised) return pending;
    SteamNetConnectionRealTimeStatus_t st;
    if (SteamNetworkingSockets()->GetConnectionRealTimeStatus(conn, &st, 0, nullptr) != k_EResultOK) return pending;
    return pending + (uint32_t)std::max(0, st.m_cbPendingReliable + st.m_cbPendingUnreliable + st.m_cbSentUnackedReliable);
}

bool Net_Send(uint32_t conn, const void* data, uint32_t size, bool reliable) {
    if (!initialised) return false;
    { NetTraffic& t = traffic[conn]; t.sentBytes += size; ++t.sentMsgs; } // (#348)
    if (Net_Condition(conn, data, size, reliable)) return true; // held or dropped (#347)
    return SendNow(conn, data, size, reliable);
}

int Net_Receive(std::vector<NetMessage>& out, int max) {
    if (!initialised) return 0;
    std::vector<ISteamNetworkingMessage*> msgs((size_t)std::max(1, max));
    const int n = SteamNetworkingSockets()->ReceiveMessagesOnPollGroup(pollGroup, msgs.data(), (int)msgs.size());
    for (int i = 0; i < n; ++i) {
        NetMessage m;
        m.conn = msgs[i]->m_conn;
        { NetTraffic& t = traffic[m.conn]; t.recvBytes += (uint64_t)msgs[i]->m_cbSize; ++t.recvMsgs; } // (#348)
        m.data.assign((const uint8_t*)msgs[i]->m_pData, (const uint8_t*)msgs[i]->m_pData + msgs[i]->m_cbSize);
        out.push_back(std::move(m));
        msgs[i]->Release();
    }
    return std::max(0, n);
}


void Net_Close(uint32_t conn) {
    if (!initialised) return;
    for (auto it = held.begin(); it != held.end();) { // what is still held goes out now, ahead of the close (#347)
        if (it->conn == conn) { SendNow(conn, it->data.data(), (uint32_t)it->data.size(), it->reliable); it = held.erase(it); } else ++it;
    }
    linkCond.erase(conn); lastReliableDue.erase(conn);
    connected.erase(conn);
    SteamNetworkingSockets()->CloseConnection(conn, 0, "bye", true);
}

#else // built without ENABLE_GNS: stubs

namespace { bool SendNow(uint32_t, const void*, uint32_t, bool) { return false; } }
bool Net_Available() { return false; }
bool Net_Init(std::string& error) { error = "built without ENABLE_GNS"; return false; }
void Net_Shutdown() {}
bool Net_Host(uint16_t) { return false; }
uint32_t Net_Connect(const char*) { return 0; }
void Net_Poll() {}
int Net_ReliableSendFailures() { return 0; }
uint32_t Net_PendingBytes(uint32_t) { return 0; }
bool Net_IsConnected(uint32_t) { return false; }
int Net_ConnectionCount() { return 0; }
bool Net_Send(uint32_t, const void*, uint32_t, bool) { return false; }
int Net_Receive(std::vector<NetMessage>&, int) { return 0; }
void Net_Close(uint32_t) {}
std::string Net_RemoteAddress(uint32_t) { return ""; }

#endif
