#include "net.h"

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

void Net_Poll() { if (initialised) SteamNetworkingSockets()->RunCallbacks(); }
bool Net_IsConnected(uint32_t conn) { return connected.count(conn) > 0; }
int Net_ConnectionCount() { return (int)connected.size(); }

bool Net_Send(uint32_t conn, const void* data, uint32_t size, bool reliable) {
    if (!initialised) return false;
    return SteamNetworkingSockets()->SendMessageToConnection(conn, data, size, reliable ? k_nSteamNetworkingSend_Reliable : k_nSteamNetworkingSend_Unreliable, nullptr) == k_EResultOK;
}

int Net_Receive(std::vector<NetMessage>& out, int max) {
    if (!initialised) return 0;
    std::vector<ISteamNetworkingMessage*> msgs((size_t)std::max(1, max));
    const int n = SteamNetworkingSockets()->ReceiveMessagesOnPollGroup(pollGroup, msgs.data(), (int)msgs.size());
    for (int i = 0; i < n; ++i) {
        NetMessage m;
        m.conn = msgs[i]->m_conn;
        m.data.assign((const uint8_t*)msgs[i]->m_pData, (const uint8_t*)msgs[i]->m_pData + msgs[i]->m_cbSize);
        out.push_back(std::move(m));
        msgs[i]->Release();
    }
    return std::max(0, n);
}

void Net_SetFakeConditions(int lagMs, float lossPercent) {
    ISteamNetworkingUtils* u = SteamNetworkingUtils();
    u->SetGlobalConfigValueInt32(k_ESteamNetworkingConfig_FakePacketLag_Send, lagMs);
    u->SetGlobalConfigValueFloat(k_ESteamNetworkingConfig_FakePacketLoss_Send, lossPercent);
}

void Net_Close(uint32_t conn) {
    if (!initialised) return;
    connected.erase(conn);
    SteamNetworkingSockets()->CloseConnection(conn, 0, "bye", true);
}

#else // built without ENABLE_GNS: stubs

bool Net_Available() { return false; }
bool Net_Init(std::string& error) { error = "built without ENABLE_GNS"; return false; }
void Net_Shutdown() {}
bool Net_Host(uint16_t) { return false; }
uint32_t Net_Connect(const char*) { return 0; }
void Net_Poll() {}
bool Net_IsConnected(uint32_t) { return false; }
int Net_ConnectionCount() { return 0; }
bool Net_Send(uint32_t, const void*, uint32_t, bool) { return false; }
int Net_Receive(std::vector<NetMessage>&, int) { return 0; }
void Net_Close(uint32_t) {}
void Net_SetFakeConditions(int, float) {}

#endif
