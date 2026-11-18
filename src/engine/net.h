#ifndef NET_H
#define NET_H

// Multiplayer transport (#310): a thin layer over Valve's GameNetworkingSockets (the same ISteamNetworkingSockets API
// the Steamworks SDK exposes, so Steam relay and lobbies can replace connect-by-IP later). Built only with the CMake
// option ENABLE_GNS (GNS_ROOT); otherwise every call is a stub and Net_Available() is false.
// One poll group holds every connection (accepted by a listen socket, or opened by Net_Connect), so a listen server
// and a client read their messages the same way.

#include <cstdint>
#include <string>
#include <vector>

// This process's part in a game (#316): Offline (single-player) and Host decide world events and roll the gameplay RNG;
// a Client applies what the host sends and rolls only cosmetics.
enum class NetRole { Offline, Host, Client };
inline NetRole gNetRole = NetRole::Offline;
inline bool Net_IsAuthority() { return gNetRole != NetRole::Client; }

struct NetMessage {
    uint32_t conn = 0;             // the connection it came on
    std::vector<uint8_t> data;
};

bool Net_Available();                                  // built with ENABLE_GNS
bool Net_Init(std::string& error);                     // once; false (with a reason) if the library won't start
void Net_Shutdown();
bool Net_Host(uint16_t port);                          // listen for connections on every local address
uint32_t Net_Connect(const char* address);             // "ip:port"; a connection handle, 0 on a bad address
void Net_Poll();                                       // run callbacks: accepts, connects, disconnects
bool Net_IsConnected(uint32_t conn);
int Net_ConnectionCount();                             // connected now
bool Net_Send(uint32_t conn, const void* data, uint32_t size, bool reliable);
int Net_Receive(std::vector<NetMessage>& out, int max = 64); // appends; returns how many
void Net_Close(uint32_t conn);
// Tests (#327, #347): every message sent from now on is delayed by lagMs and lost with lossPercent (0-100) odds: the
// default NetConditions below (seed kept), so a run replays.
void Net_SetFakeConditions(int lagMs, float lossPercent);


// Simulated network conditions (#347, tests only): latency, jitter, loss and reordering applied to what this process
// sends, per connection or as the default for every one. Deterministic (its own seeded RNG), so a run replays.
// Unreliable messages are dropped at `lossPct`; reliable ones are never dropped but a "lost" one arrives a round trip
// late (a retransmit) and they stay in order. Jitter (+-jitterMs) can reorder unreliable messages; `reorderPct` holds a
// message back by an extra 1-2 x lag. Held messages go out from Net_Poll / Net_Send once due.
// Env defaults (read by Net_Init): RAYLIBGAME_NET_LAG_MS, _JITTER_MS, _LOSS (percent), _REORDER (percent), _SEED.
struct NetConditions {
    float lagMs = 0.0f, jitterMs = 0.0f, lossPct = 0.0f, reorderPct = 0.0f;
    uint32_t seed = 1;
    bool active() const { return lagMs > 0.0f || jitterMs > 0.0f || lossPct > 0.0f || reorderPct > 0.0f; }
};
struct NetLinkStats { uint64_t sent = 0, dropped = 0, delayed = 0, reordered = 0, retransmits = 0; };
void Net_SetConditions(const NetConditions& c, uint32_t conn = 0); // conn 0: the default for every connection
NetConditions Net_GetConditions(uint32_t conn = 0);
NetLinkStats Net_GetLinkStats();                                      // totals over every connection
int Net_HeldMessages();                                               // waiting for their delivery time

#endif
