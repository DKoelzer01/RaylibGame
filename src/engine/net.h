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

#endif
