#ifndef NETCLIENT_H
#define NETCLIENT_H

// A client's side of the connection lifecycle (#321): connect by IP, the handshake (Hello -> Welcome or Reject), a ping
// every second (round trip; the server's tick), and the ways it ends: a Reject, a Bye from the server (it closed, a
// kick, a timeout), our own Bye, or nothing heard for 10 s. Several clients can live in one process (tests); a real
// game has one. Only the network code paths use it: single player never creates one.
//
// When a host runs in the same process (Host_Listening: tests), the host's poll receives every message and hands a
// client's to it (Host_TakeLocalClientMessages); otherwise the client polls the transport itself.

#include "netproto.h"
#include <cstdint>
#include <string>

struct NetClient {
    enum class State { Idle, Connecting, AwaitingWelcome, Joined, Rejected, Closed };
    State state = State::Idle;
    uint32_t conn = 0;
    std::string name, password, build;   // build: "" = the game's own (kGameVersion's commit)
    uint32_t protocol = kNetProtocol;    // tests: a different one to be refused
    MsgWelcome welcome;                  // when Joined
    MsgShower shower;                    // the latest from the server
    MsgReject reject;                    // when Rejected
    NetBye bye = NetBye::Quit;           // when Closed: why
    double rttMs = -1.0;                 // the last ping's round trip
    uint32_t serverTick = 0;             // from the last pong
    int pongs = 0;
    double started = 0.0, lastHeard = 0.0, lastPing = -1.0;
    uint32_t pingSeq = 0;
};

constexpr double kClientConnectTimeout = 10.0, kClientSilentTimeout = 10.0;

bool NetClient_Connect(NetClient& c, const char* address, const std::string& name, const std::string& password = "");
void NetClient_Update(NetClient& c);     // per frame: poll (unless a host here does), handle messages, ping, time out
void NetClient_Disconnect(NetClient& c); // a clean Bye(quit)
const char* NetClient_StateName(NetClient::State s);

#endif
