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
#include "playerinput.h"
#include <cstdint>
#include <deque>
#include <map>
#include <string>

class Scene;
#include "player.h"
#include "avatar.h"
#include <vector>

// Movement (#327): the client predicts its own player from its inputs (the same Player physics as the server) and
// keeps them; each authoritative state of its own player says which input the server applied last. Where the
// prediction then was more than kReconcileError off, the player is put back to the server's state and the inputs since
// replayed. Other players are drawn kInterpDelay behind the newest state, interpolated between the two around it.
constexpr double kReconcileError = 0.05;   // m
constexpr double kInterpDelay = 0.1;       // s
constexpr double kExtrapolate = 0.25;     // s: how far another player is carried on past their newest state
struct NetPredicted { uint32_t tick = 0; PlayerInput in; DVec3 local{}; float fx = 0, fy = 0, fz = 1; }; // fx..: the view sent with it

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
    // Movement (#327)
    uint32_t tick = 0;                           // this client's sim steps since joining (the input ticks)
    std::deque<NetPredicted> history;            // inputs not yet acknowledged, with where they left the player
    bool haveOwn = false; MsgPlayerState own;    // the newest state of our player, not yet reconciled
    int states = 0, corrections = 0;             // own states received; how many moved the player back
    double worstError = 0.0, lastError = 0.0;    // m: the prediction against the server at the acknowledged input
    std::map<uint32_t, std::deque<MsgPlayerState>> others; // other sessions' recent states (newest last)
    std::map<uint32_t, double> viewTicks;        // per other player: the step of theirs drawn now (NetClient_Update)
    // Who is playing and what they hold (#750)
    std::map<uint32_t, std::string> roster;      // everyone online by session (us too), from the server's roster
    std::map<uint32_t, MsgLoadout> loadouts;     // the others' held item, headlamp, suit tier
    int rosters = 0;                             // rosters received
    MsgLoadout loadout; bool loadoutSent = false; // ours, as last sent (NetClient_SetLoadout sends a change)
    int action = 0; float actionPhase = 0.0f;    // our hands now (the game sets them; sent with every input packet)
};

// One sim step of the joined client's own player: send the input, predict with it (in its own frame), and reconcile
// with the newest authoritative state first. `self` is the predicted player (the local player in a game).
void NetClient_Step(NetClient& c, Scene& world, Player& self, const PlayerInput& in);
// Another session's player as drawn: body-local feet kInterpDelay behind the newest state (false: none known yet).
bool NetClient_OtherView(const NetClient& c, uint32_t session, uint64_t& body, DVec3& local);
// The same, the whole state: position, velocity and view interpolated; flags, posture and hands from the later state.
bool NetClient_OtherState(const NetClient& c, uint32_t session, MsgPlayerState& out);
// Every other player as an avatar to draw (#328/#750): interpolated, named from the roster, with their loadout.
void NetClient_FillAvatars(const NetClient& c, Scene& world, std::vector<AvatarState>& out);
// Our held item, headlamp and suit tier; sent (reliable) only when they change.
void NetClient_SetLoadout(NetClient& c, const std::string& held, bool headlamp, int suitTier);
// The game's client (null: not connected to a server). main.cpp draws its other players.
extern NetClient* gNetClient;
// A player's state on the wire and back (#327; server and client use the same).
MsgPlayerState NetClient_ToMsg(const Player::NetState& n);
Player::NetState NetClient_FromMsg(const MsgPlayerState& s);

constexpr double kClientConnectTimeout = 10.0, kClientSilentTimeout = 10.0;

bool NetClient_Connect(NetClient& c, const char* address, const std::string& name, const std::string& password = "");
void NetClient_Update(NetClient& c);     // per frame: poll (unless a host here does), handle messages, ping, time out
void NetClient_Disconnect(NetClient& c); // a clean Bye(quit)
const char* NetClient_StateName(NetClient::State s);

#endif
