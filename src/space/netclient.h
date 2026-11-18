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
#include <unordered_map>
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
    std::deque<MsgChat> chat;            // (#326) chat lines heard, oldest first (the game shows and pops them; at most 64 kept)
    double clockError = 0.0, worstClockError = 0.0; // #339: the server's world time (pong + half the round trip) minus ours, s
    double rttFloor = -1.0;                      // #339: the best recent round trip (ms): pongs far above it (a resend, a stall) are skipped
    double started = 0.0, lastHeard = 0.0, lastPing = -1.0;
    uint32_t pingSeq = 0;
    // Movement (#327)
    uint32_t tick = 0;                           // this client's sim steps since joining (the input ticks)
    std::deque<NetPredicted> history;            // inputs not yet acknowledged, with where they left the player
    bool haveOwn = false; MsgPlayerState own;    // the newest state of our player, not yet reconciled
    int states = 0, corrections = 0;             // own states received; how many moved the player back
    // #769 snapshots: the newest received and the 32 before it (acked back in each MsgInputs), each player's last 32 states
    // received (bases for the deltas that follow), and counters.
    uint32_t snapNewest = 0, snapMask = 0;
    std::unordered_map<uint32_t, std::deque<std::pair<uint32_t, PlayerStateQ>>> snapHistory;
    int snapshots = 0, snapUndecodable = 0;
    int animalMsgs = 0;                          // (#334) animal snapshots heard
    int damageTaken = 0; float lastDamage = 0.0f; bool lastDamageBlocked = false; // (#334) hits on our player from the server
    double worstError = 0.0, lastError = 0.0;    // m: the prediction against the server at the acknowledged input
    std::map<uint32_t, std::deque<MsgPlayerState>> others; // other sessions' recent states (newest last)
    std::map<uint32_t, double> viewTicks;        // per other player: the step of theirs drawn now (NetClient_Update)
    // Who is playing and what they hold (#750)
    std::map<uint32_t, std::string> roster;      // everyone online by session (us too), from the server's roster
    std::map<uint32_t, MsgLoadout> loadouts;     // the others' held item, headlamp, suit tier
    int rosters = 0;                             // rosters received
    MsgLoadout loadout; bool loadoutSent = false; // ours, as last sent (NetClient_SetLoadout sends a change)
    int action = 0; float actionPhase = 0.0f;    // our hands now (the game sets them; sent with every input packet)
    // Terrain edits (#329)
    class Inventory* inventory = nullptr;        // follows the server's copy (InventoryDelta); the game's: the local one
    uint32_t editSeq = 0;                        // our requests so far
    int editsSeen = 0, editsOwn = 0, editsRefused = 0, inventoryDeltas = 0;
    int placesSeen = 0;                          // edits that placed ground (the soft brush)
    NetEditRefusal lastRefusal = NetEditRefusal::Count;
    uint32_t refusedWhy = 0;                     // a bit per NetEditRefusal seen
    uint64_t editOrder = 14695981039346656037ull; // a hash of every edit's (session, seq) in the order applied
    // Building (#330)
    uint32_t partSeq = 0;                        // our part requests so far
    int plotEdits = 0, plotEditsOwn = 0, partsRefused = 0;
    uint32_t partRefusedWhy = 0;                 // a bit per NetPartRefusal seen
    uint32_t lastNewPlot = 0;                    // the id of the last plot the server made for us
    bool keepPlotLog = false; std::vector<MsgPlotEdit> plotLog; // tests: every plot edit received, to replay elsewhere
    uint64_t plotOrder = 14695981039346656037ull; // a hash of every plot edit in the order applied
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
// Ask the server for a dig or place stroke (#329; the terrain tool's gTerrainToolNetEdit while joined). Returns its
// seq (0: not joined). The edit comes back as a TerrainEdit with it, or an EditRefused.
uint32_t NetClient_RequestEdit(NetClient& c, NetEditKind kind, const class Planetoid& body, Vector3 hit, Vector3 centre, float radius, const std::string& item);
// Ask the server to place a part (#330; the build tool's gBuildToolNetEdit while joined) into `plot` (null: a new plot
// anchored at `anchor`, body-local), or to break `plot`'s part partIndex. Returns its seq (0: not sent). The change comes
// back as PlotEdits with it, or a PartRefused.
uint32_t NetClient_RequestPart(NetClient& c, const class Planetoid& body, const class Plot* plot, Vector3 anchor, const struct Part& part, bool remove, int partIndex);
// The game's client (null: not connected to a server). main.cpp draws its other players.
extern NetClient* gNetClient;
// A player's state on the wire and back (#327; server and client use the same).
MsgPlayerState NetClient_ToMsg(const Player::NetState& n);
Player::NetState NetClient_FromMsg(const MsgPlayerState& s);

constexpr double kClientConnectTimeout = 10.0, kClientSilentTimeout = 10.0;

bool NetClient_Connect(NetClient& c, const char* address, const std::string& name, const std::string& password = "");
void NetClient_Update(NetClient& c);     // per frame: poll (unless a host here does), handle messages, ping, time out
void NetClient_Disconnect(NetClient& c); // a clean Bye(quit)
void NetClient_SendChat(NetClient& c, const std::string& text); // (#326) a chat line; "/command" lines go to ops' commands on the server
void NetClient_HitAnimal(NetClient& c, uint64_t netId);         // (#334) a swing at an animal: the server checks and applies it
void NetClient_UseAnimal(NetClient& c, uint64_t netId, int action); // (#334) 1 feed, 2 take its product (the server checks and applies)
const char* NetClient_StateName(NetClient::State s);

#endif
