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
#include "net.h"
#include "playerinput.h"
#include "skills.h"
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
struct NetPredicted { // fx..: the view sent with it; fly..: a creative flyer's state after the step, sent too (#1922)
    uint32_t tick = 0; PlayerInput in; DVec3 local{}; float fx = 0, fy = 0, fz = 1;
    bool fly = false; uint64_t flyBody = 0; DVec3 flyLocal{}; Vector3 flyVel{};
};

struct NetClient {
    enum class State { Idle, Connecting, AwaitingWelcome, Joined, Rejected, Closed };
    State state = State::Idle;
    uint32_t conn = 0;
    std::string name, password, build;   // build: "" = the game's own (kGameVersion's commit)
    std::string uuid;                    // (#346) the player's id ("" a test client: one made from the name, the same every run)
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
    int creatureFxMsgs = 0;                      // (#1343, #1364) creature shot and event messages heard
    int damageTaken = 0; float lastDamage = 0.0f; bool lastDamageBlocked = false; // (#334) hits on our player from the server
    MsgVitals vitals; int vitalsMsgs = 0;        // (#337) our player's survival as the server last said
    MsgTeamState team; int teamMsgs = 0;      // (#604) our team as the server last said (members, health, an invite waiting)
    int gravesMsgs = 0, lastGraves = 0;          // (#804) Graves lists heard, and the size of the last
    int questBooks = 0, questTeam = 0;          // (#338) team quest books heard, and our team
    int teamStores = 0, teamStoreTeam = 0;      // (#1886) team progression stores heard, and the team of the last one
    std::string teamStoreText;                  // ...and its text
    int mapPinsSent = 0;                        // (#1886) map pins added or removed for the team (NetClient_SendMapPin)
    int machineViews = 0, machineStates = 0;    // (#332) machine views and working lists heard
    int plantLists = 0, plantsRemoved = 0, farmLists = 0; // (#335) plant removals and farm field lists heard
    int claimLists = 0;                          // (#343) claims lists heard
    SkillState skills;                           // (#920) this player's skills as the server keeps them (main.cpp applies gJoin's)
    int skillStates = 0;                         // (#920) SkillState messages heard
    double worstError = 0.0, lastError = 0.0;    // m: the prediction against the server at the acknowledged input
    double respawnGraceUntil = -1.0;             // (#828) dead, or just respawned: the server moved us (not a miss); errors are not counted
    std::map<uint32_t, std::deque<MsgPlayerState>> others; // other sessions' recent states (newest last)
    std::map<uint32_t, double> viewTicks;        // per other player: the step of theirs drawn now (NetClient_Update)
    // Who is playing and what they hold (#750)
    std::map<uint32_t, std::string> roster;      // everyone online by session (us too), from the server's roster
    std::map<uint32_t, MsgLoadout> loadouts;     // the others' held item, headlamp, suit tier
    int rosters = 0;                             // rosters received
    MsgLoadout loadout; bool loadoutSent = false; // ours, as last sent (NetClient_SetLoadout sends a change)
    uint32_t suitColor = 0;                      // 0xRRGGBB picked when joining (#341; 0: the default suit), sent with the loadout
    int action = 0; float actionPhase = 0.0f;    // our hands now (the game sets them; sent with every input packet)
    // Terrain edits (#329)
    class Inventory* inventory = nullptr;        // follows the server's copy (InventoryDelta); the game's: the local one
    std::vector<std::pair<std::string, double>> pendingInventory; // (#803) deltas heard before `inventory` was attached, applied then
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
    // Join sync (#776): from the Welcome until SyncDone the world's edits stream in; live edits wait (applied after, in
    // order). Entries whose hash matches what we hold (a cached client save) are not fetched.
    bool syncing = false;
    bool syncAssumeEmpty = false;                // tests: ask for everything, as a client with nothing cached
    // A joining game holds the sync's messages until its world is generated (they arrive with the Welcome, before the
    // bodies exist: installing then found no body); clear it when ready and they are handled, in order.
    bool syncHold = false;
    std::vector<NetMessage> syncHeld;
    uint32_t syncSnapshot = 0;
    int syncPages = 0, syncPagesGot = 0;
    std::vector<NetSyncEntry> syncIndex;
    std::vector<uint8_t> syncWanted;             // per entry: 1 asked for and not complete yet
    std::map<uint32_t, std::vector<uint8_t>> syncBuf;
    std::map<uint32_t, uint32_t> syncGot;        // bytes received per entry
    std::vector<NetMessage> syncQueue;           // live edits held during the sync
    int syncEntriesWanted = 0, syncEntriesDone = 0, syncHashBad = 0, syncInstalled = 0, syncCleared = 0, syncQueued = 0;
    uint64_t syncBytes = 0, syncBytesWanted = 0;
    double syncStartedAt = 0.0, syncDoneAt = -1.0;
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
// Join sync (#776): a body's plots as bytes (the save's text), and back onto a body (replacing its plots); the server
// sends the first, a client installs with the second.
std::vector<uint8_t> Plots_Blob(const class PlotSet& plots);
bool Plots_Install(class Planetoid& body, const std::vector<uint8_t>& blob);
// The terrain's fingerprint (every body's region hashes, combined): equal on a server and a synced client.
uint64_t NetSync_TerrainFingerprint(class Scene& world);
// The client's cache of the server's terrain edits, beside its client save (`save`: its path): loaded at a join when it
// is of the same world (seed, bodies fingerprint), saved when a sync is done. Return the regions loaded / written.
int NetClient_LoadCache(class Scene& world, const std::string& save, uint32_t seed, uint64_t fingerprint);
int NetClient_SaveCache(class Scene& world, const std::string& save, uint32_t seed, uint64_t fingerprint);
// Tests (#776): scatter `strokes` dig strokes (1.5-2.5 m) over 600 m around a planet's start, seeded: a world with many
// edits to sync. Returns the strokes applied (0: no planet).
int NetSync_ScatterStrokes(class Scene& world, int strokes, uint32_t seed);
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
// (#1886) a map pin placed (add) or removed, for the team's shared pins; false when not joined (it stays a local waypoint)
bool NetClient_SendMapPin(NetClient& c, bool add, const std::string& body, Vector3 dir, const std::string& name);
// (#1886) a research queue change (op: ResearchQueueOp 1 add, 2 remove, 3 up); the host's next team store shows the result
bool NetClient_SendResearchQueue(NetClient& c, uint8_t op, const std::string& node);
void NetClient_HitAnimal(NetClient& c, uint64_t netId);         // (#334) a swing at an animal: the server checks and applies it
void NetClient_UseAnimal(NetClient& c, uint64_t netId, int action); // (#334) 1 feed, 2 take its product (the server checks and applies)
void NetClient_UseItem(NetClient& c, const std::string& item);       // (#337) eat or use an item: the server applies it to its copy
void NetClient_UseBed(NetClient& c, const std::string& body, uint32_t plot, const int cell[3]); // (#800) make that bed my spawn
void NetClient_Craft(NetClient& c, int index, const std::string& output, int count, bool useChests = false, bool recursive = false); // (#338) craft on the server
void NetClient_SkillSpend(NetClient& c, const std::string& perk);                   // (#920) a perk point, as the server checks it
void NetClient_DropItems(NetClient& c, const std::string& item, int64_t centi);        // (#865) a pile at the feet, made by the server
void NetClient_CreativeGive(NetClient& c, const std::string& item, int count);         // (#866) creative items ("" clears the inventory)
void NetClient_MachineAction(NetClient& c, uint64_t bodyNetId, uint32_t plot, int part, int action, const std::string& item, int count, int slot, int which); // (#332)
void NetClient_PlantAction(NetClient& c, int action, uint64_t bodyNetId, uint64_t plant, Vector3 ground, const std::string& seed); // (#335) 1 fell, 2 sow, 3 harvest (#1085: at the body-local ground)
const char* NetClient_StateName(NetClient::State s);

#endif
