#ifndef NETPROTO_H
#define NETPROTO_H
// The multiplayer wire protocol (#320, phase 2 of #303; decisions in docs/multiplayer.md and the ADR #309).
//
// Every message is one GameNetworkingSockets message: a byte with its id (NetMsg), then its fields, written by
// NetWriter and read back by NetReader. Integers are little-endian or LEB128 varints (zigzag for signed); positions
// are body-local doubles relative to a body's NetId (ADR: never render or world space), quantised to 1 mm; angles and
// unit vectors are quantised to 16 bits. NetReader never reads past the end: a short or garbled message sets ok to
// false and every later read returns 0, so a bad packet is dropped, never a crash (fuzz-checked by --test-netproto).
//
// kNetProtocol changes whenever any message's layout does; the handshake (Hello) refuses a mismatch with a clear reason.

#include "dvec3.h"
#include <cstdint>
#include <string>
#include <functional>
#include <vector>

constexpr uint32_t kNetProtocol = 14; // 1 was the text "join" of #311; 5: world time in Welcome and Pong (#339); 6 (unused); 7: world-gen settings in Welcome (#757); 8: terrain edit requests, refusals, inventory deltas (#329); 9: chat (#326); 10: building (#330); 11 is #776's; 12: player snapshots, delta-coded against acked states (#769); 13: animals (#334); 14: animal actions, damage to a player (#334)

enum class NetMsg : uint8_t {
    Hello = 1,     // client -> server: protocol, build, player name, password
    Welcome,       // server -> client: the session, the world to regenerate (JoinInfo), server settings
    Reject,        // server -> client: why not (then the connection closes)
    Bye,           // either way: a clean disconnect, with the reason
    Ping,          // client -> server: sequence and the client's time (keepalive, round trip)
    Pong,          // server -> client: the ping echoed, with the server's sim tick
    Input,         // client -> server: one sim step's PlayerInput (#327)
    PlayerState,   // server -> client: a player's authoritative state at a tick (#327)
    TerrainEdit,   // either way: a terrain edit (carve/fill, sphere/box) on a body (#320's first mapped command)
    Shower,        // server -> client: the meteor shower now (#316; was #311's "S" text line)
    Roster,        // server -> client: who is playing (session, name), sent on every join and leave (#750)
    Loadout,       // either way: a player's held item, headlamp, suit tier, sent when it changes (#750; to clients with the session)
    EditRequest,   // client -> server: dig or place at an aimed point (#329); the server validates, applies, pays
    EditRefused,   // server -> client: that request was refused (why)
    InventoryDelta, // server -> client: an item added to or taken from that player's inventory (the server's is the truth)
    Chat,          // either way (#326): a chat line; from the server (session 0: "say", replies to ops' commands) or a player
    PartRequest,   // client -> server: place a part (into a plot, or a new one) or break one (#330)
    PartRefused,   // server -> client: that request was refused (why)
    PlotEdit,      // server -> client: a plot made, a part placed or removed, in the server's order (#330)
    PlayerSnapshot = 24, // server -> client (#769): the player states of one snapshot in one message, each delta-coded against one
                   // the client acknowledged (20-23 are #776's)
    Animals = 25,  // server -> client (#334): the animals near the client (state, health), and the ones gone from its view
    AnimalHit = 26, // client -> server (#334): a swing at an animal, or feed / take its product (the server checks reach and cooldown)
    Damage = 27,   // server -> client (#334): a hit on the client's player (an animal's bite): it applies the damage (shield, health)
    Count
};
const char* NetMsgName(NetMsg m);

class NetWriter {
public:
    std::vector<uint8_t> bytes;
    void u8(uint8_t v) { bytes.push_back(v); }
    void u16(uint16_t v) { u8((uint8_t)v); u8((uint8_t)(v >> 8)); }
    void u32(uint32_t v) { for (int i = 0; i < 4; ++i) u8((uint8_t)(v >> (8 * i))); }
    void u64(uint64_t v) { for (int i = 0; i < 8; ++i) u8((uint8_t)(v >> (8 * i))); }
    void varu(uint64_t v) { while (v >= 0x80) { u8((uint8_t)(v | 0x80)); v >>= 7; } u8((uint8_t)v); }
    void vari(int64_t v) { varu(((uint64_t)v << 1) ^ (uint64_t)(v >> 63)); } // zigzag
    void f32(float v);
    void f64(double v);
    void str(const std::string& s) { varu(s.size()); bytes.insert(bytes.end(), s.begin(), s.end()); }
    void boolean(bool b) { u8(b ? 1 : 0); }
    void pos(const DVec3& p);           // body-local, metres, quantised to 1 mm (varint per axis)
    void unit(float x, float y, float z); // a unit vector, 3 x int16
    void angle(float radians);          // uint16 over a full turn
    void speed(float v);                // m/s, quantised to 1 cm/s, int32 varint
};

class NetReader {
public:
    NetReader(const uint8_t* data, size_t size) : p(data), end(data + size) {}
    explicit NetReader(const std::vector<uint8_t>& b) : NetReader(b.data(), b.size()) {}
    bool ok = true;                     // false once anything was missing or malformed (then every read gives 0)
    bool atEnd() const { return p == end; }
    uint8_t u8() { if (!need(1)) return 0; return *p++; }
    uint16_t u16() { const uint16_t a = u8(); return (uint16_t)(a | (u8() << 8)); }
    uint32_t u32() { uint32_t v = 0; for (int i = 0; i < 4; ++i) v |= (uint32_t)u8() << (8 * i); return v; }
    uint64_t u64() { uint64_t v = 0; for (int i = 0; i < 8; ++i) v |= (uint64_t)u8() << (8 * i); return v; }
    uint64_t varu();
    int64_t vari() { const uint64_t z = varu(); return (int64_t)(z >> 1) ^ -(int64_t)(z & 1); }
    float f32();
    double f64();
    std::string str(size_t maxLen = 4096);
    bool boolean() { const uint8_t b = u8(); if (b > 1) ok = false; return b == 1; }
    DVec3 pos();
    void unit(float& x, float& y, float& z);
    float angle();
    float speed();
private:
    bool need(size_t n) { if (!ok || (size_t)(end - p) < n) { ok = false; p = end; return false; } return true; }
    const uint8_t* p;
    const uint8_t* end;
};

// ---- Messages (each: its fields, Write appends id + fields, Read parses the bytes after the id) ----

struct MsgHello {
    uint32_t protocol = kNetProtocol;
    std::string build;    // the client's build (kGameVersion's commit): a different build may simulate differently
    std::string name;     // the player's name (IP joins); a reconnect under the same name gets its session back
    std::string password; // "" none
};
struct NetServerSettings { uint8_t difficulty = 1; bool pvp = false, claims = false; uint8_t maxPlayers = 8; };
struct MsgWelcome {
    uint32_t protocol = kNetProtocol;
    uint32_t session = 0;            // the player's session id on this server (stable across reconnects by name)
    bool reconnected = false;        // the session existed: position and inventory kept
    uint32_t seed = 0;               // JoinInfo: regenerate the world from these, then compare the fingerprint
    std::vector<std::string> bodies;
    uint64_t fingerprint = 0;
    uint32_t serverTick = 0;         // the server's sim tick when it sent this (the client's clock starts here)
    double worldTime = 0.0;          // the server's WorldTime() when it sent this (#339: spin, orbits, sun, weather)
    std::string worldgen;            // #757: the world's generation settings (WorldGenParams_ToText: the sidecar's lines)
    bool sandbox = false;            // #757: --sandbox (#662): the start planet is the flat test plane
    std::string motd;                // #326: the server's message of the day ("" none), shown on join
    NetServerSettings settings;
};
enum class NetReject : uint8_t { Protocol = 1, Build, Password, Full, Banned, NameTaken, Malformed };
const char* NetRejectText(NetReject r);
struct MsgReject { NetReject reason = NetReject::Malformed; std::string detail; };
enum class NetBye : uint8_t { Quit = 1, ServerClosed, TimedOut, Kicked };
const char* NetByeText(NetBye r);
struct MsgBye { NetBye reason = NetBye::Quit; };
struct MsgPing { uint32_t seq = 0; double clientTime = 0.0; };
struct MsgPong { uint32_t seq = 0; double clientTime = 0.0; uint32_t serverTick = 0; double worldTime = 0.0; }; // worldTime: #339
struct MsgInput {                    // one sim step of a client's PlayerInput (#327)
    uint32_t tick = 0;               // the client's sim tick it was sampled for
    uint64_t held = 0, pressed = 0;  // Act bits
    float fx = 0.0f, fy = 0.0f, fz = 1.0f; // the view direction after this step's look, body-local (unit; the frame turns with the body)
};
// The input packet (NetMsg::Input): the newest few steps (consecutive ticks, newest first), sent unreliable every step,
// so a lost packet is covered by the next one instead of stalling the rest behind a resend; the server keeps the ticks
// it hasn't seen.
constexpr int kInputRedundancy = 4;
struct MsgInputs {
    std::vector<MsgInput> steps;
    uint32_t ackSnap = 0, ackMask = 0; // #769: the newest player snapshot received (its sequence) and the 32 before it (bit i: ackSnap-1-i)
    uint8_t action = 0;              // what the hands do now (avatars, #750): 0 none, 1 dig / mine, 2 chop, 3 place / use
    float actionPhase = 0.0f;        // 0..1 within the swing
};
struct MsgPlayerState {              // a player's authoritative state (#327)
    uint32_t session = 0, tick = 0;
    uint32_t ackInput = 0;           // the last input tick the server applied (the owner replays the ones after it)
    uint64_t body = 0;               // the body's NetId the position is relative to
    DVec3 local{};                   // feet, body-local
    float vx = 0, vy = 0, vz = 0;    // velocity relative to the surface under the feet, on the body's axes
    float jumpBuffer = 0, coyote = 0, jumpLock = 0; // s (sent in ms, at most 0.255)
    float jetpackFuel = 1.0f;        // 0..1
    float yaw = 0.0f, pitch = 0.0f;  // the view direction, body-local: yaw about the body's y from +z, pitch from its xz plane
    uint8_t flags = 0;               // bit 0 grounded, 1 fly mode, 2 dead, 3 in water, 4 hovering, 5-6 posture, 7 jetpack firing
    uint8_t action = 0; float actionPhase = 0.0f; // the hands (relayed from the input, #750)
};
// #769 player snapshots. A player state as its quantised wire values (exactly what a PlayerState message carries:
// mm, cm/s, 1/65536 turn, ...): deltas are taken between these, so a decoded state is bit-identical to a full message's.
struct PlayerStateQ {
    uint32_t session = 0, ackInput = 0; uint64_t body = 0;
    int64_t p[3]{}; int64_t v[3]{};
    uint16_t yaw = 0, pitch = 0;
    uint8_t flags = 0, jumpBuffer = 0, coyote = 0, jumpLock = 0, fuel = 0, action = 0, phase = 0;
    bool operator==(const PlayerStateQ& o) const;
};
PlayerStateQ NetQuantize(const MsgPlayerState& m);
MsgPlayerState NetDequantize(const PlayerStateQ& q, uint32_t tick);
// One snapshot to one client: a sequence number (per client, acked back in MsgInputs), the server tick, and the states the
// client's interest picked. Each entry is full (baseSeq 0) or the fields that changed against the state of that player the
// client acknowledged in snapshot baseSeq (the writer fills `base`; the reader looks it up in what it received).
struct MsgPlayerSnapshot {
    uint32_t seq = 0, tick = 0;
    struct Entry { PlayerStateQ s; uint32_t baseSeq = 0; PlayerStateQ base; };
    std::vector<Entry> entries;
};
// #334 animals: the server's states of the ones near a client (position body-local to 1 mm, heading, state, health), and
// the ones that left its view or died. Sent 10 a second, unreliable (the client also drops one not heard of for 3 s).
struct MsgAnimals {
    uint32_t tick = 0;
    struct Entry { uint64_t netId = 0; uint8_t species = 0; uint64_t body = 0; DVec3 pos{}; float heading = 0.0f; uint8_t state = 0; float health = 0.0f, age = 1.0f; uint8_t flags = 0; };
    std::vector<Entry> entries;
    std::vector<uint64_t> removed;
};
struct MsgAnimalHit { uint64_t netId = 0; uint8_t action = 0; }; // 0 hit, 1 feed, 2 take its product
struct MsgDamage { float amount = 0.0f; uint8_t kind = 0; bool hasSource = false; uint64_t body = 0; DVec3 local{}; float knockback = 0.0f; uint64_t sourceId = 0; };
struct MsgShower { double start = 0.0, end = 0.0; uint32_t id = 0; }; // world seconds; id 0 none (events.h)

struct MsgRoster { struct Entry { uint32_t session = 0; std::string name; }; std::vector<Entry> players; }; // everyone online
struct MsgLoadout { uint32_t session = 0; std::string held; bool headlamp = false; uint8_t suitTier = 0; }; // session: 0 from a client
struct MsgTerrainEdit {
    uint64_t body = 0;               // NetId
    DVec3 center{};                  // body-local
    float radius = 0.0f;
    uint8_t flags = 0;               // bit 0 carve, 1 fill, 2 box, 3 soft only
    float hx = 0, hy = 0, hz = 0;    // box half extents (box only), along the body axes of ax/ay/az
    float ax[3] = { 1, 0, 0 }, ay[3] = { 0, 1, 0 };
    int8_t paint = -1;
    float strength = 0.0f;           // the soft brush (neither carve nor fill): density added (> 0) or taken
    uint32_t session = 0, seq = 0;   // whose request it applied (0: the server's own) and its number (#329)
};
// Terrain edits by players (#329): the client asks; the server checks it (reach from its own player, the tool in the
// player's loadout and its dig rule, the stroke rate, the cost) and applies it in its order, credits or charges that
// player's inventory (InventoryDelta to them), and sends the TerrainEdit to everyone (the asker too: it applies there,
// so every client gets the same edits in the same order). Refused: EditRefused to the asker, nothing changes.
enum class NetEditKind : uint8_t { Dig = 0, Place = 1, Count };
struct MsgEditRequest {
    uint32_t seq = 0;                // the client's request number (echoed in the TerrainEdit or the refusal)
    NetEditKind kind = NetEditKind::Dig;
    uint64_t body = 0;               // NetId
    DVec3 hit{};                     // the aimed surface point, body-local
    DVec3 center{};                  // the sphere (dig: sunk along the view; place: out along the normal), body-local
    float radius = 0.0f;
    std::string item;                // place: the filler that pays (dig: unused; the tool comes from the loadout)
};
enum class NetEditRefusal : uint8_t { Reach, Rate, Tool, Cost, Body, Mantle, Shape, Count };
const char* NetEditRefusalText(NetEditRefusal r);
struct MsgEditRefused { uint32_t seq = 0; NetEditRefusal why = NetEditRefusal::Reach; };
struct MsgInventoryDelta { std::string item; double amount = 0.0; }; // > 0 added, < 0 taken (fractional, as the inventory)
// Building (#330): as terrain edits. The client asks to place a part (into a plot, or a new plot anchored at a body-local
// point) or to break one; the server checks reach, the part, its room and its cost on that player's inventory, applies
// it, refunds or charges (InventoryDelta), and sends PlotEdits to everyone in its order. A new plot comes first as its
// own PlotEdit with the server's id and the exact float anchor, so every client makes the same plot.
struct NetPart { uint8_t type = 0; int32_t x = 0, y = 0, z = 0; uint8_t orient = 0, material = 0; }; // a Part (plots.h)
enum class NetPartOp : uint8_t { Place = 0, Remove = 1, NewPlot = 2, Count };
struct MsgPartRequest {
    uint32_t seq = 0;
    NetPartOp op = NetPartOp::Place;  // Place or Remove
    uint64_t body = 0;                // NetId
    uint32_t plot = 0;                // 0: place into a new plot at `anchor`
    float anchor[3] = { 0, 0, 0 };    // the new plot's body-local anchor (exact floats: the plot frame is built from them)
    NetPart part;                     // Place
    uint32_t partIndex = 0;           // Remove
};
enum class NetPartRefusal : uint8_t { Reach, Cost, Blocked, NoPlot, NoPart, Body, Kind, Rate, Count };
const char* NetPartRefusalText(NetPartRefusal r);
struct MsgPartRefused { uint32_t seq = 0; NetPartRefusal why = NetPartRefusal::Reach; };
struct MsgPlotEdit {
    uint32_t session = 0, seq = 0;    // whose request (0: the server's own)
    NetPartOp op = NetPartOp::Place;
    uint64_t body = 0;
    uint32_t plot = 0;
    float anchor[3] = { 0, 0, 0 };    // NewPlot
    NetPart part;                     // Place
    uint32_t partIndex = 0;           // Remove
};
struct MsgChat { uint32_t from = 0; std::string text; }; // #326: from = the sender's session (0 the server; a client's own is ignored)

void NetWrite(NetWriter& w, const MsgHello& m);
void NetWrite(NetWriter& w, const MsgWelcome& m);
void NetWrite(NetWriter& w, const MsgReject& m);
void NetWrite(NetWriter& w, const MsgBye& m);
void NetWrite(NetWriter& w, const MsgPing& m);
void NetWrite(NetWriter& w, const MsgPong& m);
void NetWrite(NetWriter& w, const MsgInputs& m);
void NetWrite(NetWriter& w, const MsgPlayerState& m);
void NetWrite(NetWriter& w, const MsgTerrainEdit& m);
void NetWrite(NetWriter& w, const MsgShower& m);
void NetWrite(NetWriter& w, const MsgRoster& m);
void NetWrite(NetWriter& w, const MsgLoadout& m);
void NetWrite(NetWriter& w, const MsgEditRequest& m);
void NetWrite(NetWriter& w, const MsgEditRefused& m);
void NetWrite(NetWriter& w, const MsgInventoryDelta& m);
void NetWrite(NetWriter& w, const MsgPartRequest& m);
void NetWrite(NetWriter& w, const MsgPartRefused& m);
void NetWrite(NetWriter& w, const MsgPlotEdit& m);
void NetWrite(NetWriter& w, const MsgChat& m);
bool NetRead(NetReader& r, MsgHello& m);
bool NetRead(NetReader& r, MsgWelcome& m);
bool NetRead(NetReader& r, MsgReject& m);
bool NetRead(NetReader& r, MsgBye& m);
bool NetRead(NetReader& r, MsgPing& m);
bool NetRead(NetReader& r, MsgPong& m);
bool NetRead(NetReader& r, MsgInputs& m);
bool NetRead(NetReader& r, MsgPlayerState& m);
bool NetRead(NetReader& r, MsgTerrainEdit& m);
bool NetRead(NetReader& r, MsgShower& m);
bool NetRead(NetReader& r, MsgRoster& m);
bool NetRead(NetReader& r, MsgLoadout& m);
bool NetRead(NetReader& r, MsgEditRequest& m);
bool NetRead(NetReader& r, MsgEditRefused& m);
bool NetRead(NetReader& r, MsgInventoryDelta& m);
bool NetRead(NetReader& r, MsgPartRequest& m);
bool NetRead(NetReader& r, MsgPartRefused& m);
bool NetRead(NetReader& r, MsgPlotEdit& m);
bool NetRead(NetReader& r, MsgChat& m);
void NetWrite(NetWriter& w, const MsgAnimals& m);
bool NetRead(NetReader& r, MsgAnimals& m);
void NetWrite(NetWriter& w, const MsgAnimalHit& m);
bool NetRead(NetReader& r, MsgAnimalHit& m);
void NetWrite(NetWriter& w, const MsgDamage& m);
bool NetRead(NetReader& r, MsgDamage& m);
void NetWrite(NetWriter& w, const MsgPlayerSnapshot& m);
// lookup(session, baseSeq, out): the state of that player received in snapshot baseSeq (false: unknown - the entry can't be
// decoded and the read fails). Without one (validation), deltas decode against zeros.
bool NetRead(NetReader& r, MsgPlayerSnapshot& m, const std::function<bool(uint32_t, uint32_t, PlayerStateQ&)>& lookup = nullptr);
// A message's id (the first byte; NetMsg::Count if empty or unknown).
NetMsg NetPeek(const std::vector<uint8_t>& bytes);
bool NetWellFormed(const std::vector<uint8_t>& bytes); // (#349) its reader accepts it (snapshot deltas decode against zeros)

// --test-netproto: every message round-trips field for field; 20000 random and truncated buffers never crash a
// reader and never read as valid when cut short. Prints PASS/FAIL; no network needed.
bool NetProto_Test();

#endif
