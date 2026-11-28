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
#include "factory.h" // MachineView (#332)
#include <cstdint>
#include <string>
#include <functional>
#include <vector>

constexpr uint32_t kNetProtocol = 66; // 1 was the text "join" of #311; 5: world time in Welcome and Pong (#339); 6 (unused); 7: world-gen settings in Welcome (#757); 8: terrain edit requests, refusals, inventory deltas (#329); 9: chat (#326); 10: building (#330); 11 is #776's; 12: player snapshots, delta-coded against acked states (#769); 13: animals (#334); 14: animal actions, damage to a player (#334); 15 (unused); 16: server-side survival (#337); 17: the suit colour in the loadout (#341); 18: beds (#800); 19: graves (#804); 20: the player id (UUID) in the hello (#346); 21: team quest books (#338); 22: crafting on the server (#338); 23: machines on the server (#332); 24 is #604's; 25: crafting from nearby chests, cursor takes (#854); 26 is #604's; 27: plants and farming (#335); 28: claim refusals (#343); 29: MachineAction 8, insert into one chest slot, 9: take some of an item from a chest (#862, landed after 30); 30: claim banners, the claims list (#343); 31: dropped piles in Graves, DropItems, CreativeGive (#865, #866); 32: recursive flag in CraftRequest (#946); 33: SkillState, the player's skills from the server, and SkillSpend (#920); 34: edit kind Till and the till flag (#1085); 35: a campfire's fire in MachineView (#1103); 36: the research lab part types in NetPart (#1193); 37: the magic machines (research table, mana siphon, mana dynamo) in NetPart (#1252); 38: Act bits LockOn and Dodge before MagicSight (#1216, #1218); 39: CreatureFx, creature shots and events to clients (#1343, #1364); 40: TeamState members' session and last seen, the team's claim budget and allies (#828); 41 is the Ladder part's (UI/UX); 42: the alchemy part types (distillery, essence tank, cauldron) in NetPart (#1254); 43: the Ladder part type, after the alchemy ones (#1105; 41 was its number before alchemy merged first); 44-52 are the ladder branch's; 53: a creative flyer's pose in each input step (#1922); 54-56 are the ladder branch's; 57: TeamStore, the team's progression store to its members (#1886); 58 (unused, was #2047's first number); 59 is #2059's (#1735); 60: the ruin cache part type in NetPart, MachineAction 13 opens a cache (#2072); 61-63 are #2103's and the ladder's (#2071); 64 (unused, was #2072's first number); 65: MapPin and ResearchQueue, a client's map pins and research queue changes for its team (#1886); 66: a machine's alternate recipe in MachineView, MachineAction 14 chooses it (#2072)

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
    SyncIndex,     // server -> client: the join snapshot's entries (regions, plots) with hashes, in pages (#776)
    SyncWant,      // client -> server: which of them it needs
    SyncData,      // server -> client: a piece of an entry
    SyncDone,      // server -> client: all wanted entries sent; queued live edits apply now
    PlayerSnapshot = 24, // server -> client (#769): the player states of one snapshot in one message, each delta-coded against one
                   // the client acknowledged (20-23 are #776's)
    Animals = 25,  // server -> client (#334): the animals near the client (state, health), and the ones gone from its view
    AnimalHit = 26, // client -> server (#334): a swing at an animal, or feed / take its product (the server checks reach and cooldown)
    Damage = 27,   // server -> client (#334): a hit on the client's player (an animal's bite): the hit's effect (the health is in Vitals)
    Vitals = 28,   // server -> client (#337): its player's survival numbers (health, food, oxygen, power, temperature, dead), 5 a second
    UseItem = 29,  // client -> server (#337): eat or use an item (oxygen canister, power cell) from its inventory as the server keeps it
    BedUse = 30,   // client -> server (#800): make the bed at (body, plot, cell) this player's spawn
    Graves = 31,   // server -> client (#804): the world's graves (who, where, which stone), resent when they change
    QuestBook = 32, // server -> client (#338): the player's team's quest book (its save lines), when it changes and on join
    CraftRequest = 33, // client -> server (#338): craft a recipe (count times) from its inventory as the server keeps it
    MachineAction = 34, // client -> server (#332): insert, extract, take a slot, fuel; watch / unwatch a machine's screen
    MachineView = 35,   // server -> client (#332): the full view of a machine it watches, a few times a second
    MachineStates = 36, // server -> client (#332): which machines of a plot work now (animations, smoke), when it changes
    TeamState = 37,     // server -> client (#604): the player's team (name, members online with their health) and an invite waiting
    PlantAction = 38,   // client -> server (#335): fell a plant (tree, bush, loose stone), plant or harvest a farmland part
    Plants = 39,        // server -> client (#335): plants removed (id, world time: regrowth runs on the shared clock), at the join and as they go
    FarmFields = 40,    // server -> client (#335): the crops on farmland near the client (crop, growth), when they change
    Claims = 41,        // server -> client (#343): the land claims (team, colour, where, how far; ours / allowed), when they change
    DropItems = 42,     // client -> server (#865): drop some of an item at the feet (the server makes the pile; 41 is #343's Claims)
    CreativeGive = 43,  // client -> server (#866): creative items into the inventory the server keeps ("" clears it)
    SkillState = 44,    // server -> client (#920): the player's skills (XP and perks) as the save's lines, when they change
    SkillSpend = 45,    // client -> server (#920): spend a point on a perk (the server checks it and answers with SkillState)
    CreatureFx = 46,    // server -> client (#1343, #1364): creature projectiles launched and ended, and creature events (calls, telegraphs,
                        // hurt, death) near the client, for its own flights and sounds
    TeamStore = 48,     // server -> client (#1886): the player's team's progression store (research, explored map, blueprint
                        // library) as TeamStore_Text, on join, on a team change and when it changes (47 is the ladder's SetGem)
    MapPin = 49,        // client -> server (#1886): a pin placed on (or removed from) the map, for the player's team's shared pins
    ResearchQueue = 50, // client -> server (#1886): a change to the team's research queue (add, remove, up), from the research tab
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
    void blob(const uint8_t* data, size_t n) { varu(n); bytes.insert(bytes.end(), data, data + n); } // length, then the bytes
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
    std::vector<uint8_t> blob(size_t maxLen) { // length, then the bytes (longer than maxLen: malformed)
        const uint64_t n = varu();
        if (n > maxLen || !need((size_t)n)) { ok = false; return {}; }
        std::vector<uint8_t> out(p, p + n); p += n;
        return out;
    }
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
    std::string uuid;     // (#346) the player's id (a UUID): their record on this server; the name is only what others see
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
    // (#1922) Creative flight: the client's own state after this step. The server can't step a flying remote itself (the
    // fly camera is the server's world.camera), so a creative player's server copy takes it from here instead.
    bool fly = false;
    uint64_t flyBody = 0;            // NetId the pose is relative to
    DVec3 flyLocal{};                // feet, body-local (1 mm)
    float fvx = 0, fvy = 0, fvz = 0; // velocity on the body's axes, m/s (1 cm/s)
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
struct MsgVitals { float health = 100.0f, food = 100.0f, oxygen = 100.0f, power = 100.0f, temperature = 20.0f; uint8_t flags = 0; uint16_t deaths = 0; }; // flags: 1 dead, 2 breathable, 4 vacuum, 8 underwater, 16 life support, 32 sick, 64 headlamp
struct MsgUseItem { std::string item; };
struct MsgBedUse { std::string body; uint32_t plot = 0; int32_t cell[3] = { 0, 0, 0 }; }; // #800: the bed's origin cell
struct MsgQuestBook { uint32_t team = 0; std::string text; }; // #338 (at most 64 KB)
struct MsgTeamStore { // #1886: the team's TeamStore_Text
    static constexpr size_t kMaxText = 480 * 1024; // bytes: under GameNetworkingSockets' 512 KB reliable message
    uint32_t team = 0;
    std::string text;
};
struct MsgMapPin { // #1886: the team's pins come back in the store's "pins" section
    static constexpr size_t kMaxName = 48; // bytes
    bool add = true;                       // false: remove the team's pin of this name on this body
    std::string body, name;
    float dir[3] = { 0, 1, 0 };            // body-local unit direction (f32: a unit vector's int16 is ~2 m on a 60 km planet)
};
struct MsgResearchQueue { uint8_t op = 1; std::string node; }; // #1886: op as ResearchQueueOp (1 add, 2 remove, 3 up); the store answers
struct MsgCraftRequest { uint32_t recipe = 0; std::string output; uint8_t count = 1; bool useChests = false; bool recursive = false; }; // recursive (#946): the server plans and crafts the missing intermediates first (#928) // useChests (#854): what the inventory lacks, from the chests near the player // #338 (output: a check that both sides mean the same recipe)
// #332. action: 1 insert item (count), 2 insert all it takes, 3 insert fuel (item, count), 4 extract the outputs,
// 5 take one slot (slot, which: 0 inputs, 1 outputs, 2 fuel), 6 watch (its screen opened), 7 unwatch.
struct MsgMachineAction { uint64_t body = 0; uint32_t plot = 0; int32_t part = -1; uint8_t action = 0; std::string item; uint16_t count = 0; uint8_t slot = 0, which = 0; };
struct MsgMachineView { uint64_t body = 0; uint32_t plot = 0; int32_t part = -1; MachineView view; };
struct MsgMachineStates { uint64_t body = 0; uint32_t plot = 0; std::vector<int32_t> working; };
struct MsgTeamState { // #604: at most 16 members
    struct Member { std::string name; bool online = false, owner = false, you = false; float health = 0.0f;
                    uint32_t session = 0;  // (#828) online: theirs (where they are: the avatars); 0 offline
                    uint32_t lastSeen = 0; // (#828) offline: seconds since they left (0: not known)
    };
    int32_t team = 0; std::string name; std::vector<Member> members; std::string invite; // invite: the team inviting us ("" none)
    int claimBudget = 0; std::vector<std::string> allies; // (#828) claims the team may hold; the teams it lets into its claims
};
// (#335) plants and farming
// 1 fell (plant), 2 sow (x y z: the body-local ground, seed), 3 harvest (x y z: the crop's ground) (#1085: crops by position)
struct MsgPlantAction { uint8_t action = 0; uint64_t body = 0; uint64_t plant = 0; float x = 0, y = 0, z = 0; std::string seed; };
struct MsgPlants { uint64_t body = 0; uint32_t session = 0; std::vector<std::pair<uint64_t, double>> removed; }; // session: who felled them (0: the join's list)
struct NetFarmField { float x = 0, y = 0, z = 0; std::string crop; float growth = 0.0f; }; // x y z: body-local, the crop's ground (#1085)
struct MsgFarmFields { uint64_t body = 0; std::vector<NetFarmField> fields; }; // every field on the body near the client (the full list: replaces)
struct MsgDropItems { std::string item; int64_t centi = 0; };        // #865: hundredths of a unit
struct MsgCreativeGive { std::string item; uint32_t count = 0; };   // #866: "" = clear the inventory
struct NetClaim { uint32_t id = 0; int32_t team = 0; std::string teamName; uint32_t colour = 0; uint64_t body = 0; float centre[3] = {}; float radius = 0.0f; uint8_t flags = 0; }; // flags: 1 ours, 2 we may build there (ours or an ally's)
struct MsgClaims { std::vector<NetClaim> claims; }; // (#343) every claim (the full list: replaces)
struct MsgSkillState { std::string text; }; // (#920) Skills_SaveText of the player's skills (the whole state: replaces the client's)
struct MsgSkillSpend { std::string perk; }; // (#920) a perk id (skills.txt)
struct MsgGraves { // #804 (at most 256: the nearest matter, and a world with more has bigger problems)
    struct Entry { uint32_t id = 0; std::string body; DVec3 local{}; std::string owner, name; uint8_t style = 0; int64_t diedAt = 0;
                    bool dropped = false; std::vector<std::string> items; }; // (#865) a dropped pile: what it holds, up to 3 ids (its bundle's colours)
    std::vector<Entry> graves;
};
struct MsgDamage { float amount = 0.0f; uint8_t kind = 0; bool hasSource = false; uint64_t body = 0; DVec3 local{}; float knockback = 0.0f; uint64_t sourceId = 0; };
// (#1343, #1364) What a client needs to see and hear creatures: shots launched (it flies them itself: start, velocity, kind 0 stone,
// 1 spear, 2 dart), shots ended (hit or stopped), and creature events (CreatureEvent: 0 idle, 1 alert, 2 telegraph, 3 hurt, 4
// death; the creature, its species index, where). Sent each tick there is any, unreliable (lost ones only cost a sound or a flight).
struct MsgCreatureFx {
    struct Shot { uint32_t id = 0; uint64_t body = 0; DVec3 pos{}; Vector3 vel{}; uint8_t kind = 0; uint64_t owner = 0; };
    struct Event { uint64_t netId = 0; uint8_t kind = 0; uint8_t species = 0; uint64_t body = 0; DVec3 pos{}; };
    std::vector<Shot> shots;
    std::vector<uint32_t> ended;
    std::vector<Event> events;
};
struct MsgShower { double start = 0.0, end = 0.0; uint32_t id = 0; }; // world seconds; id 0 none (events.h)

struct MsgRoster { struct Entry { uint32_t session = 0; std::string name; }; std::vector<Entry> players; }; // everyone online
struct MsgLoadout { // session: 0 from a client
    uint32_t session = 0; std::string held; bool headlamp = false; uint8_t suitTier = 0;
    uint32_t suitColor = 0;          // 0xRRGGBB, picked when joining (#341); 0: the default suit
};
struct MsgTerrainEdit {
    uint64_t body = 0;               // NetId
    DVec3 center{};                  // body-local
    float radius = 0.0f;
    uint8_t flags = 0;               // bit 0 carve, 1 fill, 2 box, 3 soft only, 4 till (#1085)
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
enum class NetEditKind : uint8_t { Dig = 0, Place = 1, Till = 2, Count }; // Till (#1085): paint soil farmland (a shovel or hoe)
struct MsgEditRequest {
    uint32_t seq = 0;                // the client's request number (echoed in the TerrainEdit or the refusal)
    NetEditKind kind = NetEditKind::Dig;
    uint64_t body = 0;               // NetId
    DVec3 hit{};                     // the aimed surface point, body-local
    DVec3 center{};                  // the sphere (dig: sunk along the view; place: out along the normal), body-local
    float radius = 0.0f;
    std::string item;                // place: the filler that pays (dig: unused; the tool comes from the loadout)
};
enum class NetEditRefusal : uint8_t { Reach, Rate, Tool, Cost, Body, Mantle, Shape, Claim, Count }; // Claim: inside another team's claim (#343)
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
enum class NetPartRefusal : uint8_t { Reach, Cost, Blocked, NoPlot, NoPart, Body, Kind, Rate, Claim, Count }; // Claim (#343)
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
// Join-in-progress sync (#776): at a join the server takes a snapshot of the world's edits (terrain regions, plots) and
// sends its index in pages; the client answers with the entries it lacks (it compares each region's hash with its own:
// a rejoin with a cached client save fetches only what changed); the server streams those, nearest first, in pieces
// within a byte budget, then SyncDone. Live edits after the snapshot follow as usual (the client applies them after).
enum class NetSyncKind : uint8_t { Region = 0, Plots = 1, Count };
struct NetSyncEntry { NetSyncKind kind = NetSyncKind::Region; uint64_t body = 0, key = 0, hash = 0; uint32_t size = 0; }; // key: the region (Plots: 0)
constexpr size_t kSyncIndexPage = 256;      // entries a page
constexpr size_t kSyncPiece = 16 * 1024;    // bytes a data message
constexpr uint32_t kSyncMaxEntries = 1u << 20;
struct MsgSyncIndex { uint32_t snapshot = 0; uint16_t page = 0, pages = 1; uint32_t total = 0; std::vector<NetSyncEntry> entries; };
struct MsgSyncWant { uint32_t snapshot = 0; uint32_t total = 0; std::vector<uint8_t> bits; }; // bit i: entry i is wanted
struct MsgSyncData { uint32_t snapshot = 0, entry = 0, offset = 0, size = 0; std::vector<uint8_t> bytes; };
struct MsgSyncDone { uint32_t snapshot = 0, entries = 0; uint64_t bytes = 0; };

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
void NetWrite(NetWriter& w, const MsgSyncIndex& m);
void NetWrite(NetWriter& w, const MsgSyncWant& m);
void NetWrite(NetWriter& w, const MsgSyncData& m);
void NetWrite(NetWriter& w, const MsgSyncDone& m);
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
bool NetRead(NetReader& r, MsgSyncIndex& m);
bool NetRead(NetReader& r, MsgSyncWant& m);
bool NetRead(NetReader& r, MsgSyncData& m);
bool NetRead(NetReader& r, MsgSyncDone& m);
void NetWrite(NetWriter& w, const MsgAnimals& m);
bool NetRead(NetReader& r, MsgAnimals& m);
void NetWrite(NetWriter& w, const MsgAnimalHit& m);
bool NetRead(NetReader& r, MsgAnimalHit& m);
void NetWrite(NetWriter& w, const MsgDamage& m);
void NetWrite(NetWriter& w, const MsgCreatureFx& m);
bool NetRead(NetReader& r, MsgCreatureFx& m);
bool NetRead(NetReader& r, MsgDamage& m);
void NetWrite(NetWriter& w, const MsgVitals& m);
bool NetRead(NetReader& r, MsgVitals& m);
void NetWrite(NetWriter& w, const MsgUseItem& m);
bool NetRead(NetReader& r, MsgUseItem& m);
void NetWrite(NetWriter& w, const MsgBedUse& m);
bool NetRead(NetReader& r, MsgBedUse& m);
void NetWrite(NetWriter& w, const MsgGraves& m);
void NetWrite(NetWriter& w, const MsgQuestBook& m);
bool NetRead(NetReader& r, MsgQuestBook& m);
void NetWrite(NetWriter& w, const MsgTeamStore& m);
bool NetRead(NetReader& r, MsgTeamStore& m);
void NetWrite(NetWriter& w, const MsgMapPin& m);
bool NetRead(NetReader& r, MsgMapPin& m);
void NetWrite(NetWriter& w, const MsgResearchQueue& m);
bool NetRead(NetReader& r, MsgResearchQueue& m);
void NetWrite(NetWriter& w, const MsgCraftRequest& m);
bool NetRead(NetReader& r, MsgCraftRequest& m);
void NetWrite(NetWriter& w, const MsgMachineAction& m);
bool NetRead(NetReader& r, MsgMachineAction& m);
void NetWrite(NetWriter& w, const MsgMachineView& m);
bool NetRead(NetReader& r, MsgMachineView& m);
void NetWrite(NetWriter& w, const MsgMachineStates& m);
bool NetRead(NetReader& r, MsgMachineStates& m);
void NetWrite(NetWriter& w, const MsgTeamState& m);
bool NetRead(NetReader& r, MsgTeamState& m);
void NetWrite(NetWriter& w, const MsgPlantAction& m);
bool NetRead(NetReader& r, MsgPlantAction& m);
void NetWrite(NetWriter& w, const MsgPlants& m);
bool NetRead(NetReader& r, MsgPlants& m);
void NetWrite(NetWriter& w, const MsgFarmFields& m);
bool NetRead(NetReader& r, MsgFarmFields& m);
void NetWrite(NetWriter& w, const MsgDropItems& m);
bool NetRead(NetReader& r, MsgDropItems& m);
void NetWrite(NetWriter& w, const MsgCreativeGive& m);
bool NetRead(NetReader& r, MsgCreativeGive& m);
void NetWrite(NetWriter& w, const MsgClaims& m);
bool NetRead(NetReader& r, MsgClaims& m);
void NetWrite(NetWriter& w, const MsgSkillState& m);
bool NetRead(NetReader& r, MsgSkillState& m);
void NetWrite(NetWriter& w, const MsgSkillSpend& m);
bool NetRead(NetReader& r, MsgSkillSpend& m);
bool NetRead(NetReader& r, MsgGraves& m);
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
