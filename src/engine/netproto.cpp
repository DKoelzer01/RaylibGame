#include "netproto.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>

const char* NetMsgName(NetMsg m) {
    switch (m) {
        case NetMsg::Hello: return "hello"; case NetMsg::Welcome: return "welcome"; case NetMsg::Reject: return "reject";
        case NetMsg::Bye: return "bye"; case NetMsg::Ping: return "ping"; case NetMsg::Pong: return "pong";
        case NetMsg::Input: return "input"; case NetMsg::PlayerState: return "player state"; case NetMsg::TerrainEdit: return "terrain edit"; case NetMsg::Shower: return "shower"; case NetMsg::Roster: return "roster"; case NetMsg::Loadout: return "loadout";
        case NetMsg::PartRequest: return "part request"; case NetMsg::PartRefused: return "part refused"; case NetMsg::PlotEdit: return "plot edit";
        case NetMsg::EditRequest: return "edit request"; case NetMsg::EditRefused: return "edit refused"; case NetMsg::InventoryDelta: return "inventory delta"; case NetMsg::Chat: return "chat"; case NetMsg::PlayerSnapshot: return "player snapshot"; case NetMsg::Animals: return "animals"; case NetMsg::AnimalHit: return "animal hit"; case NetMsg::Damage: return "damage"; case NetMsg::Vitals: return "vitals"; case NetMsg::UseItem: return "use item"; case NetMsg::BedUse: return "bed use"; case NetMsg::Graves: return "graves"; case NetMsg::QuestBook: return "quest book"; case NetMsg::CraftRequest: return "craft request"; case NetMsg::MachineAction: return "machine action"; case NetMsg::MachineView: return "machine view"; case NetMsg::MachineStates: return "machine states"; case NetMsg::TeamState: return "team state"; case NetMsg::PlantAction: return "plant action"; case NetMsg::Plants: return "plants"; case NetMsg::FarmFields: return "farm fields"; case NetMsg::Claims: return "claims"; case NetMsg::DropItems: return "drop items"; case NetMsg::CreativeGive: return "creative give"; case NetMsg::SkillState: return "skill state"; case NetMsg::SkillSpend: return "skill spend";
        case NetMsg::SyncIndex: return "sync index"; case NetMsg::SyncWant: return "sync want"; case NetMsg::SyncData: return "sync data"; case NetMsg::SyncDone: return "sync done";
        default: return "?";
    }
}
const char* NetRejectText(NetReject r) {
    switch (r) {
        case NetReject::Protocol: return "a different game version (network protocol)";
        case NetReject::Build: return "a different game build";
        case NetReject::Password: return "wrong password";
        case NetReject::Full: return "the server is full";
        case NetReject::Banned: return "banned from this server";
        case NetReject::NameTaken: return "that name is taken"; // (#346: playing now, or another player's on this server; the detail says which)
        case NetReject::Malformed: return "a garbled join request";
    }
    return "?";
}
const char* NetByeText(NetBye r) {
    switch (r) {
        case NetBye::Quit: return "left"; case NetBye::ServerClosed: return "the server closed";
        case NetBye::TimedOut: return "timed out"; case NetBye::Kicked: return "kicked";
    }
    return "?";
}

// ---- Writer / reader ----

void NetWriter::f32(float v) { uint32_t b; memcpy(&b, &v, 4); u32(b); }
void NetWriter::f64(double v) { uint64_t b; memcpy(&b, &v, 8); u64(b); }
void NetWriter::pos(const DVec3& p) { vari(llround(p.x * 1000.0)); vari(llround(p.y * 1000.0)); vari(llround(p.z * 1000.0)); }
static int16_t Q16(float v) { return (int16_t)lroundf(std::fmax(-1.0f, std::fmin(1.0f, v)) * 32767.0f); }
void NetWriter::unit(float x, float y, float z) { u16((uint16_t)Q16(x)); u16((uint16_t)Q16(y)); u16((uint16_t)Q16(z)); }
void NetWriter::angle(float a) {
    const float turn = 2.0f * 3.14159265358979f;
    float t = std::fmod(a, turn); if (t < 0.0f) t += turn;
    u16((uint16_t)((uint32_t)lroundf(t / turn * 65536.0f) & 0xFFFF));
}
void NetWriter::speed(float v) { vari(lroundf(v * 100.0f)); }

uint64_t NetReader::varu() {
    uint64_t v = 0;
    for (int shift = 0; shift < 64; shift += 7) {
        const uint8_t b = u8();
        if (!ok) return 0;
        v |= (uint64_t)(b & 0x7F) << shift;
        if (!(b & 0x80)) return v;
    }
    ok = false; // more than 10 bytes: not a varint
    return 0;
}
float NetReader::f32() { const uint32_t b = u32(); float v; memcpy(&v, &b, 4); return ok && std::isfinite(v) ? v : (ok = false, 0.0f); }
double NetReader::f64() { const uint64_t b = u64(); double v; memcpy(&v, &b, 8); return ok && std::isfinite(v) ? v : (ok = false, 0.0); }
std::string NetReader::str(size_t maxLen) {
    const uint64_t n = varu();
    if (!ok || n > maxLen || !need((size_t)n)) { ok = false; return {}; }
    std::string s((const char*)p, (size_t)n);
    p += n;
    return s;
}
DVec3 NetReader::pos() { const int64_t x = vari(), y = vari(), z = vari(); return ok ? DVec3{ x / 1000.0, y / 1000.0, z / 1000.0 } : DVec3{}; }
void NetReader::unit(float& x, float& y, float& z) {
    x = (int16_t)u16() / 32767.0f; y = (int16_t)u16() / 32767.0f; z = (int16_t)u16() / 32767.0f;
    if (!ok) x = y = z = 0.0f;
}
float NetReader::angle() { return u16() / 65536.0f * 2.0f * 3.14159265358979f; }
float NetReader::speed() { return vari() / 100.0f; }

// ---- Messages ----

static void Id(NetWriter& w, NetMsg m) { w.u8((uint8_t)m); }
// After the fields: the message must end exactly there (trailing bytes mean a layout mismatch).
static bool Done(NetReader& r) { if (!r.atEnd()) r.ok = false; return r.ok; }

void NetWrite(NetWriter& w, const MsgHello& m) { Id(w, NetMsg::Hello); w.u32(m.protocol); w.str(m.build); w.str(m.name); w.str(m.password); w.str(m.uuid); }
bool NetRead(NetReader& r, MsgHello& m) {
    m.protocol = r.u32();
    // An older or newer protocol: read nothing more (its layout may differ), the server answers Reject::Protocol.
    if (r.ok && m.protocol != kNetProtocol) return true;
    m.build = r.str(64); m.name = r.str(32); m.password = r.str(64); m.uuid = r.str(40);
    return Done(r) && !m.name.empty();
}

void NetWrite(NetWriter& w, const MsgWelcome& m) {
    Id(w, NetMsg::Welcome);
    w.u32(m.protocol); w.u32(m.session); w.boolean(m.reconnected); w.u32(m.seed);
    w.varu(m.bodies.size()); for (const std::string& b : m.bodies) w.str(b);
    w.u64(m.fingerprint); w.u32(m.serverTick); w.f64(m.worldTime); w.str(m.worldgen); w.boolean(m.sandbox); w.str(m.motd.substr(0, 200));
    w.u8(m.settings.difficulty); w.boolean(m.settings.pvp); w.boolean(m.settings.claims); w.u8(m.settings.maxPlayers);
}
bool NetRead(NetReader& r, MsgWelcome& m) {
    m.protocol = r.u32(); m.session = r.u32(); m.reconnected = r.boolean(); m.seed = r.u32();
    const uint64_t n = r.varu();
    if (n > 256) r.ok = false;
    m.bodies.clear();
    for (uint64_t i = 0; i < n && r.ok; ++i) m.bodies.push_back(r.str(1024));
    m.fingerprint = r.u64(); m.serverTick = r.u32(); m.worldTime = r.f64(); m.worldgen = r.str(4096); m.sandbox = r.boolean(); m.motd = r.str(200);
    m.settings.difficulty = r.u8(); m.settings.pvp = r.boolean(); m.settings.claims = r.boolean(); m.settings.maxPlayers = r.u8();
    if (m.settings.difficulty > 2) r.ok = false;
    return Done(r);
}

void NetWrite(NetWriter& w, const MsgReject& m) { Id(w, NetMsg::Reject); w.u8((uint8_t)m.reason); w.str(m.detail); }
bool NetRead(NetReader& r, MsgReject& m) {
    const uint8_t v = r.u8();
    if (v < (uint8_t)NetReject::Protocol || v > (uint8_t)NetReject::Malformed) r.ok = false;
    m.reason = (NetReject)v; m.detail = r.str(256);
    return Done(r);
}

void NetWrite(NetWriter& w, const MsgBye& m) { Id(w, NetMsg::Bye); w.u8((uint8_t)m.reason); }
bool NetRead(NetReader& r, MsgBye& m) {
    const uint8_t v = r.u8();
    if (v < (uint8_t)NetBye::Quit || v > (uint8_t)NetBye::Kicked) r.ok = false;
    m.reason = (NetBye)v;
    return Done(r);
}

void NetWrite(NetWriter& w, const MsgPing& m) { Id(w, NetMsg::Ping); w.u32(m.seq); w.f64(m.clientTime); }
bool NetRead(NetReader& r, MsgPing& m) { m.seq = r.u32(); m.clientTime = r.f64(); return Done(r); }
void NetWrite(NetWriter& w, const MsgPong& m) { Id(w, NetMsg::Pong); w.u32(m.seq); w.f64(m.clientTime); w.u32(m.serverTick); w.f64(m.worldTime); }
bool NetRead(NetReader& r, MsgPong& m) { m.seq = r.u32(); m.clientTime = r.f64(); m.serverTick = r.u32(); m.worldTime = r.f64(); return Done(r); }

void NetWrite(NetWriter& w, const MsgInputs& m) { // newest tick, then the steps newest first (ticks implicit)
    Id(w, NetMsg::Input);
    const size_t n = std::min<size_t>(m.steps.size(), 8);
    w.u32(n ? m.steps[0].tick : 0); w.u8((uint8_t)n);
    for (size_t i = 0; i < n; ++i) { const MsgInput& s = m.steps[i]; w.varu(s.held); w.varu(s.pressed); w.unit(s.fx, s.fy, s.fz); }
    w.u8(m.action); w.u8((uint8_t)std::clamp((int)std::lround(m.actionPhase * 255.0f), 0, 255));
    w.u32(m.ackSnap); w.u32(m.ackMask); // (#769)
}
bool NetRead(NetReader& r, MsgInputs& m) {
    const uint32_t newest = r.u32();
    const uint8_t n = r.u8();
    if (n == 0 || n > 8 || n > newest) r.ok = false;
    m.steps.clear();
    for (uint32_t i = 0; i < n && r.ok; ++i) {
        MsgInput s; s.tick = newest - i;
        s.held = r.varu(); s.pressed = r.varu(); r.unit(s.fx, s.fy, s.fz);
        m.steps.push_back(s);
    }
    m.action = r.u8(); m.actionPhase = r.u8() / 255.0f;
    m.ackSnap = r.u32(); m.ackMask = r.u32(); // (#769)
    if (m.action > 3) r.ok = false;
    return Done(r);
}

void NetWrite(NetWriter& w, const MsgPlayerState& m) {
    Id(w, NetMsg::PlayerState);
    w.varu(m.session); w.u32(m.tick); w.u32(m.ackInput); w.varu(m.body); w.pos(m.local);
    w.speed(m.vx); w.speed(m.vy); w.speed(m.vz); w.angle(m.yaw); w.u16((uint16_t)Q16(m.pitch / 1.6f)); w.u8(m.flags);
    auto ms = [&](float s) { w.u8((uint8_t)std::clamp((int)std::lround(s * 1000.0f), 0, 255)); };
    ms(m.jumpBuffer); ms(m.coyote); ms(m.jumpLock);
    w.u8((uint8_t)std::clamp((int)std::lround(m.jetpackFuel * 255.0f), 0, 255));
    w.u8(m.action); w.u8((uint8_t)std::clamp((int)std::lround(m.actionPhase * 255.0f), 0, 255));
}
bool NetRead(NetReader& r, MsgPlayerState& m) {
    m.session = (uint32_t)r.varu(); m.tick = r.u32(); m.ackInput = r.u32(); m.body = r.varu(); m.local = r.pos();
    m.vx = r.speed(); m.vy = r.speed(); m.vz = r.speed(); m.yaw = r.angle(); m.pitch = (int16_t)r.u16() / 32767.0f * 1.6f; m.flags = r.u8();
    m.jumpBuffer = r.u8() / 1000.0f; m.coyote = r.u8() / 1000.0f; m.jumpLock = r.u8() / 1000.0f; m.jetpackFuel = r.u8() / 255.0f;
    m.action = r.u8(); m.actionPhase = r.u8() / 255.0f;
    if (m.action > 3) r.ok = false;
    return Done(r);
}

// ---- #769 player snapshots ----
bool PlayerStateQ::operator==(const PlayerStateQ& o) const {
    return session == o.session && ackInput == o.ackInput && body == o.body && p[0] == o.p[0] && p[1] == o.p[1] && p[2] == o.p[2] &&
           v[0] == o.v[0] && v[1] == o.v[1] && v[2] == o.v[2] && yaw == o.yaw && pitch == o.pitch && flags == o.flags && jumpBuffer == o.jumpBuffer &&
           coyote == o.coyote && jumpLock == o.jumpLock && fuel == o.fuel && action == o.action && phase == o.phase;
}
// Through the PlayerState message itself: the quantised values are exactly what it carries, and back exactly what it reads.
PlayerStateQ NetQuantize(const MsgPlayerState& m) {
    NetWriter w; NetWrite(w, m);
    NetReader r(w.bytes); r.u8(); // (the id)
    PlayerStateQ q;
    q.session = (uint32_t)r.varu(); r.u32(); q.ackInput = r.u32(); q.body = r.varu();
    for (int64_t& x : q.p) x = r.vari();
    for (int64_t& x : q.v) x = r.vari();
    q.yaw = r.u16(); q.pitch = r.u16(); q.flags = r.u8();
    q.jumpBuffer = r.u8(); q.coyote = r.u8(); q.jumpLock = r.u8(); q.fuel = r.u8(); q.action = r.u8(); q.phase = r.u8();
    return q;
}
MsgPlayerState NetDequantize(const PlayerStateQ& q, uint32_t tick) {
    NetWriter w; Id(w, NetMsg::PlayerState);
    w.varu(q.session); w.u32(tick); w.u32(q.ackInput); w.varu(q.body);
    for (int64_t x : q.p) w.vari(x);
    for (int64_t x : q.v) w.vari(x);
    w.u16(q.yaw); w.u16(q.pitch); w.u8(q.flags);
    w.u8(q.jumpBuffer); w.u8(q.coyote); w.u8(q.jumpLock); w.u8(q.fuel); w.u8(q.action); w.u8(q.phase);
    NetReader r(w.bytes); r.u8();
    MsgPlayerState m; NetRead(r, m);
    return m;
}
// Entry: session, the base's age in snapshots (0 = full), a mask of the fields present, then those fields - against the base:
// differences for the numbers that move (input tick, position, velocity), new values for the rest.
namespace {
enum : uint16_t { kQAck = 1, kQBody = 2, kQPos = 4, kQVel = 8, kQYaw = 16, kQPitch = 32, kQFlags = 64, kQTimers = 128, kQFuel = 256, kQHands = 512, kQAll = 1023 };
}
void NetWrite(NetWriter& w, const MsgPlayerSnapshot& m) {
    Id(w, NetMsg::PlayerSnapshot);
    w.u32(m.seq); w.u32(m.tick); w.varu(m.entries.size());
    for (const MsgPlayerSnapshot::Entry& e : m.entries) {
        const PlayerStateQ& s = e.s; const PlayerStateQ b = e.baseSeq ? e.base : PlayerStateQ{};
        uint16_t mask = kQAll;
        if (e.baseSeq) {
            mask = 0;
            if (s.ackInput != b.ackInput) mask |= kQAck;
            if (s.body != b.body) mask |= kQBody;
            if (s.p[0] != b.p[0] || s.p[1] != b.p[1] || s.p[2] != b.p[2]) mask |= kQPos;
            if (s.v[0] != b.v[0] || s.v[1] != b.v[1] || s.v[2] != b.v[2]) mask |= kQVel;
            if (s.yaw != b.yaw) mask |= kQYaw;
            if (s.pitch != b.pitch) mask |= kQPitch;
            if (s.flags != b.flags) mask |= kQFlags;
            if (s.jumpBuffer != b.jumpBuffer || s.coyote != b.coyote || s.jumpLock != b.jumpLock) mask |= kQTimers;
            if (s.fuel != b.fuel) mask |= kQFuel;
            if (s.action != b.action || s.phase != b.phase) mask |= kQHands;
        }
        w.varu(s.session); w.varu(e.baseSeq ? m.seq - e.baseSeq : 0); w.varu(mask);
        if (mask & kQAck) w.vari((int64_t)s.ackInput - (int64_t)b.ackInput);
        if (mask & kQBody) w.varu(s.body);
        if (mask & kQPos) for (int i = 0; i < 3; ++i) w.vari(s.p[i] - b.p[i]);
        if (mask & kQVel) for (int i = 0; i < 3; ++i) w.vari(s.v[i] - b.v[i]);
        if (mask & kQYaw) w.u16(s.yaw);
        if (mask & kQPitch) w.u16(s.pitch);
        if (mask & kQFlags) w.u8(s.flags);
        if (mask & kQTimers) { w.u8(s.jumpBuffer); w.u8(s.coyote); w.u8(s.jumpLock); }
        if (mask & kQFuel) w.u8(s.fuel);
        if (mask & kQHands) { w.u8(s.action); w.u8(s.phase); }
    }
}
bool NetRead(NetReader& r, MsgPlayerSnapshot& m, const std::function<bool(uint32_t, uint32_t, PlayerStateQ&)>& lookup) {
    m.seq = r.u32(); m.tick = r.u32();
    const uint64_t n = r.varu();
    if (n > 256) r.ok = false;
    m.entries.clear();
    for (uint64_t i = 0; i < n && r.ok; ++i) {
        MsgPlayerSnapshot::Entry e;
        const uint32_t session = (uint32_t)r.varu();
        const uint64_t age = r.varu();
        const uint64_t mask = r.varu();
        if (mask > kQAll || age > m.seq || (age == 0 && mask != kQAll)) { r.ok = false; break; }
        PlayerStateQ b;
        if (age) { e.baseSeq = m.seq - (uint32_t)age; if (lookup && !lookup(session, e.baseSeq, b)) { r.ok = false; break; } e.base = b; }
        PlayerStateQ s = b; s.session = session;
        if (mask & kQAck) s.ackInput = (uint32_t)((int64_t)b.ackInput + r.vari());
        if (mask & kQBody) s.body = r.varu();
        if (mask & kQPos) for (int k = 0; k < 3; ++k) s.p[k] = b.p[k] + r.vari();
        if (mask & kQVel) for (int k = 0; k < 3; ++k) s.v[k] = b.v[k] + r.vari();
        if (mask & kQYaw) s.yaw = r.u16();
        if (mask & kQPitch) s.pitch = r.u16();
        if (mask & kQFlags) s.flags = r.u8();
        if (mask & kQTimers) { s.jumpBuffer = r.u8(); s.coyote = r.u8(); s.jumpLock = r.u8(); }
        if (mask & kQFuel) s.fuel = r.u8();
        if (mask & kQHands) { s.action = r.u8(); s.phase = r.u8(); if (s.action > 3) r.ok = false; }
        e.s = s;
        m.entries.push_back(e);
    }
    return Done(r);
}

void NetWrite(NetWriter& w, const MsgTerrainEdit& m) {
    Id(w, NetMsg::TerrainEdit);
    w.varu(m.body); w.pos(m.center); w.f32(m.radius); w.u8(m.flags); w.u8((uint8_t)m.paint);
    if (m.flags & 4) { w.f32(m.hx); w.f32(m.hy); w.f32(m.hz); w.unit(m.ax[0], m.ax[1], m.ax[2]); w.unit(m.ay[0], m.ay[1], m.ay[2]); }
    w.f32(m.strength); w.varu(m.session); w.u32(m.seq);
}
bool NetRead(NetReader& r, MsgTerrainEdit& m) {
    m.body = r.varu(); m.center = r.pos(); m.radius = r.f32(); m.flags = r.u8(); m.paint = (int8_t)r.u8();
    if (m.flags & 4) { m.hx = r.f32(); m.hy = r.f32(); m.hz = r.f32(); r.unit(m.ax[0], m.ax[1], m.ax[2]); r.unit(m.ay[0], m.ay[1], m.ay[2]); }
    m.strength = r.f32();
    const uint64_t s = r.varu();
    m.session = (uint32_t)s; m.seq = r.u32();
    if (m.radius < 0.0f || m.radius > 64.0f || (m.flags & ~0x1F) /* 16: a till (#1085, #1365) */ || !(std::fabs(m.strength) <= 4.0f) || s > 0xFFFFFFFFull) r.ok = false; // a dig is metres, not a planet
    return Done(r);
}

const char* NetEditRefusalText(NetEditRefusal r) {
    switch (r) {
        case NetEditRefusal::Reach: return "out of reach"; case NetEditRefusal::Rate: return "too fast";
        case NetEditRefusal::Tool: return "the tool in hand can't dig that"; case NetEditRefusal::Cost: return "not enough to place";
        case NetEditRefusal::Body: return "no such body"; case NetEditRefusal::Mantle: return "the mantle can't be dug";
        case NetEditRefusal::Claim: return "another team has claimed this land";
        case NetEditRefusal::Shape: return "not a stroke the tool makes";
        default: return "?";
    }
}
void NetWrite(NetWriter& w, const MsgEditRequest& m) {
    Id(w, NetMsg::EditRequest);
    w.u32(m.seq); w.u8((uint8_t)m.kind); w.varu(m.body); w.pos(m.hit); w.pos(m.center); w.f32(m.radius); w.str(m.item.substr(0, 48));
}
bool NetRead(NetReader& r, MsgEditRequest& m) {
    m.seq = r.u32();
    const uint8_t k = r.u8();
    if (k >= (uint8_t)NetEditKind::Count) r.ok = false;
    m.kind = (NetEditKind)k; m.body = r.varu(); m.hit = r.pos(); m.center = r.pos(); m.radius = r.f32(); m.item = r.str(48);
    if (!(m.radius > 0.0f && m.radius <= 16.0f)) r.ok = false;
    return Done(r);
}
void NetWrite(NetWriter& w, const MsgEditRefused& m) { Id(w, NetMsg::EditRefused); w.u32(m.seq); w.u8((uint8_t)m.why); }
bool NetRead(NetReader& r, MsgEditRefused& m) {
    m.seq = r.u32();
    const uint8_t why = r.u8();
    if (why >= (uint8_t)NetEditRefusal::Count) r.ok = false;
    m.why = (NetEditRefusal)why;
    return Done(r);
}
void NetWrite(NetWriter& w, const MsgInventoryDelta& m) { Id(w, NetMsg::InventoryDelta); w.str(m.item.substr(0, 48)); w.f64(m.amount); }
void NetWrite(NetWriter& w, const MsgChat& m) { Id(w, NetMsg::Chat); w.u32(m.from); w.str(m.text.substr(0, 256)); }
bool NetRead(NetReader& r, MsgChat& m) { m.from = r.u32(); m.text = r.str(256); return Done(r); }

// ---- Join sync (#776) ----
void NetWrite(NetWriter& w, const MsgSyncIndex& m) {
    Id(w, NetMsg::SyncIndex);
    w.u32(m.snapshot); w.u16(m.page); w.u16(m.pages); w.varu(m.total);
    const size_t n = std::min(m.entries.size(), kSyncIndexPage);
    w.varu(n);
    for (size_t i = 0; i < n; ++i) { const NetSyncEntry& e = m.entries[i]; w.u8((uint8_t)e.kind); w.varu(e.body); w.u64(e.key); w.u64(e.hash); w.varu(e.size); }
}
bool NetRead(NetReader& r, MsgSyncIndex& m) {
    m.snapshot = r.u32(); m.page = r.u16(); m.pages = r.u16();
    const uint64_t total = r.varu(), n = r.varu();
    if (total > kSyncMaxEntries || n > kSyncIndexPage || m.pages == 0 || m.page >= m.pages) r.ok = false;
    m.total = (uint32_t)total;
    m.entries.clear();
    for (uint64_t i = 0; i < n && r.ok; ++i) {
        NetSyncEntry e;
        const uint8_t kind = r.u8();
        if (kind >= (uint8_t)NetSyncKind::Count) r.ok = false;
        e.kind = (NetSyncKind)kind; e.body = r.varu(); e.key = r.u64(); e.hash = r.u64();
        const uint64_t size = r.varu();
        if (size > (64u << 20)) r.ok = false; // a region or a plot set: megabytes at most
        e.size = (uint32_t)size;
        m.entries.push_back(e);
    }
    return Done(r);
}
void NetWrite(NetWriter& w, const MsgSyncWant& m) { Id(w, NetMsg::SyncWant); w.u32(m.snapshot); w.varu(m.total); w.blob(m.bits.data(), m.bits.size()); }
bool NetRead(NetReader& r, MsgSyncWant& m) {
    m.snapshot = r.u32();
    const uint64_t total = r.varu();
    if (total > kSyncMaxEntries) r.ok = false;
    m.total = (uint32_t)total;
    m.bits = r.blob(kSyncMaxEntries / 8);
    if (m.bits.size() != (size_t)((m.total + 7) / 8)) r.ok = false;
    return Done(r);
}
void NetWrite(NetWriter& w, const MsgSyncData& m) {
    Id(w, NetMsg::SyncData);
    w.u32(m.snapshot); w.varu(m.entry); w.varu(m.offset); w.varu(m.size); w.blob(m.bytes.data(), m.bytes.size());
}
bool NetRead(NetReader& r, MsgSyncData& m) {
    m.snapshot = r.u32();
    const uint64_t entry = r.varu(), offset = r.varu(), size = r.varu();
    if (entry >= kSyncMaxEntries || size > (64u << 20) || offset > size) r.ok = false;
    m.entry = (uint32_t)entry; m.offset = (uint32_t)offset; m.size = (uint32_t)size;
    m.bytes = r.blob(kSyncPiece);
    if ((uint64_t)m.offset + m.bytes.size() > m.size) r.ok = false;
    return Done(r);
}
void NetWrite(NetWriter& w, const MsgSyncDone& m) { Id(w, NetMsg::SyncDone); w.u32(m.snapshot); w.varu(m.entries); w.varu(m.bytes); }
bool NetRead(NetReader& r, MsgSyncDone& m) {
    m.snapshot = r.u32();
    const uint64_t e = r.varu();
    if (e > kSyncMaxEntries) r.ok = false;
    m.entries = (uint32_t)e; m.bytes = r.varu();
    return Done(r);
}
bool NetRead(NetReader& r, MsgInventoryDelta& m) {
    m.item = r.str(48); m.amount = r.f64();
    if (m.item.empty() || !(std::fabs(m.amount) < 1e6)) r.ok = false;
    return Done(r);
}

// ---- Building (#330) ----
const char* NetPartRefusalText(NetPartRefusal r) {
    switch (r) {
        case NetPartRefusal::Reach: return "out of reach"; case NetPartRefusal::Cost: return "not enough to build it";
        case NetPartRefusal::Blocked: return "something is in the way"; case NetPartRefusal::NoPlot: return "no such plot";
        case NetPartRefusal::NoPart: return "no such part"; case NetPartRefusal::Body: return "no such body";
        case NetPartRefusal::Claim: return "another team has claimed this land";
        case NetPartRefusal::Kind: return "not a part players build"; case NetPartRefusal::Rate: return "too fast";
        default: return "?";
    }
}
static void WritePart(NetWriter& w, const NetPart& p) { w.u8(p.type); w.vari(p.x); w.vari(p.y); w.vari(p.z); w.u8(p.orient); w.u8(p.material); }
static void ReadPart(NetReader& r, NetPart& p) {
    p.type = r.u8();
    const int64_t x = r.vari(), y = r.vari(), z = r.vari();
    if (std::llabs(x) > 1000000 || std::llabs(y) > 1000000 || std::llabs(z) > 1000000) r.ok = false; // cells: half a metre each
    p.x = (int32_t)x; p.y = (int32_t)y; p.z = (int32_t)z;
    p.orient = r.u8(); p.material = r.u8();
    if (p.orient >= 24) r.ok = false;
}
static void ReadAnchor(NetReader& r, float a[3]) { a[0] = r.f32(); a[1] = r.f32(); a[2] = r.f32(); for (int i = 0; i < 3; ++i) if (!(std::fabs(a[i]) < 1e8f)) r.ok = false; }
void NetWrite(NetWriter& w, const MsgPartRequest& m) {
    Id(w, NetMsg::PartRequest);
    w.u32(m.seq); w.u8((uint8_t)m.op); w.varu(m.body); w.varu(m.plot);
    if (m.op == NetPartOp::Place) { if (m.plot == 0) { w.f32(m.anchor[0]); w.f32(m.anchor[1]); w.f32(m.anchor[2]); } WritePart(w, m.part); }
    else w.varu(m.partIndex);
}
bool NetRead(NetReader& r, MsgPartRequest& m) {
    m.seq = r.u32();
    const uint8_t op = r.u8();
    if (op > (uint8_t)NetPartOp::Remove) r.ok = false;
    m.op = (NetPartOp)op; m.body = r.varu();
    const uint64_t plot = r.varu();
    if (plot > 0xFFFFFFFFull) r.ok = false;
    m.plot = (uint32_t)plot;
    if (m.op == NetPartOp::Place) { if (m.plot == 0) ReadAnchor(r, m.anchor); ReadPart(r, m.part); }
    else {
        const uint64_t i = r.varu();
        if (i > 0xFFFFFFFFull || m.plot == 0) r.ok = false;
        m.partIndex = (uint32_t)i;
    }
    return Done(r);
}
void NetWrite(NetWriter& w, const MsgPartRefused& m) { Id(w, NetMsg::PartRefused); w.u32(m.seq); w.u8((uint8_t)m.why); }
bool NetRead(NetReader& r, MsgPartRefused& m) {
    m.seq = r.u32();
    const uint8_t why = r.u8();
    if (why >= (uint8_t)NetPartRefusal::Count) r.ok = false;
    m.why = (NetPartRefusal)why;
    return Done(r);
}
void NetWrite(NetWriter& w, const MsgPlotEdit& m) {
    Id(w, NetMsg::PlotEdit);
    w.varu(m.session); w.u32(m.seq); w.u8((uint8_t)m.op); w.varu(m.body); w.varu(m.plot);
    if (m.op == NetPartOp::NewPlot) { w.f32(m.anchor[0]); w.f32(m.anchor[1]); w.f32(m.anchor[2]); }
    else if (m.op == NetPartOp::Place) WritePart(w, m.part);
    else w.varu(m.partIndex);
}
bool NetRead(NetReader& r, MsgPlotEdit& m) {
    const uint64_t s = r.varu();
    m.session = (uint32_t)s; m.seq = r.u32();
    const uint8_t op = r.u8();
    if (op >= (uint8_t)NetPartOp::Count || s > 0xFFFFFFFFull) r.ok = false;
    m.op = (NetPartOp)op; m.body = r.varu();
    const uint64_t plot = r.varu();
    if (plot == 0 || plot > 0xFFFFFFFFull) r.ok = false;
    m.plot = (uint32_t)plot;
    if (m.op == NetPartOp::NewPlot) ReadAnchor(r, m.anchor);
    else if (m.op == NetPartOp::Place) ReadPart(r, m.part);
    else { const uint64_t i = r.varu(); if (i > 0xFFFFFFFFull) r.ok = false; m.partIndex = (uint32_t)i; }
    return Done(r);
}

void NetWrite(NetWriter& w, const MsgShower& m) { Id(w, NetMsg::Shower); w.f64(m.start); w.f64(m.end); w.u32(m.id); }
bool NetRead(NetReader& r, MsgShower& m) { m.start = r.f64(); m.end = r.f64(); m.id = r.u32(); return Done(r); }

constexpr size_t kMaxRoster = 64, kMaxNameLen = 32, kMaxItemLen = 48;
void NetWrite(NetWriter& w, const MsgRoster& m) {
    Id(w, NetMsg::Roster);
    const size_t n = std::min(m.players.size(), kMaxRoster);
    w.varu(n);
    for (size_t i = 0; i < n; ++i) { w.varu(m.players[i].session); w.str(m.players[i].name.substr(0, kMaxNameLen)); }
}
bool NetRead(NetReader& r, MsgRoster& m) {
    const uint64_t n = r.varu();
    if (n > kMaxRoster) r.ok = false;
    m.players.clear();
    for (uint64_t i = 0; i < n && r.ok; ++i) {
        MsgRoster::Entry e;
        const uint64_t s = r.varu();
        if (s == 0 || s > 0xFFFFFFFFull) r.ok = false;
        e.session = (uint32_t)s; e.name = r.str(kMaxNameLen);
        m.players.push_back(e);
    }
    return Done(r);
}
void NetWrite(NetWriter& w, const MsgLoadout& m) {
    Id(w, NetMsg::Loadout);
    w.varu(m.session); w.str(m.held.substr(0, kMaxItemLen)); w.boolean(m.headlamp); w.u8(m.suitTier);
    w.u8((uint8_t)(m.suitColor >> 16)); w.u8((uint8_t)(m.suitColor >> 8)); w.u8((uint8_t)m.suitColor); // (#341, last)
}
bool NetRead(NetReader& r, MsgLoadout& m) {
    const uint64_t s = r.varu();
    if (s > 0xFFFFFFFFull) r.ok = false;
    m.session = (uint32_t)s; m.held = r.str(kMaxItemLen); m.headlamp = r.boolean(); m.suitTier = r.u8();
    const uint32_t cr = r.u8(), cg = r.u8(), cb = r.u8();
    m.suitColor = (cr << 16) | (cg << 8) | cb;
    return Done(r);
}

NetMsg NetPeek(const std::vector<uint8_t>& b) {
    if (b.empty() || b[0] == 0 || b[0] >= (uint8_t)NetMsg::Count) return NetMsg::Count;
    return (NetMsg)b[0];
}

// ---- Test ----

namespace {
template <class M> bool RoundTrip(const M& m, M& out, std::vector<uint8_t>* bytes = nullptr) {
    NetWriter w; NetWrite(w, m);
    if (bytes) *bytes = w.bytes;
    NetReader r(w.bytes.data() + 1, w.bytes.size() - 1);
    return NetRead(r, out);
}
// Parses any message with its reader (as the receive loop will).
bool ParseAny(const std::vector<uint8_t>& b) {
    const NetMsg id = NetPeek(b);
    if (id == NetMsg::Count) return false;
    NetReader r(b.data() + 1, b.size() - 1);
    switch (id) {
        case NetMsg::Hello: { MsgHello m; return NetRead(r, m); }
        case NetMsg::Welcome: { MsgWelcome m; return NetRead(r, m); }
        case NetMsg::Reject: { MsgReject m; return NetRead(r, m); }
        case NetMsg::Bye: { MsgBye m; return NetRead(r, m); }
        case NetMsg::Ping: { MsgPing m; return NetRead(r, m); }
        case NetMsg::Pong: { MsgPong m; return NetRead(r, m); }
        case NetMsg::Input: { MsgInputs m; return NetRead(r, m); }
        case NetMsg::PlayerState: { MsgPlayerState m; return NetRead(r, m); }
        case NetMsg::TerrainEdit: { MsgTerrainEdit m; return NetRead(r, m); }
        case NetMsg::Shower: { MsgShower m; return NetRead(r, m); }
        case NetMsg::Roster: { MsgRoster m; return NetRead(r, m); }
        case NetMsg::Loadout: { MsgLoadout m; return NetRead(r, m); }
        case NetMsg::EditRequest: { MsgEditRequest m; return NetRead(r, m); }
        case NetMsg::EditRefused: { MsgEditRefused m; return NetRead(r, m); }
        case NetMsg::InventoryDelta: { MsgInventoryDelta m; return NetRead(r, m); }
        case NetMsg::PartRequest: { MsgPartRequest m; return NetRead(r, m); }
        case NetMsg::PartRefused: { MsgPartRefused m; return NetRead(r, m); }
        case NetMsg::PlotEdit: { MsgPlotEdit m; return NetRead(r, m); }
        case NetMsg::Chat: { MsgChat m; return NetRead(r, m); }
        case NetMsg::SyncIndex: { MsgSyncIndex m; return NetRead(r, m); }
        case NetMsg::SyncWant: { MsgSyncWant m; return NetRead(r, m); }
        case NetMsg::SyncData: { MsgSyncData m; return NetRead(r, m); }
        case NetMsg::SyncDone: { MsgSyncDone m; return NetRead(r, m); }
        case NetMsg::PlayerSnapshot: { MsgPlayerSnapshot m; return NetRead(r, m); }
        case NetMsg::Animals: { MsgAnimals m; return NetRead(r, m); }
        case NetMsg::AnimalHit: { MsgAnimalHit m; return NetRead(r, m); }
        case NetMsg::Damage: { MsgDamage m; return NetRead(r, m); }
        case NetMsg::Vitals: { MsgVitals m; return NetRead(r, m); }
        case NetMsg::UseItem: { MsgUseItem m; return NetRead(r, m); }
        case NetMsg::BedUse: { MsgBedUse m; return NetRead(r, m); }
        case NetMsg::Graves: { MsgGraves m; return NetRead(r, m); }
        case NetMsg::QuestBook: { MsgQuestBook m; return NetRead(r, m); }
        case NetMsg::CraftRequest: { MsgCraftRequest m; return NetRead(r, m); }
        case NetMsg::MachineAction: { MsgMachineAction m; return NetRead(r, m); }
        case NetMsg::MachineView: { MsgMachineView m; return NetRead(r, m); }
        case NetMsg::MachineStates: { MsgMachineStates m; return NetRead(r, m); }
        case NetMsg::TeamState: { MsgTeamState m; return NetRead(r, m); }
        case NetMsg::PlantAction: { MsgPlantAction m; return NetRead(r, m); }
        case NetMsg::Plants: { MsgPlants m; return NetRead(r, m); }
        case NetMsg::FarmFields: { MsgFarmFields m; return NetRead(r, m); }
        case NetMsg::DropItems: { MsgDropItems m; return NetRead(r, m); }
        case NetMsg::CreativeGive: { MsgCreativeGive m; return NetRead(r, m); }
        case NetMsg::Claims: { MsgClaims m; return NetRead(r, m); }
        case NetMsg::SkillState: { MsgSkillState m; return NetRead(r, m); }
        case NetMsg::SkillSpend: { MsgSkillSpend m; return NetRead(r, m); }
        default: return false;
    }
}
}

bool NetProto_Test() {
    bool ok = true;
    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: netproto: %s: NO\n", what); ok = false; } };
    std::vector<std::vector<uint8_t>> samples; // each message's bytes, for the truncation and fuzz passes
    {
        MsgHello a{ kNetProtocol, "m84-integration abc1234", "Dewey", "hunter2", "0f0f0f0f-1234-4abc-8def-0123456789ab" }, b; std::vector<uint8_t> s;
        check(RoundTrip(a, b, &s) && b.protocol == a.protocol && b.build == a.build && b.name == a.name && b.password == a.password && b.uuid == a.uuid, "hello");
        samples.push_back(s);
        MsgHello old{ 1, "x", "y", "" }, ob;
        check(RoundTrip(old, ob) && ob.protocol == 1, "an old protocol's hello still reads its version (to reject it clearly)");
    }
    {
        MsgWelcome a; a.session = 7; a.reconnected = true; a.seed = 42; a.bodies = { "planet Terra 0 0 0 60000", "moon Luna 1 2 3 9000" };
        a.fingerprint = 0x1234567890abcdefull; a.serverTick = 123456; a.worldTime = 4321.5; a.worldgen = "ocean 0.5\nkarst 0\n"; a.sandbox = true; a.motd = "be nice"; a.settings = { 2, true, false, 16 };
        MsgWelcome b; std::vector<uint8_t> s;
        check(RoundTrip(a, b, &s) && b.session == 7 && b.reconnected && b.seed == 42 && b.bodies == a.bodies && b.fingerprint == a.fingerprint &&
              b.serverTick == 123456 && b.worldTime == 4321.5 && b.worldgen == a.worldgen && b.sandbox && b.motd == "be nice" && b.settings.difficulty == 2 && b.settings.pvp && !b.settings.claims && b.settings.maxPlayers == 16, "welcome");
        samples.push_back(s);
    }
    { MsgReject a{ NetReject::Password, "try again" }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.reason == a.reason && b.detail == a.detail, "reject"); samples.push_back(s); }
    { MsgBye a{ NetBye::ServerClosed }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.reason == a.reason, "bye"); samples.push_back(s); }
    { MsgPing a{ 99, 12.5 }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.seq == 99 && b.clientTime == 12.5, "ping"); samples.push_back(s); }
    { MsgPong a{ 99, 12.5, 777, 86400.25 }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.seq == 99 && b.clientTime == 12.5 && b.serverTick == 777 && b.worldTime == 86400.25, "pong"); samples.push_back(s); }
    {
        MsgInputs a, b; std::vector<uint8_t> s;
        a.steps = { { 5000, 0x1234ull, 0x10ull, 0.6f, 0.0f, -0.8f }, { 4999, 0x1234ull, 0ull, 0.0f, 1.0f, 0.0f }, { 4998, 0x1ull, 0ull, 1.0f, 0.0f, 0.0f } };

        a.action = 2; a.actionPhase = 0.5f;
        const bool rt = RoundTrip(a, b, &s);
        check(rt && b.steps.size() == 3 && b.steps[0].tick == 5000 && b.steps[2].tick == 4998 && b.steps[0].held == 0x1234ull && b.action == 2 && fabsf(b.actionPhase - 0.5f) < 3e-3f && b.steps[0].pressed == 0x10ull &&
              b.steps[2].held == 1ull && fabsf(b.steps[0].fx - 0.6f) < 1e-3f && fabsf(b.steps[0].fz + 0.8f) < 1e-3f && fabsf(b.steps[1].fy - 1.0f) < 1e-3f, "inputs (3 steps)");
        samples.push_back(s);
    }
    {
        MsgPlayerState a; a.session = 3; a.tick = 9001; a.ackInput = 8999; a.body = 0xABCDull; a.local = { 60012.3456, -123.4567, 0.0005 };
        a.vx = 4.21f; a.vy = -9.81f; a.vz = 0.0f; a.yaw = 6.0f; a.pitch = 1.2f; a.flags = 0b1011; a.jumpBuffer = 0.15f; a.jumpLock = 0.08f; a.jetpackFuel = 0.5f; a.action = 1; a.actionPhase = 0.25f;
        MsgPlayerState b; std::vector<uint8_t> s;
        const bool rt = RoundTrip(a, b, &s);
        const double dp = std::fmax(std::fabs(b.local.x - a.local.x), std::fmax(std::fabs(b.local.y - a.local.y), std::fabs(b.local.z - a.local.z)));
        check(rt && b.session == 3 && b.tick == 9001 && b.ackInput == 8999 && b.body == a.body && dp <= 0.0005 + 1e-9 &&
              fabsf(b.vx - a.vx) < 0.006f && fabsf(b.vy - a.vy) < 0.006f && fabsf(b.yaw - a.yaw) < 1e-3f && fabsf(b.pitch - a.pitch) < 1e-3f && b.flags == a.flags && fabsf(b.jumpBuffer - 0.15f) < 1e-3f && fabsf(b.jumpLock - 0.08f) < 1e-3f && b.coyote == 0.0f && fabsf(b.jetpackFuel - 0.5f) < 3e-3f && b.action == 1 && fabsf(b.actionPhase - 0.25f) < 3e-3f,
              "player state (1 mm, 1 cm/s, 1e-3 rad at 60 km)");
        printf("bench: netproto: a player state is %zu bytes (position error %.4f mm at 60 km)\n", s.size(), dp * 1000.0);
        samples.push_back(s);
    }
    {
        MsgTerrainEdit a; a.body = 1; a.center = { 100.25, 59950.5, -3.75 }; a.radius = 1.5f; a.flags = 1 | 4; a.hx = 0.25f; a.hy = 0.5f; a.hz = 0.25f;
        a.ay[0] = 0.0f; a.ay[1] = 0.6f; a.ay[2] = 0.8f; a.paint = 3; a.strength = 0.6f; a.session = 4; a.seq = 77;
        MsgTerrainEdit b; std::vector<uint8_t> s;
        check(RoundTrip(a, b, &s) && b.body == 1 && std::fabs(b.center.y - a.center.y) < 0.001 && b.radius == a.radius && b.flags == a.flags &&
              b.hy == a.hy && fabsf(b.ay[2] - 0.8f) < 1e-3f && b.paint == 3 && b.strength == 0.6f && b.session == 4 && b.seq == 77, "terrain edit (box)");
        samples.push_back(s);
    }
    { MsgShower a{ 3600.5, 3900.25, 17 }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.start == a.start && b.end == a.end && b.id == 17, "shower"); samples.push_back(s); }
    {
        MsgEditRequest a, b; a.seq = 9; a.kind = NetEditKind::Place; a.body = 2; a.hit = { 10.5, 60000.25, -2.0 }; a.center = { 10.5, 60000.5, -2.0 }; a.radius = 0.65f; a.item = "dirt";
        std::vector<uint8_t> s;
        check(RoundTrip(a, b, &s) && b.seq == 9 && b.kind == NetEditKind::Place && b.body == 2 && std::fabs(b.center.y - a.center.y) < 0.001 && b.radius == a.radius && b.item == "dirt", "edit request");
        samples.push_back(s);
    }
    { MsgEditRefused a{ 9, NetEditRefusal::Cost }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.seq == 9 && b.why == NetEditRefusal::Cost, "edit refused"); samples.push_back(s); }
    { MsgInventoryDelta a{ "gravel", -0.375 }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.item == "gravel" && b.amount == -0.375, "inventory delta"); samples.push_back(s); }
    {
        MsgPartRequest a, b; a.seq = 4; a.op = NetPartOp::Place; a.body = 2; a.plot = 0; a.anchor[0] = 12.5f; a.anchor[1] = 59991.123f; a.anchor[2] = -3.0f;
        a.part = { 1, -3, 7, 120, 5, 2 };
        std::vector<uint8_t> s;
        check(RoundTrip(a, b, &s) && b.seq == 4 && b.op == NetPartOp::Place && b.plot == 0 && b.anchor[1] == a.anchor[1] && b.part.type == 1 && b.part.x == -3 && b.part.z == 120 &&
              b.part.orient == 5 && b.part.material == 2, "part request (a new plot)");
        samples.push_back(s);
        MsgPartRequest c, d; c.seq = 5; c.op = NetPartOp::Remove; c.body = 2; c.plot = 3; c.partIndex = 17;
        check(RoundTrip(c, d, &s) && d.op == NetPartOp::Remove && d.plot == 3 && d.partIndex == 17, "part request (remove)");
        samples.push_back(s);
    }
    { MsgPartRefused a{ 6, NetPartRefusal::Blocked }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.seq == 6 && b.why == NetPartRefusal::Blocked, "part refused"); samples.push_back(s); }
    {
        MsgPlotEdit a, b; a.session = 2; a.seq = 4; a.op = NetPartOp::NewPlot; a.body = 2; a.plot = 9; a.anchor[0] = 1.25f; a.anchor[1] = 60000.0625f; a.anchor[2] = 0.5f;
        std::vector<uint8_t> s;
        check(RoundTrip(a, b, &s) && b.session == 2 && b.op == NetPartOp::NewPlot && b.plot == 9 && b.anchor[1] == a.anchor[1], "plot edit (new plot)");
        samples.push_back(s);
        MsgPlotEdit c, d; c.session = 2; c.seq = 4; c.op = NetPartOp::Place; c.body = 2; c.plot = 9; c.part = { 9, 0, 1, 0, 3, 0 };
        check(RoundTrip(c, d, &s) && d.op == NetPartOp::Place && d.part.type == 9 && d.part.y == 1 && d.part.orient == 3, "plot edit (place)");
        samples.push_back(s);
    }
    { MsgChat a{ 3, "hello, world" }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.from == 3 && b.text == "hello, world", "chat"); samples.push_back(s); }
    {
        MsgSyncIndex a, b; a.snapshot = 7; a.page = 1; a.pages = 3; a.total = 600;
        a.entries = { { NetSyncKind::Region, 2, 0x123456789ABCull, 0xDEADBEEFCAFEull, 40000 }, { NetSyncKind::Plots, 2, 0, 42, 900 } };
        std::vector<uint8_t> s;
        check(RoundTrip(a, b, &s) && b.snapshot == 7 && b.page == 1 && b.pages == 3 && b.total == 600 && b.entries.size() == 2 && b.entries[0].key == a.entries[0].key &&
              b.entries[0].hash == a.entries[0].hash && b.entries[0].size == 40000 && b.entries[1].kind == NetSyncKind::Plots, "sync index");
        samples.push_back(s);
    }
    {
        MsgSyncWant a, b; a.snapshot = 7; a.total = 10; a.bits = { 0x05, 0x02 }; std::vector<uint8_t> s;
        check(RoundTrip(a, b, &s) && b.total == 10 && b.bits == a.bits, "sync want"); samples.push_back(s);
    }
    {
        MsgSyncData a, b; a.snapshot = 7; a.entry = 3; a.offset = 16384; a.size = 40000; a.bytes.assign(1000, 0xAB); std::vector<uint8_t> s;
        check(RoundTrip(a, b, &s) && b.entry == 3 && b.offset == 16384 && b.size == 40000 && b.bytes == a.bytes, "sync data"); samples.push_back(s);
    }
    { MsgSyncDone a{ 7, 12, 123456789 }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.entries == 12 && b.bytes == 123456789, "sync done"); samples.push_back(s); }
    { // #769: quantised states round-trip exactly; a snapshot with a full and a delta entry decodes to the same states
        MsgPlayerState ps; ps.session = 4; ps.tick = 900; ps.ackInput = 812; ps.body = 1; ps.local = { 12345.678, -60001.25, 7.5 };
        ps.vx = 3.21f; ps.vy = -0.5f; ps.vz = 6.0f; ps.yaw = 2.0f; ps.pitch = -0.3f; ps.flags = 1 | 32; ps.jetpackFuel = 0.5f; ps.action = 1; ps.actionPhase = 0.25f;
        const PlayerStateQ q0 = NetQuantize(ps);
        const PlayerStateQ q0b = NetQuantize(NetDequantize(q0, ps.tick));
        MsgPlayerState moved = NetDequantize(q0, ps.tick); moved.local.x += 0.2; moved.local.z -= 0.013; moved.ackInput += 2; moved.yaw += 0.1f;
        const PlayerStateQ q1 = NetQuantize(moved);
        MsgPlayerSnapshot a; a.seq = 50; a.tick = 902;
        a.entries.push_back({ q1, 0, {} });              // full
        a.entries.push_back({ q1, 47, q0 });             // against the state acknowledged in snapshot 47
        a.entries[1].s.session = 4;
        NetWriter w; NetWrite(w, a);
        NetReader r(w.bytes); r.u8();
        MsgPlayerSnapshot b;
        const bool ok = NetRead(r, b, [&](uint32_t s, uint32_t seq, PlayerStateQ& out) { if (s != 4 || seq != 47) return false; out = q0; return true; });
        size_t deltaBytes = 0; { NetWriter one; MsgPlayerSnapshot d = a; d.entries.erase(d.entries.begin()); NetWrite(one, d); deltaBytes = one.bytes.size(); }
        size_t fullBytes = 0; { NetWriter one; MsgPlayerSnapshot d = a; d.entries.pop_back(); NetWrite(one, d); fullBytes = one.bytes.size(); }
        check(q0b == q0 && ok && b.entries.size() == 2 && b.entries[0].s == q1 && b.entries[1].s == q1 && deltaBytes + 10 < fullBytes, "player snapshot (delta)");
        printf("bench: netproto: a walking player's state: %zu bytes full, %zu as a delta (in a snapshot message, header included)\n", fullBytes, deltaBytes);
        samples.push_back(w.bytes);
    }
    { // #334 animals
        MsgAnimals a; a.tick = 77;
        a.entries.push_back({ 0x2000000000000005ull, 3, 1, DVec3{ 100.25, -59980.5, 12.0 }, 1.25f, 2, 17.5f, 0.6f, 3 });
        a.removed = { 0x2000000000000009ull };
        MsgAnimals b; std::vector<uint8_t> s;
        check(RoundTrip(a, b, &s) && b.tick == 77 && b.entries.size() == 1 && b.entries[0].netId == a.entries[0].netId && b.entries[0].species == 3 &&
              fabs(b.entries[0].pos.y - a.entries[0].pos.y) < 0.001 && fabs(b.entries[0].health - 17.5f) < 0.06f && b.entries[0].flags == 3 && b.removed == a.removed, "animals");
        samples.push_back(s);
        MsgAnimalHit h{ 0x2000000000000005ull, 2 }, hb; std::vector<uint8_t> s2;
        check(RoundTrip(h, hb, &s2) && hb.netId == h.netId && hb.action == 2, "animal hit"); samples.push_back(s2);
        MsgDamage dm{ 8.0f, 2, true, 1, DVec3{ 10.0, 60000.0, -3.0 }, 3.0f, 0x2000000000000005ull }, db; std::vector<uint8_t> s3;
        check(RoundTrip(dm, db, &s3) && db.amount == 8.0f && db.kind == 2 && db.hasSource && db.knockback == 3.0f && db.sourceId == dm.sourceId && fabs(db.local.y - 60000.0) < 0.001, "damage"); samples.push_back(s3);
        MsgVitals va{ 72.5f, 40.25f, 99.0f, 0.0f, -3.5f, 1 | 8, 3 }, vb; std::vector<uint8_t> s4;
        check(RoundTrip(va, vb, &s4) && fabs(vb.health - 72.5f) < 0.006f && fabs(vb.food - 40.25f) < 0.006f && vb.temperature == -3.5f && vb.flags == 9 && vb.deaths == 3, "vitals"); samples.push_back(s4);
        MsgUseItem ua{ "berries" }, ub; std::vector<uint8_t> s5;
        check(RoundTrip(ua, ub, &s5) && ub.item == "berries", "use item"); samples.push_back(s5);
        MsgBedUse ba{ "terra", 7, { -3, 12, 40 } }, bb; std::vector<uint8_t> s6;
        check(RoundTrip(ba, bb, &s6) && bb.body == "terra" && bb.plot == 7 && bb.cell[0] == -3 && bb.cell[1] == 12 && bb.cell[2] == 40, "bed use"); samples.push_back(s6);
        MsgGraves ga, gb; ga.graves.push_back({ 7, "terra", DVec3{ 1.5, 60000.25, -3.0 }, "victim", "victim's grave 7", 2, 1795000000 }); std::vector<uint8_t> s7;
        check(RoundTrip(ga, gb, &s7) && gb.graves.size() == 1 && gb.graves[0].id == 7 && gb.graves[0].owner == "victim" && gb.graves[0].style == 2 && gb.graves[0].diedAt == 1795000000 && fabs(gb.graves[0].local.y - 60000.25) < 1e-6, "graves"); samples.push_back(s7);
        { // #865, #866: a dropped pile in the graves list; drops and creative items
            MsgGraves a, b; std::vector<uint8_t> s;
            MsgGraves::Entry e; e.id = 9; e.body = "terra"; e.dropped = true; e.items = { "stone", "oak_log" }; a.graves.push_back(e);
            check(RoundTrip(a, b, &s) && b.graves.size() == 1 && b.graves[0].dropped && b.graves[0].items.size() == 2 && b.graves[0].items[1] == "oak_log", "graves: a dropped pile"); samples.push_back(s);
            MsgDropItems da{ "stone", 250 }, db; std::vector<uint8_t> s2;
            check(RoundTrip(da, db, &s2) && db.item == "stone" && db.centi == 250, "drop items"); samples.push_back(s2);
            MsgCreativeGive ca{ "diamond", 64 }, cb; std::vector<uint8_t> s3;
            check(RoundTrip(ca, cb, &s3) && cb.item == "diamond" && cb.count == 64, "creative give"); samples.push_back(s3);
        }
        MsgQuestBook qa{ 3, "quest 2\nside_done e0f\nquest_count eat 1.00\n" }, qb; std::vector<uint8_t> s8;
        check(RoundTrip(qa, qb, &s8) && qb.team == 3 && qb.text == qa.text, "quest book"); samples.push_back(s8);
        { // #604
            MsgTeamState a, b; std::vector<uint8_t> s;
            a.team = 3; a.name = "alice's_team"; a.invite = "bob's_team";
            a.members = { { "alice", true, true, true, 87.5f, 7, 0 }, { "carol", false, false, false, 0.0f, 0, 172800 } };
            a.claimBudget = 5; a.allies = { "dave's_team" }; // (#828)
            check(RoundTrip(a, b, &s) && b.team == 3 && b.name == a.name && b.invite == a.invite && b.members.size() == 2 && b.members[0].owner && b.members[0].you &&
                  b.members[0].online && fabsf(b.members[0].health - 87.5f) < 0.01f && b.members[1].name == "carol" && !b.members[1].online &&
                  b.members[0].session == 7 && b.members[1].lastSeen == 172800 && b.claimBudget == 5 && b.allies.size() == 1 && b.allies[0] == "dave's_team", "team state");
            samples.push_back(s);
        }
        MsgCraftRequest ca{ 12, "oak_planks", 5, true, true }, cb; std::vector<uint8_t> s9;
        check(RoundTrip(ca, cb, &s9) && cb.recipe == 12 && cb.output == "oak_planks" && cb.count == 5 && cb.useChests && cb.recursive, "craft request (recursive: #946)"); samples.push_back(s9);
        MsgMachineAction ma{ 0x2000000000000003ull, 4, 17, 5, "", 0, 2, 1 }, mb; std::vector<uint8_t> s10;
        check(RoundTrip(ma, mb, &s10) && mb.body == ma.body && mb.plot == 4 && mb.part == 17 && mb.action == 5 && mb.slot == 2 && mb.which == 1, "machine action"); samples.push_back(s10);
        { // #862: into one chest slot
            MsgMachineAction a{ 7, 4, 17, 8, "stone", 12, 31, 0 }, b; std::vector<uint8_t> s;
            check(RoundTrip(a, b, &s) && b.action == 8 && b.item == "stone" && b.count == 12 && b.slot == 31, "machine action: into a slot"); samples.push_back(s);
        }
        MsgMachineView va2, vb2; va2.body = 9; va2.plot = 2; va2.part = 3; va2.view.status = MachineStatus::Working; va2.view.progress = 40; va2.view.duration = 120;
        va2.view.inSlots = { { "iron_ore", 3 }, { "", 0 } }; va2.view.outSlots = { { "iron_ingot", 2 } }; va2.view.hasFuelSlot = true; va2.view.fuelSlot = { "coal", 5 }; va2.view.hotOutputs = true; va2.view.canLight = true; va2.view.fire = kFireOut; va2.view.valid = true;
        std::vector<uint8_t> s11;
        check(RoundTrip(va2, vb2, &s11) && vb2.view.status == MachineStatus::Working && vb2.view.progress == 40 && vb2.view.inSlots.size() == 2 && vb2.view.inSlots[0].first == "iron_ore"
              && vb2.view.outSlots[0].second == 2 && vb2.view.fuelSlot.first == "coal" && vb2.view.hotOutputs && vb2.view.canLight && vb2.view.fire == kFireOut && vb2.view.valid, "machine view"); samples.push_back(s11);
        MsgMachineStates sa{ 9, 2, { 3, 5, 8 } }, sb; std::vector<uint8_t> s12;
        check(RoundTrip(sa, sb, &s12) && sb.working.size() == 3 && sb.working[2] == 8, "machine states"); samples.push_back(s12);
        MsgPlantAction pa{ 2, 0x2000000000000003ull, 0x0102030405060708ull, 12.25f, 60001.5f, -3.75f, "grain_seeds" }, pb; std::vector<uint8_t> s13;
        check(RoundTrip(pa, pb, &s13) && pb.action == 2 && pb.body == pa.body && pb.plant == pa.plant && pb.x == 12.25f && pb.y == 60001.5f && pb.z == -3.75f && pb.seed == "grain_seeds", "plant action"); samples.push_back(s13);
        MsgPlants la{ 5, 3, { { 0xFFEEDDCCBBAA9988ull, 12345.25 }, { 7, 0.5 } } }, lb; std::vector<uint8_t> s14;
        check(RoundTrip(la, lb, &s14) && lb.body == 5 && lb.session == 3 && lb.removed.size() == 2 && lb.removed[0].first == 0xFFEEDDCCBBAA9988ull && lb.removed[0].second == 12345.25, "plants"); samples.push_back(s14);
        MsgClaims ka, kb; NetClaim kc; kc.id = 7; kc.team = 3; kc.teamName = "rangers"; kc.colour = 0x3C8CE0u; kc.body = 0x2000000000000001ull; kc.centre[0] = 1.5f; kc.centre[1] = 60000.25f; kc.radius = 400.0f; kc.flags = 2; ka.claims.push_back(kc); std::vector<uint8_t> s16;
        check(RoundTrip(ka, kb, &s16) && kb.claims.size() == 1 && kb.claims[0].teamName == "rangers" && kb.claims[0].colour == 0x3C8CE0u && kb.claims[0].body == kc.body && kb.claims[0].centre[1] == 60000.25f && kb.claims[0].radius == 400.0f && kb.claims[0].flags == 2, "claims"); samples.push_back(s16);
        MsgFarmFields fa, fb; fa.body = 5; fa.fields = { { 1.5f, 60000.25f, -2.0f, "grain", 0.5f }, { 2.0f, 60000.3f, -2.0f, "roots", 1.0f } }; std::vector<uint8_t> s15;
        check(RoundTrip(fa, fb, &s15) && fb.fields.size() == 2 && fb.fields[1].crop == "roots" && fb.fields[1].growth == 1.0f && fabsf(fb.fields[0].growth - 0.5f) < 1e-4f && fb.fields[0].y == 60000.25f, "farm fields"); samples.push_back(s15);
    }
    { MsgRoster a, b; a.players = { { 1, "alice" }, { 7, "bob" } }; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.players.size() == 2 && b.players[1].session == 7 && b.players[1].name == "bob", "roster"); samples.push_back(s); }
    { MsgSkillState a{ "skill_xp mining 2478.5\nperk quick_pick 2\n" }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.text == a.text, "skill state (#920)"); samples.push_back(s); }
    { MsgSkillSpend a{ "quick_pick" }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.perk == "quick_pick", "skill spend (#920)"); samples.push_back(s); }
    { MsgLoadout a{ 3, "stone_pickaxe", true, 2, 0xD04A2Fu }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.session == 3 && b.held == "stone_pickaxe" && b.headlamp && b.suitTier == 2 && b.suitColor == 0xD04A2Fu, "loadout"); samples.push_back(s); }
    // Every message cut short anywhere must not read as valid; garbage must never crash a reader.
    int truncValid = 0, truncTotal = 0;
    for (const auto& s : samples)
        for (size_t n = 0; n < s.size(); ++n) { ++truncTotal; if (ParseAny(std::vector<uint8_t>(s.begin(), s.begin() + n))) ++truncValid; }
    char what[96];
    snprintf(what, sizeof what, "a message cut short reads as valid (%d of %d cuts)", truncValid, truncTotal);
    check(truncValid == 0, what);
    std::mt19937 rng(1234);
    int fuzzValid = 0;
    for (int i = 0; i < 20000; ++i) {
        std::vector<uint8_t> b(rng() % 64);
        for (auto& x : b) x = (uint8_t)rng();
        if (!b.empty() && (i & 1)) b[0] = (uint8_t)(1 + rng() % ((uint8_t)NetMsg::Count - 1)); // half with a real id
        if (ParseAny(b)) ++fuzzValid;
        // and mutate a real message: a flipped byte must parse or fail cleanly
        std::vector<uint8_t> m = samples[i % samples.size()];
        m[rng() % m.size()] ^= (uint8_t)(1u << (rng() % 8));
        ParseAny(m);
    }
    printf("bench: netproto: %zu message kinds round-trip; %d truncations, %d read as valid; 20000 random buffers (%d happened to parse) and 20000 bit flips read without a crash: %s\n",
           samples.size(), truncTotal, truncValid, fuzzValid, ok ? "PASS" : "FAIL");
    return ok;
}

// (#349) whether a message parses with its reader - the server's guard counts the ones that don't
bool NetWellFormed(const std::vector<uint8_t>& bytes) { return ParseAny(bytes); }

// ---- #334 animals ----
void NetWrite(NetWriter& w, const MsgAnimals& m) {
    Id(w, NetMsg::Animals);
    w.u32(m.tick);
    const size_t n = std::min<size_t>(m.entries.size(), 256), k = std::min<size_t>(m.removed.size(), 256);
    w.varu(n);
    for (size_t i = 0; i < n; ++i) {
        const MsgAnimals::Entry& e = m.entries[i];
        w.varu(e.netId); w.u8(e.species); w.varu(e.body); w.pos(e.pos); w.angle(e.heading); w.u8(e.state);
        w.u16((uint16_t)std::clamp((int)std::lround(e.health * 10.0f), 0, 65535)); w.u8((uint8_t)std::clamp((int)std::lround(e.age * 255.0f), 0, 255)); w.u8(e.flags);
    }
    w.varu(k);
    for (size_t i = 0; i < k; ++i) w.varu(m.removed[i]);
}
bool NetRead(NetReader& r, MsgAnimals& m) {
    m.tick = r.u32();
    const uint64_t n = r.varu();
    if (n > 256) r.ok = false;
    m.entries.clear();
    for (uint64_t i = 0; i < n && r.ok; ++i) {
        MsgAnimals::Entry e;
        e.netId = r.varu(); e.species = r.u8(); e.body = r.varu(); e.pos = r.pos(); e.heading = r.angle(); e.state = r.u8();
        e.health = r.u16() / 10.0f; e.age = r.u8() / 255.0f; e.flags = r.u8();
        if (e.state > 5) r.ok = false;
        m.entries.push_back(e);
    }
    const uint64_t k = r.varu();
    if (k > 256) r.ok = false;
    m.removed.clear();
    for (uint64_t i = 0; i < k && r.ok; ++i) m.removed.push_back(r.varu());
    return Done(r);
}
void NetWrite(NetWriter& w, const MsgAnimalHit& m) { Id(w, NetMsg::AnimalHit); w.varu(m.netId); w.u8(m.action); }
bool NetRead(NetReader& r, MsgAnimalHit& m) { m.netId = r.varu(); m.action = r.u8(); if (m.action > 2) r.ok = false; return Done(r); }
void NetWrite(NetWriter& w, const MsgDamage& m) {
    Id(w, NetMsg::Damage);
    w.f32(m.amount); w.u8(m.kind); w.boolean(m.hasSource); w.varu(m.body); w.pos(m.local); w.f32(m.knockback); w.varu(m.sourceId);
}
bool NetRead(NetReader& r, MsgDamage& m) {
    m.amount = r.f32(); m.kind = r.u8(); m.hasSource = r.boolean(); m.body = r.varu(); m.local = r.pos(); m.knockback = r.f32(); m.sourceId = r.varu();
    if (!(m.amount >= 0.0f && m.amount < 1000.0f) || !(m.knockback >= 0.0f && m.knockback < 100.0f)) r.ok = false;
    return Done(r);
}
// ---- #337 survival ----
void NetWrite(NetWriter& w, const MsgVitals& m) {
    Id(w, NetMsg::Vitals);
    auto q = [&](float v) { w.u16((uint16_t)std::clamp((int)std::lround(v * 100.0f), 0, 65535)); }; // 0..100 to 0.01
    q(m.health); q(m.food); q(m.oxygen); q(m.power); w.f32(m.temperature); w.u8(m.flags); w.u16(m.deaths);
}
bool NetRead(NetReader& r, MsgVitals& m) {
    m.health = r.u16() / 100.0f; m.food = r.u16() / 100.0f; m.oxygen = r.u16() / 100.0f; m.power = r.u16() / 100.0f;
    m.temperature = r.f32(); m.flags = r.u8(); m.deaths = r.u16();
    if (m.health > 100.01f || m.food > 100.01f || m.oxygen > 100.01f || m.power > 100.01f) r.ok = false;
    return Done(r);
}
void NetWrite(NetWriter& w, const MsgUseItem& m) { Id(w, NetMsg::UseItem); w.str(m.item.substr(0, 48)); }
bool NetRead(NetReader& r, MsgUseItem& m) { m.item = r.str(48); return Done(r) && !m.item.empty(); }
void NetWrite(NetWriter& w, const MsgBedUse& m) { Id(w, NetMsg::BedUse); w.str(m.body.substr(0, 64)); w.u32(m.plot); for (int i = 0; i < 3; ++i) w.vari(m.cell[i]); }
bool NetRead(NetReader& r, MsgBedUse& m) {
    m.body = r.str(64); m.plot = r.u32();
    for (int i = 0; i < 3; ++i) m.cell[i] = (int32_t)std::clamp<int64_t>(r.vari(), INT32_MIN, INT32_MAX);
    return Done(r) && !m.body.empty();
}
void NetWrite(NetWriter& w, const MsgGraves& m) {
    Id(w, NetMsg::Graves);
    const size_t n = std::min<size_t>(m.graves.size(), 256);
    w.varu(n);
    for (size_t i = 0; i < n; ++i) {
        const auto& g = m.graves[i];
        w.u32(g.id); w.str(g.body.substr(0, 64)); w.f64(g.local.x); w.f64(g.local.y); w.f64(g.local.z);
        w.str(g.owner.substr(0, 32)); w.str(g.name.substr(0, 64)); w.u8(g.style); w.vari(g.diedAt);
        const size_t k = std::min<size_t>(g.items.size(), 3); // (#865) a dropped pile and what it holds (its bundle's colours)
        w.boolean(g.dropped); w.u8((uint8_t)k);
        for (size_t j = 0; j < k; ++j) w.str(g.items[j].substr(0, 64));
    }
}
bool NetRead(NetReader& r, MsgGraves& m) {
    const uint64_t n = r.varu();
    if (n > 256) { r.ok = false; return false; }
    m.graves.clear();
    for (uint64_t i = 0; i < n && r.ok; ++i) {
        MsgGraves::Entry g;
        g.id = r.u32(); g.body = r.str(64); g.local.x = r.f64(); g.local.y = r.f64(); g.local.z = r.f64();
        g.owner = r.str(32); g.name = r.str(64); g.style = r.u8(); g.diedAt = r.vari();
        g.dropped = r.boolean();
        const uint8_t k = r.u8();
        if (k > 3) r.ok = false;
        for (uint8_t j = 0; j < k && r.ok; ++j) g.items.push_back(r.str(64));
        if (!(std::fabs(g.local.x) < 1e9 && std::fabs(g.local.y) < 1e9 && std::fabs(g.local.z) < 1e9) || g.style > 2) r.ok = false;
        m.graves.push_back(std::move(g));
    }
    return Done(r);
}
void NetWrite(NetWriter& w, const MsgQuestBook& m) { Id(w, NetMsg::QuestBook); w.u32(m.team); w.str(m.text.substr(0, 65536)); }
bool NetRead(NetReader& r, MsgQuestBook& m) { m.team = r.u32(); m.text = r.str(65536); return Done(r); }
void NetWrite(NetWriter& w, const MsgCraftRequest& m) { Id(w, NetMsg::CraftRequest); w.u32(m.recipe); w.str(m.output.substr(0, 64)); w.u8(m.count); w.boolean(m.useChests); w.boolean(m.recursive); }
bool NetRead(NetReader& r, MsgCraftRequest& m) { m.recipe = r.u32(); m.output = r.str(64); m.count = r.u8(); m.useChests = r.boolean(); m.recursive = r.boolean(); return Done(r) && !m.output.empty() && m.count >= 1 && m.count <= 64; }
// (#865, #866) drops and creative items
void NetWrite(NetWriter& w, const MsgDropItems& m) { Id(w, NetMsg::DropItems); w.str(m.item.substr(0, 64)); w.vari(m.centi); }
bool NetRead(NetReader& r, MsgDropItems& m) { m.item = r.str(64); m.centi = r.vari(); return Done(r) && !m.item.empty() && m.centi > 0 && m.centi <= 100000000; }
void NetWrite(NetWriter& w, const MsgCreativeGive& m) { Id(w, NetMsg::CreativeGive); w.str(m.item.substr(0, 64)); w.varu(m.count); }
bool NetRead(NetReader& r, MsgCreativeGive& m) { m.item = r.str(64); m.count = (uint32_t)std::min<uint64_t>(r.varu(), 1000000); return Done(r); }

// ---- Machines on the server (#332) ----
void NetWrite(NetWriter& w, const MsgMachineAction& m) {
    Id(w, NetMsg::MachineAction); w.varu(m.body); w.u32(m.plot); w.vari(m.part); w.u8(m.action); w.str(m.item.substr(0, 64)); w.u16(m.count); w.u8(m.slot); w.u8(m.which);
}
bool NetRead(NetReader& r, MsgMachineAction& m) {
    m.body = r.varu(); m.plot = r.u32(); const int64_t part = r.vari(); m.action = r.u8(); m.item = r.str(64); m.count = r.u16(); m.slot = r.u8(); m.which = r.u8();
    if (part < 0 || part > 1000000 || m.action < 1 || m.action > 10 || m.which > 2) r.ok = false; // (8, 9: #862; 10: a campfire, #1103 - widen with each new action, #1365)
    m.part = (int32_t)part;
    return Done(r);
}
namespace {
    constexpr size_t kMaxViewSlots = 64; // a chest has 45; machines a handful
    void WriteSlots(NetWriter& w, const std::vector<std::pair<std::string, int>>& s) {
        const size_t n = std::min(s.size(), kMaxViewSlots);
        w.varu(n);
        for (size_t i = 0; i < n; ++i) { w.str(s[i].first.substr(0, 64)); w.vari(s[i].second); }
    }
    void ReadSlots(NetReader& r, std::vector<std::pair<std::string, int>>& s) {
        const uint64_t n = r.varu();
        if (n > kMaxViewSlots) { r.ok = false; return; }
        s.clear();
        for (uint64_t i = 0; i < n && r.ok; ++i) { std::string id = r.str(64); const int64_t c = r.vari(); if (c < 0 || c > 100000000) r.ok = false; s.push_back({ id, (int)c }); }
    }
}
void NetWrite(NetWriter& w, const MsgMachineView& m) {
    Id(w, NetMsg::MachineView); w.varu(m.body); w.u32(m.plot); w.vari(m.part);
    const MachineView& v = m.view;
    w.u8((uint8_t)v.status); w.vari(v.progress); w.vari(v.duration); w.vari(v.energy); w.vari(v.energyMax);
    w.str(v.input.substr(0, 128)); w.str(v.output.substr(0, 128));
    WriteSlots(w, v.inSlots); WriteSlots(w, v.outSlots);
    w.boolean(v.hasFuelSlot); w.str(v.fuelSlot.first.substr(0, 64)); w.vari(v.fuelSlot.second);
    w.vari(v.water); w.vari(v.waterMax); w.vari(v.heat); w.vari(v.heatMax); w.boolean(v.hotOutputs); w.boolean(v.canLight); w.u8(v.fire); w.boolean(v.valid);
}
bool NetRead(NetReader& r, MsgMachineView& m) {
    m.body = r.varu(); m.plot = r.u32(); const int64_t part = r.vari();
    MachineView& v = m.view;
    const uint8_t st = r.u8();
    if (st > (uint8_t)MachineStatus::NoWater || part < 0 || part > 1000000) r.ok = false;
    v.status = (MachineStatus)st; m.part = (int32_t)part;
    v.progress = (int)r.vari(); v.duration = (int)r.vari(); v.energy = (int)r.vari(); v.energyMax = (int)r.vari();
    v.input = r.str(128); v.output = r.str(128);
    ReadSlots(r, v.inSlots); ReadSlots(r, v.outSlots);
    v.hasFuelSlot = r.boolean(); v.fuelSlot.first = r.str(64); v.fuelSlot.second = (int)r.vari();
    v.water = (int)r.vari(); v.waterMax = (int)r.vari(); v.heat = (int)r.vari(); v.heatMax = (int)r.vari(); v.hotOutputs = r.boolean(); v.canLight = r.boolean(); v.fire = r.u8(); v.valid = r.boolean();
    return Done(r);
}
void NetWrite(NetWriter& w, const MsgMachineStates& m) {
    Id(w, NetMsg::MachineStates); w.varu(m.body); w.u32(m.plot);
    const size_t n = std::min<size_t>(m.working.size(), 4096);
    w.varu(n);
    for (size_t i = 0; i < n; ++i) w.vari(m.working[i]);
}
bool NetRead(NetReader& r, MsgMachineStates& m) {
    m.body = r.varu(); m.plot = r.u32();
    const uint64_t n = r.varu();
    if (n > 4096) { r.ok = false; return false; }
    m.working.clear();
    for (uint64_t i = 0; i < n && r.ok; ++i) { const int64_t p = r.vari(); if (p < 0 || p > 1000000) r.ok = false; m.working.push_back((int32_t)p); }
    return Done(r);
}

// ---- Plants and farming (#335) ----
void NetWrite(NetWriter& w, const MsgPlantAction& m) {
    Id(w, NetMsg::PlantAction); w.u8(m.action); w.varu(m.body); w.u64(m.plant); w.f32(m.x); w.f32(m.y); w.f32(m.z); w.str(m.seed.substr(0, 64));
}
bool NetRead(NetReader& r, MsgPlantAction& m) {
    m.action = r.u8(); m.body = r.varu(); m.plant = r.u64(); m.x = r.f32(); m.y = r.f32(); m.z = r.f32(); m.seed = r.str(64);
    if (m.action < 1 || m.action > 3 || !std::isfinite(m.x) || !std::isfinite(m.y) || !std::isfinite(m.z)) r.ok = false;
    return Done(r);
}
namespace { constexpr size_t kMaxPlantsPerMsg = 4096, kMaxFarmFields = 4096; }
void NetWrite(NetWriter& w, const MsgPlants& m) {
    Id(w, NetMsg::Plants); w.varu(m.body); w.u32(m.session);
    const size_t n = std::min(m.removed.size(), kMaxPlantsPerMsg);
    w.varu(n);
    for (size_t i = 0; i < n; ++i) { w.u64(m.removed[i].first); w.f64(m.removed[i].second); }
}
bool NetRead(NetReader& r, MsgPlants& m) {
    m.body = r.varu(); m.session = r.u32();
    const uint64_t n = r.varu();
    if (n > kMaxPlantsPerMsg) { r.ok = false; return false; }
    m.removed.clear();
    for (uint64_t i = 0; i < n && r.ok; ++i) { const uint64_t id = r.u64(); const double t = r.f64(); m.removed.push_back({ id, t }); }
    return Done(r);
}
void NetWrite(NetWriter& w, const MsgFarmFields& m) {
    Id(w, NetMsg::FarmFields); w.varu(m.body);
    const size_t n = std::min(m.fields.size(), kMaxFarmFields);
    w.varu(n);
    for (size_t i = 0; i < n; ++i) {
        const NetFarmField& f = m.fields[i];
        w.f32(f.x); w.f32(f.y); w.f32(f.z); w.str(f.crop.substr(0, 64)); w.u16((uint16_t)lroundf(std::clamp(f.growth, 0.0f, 1.0f) * 65535.0f));
    }
}
bool NetRead(NetReader& r, MsgFarmFields& m) {
    m.body = r.varu();
    const uint64_t n = r.varu();
    if (n > kMaxFarmFields) { r.ok = false; return false; }
    m.fields.clear();
    for (uint64_t i = 0; i < n && r.ok; ++i) {
        NetFarmField f; f.x = r.f32(); f.y = r.f32(); f.z = r.f32(); f.crop = r.str(64); f.growth = r.u16() / 65535.0f;
        if (!std::isfinite(f.x) || !std::isfinite(f.y) || !std::isfinite(f.z) || f.crop.empty()) r.ok = false;
        m.fields.push_back(std::move(f));
    }
    return Done(r);
}
// ---- Teams (#604) ----
void NetWrite(NetWriter& w, const MsgTeamState& m) {
    Id(w, NetMsg::TeamState); w.u32((uint32_t)m.team); w.str(m.name.substr(0, 40)); w.str(m.invite.substr(0, 40));
    const size_t n = std::min<size_t>(m.members.size(), 16);
    w.u8((uint8_t)n);
    for (size_t i = 0; i < n; ++i) {
        const MsgTeamState::Member& x = m.members[i];
        w.str(x.name.substr(0, 32)); w.u8((uint8_t)((x.online ? 1 : 0) | (x.owner ? 2 : 0) | (x.you ? 4 : 0))); w.f32(x.health);
        w.varu(x.session); w.varu(x.lastSeen); // (#828)
    }
    w.u8((uint8_t)std::clamp(m.claimBudget, 0, 255)); // (#828)
    const size_t a = std::min<size_t>(m.allies.size(), 16);
    w.u8((uint8_t)a);
    for (size_t i = 0; i < a; ++i) w.str(m.allies[i].substr(0, 40));
}
bool NetRead(NetReader& r, MsgTeamState& m) {
    m.team = (int32_t)r.u32(); m.name = r.str(40); m.invite = r.str(40);
    const int n = r.u8();
    if (n > 16) return false;
    m.members.clear();
    for (int i = 0; i < n && r.ok; ++i) {
        MsgTeamState::Member x; x.name = r.str(32);
        const uint8_t f = r.u8(); x.online = f & 1; x.owner = f & 2; x.you = f & 4; x.health = r.f32();
        const uint64_t ses = r.varu(), seen = r.varu(); // (#828)
        if (ses > 0xFFFFFFFFull || seen > 0xFFFFFFFFull) r.ok = false;
        x.session = (uint32_t)ses; x.lastSeen = (uint32_t)seen;
        m.members.push_back(x);
    }
    m.claimBudget = r.u8();
    const int a = r.u8();
    if (a > 16) return false;
    m.allies.clear();
    for (int i = 0; i < a && r.ok; ++i) m.allies.push_back(r.str(40));
    return Done(r);
}

// ---- Land claims (#343) ----
namespace { constexpr size_t kMaxClaimsPerMsg = 4096; }
// (#920) the player's skills: the save's text lines (a few dozen at most)
void NetWrite(NetWriter& w, const MsgSkillState& m) { Id(w, NetMsg::SkillState); w.str(m.text.substr(0, 8192)); }
bool NetRead(NetReader& r, MsgSkillState& m) { m.text = r.str(8192); return Done(r); }
void NetWrite(NetWriter& w, const MsgSkillSpend& m) { Id(w, NetMsg::SkillSpend); w.str(m.perk.substr(0, 64)); }
bool NetRead(NetReader& r, MsgSkillSpend& m) { m.perk = r.str(64); return Done(r) && !m.perk.empty(); }

void NetWrite(NetWriter& w, const MsgClaims& m) {
    Id(w, NetMsg::Claims);
    const size_t n = std::min(m.claims.size(), kMaxClaimsPerMsg);
    w.varu(n);
    for (size_t i = 0; i < n; ++i) {
        const NetClaim& c = m.claims[i];
        w.u32(c.id); w.vari(c.team); w.str(c.teamName.substr(0, 40)); w.u32(c.colour); w.varu(c.body);
        w.f32(c.centre[0]); w.f32(c.centre[1]); w.f32(c.centre[2]); w.f32(c.radius); w.u8(c.flags);
    }
}
bool NetRead(NetReader& r, MsgClaims& m) {
    const uint64_t n = r.varu();
    if (n > kMaxClaimsPerMsg) { r.ok = false; return false; }
    m.claims.clear();
    for (uint64_t i = 0; i < n && r.ok; ++i) {
        NetClaim c; c.id = r.u32(); c.team = (int32_t)r.vari(); c.teamName = r.str(40); c.colour = r.u32(); c.body = r.varu();
        c.centre[0] = r.f32(); c.centre[1] = r.f32(); c.centre[2] = r.f32(); c.radius = r.f32(); c.flags = r.u8();
        if (c.radius < 0.0f || c.radius > 100000.0f || c.flags > 3) r.ok = false;
        m.claims.push_back(std::move(c));
    }
    return Done(r);
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#869): every message round-trips; truncated or random bytes never parse as valid or crash ----
#include "unit_helpers.h"

// The whole --test-netproto self-test (every message kind, all truncations, 20000 random buffers and bit flips). It stays
// the one body of checks while the flag lives (step 4 of #868 moves the flag's users here and removes it).
TEST(NetProto, SelfTestEveryMessageRoundTrips) {
    testing::internal::CaptureStdout();
    const bool ok = NetProto_Test();
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_TRUE(ok) << out;
}

TEST(NetProto, UnknownIdsAreRejected) {
    EXPECT_EQ(NetPeek({}), NetMsg::Count);
    EXPECT_EQ(NetPeek({ (uint8_t)NetMsg::Count }), NetMsg::Count);
    EXPECT_EQ(NetPeek({ 0xFF, 1, 2, 3 }), NetMsg::Count);
    EXPECT_FALSE(ParseAny({ 0xFF, 1, 2, 3 }));
}

TEST(NetProto, EveryIdAloneParsesWithoutACrash) {
    for (int id = 0; id < (int)NetMsg::Count; ++id) {
        std::vector<uint8_t> b{ (uint8_t)id };
        ParseAny(b); // an id with no body: must fail or pass cleanly, never read past the end
        b.resize(64, 0xFF);
        ParseAny(b);
    }
    SUCCEED();
}

// A player state holds its position to 1 mm anywhere within 100 km of a body's centre, whatever the sign (#769).
TEST(NetProto, PlayerStatePositionPrecisionEverywhere) {
    std::mt19937 rng(856);
    std::uniform_real_distribution<double> d(-100000.0, 100000.0);
    double worst = 0.0;
    for (int i = 0; i < 2000; ++i) {
        MsgPlayerState a; a.session = 1; a.tick = 100; a.body = 7; a.local = { d(rng), d(rng), d(rng) };
        MsgPlayerState b;
        ASSERT_TRUE(RoundTrip(a, b));
        worst = std::fmax(worst, std::fmax(std::fabs(b.local.x - a.local.x), std::fmax(std::fabs(b.local.y - a.local.y), std::fabs(b.local.z - a.local.z))));
    }
    EXPECT_LE(worst, 0.0005 + 1e-9);
}

// Quantising a dequantised state changes nothing (the server and client compare quantised states, #769).
TEST(NetProto, QuantizeIsIdempotent) {
    std::mt19937 rng(769);
    std::uniform_real_distribution<float> u(-1.0f, 1.0f);
    for (int i = 0; i < 1000; ++i) {
        MsgPlayerState s; s.session = 2; s.tick = 1000 + i; s.ackInput = 990 + i; s.body = 3;
        s.local = { 60000.0 * u(rng), 60000.0 * u(rng), 60000.0 * u(rng) };
        s.vx = 20.0f * u(rng); s.vy = 20.0f * u(rng); s.vz = 20.0f * u(rng);
        s.yaw = 3.0f * u(rng) + 3.0f; s.pitch = 1.5f * u(rng); s.flags = (uint8_t)(i & 63);
        s.jumpBuffer = 0.1f + 0.1f * u(rng); s.jetpackFuel = 0.5f + 0.5f * u(rng); s.action = (uint8_t)(i % 3); s.actionPhase = 0.5f + 0.5f * u(rng);
        const PlayerStateQ q = NetQuantize(s);
        ASSERT_TRUE(NetQuantize(NetDequantize(q, s.tick)) == q) << "state " << i;
    }
}

// Mutated real messages (one to eight flipped bits) parse or fail cleanly; a different seed from the self-test's.
TEST(NetProto, BitFlipsNeverCrash) {
    std::vector<std::vector<uint8_t>> samples;
    { MsgHello m{ kNetProtocol, "build", "name", "pw", "0f0f0f0f-1234-4abc-8def-0123456789ab" }; NetWriter w; NetWrite(w, m); samples.push_back(w.bytes); }
    { MsgChat m{ 3, "hello" }; NetWriter w; NetWrite(w, m); samples.push_back(w.bytes); }
    { MsgSyncData m; m.snapshot = 1; m.entry = 2; m.offset = 0; m.size = 300; m.bytes.assign(300, 7); NetWriter w; NetWrite(w, m); samples.push_back(w.bytes); }
    { MsgGraves m; m.graves.push_back({ 1, "terra", DVec3{ 1, 2, 3 }, "a", "b", 1, 5 }); NetWriter w; NetWrite(w, m); samples.push_back(w.bytes); }
    { MsgTeamState m; m.team = 1; m.name = "t"; m.members = { { "x", true, false, true, 50.0f } }; NetWriter w; NetWrite(w, m); samples.push_back(w.bytes); }
    std::mt19937 rng(20261118);
    for (int i = 0; i < 20000; ++i) {
        std::vector<uint8_t> m = samples[i % samples.size()];
        const int flips = 1 + (int)(rng() % 8);
        for (int k = 0; k < flips; ++k) m[rng() % m.size()] ^= (uint8_t)(1u << (rng() % 8));
        ParseAny(m);
    }
    SUCCEED();
}

// (#1365) The largest value each sender writes must read back: a read range left at an older maximum drops the message on
// the receiving side without a trace (a till's flag 16 on every client; a campfire's machine action 10 on the server).
// Widen a range here and in its NetRead together when a new flag, kind or action is added.
TEST(NetProto, TheLargestValueEachSenderWritesReadsBack) {
    MsgTerrainEdit till; till.body = 1; till.center = { 10.0, 60000.0, -2.0 }; till.radius = 0.6f; till.flags = 16; till.session = 2; till.seq = 1;
    MsgTerrainEdit tillBack;
    EXPECT_TRUE(RoundTrip(till, tillBack)) << "a till (flag 16, #1085) reaches the clients";
    EXPECT_EQ(tillBack.flags, 16);
    MsgTerrainEdit all = till; all.flags = 1 | 2 | 4 | 8 | 16; all.hx = all.hy = all.hz = 0.5f; all.ax[0] = 1.0f; all.ay[1] = 1.0f;
    EXPECT_TRUE(RoundTrip(all, tillBack)) << "every edit flag at once";
    MsgTerrainEdit beyond = till; beyond.flags = 32;
    EXPECT_FALSE(RoundTrip(beyond, tillBack)) << "a flag nobody sends is still refused";
    MsgMachineAction fire; fire.body = 1; fire.plot = 3; fire.part = 0; fire.action = 10; fire.count = 1;
    MsgMachineAction fireBack;
    EXPECT_TRUE(RoundTrip(fire, fireBack)) << "a campfire's light / put out (action 10, #1103) reaches the server";
    EXPECT_EQ(fireBack.action, 10);
    MsgMachineAction past = fire; past.action = 11;
    EXPECT_FALSE(RoundTrip(past, fireBack));
    MsgEditRequest q; q.seq = 1; q.kind = NetEditKind::Till; q.body = 1; q.hit = q.center = { 10.0, 60000.0, -2.0 }; q.radius = 0.6f;
    MsgEditRequest qBack;
    EXPECT_TRUE(RoundTrip(q, qBack)) << "the till request (kind Till) reaches the server";
}
#endif
