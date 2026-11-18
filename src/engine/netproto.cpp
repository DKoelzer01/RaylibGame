#include "netproto.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>

const char* NetMsgName(NetMsg m) {
    switch (m) {
        case NetMsg::Hello: return "hello"; case NetMsg::Welcome: return "welcome"; case NetMsg::Reject: return "reject";
        case NetMsg::Bye: return "bye"; case NetMsg::Ping: return "ping"; case NetMsg::Pong: return "pong";
        case NetMsg::Input: return "input"; case NetMsg::PlayerState: return "player state"; case NetMsg::TerrainEdit: return "terrain edit";
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
        case NetReject::NameTaken: return "that name is already playing";
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

void NetWrite(NetWriter& w, const MsgHello& m) { Id(w, NetMsg::Hello); w.u32(m.protocol); w.str(m.build); w.str(m.name); w.str(m.password); }
bool NetRead(NetReader& r, MsgHello& m) {
    m.protocol = r.u32();
    // An older or newer protocol: read nothing more (its layout may differ), the server answers Reject::Protocol.
    if (r.ok && m.protocol != kNetProtocol) return true;
    m.build = r.str(64); m.name = r.str(32); m.password = r.str(64);
    return Done(r) && !m.name.empty();
}

void NetWrite(NetWriter& w, const MsgWelcome& m) {
    Id(w, NetMsg::Welcome);
    w.u32(m.protocol); w.u32(m.session); w.boolean(m.reconnected); w.u32(m.seed);
    w.varu(m.bodies.size()); for (const std::string& b : m.bodies) w.str(b);
    w.u64(m.fingerprint); w.u32(m.serverTick);
    w.u8(m.settings.difficulty); w.boolean(m.settings.pvp); w.boolean(m.settings.claims); w.u8(m.settings.maxPlayers);
}
bool NetRead(NetReader& r, MsgWelcome& m) {
    m.protocol = r.u32(); m.session = r.u32(); m.reconnected = r.boolean(); m.seed = r.u32();
    const uint64_t n = r.varu();
    if (n > 256) r.ok = false;
    m.bodies.clear();
    for (uint64_t i = 0; i < n && r.ok; ++i) m.bodies.push_back(r.str(1024));
    m.fingerprint = r.u64(); m.serverTick = r.u32();
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
void NetWrite(NetWriter& w, const MsgPong& m) { Id(w, NetMsg::Pong); w.u32(m.seq); w.f64(m.clientTime); w.u32(m.serverTick); }
bool NetRead(NetReader& r, MsgPong& m) { m.seq = r.u32(); m.clientTime = r.f64(); m.serverTick = r.u32(); return Done(r); }

void NetWrite(NetWriter& w, const MsgInput& m) { Id(w, NetMsg::Input); w.u32(m.tick); w.varu(m.held); w.varu(m.pressed); w.angle(m.yaw); w.u16((uint16_t)Q16(m.pitch / 1.6f)); }
bool NetRead(NetReader& r, MsgInput& m) {
    m.tick = r.u32(); m.held = r.varu(); m.pressed = r.varu(); m.yaw = r.angle(); m.pitch = (int16_t)r.u16() / 32767.0f * 1.6f;
    return Done(r);
}

void NetWrite(NetWriter& w, const MsgPlayerState& m) {
    Id(w, NetMsg::PlayerState);
    w.varu(m.session); w.u32(m.tick); w.u32(m.ackInput); w.varu(m.body); w.pos(m.local);
    w.speed(m.vx); w.speed(m.vy); w.speed(m.vz); w.angle(m.yaw); w.u16((uint16_t)Q16(m.pitch / 1.6f)); w.u8(m.flags);
}
bool NetRead(NetReader& r, MsgPlayerState& m) {
    m.session = (uint32_t)r.varu(); m.tick = r.u32(); m.ackInput = r.u32(); m.body = r.varu(); m.local = r.pos();
    m.vx = r.speed(); m.vy = r.speed(); m.vz = r.speed(); m.yaw = r.angle(); m.pitch = (int16_t)r.u16() / 32767.0f * 1.6f; m.flags = r.u8();
    return Done(r);
}

void NetWrite(NetWriter& w, const MsgTerrainEdit& m) {
    Id(w, NetMsg::TerrainEdit);
    w.varu(m.body); w.pos(m.center); w.f32(m.radius); w.u8(m.flags); w.u8((uint8_t)m.paint);
    if (m.flags & 4) { w.f32(m.hx); w.f32(m.hy); w.f32(m.hz); w.unit(m.ax[0], m.ax[1], m.ax[2]); w.unit(m.ay[0], m.ay[1], m.ay[2]); }
}
bool NetRead(NetReader& r, MsgTerrainEdit& m) {
    m.body = r.varu(); m.center = r.pos(); m.radius = r.f32(); m.flags = r.u8(); m.paint = (int8_t)r.u8();
    if (m.flags & 4) { m.hx = r.f32(); m.hy = r.f32(); m.hz = r.f32(); r.unit(m.ax[0], m.ax[1], m.ax[2]); r.unit(m.ay[0], m.ay[1], m.ay[2]); }
    if (m.radius < 0.0f || m.radius > 64.0f || (m.flags & ~0x0F)) r.ok = false; // a dig is metres, not a planet
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
        case NetMsg::Input: { MsgInput m; return NetRead(r, m); }
        case NetMsg::PlayerState: { MsgPlayerState m; return NetRead(r, m); }
        case NetMsg::TerrainEdit: { MsgTerrainEdit m; return NetRead(r, m); }
        default: return false;
    }
}
}

bool NetProto_Test() {
    bool ok = true;
    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: netproto: %s: NO\n", what); ok = false; } };
    std::vector<std::vector<uint8_t>> samples; // each message's bytes, for the truncation and fuzz passes
    {
        MsgHello a{ kNetProtocol, "m84-integration abc1234", "Dewey", "hunter2" }, b; std::vector<uint8_t> s;
        check(RoundTrip(a, b, &s) && b.protocol == a.protocol && b.build == a.build && b.name == a.name && b.password == a.password, "hello");
        samples.push_back(s);
        MsgHello old{ 1, "x", "y", "" }, ob;
        check(RoundTrip(old, ob) && ob.protocol == 1, "an old protocol's hello still reads its version (to reject it clearly)");
    }
    {
        MsgWelcome a; a.session = 7; a.reconnected = true; a.seed = 42; a.bodies = { "planet Terra 0 0 0 60000", "moon Luna 1 2 3 9000" };
        a.fingerprint = 0x1234567890abcdefull; a.serverTick = 123456; a.settings = { 2, true, false, 16 };
        MsgWelcome b; std::vector<uint8_t> s;
        check(RoundTrip(a, b, &s) && b.session == 7 && b.reconnected && b.seed == 42 && b.bodies == a.bodies && b.fingerprint == a.fingerprint &&
              b.serverTick == 123456 && b.settings.difficulty == 2 && b.settings.pvp && !b.settings.claims && b.settings.maxPlayers == 16, "welcome");
        samples.push_back(s);
    }
    { MsgReject a{ NetReject::Password, "try again" }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.reason == a.reason && b.detail == a.detail, "reject"); samples.push_back(s); }
    { MsgBye a{ NetBye::ServerClosed }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.reason == a.reason, "bye"); samples.push_back(s); }
    { MsgPing a{ 99, 12.5 }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.seq == 99 && b.clientTime == 12.5, "ping"); samples.push_back(s); }
    { MsgPong a{ 99, 12.5, 777 }, b; std::vector<uint8_t> s; check(RoundTrip(a, b, &s) && b.seq == 99 && b.clientTime == 12.5 && b.serverTick == 777, "pong"); samples.push_back(s); }
    {
        MsgInput a{ 5000, 0x1234ull, 0x10ull, 1.25f, -0.4f }, b; std::vector<uint8_t> s;
        check(RoundTrip(a, b, &s) && b.tick == 5000 && b.held == a.held && b.pressed == a.pressed && fabsf(b.yaw - a.yaw) < 1e-3f && fabsf(b.pitch - a.pitch) < 1e-3f, "input");
        samples.push_back(s);
    }
    {
        MsgPlayerState a; a.session = 3; a.tick = 9001; a.ackInput = 8999; a.body = 0xABCDull; a.local = { 60012.3456, -123.4567, 0.0005 };
        a.vx = 4.21f; a.vy = -9.81f; a.vz = 0.0f; a.yaw = 6.0f; a.pitch = 1.2f; a.flags = 0b1011;
        MsgPlayerState b; std::vector<uint8_t> s;
        const bool rt = RoundTrip(a, b, &s);
        const double dp = std::fmax(std::fabs(b.local.x - a.local.x), std::fmax(std::fabs(b.local.y - a.local.y), std::fabs(b.local.z - a.local.z)));
        check(rt && b.session == 3 && b.tick == 9001 && b.ackInput == 8999 && b.body == a.body && dp <= 0.0005 + 1e-9 &&
              fabsf(b.vx - a.vx) < 0.006f && fabsf(b.vy - a.vy) < 0.006f && fabsf(b.yaw - a.yaw) < 1e-3f && fabsf(b.pitch - a.pitch) < 1e-3f && b.flags == a.flags,
              "player state (1 mm, 1 cm/s, 1e-3 rad at 60 km)");
        printf("bench: netproto: a player state is %zu bytes (position error %.4f mm at 60 km)\n", s.size(), dp * 1000.0);
        samples.push_back(s);
    }
    {
        MsgTerrainEdit a; a.body = 1; a.center = { 100.25, 59950.5, -3.75 }; a.radius = 1.5f; a.flags = 1 | 4; a.hx = 0.25f; a.hy = 0.5f; a.hz = 0.25f;
        a.ay[0] = 0.0f; a.ay[1] = 0.6f; a.ay[2] = 0.8f; a.paint = 3;
        MsgTerrainEdit b; std::vector<uint8_t> s;
        check(RoundTrip(a, b, &s) && b.body == 1 && std::fabs(b.center.y - a.center.y) < 0.001 && b.radius == a.radius && b.flags == a.flags &&
              b.hy == a.hy && fabsf(b.ay[2] - 0.8f) < 1e-3f && b.paint == 3, "terrain edit (box)");
        samples.push_back(s);
    }
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
