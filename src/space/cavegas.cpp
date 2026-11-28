#include "cavegas.h"
#include "audio.h"
#include "camera.h"
#include "hotbar.h"
#include "inventory.h"
#include "net.h"
#include "particles.h"
#include "player.h"
#include "playersession.h"
#include "plots.h"
#include "scene.h"
#include "simclock.h"
#include "suit.h"
#include "survival.h"
#include "voxel.h"
#include "world.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

namespace {
    // Chances that a cave span holds a pocket, by depth under sea level (GUESS; the depth table #1952 takes them over).
    constexpr float kChanceShallow = 0.04f;   // 25-150 m
    constexpr float kChanceMiddle = 0.08f;    // 150-500 m
    constexpr float kChanceDeep = 0.14f;      // deeper
    constexpr float kChamberFactor = 2.0f;    // a chamber's air is still: twice as likely
    constexpr float kMethaneShare = 0.4f;     // of the pockets; the rest are CO2
    constexpr float kLayer = 0.3f;            // CO2 lies below centre + 0.3 r, methane above centre - 0.3 r
    constexpr float kWarnBasic = 3.0f;        // m from a pocket's edge: any suit warns
    constexpr float kWarnSensor = 15.0f;      // m: with a gas sensor module (and a direction)
    constexpr float kHissReach = 6.0f;        // m from the edge: the hiss is heard
    constexpr float kHazeReach = 30.0f;       // m from the player: pockets that show their haze
    constexpr float kNearRange = 40.0f;       // m: the pockets kept round the player
    constexpr float kBlastDamage = 40.0f;     // hp at the centre (GUESS), falling to 0 at r + kBlastReach
    constexpr float kBlastReach = 4.0f;       // m past the pocket's edge
    constexpr double kBurnTime = 600.0;       // s a burnt methane pocket stays clear (GUESS)
    constexpr double kBlastWindow = 0.5;      // s a blast is there for the survival steps to take

    struct Blast { uint32_t id; const Planetoid* body; Vector3 centre; float radius; double at; };
    struct CacheEntry { const Planetoid* body = nullptr; Vector3 centre{}; double at = -1e9; std::vector<CaveGasPocket> pockets; };

    std::map<uint32_t, double> burntUntil;            // methane pockets set off: id -> sim time they come back
    std::vector<Blast> blasts;                        // recent ones (kBlastWindow)
    uint32_t nextBlastId = 1;
    std::vector<std::pair<const Planetoid*, CaveGasPocket>> testPockets;
    CacheEntry cache[8];                              // CaveGas_At's pockets near a point (the local player, remote players)
    std::vector<CaveGasPocket> near;                  // round the local player (Update)
    std::vector<Vector3> flames;                      // body-local: placed torches and campfires round the player
    const Planetoid* flamesBody = nullptr;
    double nextNear = -1.0, nextFlames = -1.0, nextHaze = 0.0;
    int hiss = -2;                                    // the hiss emitter (-2 not made yet, -1 none: audio off)
    CaveGasStats stats;
    Vector3 warnDir{};                                // body-local, head to the pocket (the sensor's direction)

    uint64_t Mix(uint64_t h) { h ^= h >> 33; h *= 0xff51afd7ed558ccdULL; h ^= h >> 33; h *= 0xc4ceb9fe1a85ec53ULL; h ^= h >> 33; return h; }
    uint64_t Q(float v) { return (uint64_t)(int64_t)floorf(v); }  // a span end, to the metre (cave segments are bit-exact per seed)
    float U(uint64_t& h) { h = Mix(h + 0x9E3779B97F4A7C15ULL); return (float)(h >> 40) / 16777216.0f; } // [0, 1)

    float ChanceAt(float depth) { return depth < 150.0f ? kChanceShallow : depth < 500.0f ? kChanceMiddle : kChanceDeep; }

    bool Burnt(uint32_t id) { auto it = burntUntil.find(id); return it != burntUntil.end() && SimTime() < it->second; }

    const PlanetField* FieldOf(const Planetoid* b) { return b ? dynamic_cast<const PlanetField*>(b->field.get()) : nullptr; }

    void Near(const Planetoid* body, Vector3 local, float range, std::vector<CaveGasPocket>& out) {
        out.clear();
        if (const PlanetField* pf = FieldOf(body)) CaveGas_PocketsNear(*pf, local, range, out);
        for (const auto& [b, p] : testPockets)
            if (b == body && Vector3Distance(p.centre, local) < range + p.radius) out.push_back(p);
    }

    // The gas at a body-local point among these pockets (a burnt methane pocket is clear).
    CaveGasKind GasIn(const std::vector<CaveGasPocket>& ps, Vector3 local) {
        CaveGasKind k = CaveGasKind::None;
        for (const CaveGasPocket& p : ps) {
            if (!CaveGas_InPocket(p, local) || (p.kind == CaveGasKind::Methane && Burnt(p.id))) continue;
            if (p.kind == CaveGasKind::CO2) return CaveGasKind::CO2; // the worse one
            k = p.kind;
        }
        return k;
    }

    bool LocalTorchHeld() {
        const PlayerSession& s = LocalSession();
        const std::string& held = s.inventory.grid()[Inventory::kGrid + s.hotbar.selected].id;
        return !held.empty() && Hotbar_IsLight(held) && s.inventory.count(held) > 0;
    }

    void ScanFlames(const Planetoid* body, Vector3 me) {
        flames.clear();
        flamesBody = body;
        if (!body || !body->plots) return;
        for (const auto& plot : body->plots->plots) {
            if (Vector3Distance(plot->toBody({ 0, 0, 0 }), me) > 120.0f) continue;
            for (const Part& part : plot->parts)
                if (part.type == PartType::Torch || part.type == PartType::Campfire)
                    flames.push_back(plot->toBody({ (part.origin.x + 0.5f) * kCellSize, (part.origin.y + 0.5f) * kCellSize, (part.origin.z + 0.5f) * kCellSize }));
        }
    }

    void SetOff(Scene& world, const Planetoid& body, const CaveGasPocket& p) {
        const double now = SimTime();
        burntUntil[p.id] = now + kBurnTime;
        blasts.push_back({ nextBlastId++, &body, p.centre, p.radius, now });
        ++stats.blasts;
        ParticleBurst f; f.kind = ParticleKind::Fire; f.color = { 255, 170, 60, 255 };
        f.count = 60; f.size = 0.9f; f.life = 0.8f; f.speed = p.radius * 1.5f; f.spread = 1.0f;
        Particles_Emit(body, p.centre, f);
        ParticleBurst s; s.kind = ParticleKind::Smoke; s.color = { 70, 64, 58, 160 };
        s.count = 30; s.size = 1.2f; s.life = 4.0f; s.speed = 1.5f; s.spread = 1.0f;
        Particles_Emit(body, p.centre, s);
        const DVec3 at = body.toWorldD(DVec3(p.centre)); // (absolute world, doubles: no render origin in the sim)
        if (const int boom = Audio_Sound("thunder_near"); boom >= 0) Audio_PlaySoundAt(boom, at, AudioBus::Effects, AudioPriority::High, 1.0f, 1.6f);
        else Audio_PlayAt(Sfx::Dig, at, 1.0f, 0.4f);
    }
}

void CaveGas_PocketsNear(const PlanetField& f, Vector3 local, float range, std::vector<CaveGasPocket>& out) {
    constexpr float kMaxR = 12.0f;
    std::vector<PlanetField::CaveSpan> spans;
    f.caveSpansNear(Vector3SubtractValue(local, range + kMaxR), Vector3AddValue(local, range + kMaxR), spans);
    const uint64_t seed = f.terrainSeed();
    for (const PlanetField::CaveSpan& sp : spans) {
        if (sp.kind != 0 && sp.kind != 3) continue; // passages and chambers (not mouths: they breathe; not rubble or slots)
        uint64_t h = Mix(seed * 0x9E3779B97F4A7C15ULL ^ 0x6A5C0DE5ULL);
        h = Mix(h ^ Q(sp.a.x)); h = Mix(h ^ Q(sp.a.y)); h = Mix(h ^ Q(sp.a.z));
        h = Mix(h ^ Q(sp.b.x)); h = Mix(h ^ Q(sp.b.y)); h = Mix(h ^ Q(sp.b.z));
        const uint32_t id = (uint32_t)(h >> 32) & 0x7FFFFFFFu;
        const float t = U(h), roll = U(h), size = U(h), kind = U(h);
        const Vector3 c = Vector3Lerp(sp.a, sp.b, t);
        const float depth = f.radius - Vector3Length(c);
        if (roll >= ChanceAt(depth) * (sp.kind == 3 ? kChamberFactor : 1.0f)) continue;
        const float spanR = sp.ra + (sp.rb - sp.ra) * t;
        CaveGasPocket p;
        p.id = id | 1u;
        p.centre = c;
        p.up = Vector3Normalize(sp.up);
        p.radius = std::min(kMaxR, std::max(4.0f + 6.0f * size, spanR * 1.2f)); // it fills the passage across
        p.kind = kind < kMethaneShare ? CaveGasKind::Methane : CaveGasKind::CO2;
        if (Vector3Distance(c, local) < range + p.radius) out.push_back(p);
    }
}

bool CaveGas_InPocket(const CaveGasPocket& p, Vector3 local) {
    const Vector3 d = Vector3Subtract(local, p.centre);
    if (Vector3LengthSqr(d) >= p.radius * p.radius) return false;
    const float h = Vector3DotProduct(d, p.up);
    return p.kind == CaveGasKind::CO2 ? h < kLayer * p.radius : h > -kLayer * p.radius;
}

CaveGasKind CaveGas_At(const Scene& world, const Planetoid* body, Vector3 renderPos) {
    if (!body || (!FieldOf(body) && testPockets.empty())) return CaveGasKind::None;
    (void)world;
    const Vector3 local = body->toLocal(renderPos);
    const double now = SimTime();
    // The pockets within 30 m of a point, kept a second (or until the point moves 15 m): one cave query per player a second.
    CacheEntry* e = nullptr;
    for (CacheEntry& c : cache)
        if (c.body == body && now - c.at < 1.0 && now >= c.at && Vector3Distance(c.centre, local) < 15.0f) { e = &c; break; }
    if (!e) {
        e = &cache[0];
        for (CacheEntry& c : cache) if (c.at < e->at) e = &c;
        e->body = body; e->centre = local; e->at = now;
        Near(body, local, 30.0f, e->pockets);
    }
    return GasIn(e->pockets, local);
}

float CaveGas_BlastDamage(Vector3 renderPos, uint32_t& lastSeen) {
    float dmg = 0.0f;
    const double now = SimTime();
    for (const Blast& b : blasts) {
        if (b.id <= lastSeen || now - b.at > kBlastWindow) continue;
        lastSeen = b.id;
        const float d = Vector3Distance(b.body->toLocal(renderPos), b.centre);
        const float reach = b.radius + kBlastReach;
        if (d < reach) dmg += kBlastDamage * (1.0f - d / reach);
    }
    return dmg;
}

void CaveGas_Update(Scene& world) {
    const double now = SimTime();
    blasts.erase(std::remove_if(blasts.begin(), blasts.end(), [&](const Blast& b) { return now - b.at > kBlastWindow; }), blasts.end());
    const Planetoid* body = player.nearest;
    stats.warnKind = 0; stats.warnDist = -1.0f; stats.edgeDist = -1.0f; stats.hissing = false;
    if (!body || (!FieldOf(body) && testPockets.empty()) || LocalSession().survival.dead) {
        near.clear(); stats.pockets = 0;
        if (hiss >= 0) Audio_EmitterSetPlaying(hiss, false);
        return;
    }
    const Vector3 feet = body->toLocalD(player.position).toFloat();
    const Vector3 head = Vector3Add(feet, Vector3Scale(body->dirToLocal(player.up), player.eyeHeight));
    if (now >= nextNear || now < nextNear - 1.0) { nextNear = now + 0.25; Near(body, head, kNearRange, near); }
    if (now >= nextFlames || now < nextFlames - 1.0 || flamesBody != body) { nextFlames = now + 0.5; ScanFlames(body, feet); }
    stats.pockets = (int)near.size();

    // Flames in methane: the authority sets the pocket off (a placed torch or campfire: the same on every peer; the local
    // player's held torch). A client's held torch waits for a message to the host (later).
    if (Net_IsAuthority()) {
        const bool torch = LocalTorchHeld();
        for (const CaveGasPocket& p : near) {
            if (p.kind != CaveGasKind::Methane || Burnt(p.id)) continue;
            bool lit = torch && CaveGas_InPocket(p, head);
            for (const Vector3& fl : flames) lit = lit || CaveGas_InPocket(p, fl);
            if (lit) SetOff(world, *body, p);
        }
    }

    // The nearest live pocket's edge: the warning and the hiss.
    const bool sensor = Suit_Stats().gasSensor;
    float best = 1e9f; const CaveGasPocket* bp = nullptr;
    for (const CaveGasPocket& p : near) {
        if (p.kind == CaveGasKind::Methane && Burnt(p.id)) continue;
        const float d = CaveGas_InPocket(p, head) ? 0.0f : std::max(0.0f, Vector3Distance(head, p.centre) - p.radius);
        if (d < best) { best = d; bp = &p; }
    }
    stats.edgeDist = bp ? best : -1.0f;
    if (bp && best < (sensor ? kWarnSensor : kWarnBasic)) {
        stats.warnKind = (int)bp->kind; stats.warnDist = best;
        warnDir = Vector3Subtract(bp->centre, head);
    }
    if (hiss == -2) {
        AudioAttenuation att; att.refDistance = 1.5f; att.maxDistance = 14.0f;
        const int snd = Audio_Sound("gas_hiss");
        hiss = snd >= 0 ? Audio_EmitterCreate(snd, AudioBus::Ambience, AudioPriority::High, att) : -1;
    }
    if (bp && best < kHissReach) {
        stats.hissing = true;
        if (hiss >= 0) {
            const Vector3 toHead = Vector3Subtract(head, bp->centre);
            const float len = Vector3Length(toHead);
            const Vector3 edge = Vector3Add(bp->centre, Vector3Scale(len > 1e-3f ? Vector3Scale(toHead, 1.0f / len) : bp->up, std::min(len, bp->radius)));
            Audio_EmitterSetPosition(hiss, body->toWorldD(DVec3(edge)));
            Audio_EmitterSetVolume(hiss, bp->kind == CaveGasKind::CO2 ? 0.55f : 0.45f);
            Audio_EmitterSetPitch(hiss, bp->kind == CaveGasKind::CO2 ? 0.9f : 1.2f);
            Audio_EmitterSetPlaying(hiss, true);
        }
    } else if (hiss >= 0) Audio_EmitterSetPlaying(hiss, false);

    // The haze: a few slow puffs at random air points of each live pocket within 30 m (CO2 low and grey-green, methane a
    // faint yellow shimmer).
    if (now >= nextHaze || now < nextHaze - 1.0) {
        nextHaze = now + 0.2;
        static uint64_t r = 0x51ED2701ULL;
        for (const CaveGasPocket& p : near) {
            if (p.kind == CaveGasKind::Methane && Burnt(p.id)) continue;
            if (Vector3Distance(p.centre, head) > kHazeReach + p.radius) continue;
            for (int k = 0, shown = 0; k < 32 && shown < 10; ++k) {
                // Half round the pocket, half within 6 m of the head (where the player looks: a pocket's sphere often
                // spans several passages, and puffs behind the rock show nothing).
                const Vector3 around = k % 2 ? head : p.centre;
                const float reach = k % 2 ? 6.0f : p.radius * 0.9f;
                const Vector3 q = Vector3Add(around, Vector3Scale({ U(r) * 2 - 1, U(r) * 2 - 1, U(r) * 2 - 1 }, reach));
                if (!CaveGas_InPocket(p, q) || body->densityAt(q) >= ISO_LEVEL) continue;
                ++shown; ++stats.puffs;
                ParticleBurst s;
                s.count = 1; s.spread = 1.0f;
                // CO2 olive, darker than lamp-lit rock (GUESS)
                if (p.kind == CaveGasKind::CO2) { s.kind = ParticleKind::Dust; s.color = { 96, 128, 48, 235 }; s.size = 1.8f; s.life = 5.0f; s.speed = 0.06f; s.dir = Vector3Negate(p.up); }
                else { s.kind = ParticleKind::Smoke; s.color = { 222, 206, 128, 70 }; s.size = 1.1f; s.life = 4.0f; s.speed = 0.10f; s.dir = p.up; }
                Particles_Emit(*body, q, s);
            }
        }
    }
}

void CaveGas_DrawHUD(const Scene& world) {
    if (stats.warnKind == 0) return;
    const SurvivalState& sv = LocalSession().survival;
    const bool co2 = stats.warnKind == (int)CaveGasKind::CO2;
    const bool in = sv.caveGas == stats.warnKind;
    const char* what = co2 ? (in ? "IN CO2: unbreathable air, oxygen draining x3" : "WARNING: CO2 ahead - unbreathable air")
                           : (in ? "IN METHANE: no flames! oxygen draining" : "WARNING: methane ahead - no flames");
    std::string line = what;
    if (!in && Suit_Stats().gasSensor && stats.warnDist >= 0.0f) {
        // Where it is, for the sensor: against the camera's view.
        const Planetoid* body = player.nearest;
        const Vector3 d = body ? body->dirToWorld(warnDir) : warnDir;
        const Vector3 fwd = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
        const Vector3 up = Vector3Normalize(world.camera.up), right = Vector3CrossProduct(fwd, up);
        const float f = Vector3DotProduct(d, fwd), rt = Vector3DotProduct(d, right), u = Vector3DotProduct(d, up);
        const char* where = fabsf(u) > std::max(fabsf(f), fabsf(rt)) ? (u > 0 ? "above" : "below")
                          : fabsf(rt) > fabsf(f) ? (rt > 0 ? "right" : "left") : (f > 0 ? "ahead" : "behind");
        line += TextFormat("  (%.0f m, %s)", stats.warnDist, where);
    }
    const int fs = 20;
    const int w = MeasureText(line.c_str(), fs);
    const int x = (GetScreenWidth() - w) / 2, y = (int)(GetScreenHeight() * 0.26f);
    const float pulse = 0.75f + 0.25f * sinf((float)GetTime() * 6.0f);
    DrawRectangle(x - 10, y - 6, w + 20, fs + 12, Fade(BLACK, 0.45f));
    DrawText(line.c_str(), x, y, fs, Fade(co2 ? Color{ 170, 230, 120, 255 } : Color{ 255, 200, 90, 255 }, pulse));
}

void CaveGas_Reset() {
    burntUntil.clear(); blasts.clear(); testPockets.clear(); near.clear(); flames.clear(); flamesBody = nullptr;
    for (CacheEntry& c : cache) c = CacheEntry{};
    nextNear = nextFlames = -1.0; stats = {};
}

void CaveGas_TestAdd(const Planetoid& body, const CaveGasPocket& p) {
    CaveGasPocket q = p;
    q.id = 0xF0000000u + (uint32_t)testPockets.size();
    testPockets.push_back({ &body, q });
    for (CacheEntry& c : cache) c = CacheEntry{};
    nextNear = -1.0;
}

const CaveGasStats& CaveGas_Stats() { return stats; }

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1472): pockets are a pure function of the seed, sit in real cave spans, and layer by kind ----
#include "unit_helpers.h"

TEST(CaveGas, PocketsComeFromTheSeedAndSitInCaveSpans) {
    PlanetField::tIgnoreBake = true;
    PlanetField f(42u, 60000.0f, false);
    const Vector3 dir = f.findSpawnDir(42u);
    const Vector3 at = Vector3Scale(dir, f.radius - 90.0f); // the passage band, 90 m under sea level
    constexpr float kRange = 600.0f;
    std::vector<CaveGasPocket> a, b;
    CaveGas_PocketsNear(f, at, kRange, a);
    CaveGas_PocketsNear(f, at, kRange, b);
    std::vector<PlanetField::CaveSpan> spans;
    f.caveSpansNear(Vector3SubtractValue(at, kRange), Vector3AddValue(at, kRange), spans);
    int open = 0; for (const auto& s : spans) open += (s.kind == 0 || s.kind == 3) ? 1 : 0;
    ASSERT_GT(open, 20) << "no cave passages near the seed-42 spawn to put pockets in";
    ASSERT_GT((int)a.size(), 0) << "no pockets among " << open << " cave spans";
    EXPECT_LT((float)a.size() / open, 0.35f) << "pockets are rare: " << a.size() << " of " << open << " spans";
    ASSERT_EQ(a.size(), b.size());
    int co2 = 0, methane = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        EXPECT_EQ(a[i].id, b[i].id);
        EXPECT_EQ(a[i].kind, b[i].kind);
        EXPECT_GE(a[i].radius, 4.0f); EXPECT_LE(a[i].radius, 12.0f);
        EXPECT_LT(a[i].id, 0xF0000000u) << "below the test pockets' ids";
        (a[i].kind == CaveGasKind::CO2 ? co2 : methane)++;
        // Its centre lies on a cave span's axis (within the span's radius of it): in a cave, not in rock.
        bool onSpan = false;
        for (const auto& s : spans) {
            const Vector3 ab = Vector3Subtract(s.b, s.a);
            const float t = std::clamp(Vector3DotProduct(Vector3Subtract(a[i].centre, s.a), ab) / std::max(1e-6f, Vector3LengthSqr(ab)), 0.0f, 1.0f);
            onSpan = onSpan || Vector3Distance(a[i].centre, Vector3Add(s.a, Vector3Scale(ab, t))) < 0.01f;
        }
        EXPECT_TRUE(onSpan) << "pocket " << a[i].id << " off every cave axis";
    }
    EXPECT_GT(co2, 0); EXPECT_GT(methane, 0);
    printf("cave gas (seed 42, %.0f m round the spawn, 90 m down): %d pockets (%d CO2, %d methane) in %d cave spans\n", kRange, (int)a.size(), co2, methane, open);
    // Another seed: other pockets.
    PlanetField g(43u, 60000.0f, false);
    std::vector<CaveGasPocket> c;
    CaveGas_PocketsNear(g, Vector3Scale(g.findSpawnDir(43u), g.radius - 90.0f), kRange, c);
    bool same = c.size() == a.size();
    for (size_t i = 0; same && i < a.size(); ++i) same = a[i].id == c[i].id;
    EXPECT_FALSE(same);
    PlanetField::tIgnoreBake = false;
}

TEST(CaveGas, CarbonDioxideLiesLowAndMethaneRises) {
    CaveGasPocket p;
    p.centre = { 0, 0, 0 }; p.up = { 0, 1, 0 }; p.radius = 6.0f;
    p.kind = CaveGasKind::CO2;
    EXPECT_TRUE(CaveGas_InPocket(p, { 0, -4, 0 }));   // the floor of it
    EXPECT_TRUE(CaveGas_InPocket(p, { 0, 1.5f, 0 })); // up to 0.3 r over the centre
    EXPECT_FALSE(CaveGas_InPocket(p, { 0, 2.5f, 0 }));
    EXPECT_FALSE(CaveGas_InPocket(p, { 6.5f, -1, 0 })); // outside the sphere
    p.kind = CaveGasKind::Methane;
    EXPECT_TRUE(CaveGas_InPocket(p, { 0, 4, 0 }));
    EXPECT_TRUE(CaveGas_InPocket(p, { 0, -1.5f, 0 }));
    EXPECT_FALSE(CaveGas_InPocket(p, { 0, -2.5f, 0 }));
}
#endif
