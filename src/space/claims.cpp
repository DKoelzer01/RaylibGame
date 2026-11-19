#include "claims.h"
#include "worldtypes.h"
#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <set>

namespace {
    std::vector<Claim> claims;
    std::set<std::pair<int, int>> allies; // (team, the team it lets in)
    uint32_t nextId = 1, version = 1;
    std::vector<ClaimView> mirror;
}

void Claims_Reset() { claims.clear(); allies.clear(); mirror.clear(); nextId = 1; ++version; }

int Claims_Budget(int bodiesReached) { return 1 + 2 * std::max(0, bodiesReached); }

bool Claims_Contains(const Claim& c, Vector3 p) {
    const Vector3 d = Vector3Subtract(p, c.centre);
    const float reach = c.radius + kClaimHalfHeight;
    if (Vector3DotProduct(d, d) > reach * reach) return false; // (the cheap test first: most claims are far)
    const float r = Vector3Length(c.centre);
    if (r < 1.0f) return false;
    const Vector3 up = Vector3Scale(c.centre, 1.0f / r);
    const float h = Vector3DotProduct(d, up);
    if (fabsf(h) > kClaimHalfHeight) return false;
    return Vector3Length(Vector3Subtract(d, Vector3Scale(up, h))) <= c.radius;
}

const Claim* Claims_At(const std::string& body, Vector3 local) {
    for (const Claim& c : claims) if (Claims_Contains(c, local) && c.body == body) return &c;
    return nullptr;
}

std::vector<Claim> Claims_OfTeam(int team) {
    std::vector<Claim> out;
    for (const Claim& c : claims) if (c.team == team) out.push_back(c);
    return out;
}

uint32_t Claims_Make(int team, int budget, const std::string& body, Vector3 local, float radius, int64_t now, std::string& why, uint32_t bannerPlot, int bannerPart) {
    if (team <= 0) { why = "You need a team to claim land"; return 0; }
    const int have = (int)std::count_if(claims.begin(), claims.end(), [team](const Claim& c) { return c.team == team; });
    if (have >= budget) { why = "Your team holds " + std::to_string(have) + " of its " + std::to_string(budget) + " claims (each planet or moon your team reaches adds 2)"; return 0; }
    for (const Claim& c : claims) {
        if (c.body != body) continue;
        const float apart = Vector3Distance(c.centre, local);
        if (c.team != team && apart < c.radius + radius) { why = TextFormat("Too close to another team's claim (%.0f m between the centres, %.0f needed)", apart, c.radius + radius); return 0; }
    }
    Claim c;
    c.id = nextId++; c.team = team; c.body = body; c.centre = local; c.radius = radius; c.bannerPlot = bannerPlot; c.bannerPart = bannerPart;
    c.created = c.lastSeen = now;
    claims.push_back(c);
    ++version;
    return c.id;
}

bool Claims_Remove(uint32_t id) {
    const size_t n = claims.size();
    claims.erase(std::remove_if(claims.begin(), claims.end(), [id](const Claim& c) { return c.id == id; }), claims.end());
    if (claims.size() == n) return false;
    ++version;
    return true;
}

uint32_t Claims_OfBanner(const std::string& body, uint32_t plot, int part) {
    for (const Claim& c : claims) if (c.bannerPlot == plot && c.bannerPart == part && plot != 0 && c.body == body) return c.id;
    return 0;
}

void Claims_MoveTeam(int from, int to) {
    if (!from || !to) return;
    for (Claim& c : claims) if (c.team == from) c.team = to;
    std::set<std::pair<int, int>> moved;
    for (const auto& [a, b] : allies) { const int x = a == from ? to : a, y = b == from ? to : b; if (x != y) moved.insert({ x, y }); }
    allies.swap(moved);
    ++version;
}

void Claims_SetAlly(int team, int other, bool on) {
    if (!team || !other || team == other) return;
    if (on) allies.insert({ team, other }); else allies.erase({ team, other });
    ++version;
}
bool Claims_Allows(int claimTeam, int team) { return claimTeam == team || (team && allies.count({ claimTeam, team })); }
std::vector<int> Claims_AlliesOf(int team) {
    std::vector<int> out;
    for (const auto& [a, b] : allies) if (a == team) out.push_back(b);
    return out;
}

int Claims_Update(const std::vector<int>& teamsOnline, int64_t now) {
    int removed = 0;
    for (size_t i = 0; i < claims.size();) {
        Claim& c = claims[i];
        if (std::find(teamsOnline.begin(), teamsOnline.end(), c.team) != teamsOnline.end()) c.lastSeen = now;
        if (now - c.lastSeen > kClaimDecaySeconds) { claims.erase(claims.begin() + i); ++removed; continue; }
        ++i;
    }
    if (removed) ++version;
    return removed;
}

const std::vector<Claim>& Claims_All() { return claims; }
uint32_t Claims_Version() { return version; }

void Claims_Save(FILE* f) {
    for (const Claim& c : claims)
        fprintf(f, "claim %u %d %s %.3f %.3f %.3f %.2f %lld %lld %u %d\n", c.id, c.team, c.body.c_str(), c.centre.x, c.centre.y, c.centre.z, c.radius,
                (long long)c.created, (long long)c.lastSeen, c.bannerPlot, c.bannerPart);
    for (const auto& [a, b] : allies) fprintf(f, "claim_ally %d %d\n", a, b);
}

bool Claims_LoadLine(const char* tok, FILE* f) {
    if (!strcmp(tok, "claim_ally")) { int a = 0, b = 0; if (fscanf(f, " %d %d", &a, &b) == 2) allies.insert({ a, b }); return true; }
    if (strcmp(tok, "claim") != 0) return false;
    Claim c; char body[128] = {}; long long created = 0, seen = 0;
    char line[512] = {};
    if (!fgets(line, sizeof line, f)) return true; // (the rest of the line: #864's lines lack the banner fields)
    const int got = sscanf(line, " %u %d %127s %f %f %f %f %lld %lld %u %d", &c.id, &c.team, body, &c.centre.x, &c.centre.y, &c.centre.z, &c.radius, &created, &seen, &c.bannerPlot, &c.bannerPart);
    if (got == 9 || got == 11) {
        c.body = body; c.created = created; c.lastSeen = seen;
        claims.push_back(c);
        nextId = std::max(nextId, c.id + 1);
        ++version;
    }
    return true;
}

void Claims_SetMirror(const std::vector<ClaimView>& list) { mirror = list; }
const std::vector<ClaimView>& Claims_Mirror() { return mirror; }

uint32_t Claims_TeamColour(int team, uint32_t set) {
    if (set) return set;
    static const uint32_t kPalette[8] = { 0xE0533Cu, 0x3C8CE0u, 0x52B848u, 0xE0B23Cu, 0xA05AD8u, 0x38C2B8u, 0xE07AB0u, 0x8C9C3Cu };
    return kPalette[(unsigned)team % 8];
}

bool Claims_SelfTest() {
    bool ok = true;
    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: claims: %s: NO\n", what); ok = false; } };
    const std::vector<Claim> saved = claims; const auto savedAllies = allies; const uint32_t savedNext = nextId;
    Claims_Reset();
    const float R = 60000.0f;
    auto on = [&](float east, float north, float up) { return Vector3{ east, R + up, north }; }; // near the pole: flat enough
    std::string why;
    // budget: 1 + 2 x bodies reached
    check(Claims_Budget(0) == 1 && Claims_Budget(1) == 3 && Claims_Budget(2) == 5, "budget 1 + 2 x bodies reached");
    const uint32_t a = Claims_Make(1, Claims_Budget(1), "terra", on(0, 0, 0), 50.0f, 1000, why, 7, 3);
    check(a != 0, "a team claims 100 m across with a banner");
    check(Claims_OfBanner("terra", 7, 3) == a && Claims_OfBanner("terra", 7, 4) == 0, "the claim is found by its banner");
    check(Claims_At("terra", on(45, 0, 20)) != nullptr, "45 m away, 20 m up: inside");
    check(Claims_At("terra", on(0, 0, -100)) != nullptr, "100 m down (a mine): inside");
    check(Claims_At("terra", on(55, 0, 0)) == nullptr, "55 m away: outside");
    check(Claims_At("terra", on(0, 0, 140)) == nullptr, "140 m up: outside");
    check(Claims_At("luna", on(0, 0, 0)) == nullptr, "another body: outside");
    check(Claims_Make(2, 3, "terra", on(400, 0, 0), 400.0f, 1000, why) == 0, "an 800 m claim 400 m off overlaps the 100 m one");
    const uint32_t big = Claims_Make(2, 3, "terra", on(460, 0, 0), 400.0f, 1000, why);
    check(big != 0, "800 m across at 460 m off: clear of it");
    check(Claims_At("terra", on(850, 0, 0)) != nullptr && Claims_At("terra", on(870, 0, 0)) == nullptr, "the 800 m claim reaches 400 m");
    check(Claims_Make(1, 3, "terra", on(0, 0, 200), 50.0f, 1000, why) != 0 && Claims_Make(1, 3, "terra", on(0, 0, 400), 50.0f, 1000, why) != 0, "the team's 2nd and 3rd claims");
    check(Claims_Make(1, 3, "terra", on(0, 0, 600), 50.0f, 1000, why) == 0, "no 4th with one body reached");
    // allies
    check(!Claims_Allows(1, 2), "another team isn't allowed");
    Claims_SetAlly(1, 2, true);
    check(Claims_Allows(1, 2) && !Claims_Allows(2, 1), "an ally is allowed (one way)");
    Claims_MoveTeam(2, 4);
    check(Claims_Allows(1, 4) && Claims_OfTeam(4).size() == 1, "a merged team keeps its claims and alliances");
    Claims_SetAlly(1, 4, false);
    check(!Claims_Allows(1, 4), "an ally revoked");
    check(Claims_Update({ 1 }, 1000 + kClaimDecaySeconds + 10) == 1 && Claims_OfTeam(4).empty() && Claims_OfTeam(1).size() == 3, "14 days away: team 4's claim decays, team 1 (online) keeps theirs");
    check(Claims_Remove(a) && Claims_OfTeam(1).size() == 2, "a banner broken: unclaimed");
    { // save and load, with an ally and a #864 line (no banner fields) before the next line
        Claims_SetAlly(1, 5, true);
        const std::string path = TextFormat("claims_test_%d.tmp", (int)(std::chrono::steady_clock::now().time_since_epoch().count() % 100000));
        if (FILE* f = fopen(path.c_str(), "w")) { Claims_Save(f); fprintf(f, "claim 40 6 terra 1.0 60000.0 2.0 50.00 900 900\nclaim_ally 6 1\n"); fclose(f); }
        Claims_Reset();
        if (FILE* f = fopen(path.c_str(), "r")) { char tok[64]; while (fscanf(f, " %63s", tok) == 1) if (!Claims_LoadLine(tok, f)) break; fclose(f); }
        remove(path.c_str());
        check(Claims_OfTeam(1).size() == 2 && Claims_OfTeam(6).size() == 1 && Claims_Allows(1, 5) && Claims_Allows(6, 1), "saved and loaded: claims, allies, and a line without banner fields");
    }

    // the lookup's cost: 1000 claims of 800 m across spread over a 60 km planet, a million points looked up
    Claims_Reset();
    uint32_t rng = 12345; auto rnd = [&rng] { rng = rng * 1664525u + 1013904223u; return (rng >> 8) / 16777216.0f; };
    for (int i = 0; i < 1000; ++i) {
        const Vector3 d = Vector3Normalize({ rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1 });
        Claim c; c.id = (uint32_t)i + 1; c.team = i + 1; c.body = "terra"; c.centre = Vector3Scale(d, R); c.radius = 400.0f; claims.push_back(c);
    }
    int hits = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < 1000000; ++i) { const Vector3 d = Vector3Normalize({ rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1 }); hits += Claims_At("terra", Vector3Scale(d, R)) != nullptr; }
    const double us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count() / 1e6;
    // 1000 x pi 400^2 = 5.0e8 m^2 of 4.5e10: ~1.1% of points inside (overlaps aside)
    check(hits > 5000 && hits < 20000, "about 1% of random points inside 1000 claims of 800 m");
    printf("bench: claims: lookup with 1000 claims of 800 m across: %.2f us a point (%d of 1000000 inside)\n", us, hits);
    printf("bench: claims: containment, budget, overlap, banners, allies, decay: %s\n", ok ? "PASS" : "FAIL");
    claims = saved; allies = savedAllies; nextId = savedNext; ++version;
    return ok;
}

void Claims_SetBanner(uint32_t id, uint32_t plot, int part) {
    for (Claim& c : claims) if (c.id == id) { c.bannerPlot = plot; c.bannerPart = part; ++version; }
}

// ---- The border near another team's claim (#343): a faint wall in its colour, only within kBorderShow of its edge ----
#include "world.h"
#include <rlgl.h>
void Claims_DrawBorder(Scene& world) {
    constexpr float kBorderShow = 30.0f, kArc = 40.0f;
    for (const ClaimView& c : mirror) {
        if (c.allowed) continue;
        const Planetoid* bp = nullptr;
        for (Planetoid* p : Planetoid::All()) if (p->name == c.body) { bp = p; break; }
        if (!bp) continue;
        const Planetoid& body = *bp;
        const Vector3 me = body.toLocal(world.camera.position); // (the camera's: what is seen; render space -> body-local)
        const float rc = Vector3Length(c.centre);
        if (rc < 1.0f) continue;
        const Vector3 up = Vector3Scale(c.centre, 1.0f / rc);
        const Vector3 d = Vector3Subtract(me, c.centre);
        const float h = Vector3DotProduct(d, up);
        const Vector3 horiz = Vector3Subtract(d, Vector3Scale(up, h));
        const float dist = Vector3Length(horiz);
        const float off = fabsf(dist - c.radius);
        if (dist < 1.0f || off > kBorderShow || fabsf(h) > kClaimHalfHeight) continue;
        const Vector3 e1 = Vector3Scale(horiz, 1.0f / dist), e2 = Vector3CrossProduct(up, e1);
        const float span = std::min(kArc / c.radius, PI);
        auto edge = [&](float a, float z) { return Vector3Add(Vector3Add(c.centre, Vector3Scale(up, h + z)), Vector3Add(Vector3Scale(e1, cosf(a) * c.radius), Vector3Scale(e2, sinf(a) * c.radius))); };
        const Vector3 anchor = edge(0.0f, 0.0f);
        const Matrix m = body.localModel(DVec3{ anchor.x, anchor.y, anchor.z }, world.renderOrigin);
        const float alpha = 0.3f * (1.0f - off / kBorderShow);
        const Color col = { (unsigned char)(c.colour >> 16), (unsigned char)(c.colour >> 8), (unsigned char)c.colour, (unsigned char)(255 * alpha) };
        rlDrawRenderBatchActive();
        rlDisableBackfaceCulling();
        rlDisableDepthMask();
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(m));
        constexpr int kSegs = 24;
        for (int i = 0; i < kSegs; ++i) {
            const float a0 = -span + 2.0f * span * i / kSegs, a1 = -span + 2.0f * span * (i + 1) / kSegs;
            const Vector3 p0 = Vector3Subtract(edge(a0, -2.5f), anchor), p1 = Vector3Subtract(edge(a1, -2.5f), anchor);
            const Vector3 q0 = Vector3Subtract(edge(a0, 0.6f), anchor), q1 = Vector3Subtract(edge(a1, 0.6f), anchor);
            DrawTriangle3D(p0, p1, q1, col); DrawTriangle3D(p0, q1, q0, col);
            const Color line = { col.r, col.g, col.b, (unsigned char)std::min(255.0f, 255 * alpha * 4.5f) };
            DrawLine3D(Vector3Subtract(edge(a0, 0.6f), anchor), Vector3Subtract(edge(a1, 0.6f), anchor), line); // its top edge: ~2.5 m over the ground (the camera is ~1.9 m over it)
            DrawLine3D(Vector3Subtract(edge(a0, -1.5f), anchor), Vector3Subtract(edge(a1, -1.5f), anchor), line); // near the ground
            DrawLine3D(Vector3Subtract(edge(a0, -2.5f), anchor), Vector3Subtract(edge(a0, 0.6f), anchor), line); // a post every segment
        }
        rlDrawRenderBatchActive();
        rlPopMatrix();
        rlEnableDepthMask();
        rlEnableBackfaceCulling();
    }
}
