#include "claims.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <cstring>

namespace {
    std::vector<Claim> claims;
    uint32_t nextId = 1;
}

void Claims_Reset() { claims.clear(); nextId = 1; }

int Claims_Budget(int members) { return std::clamp(1 + members, 2, 9); }

bool Claims_Contains(const Claim& c, Vector3 p) {
    const float r = Vector3Length(c.centre);
    if (r < 1.0f) return false;
    const Vector3 up = Vector3Scale(c.centre, 1.0f / r);
    const Vector3 d = Vector3Subtract(p, c.centre);
    const float h = Vector3DotProduct(d, up);
    if (fabsf(h) > kClaimHalfHeight) return false;
    return Vector3Length(Vector3Subtract(d, Vector3Scale(up, h))) <= c.radius;
}

const Claim* Claims_At(const std::string& body, Vector3 local) {
    for (const Claim& c : claims) if (c.body == body && Claims_Contains(c, local)) return &c;
    return nullptr;
}

std::vector<Claim> Claims_OfTeam(int team) {
    std::vector<Claim> out;
    for (const Claim& c : claims) if (c.team == team) out.push_back(c);
    return out;
}

uint32_t Claims_Make(int team, int members, const std::string& body, Vector3 local, int64_t now, std::string& why) {
    if (team <= 0) { why = "You need a team to claim land"; return 0; }
    const int have = (int)std::count_if(claims.begin(), claims.end(), [team](const Claim& c) { return c.team == team; });
    const int budget = Claims_Budget(members);
    if (have >= budget) { why = "Your team holds " + std::to_string(have) + " of its " + std::to_string(budget) + " claims (more members: more claims; /unclaim one)"; return 0; }
    for (const Claim& c : claims) {
        if (c.body != body) continue;
        const float apart = Vector3Distance(c.centre, local);
        if (c.team != team && apart < c.radius + kClaimRadius) { why = "Too close to another team's claim"; return 0; }
        if (c.team == team && apart < kClaimRadius) { why = "Your team has claimed this already"; return 0; }
    }
    Claim c;
    c.id = nextId++; c.team = team; c.body = body; c.centre = local; c.created = c.lastSeen = now;
    claims.push_back(c);
    return c.id;
}

bool Claims_Remove(uint32_t id) {
    const size_t n = claims.size();
    claims.erase(std::remove_if(claims.begin(), claims.end(), [id](const Claim& c) { return c.id == id; }), claims.end());
    return claims.size() != n;
}

void Claims_MoveTeam(int from, int to) { if (from && to) for (Claim& c : claims) if (c.team == from) c.team = to; }

int Claims_Update(const std::vector<int>& teamsOnline, int64_t now) {
    int removed = 0;
    for (size_t i = 0; i < claims.size();) {
        Claim& c = claims[i];
        if (std::find(teamsOnline.begin(), teamsOnline.end(), c.team) != teamsOnline.end()) c.lastSeen = now;
        if (now - c.lastSeen > kClaimDecaySeconds) { claims.erase(claims.begin() + i); ++removed; continue; }
        ++i;
    }
    return removed;
}

const std::vector<Claim>& Claims_All() { return claims; }

void Claims_Save(FILE* f) {
    for (const Claim& c : claims)
        fprintf(f, "claim %u %d %s %.3f %.3f %.3f %.2f %lld %lld\n", c.id, c.team, c.body.c_str(), c.centre.x, c.centre.y, c.centre.z, c.radius, (long long)c.created, (long long)c.lastSeen);
}

bool Claims_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "claim") != 0) return false;
    Claim c; char body[128] = {}; long long created = 0, seen = 0;
    if (fscanf(f, " %u %d %127s %f %f %f %f %lld %lld", &c.id, &c.team, body, &c.centre.x, &c.centre.y, &c.centre.z, &c.radius, &created, &seen) == 9) {
        c.body = body; c.created = created; c.lastSeen = seen;
        claims.push_back(c);
        nextId = std::max(nextId, c.id + 1);
    }
    return true;
}

bool Claims_SelfTest() {
    bool ok = true;
    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: claims: %s: NO\n", what); ok = false; } };
    const std::vector<Claim> saved = claims; const uint32_t savedNext = nextId;
    Claims_Reset();
    const Vector3 at{ 0, 60000, 0 }; std::string why;
    const uint32_t a = Claims_Make(1, 1, "terra", at, 1000, why);
    check(a != 0, "a team of one claims");
    check(Claims_At("terra", { 20, 60010, 0 }) != nullptr, "a point 20 m away, 10 m up, is inside");
    check(Claims_At("terra", { 30, 60000, 0 }) == nullptr, "30 m away is outside");
    check(Claims_At("terra", { 0, 60060, 0 }) == nullptr, "60 m up is outside");
    check(Claims_At("luna", { 0, 60000, 0 }) == nullptr, "another body is outside");
    check(Claims_Make(2, 1, "terra", { 40, 60000, 0 }, 1000, why) == 0, "another team can't claim 40 m away (overlap)");
    check(Claims_Make(2, 1, "terra", { 60, 60000, 0 }, 1000, why) != 0, "another team claims 60 m away");
    check(Claims_Make(1, 1, "terra", { 0, 60000, 100 }, 1000, why) != 0, "the team's second claim (budget 2)");
    check(Claims_Make(1, 1, "terra", { 0, 60000, 200 }, 1000, why) == 0, "no third for a team of one");
    check(Claims_Make(1, 2, "terra", { 0, 60000, 200 }, 1000, why) != 0, "a third with a second member");
    Claims_MoveTeam(1, 3);
    check(Claims_OfTeam(3).size() == 3 && Claims_OfTeam(1).empty(), "an emptied team's claims go with its members");
    check(Claims_Update({ 3 }, 1000 + kClaimDecaySeconds + 10) == 1 && Claims_OfTeam(2).empty() && Claims_OfTeam(3).size() == 3, "14 days away: team 2's claim decays, team 3 (online) keeps theirs");
    check(Claims_Remove(a) && Claims_OfTeam(3).size() == 2, "unclaim");
    printf("bench: claims: containment, budget, overlap, team moves, decay: %s\n", ok ? "PASS" : "FAIL");
    claims = saved; nextId = savedNext;
    return ok;
}
