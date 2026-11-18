#pragma once

// Land claims (#343; the user's rules on #343, 2026-11-16): claims belong to TEAMS (#604). A team's claim budget grows with
// its members; a team away (no member online) for kClaimDecaySeconds (14 days: the graves' clock, #804) loses its claims.
// Inside another team's claim a player may not edit terrain, build or remove parts, open chests or machines, fell plants or
// work farmland (the server refuses each request); ops (serveradmin) may. Graves inside their owner's team claim stay closed
// to others (combat.h gInOwnerTeamClaim).
//
// A claim is an upright cylinder on a body: kClaimRadius around its centre (body-local) across the ground, kClaimHalfHeight
// above and below it along the body's up there. Made by a team member where they stand (/claim), removed by one (/unclaim).
// The server keeps them in the world save ("claim" lines); clients hear about claims only by chat for now (design on #343).

#include <cstdint>
#include <string>
#include <vector>
#include <cstdio>
#include <raylib.h>

constexpr float kClaimRadius = 24.0f;       // m across the ground
constexpr float kClaimHalfHeight = 48.0f;   // m up and down
constexpr int64_t kClaimDecaySeconds = 14 * 24 * 3600;

struct Claim {
    uint32_t id = 0;
    int team = 0;
    std::string body;          // the body's name
    Vector3 centre{};          // body-local
    float radius = kClaimRadius;
    int64_t created = 0, lastSeen = 0; // wall-clock seconds; lastSeen: a member of the team was online
};

void Claims_Reset();
int Claims_Budget(int members);                 // claims a team of that many may hold (1 + members, at most 9)
bool Claims_Contains(const Claim& c, Vector3 local);
const Claim* Claims_At(const std::string& body, Vector3 local); // the claim covering that point (nullptr: none)
std::vector<Claim> Claims_OfTeam(int team);
// A new claim for `team` centred at `local`: refused (0, with why) over the budget or overlapping another team's claim.
uint32_t Claims_Make(int team, int members, const std::string& body, Vector3 local, int64_t now, std::string& why);
bool Claims_Remove(uint32_t id);
void Claims_MoveTeam(int from, int to);          // a team emptied into another (#604 Teams_Move): its claims go with them
// Decay: the teams with a member online now keep their claims fresh; claims of teams away kClaimDecaySeconds go. Returns
// the claims removed.
int Claims_Update(const std::vector<int>& teamsOnline, int64_t now);
const std::vector<Claim>& Claims_All();

void Claims_Save(FILE* f);
bool Claims_LoadLine(const char* tok, FILE* f);

// (#343 tests) the pure rules: containment, budget, overlap, decay. Prints "bench: claims: ..." and PASS/FAIL.
bool Claims_SelfTest();
