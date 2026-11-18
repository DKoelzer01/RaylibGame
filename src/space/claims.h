#pragma once

// Land claims (#343; the user's rules on #343, 2026-11-16 and 2026-11-18): claims belong to TEAMS (#604). A team claims
// land by placing a CLAIM BANNER (four tiers through the quest tiers: 100 / 200 / 400 / 800 m across); breaking its banner
// unclaims (only the team can: the banner stands in its own claim). Ops keep the /claim commands. A team may hold
// 1 + 2 x (planets and moons it has reached) claims; a team away (no member online) for kClaimDecaySeconds (14 days: the
// graves' clock, #804) loses its claims. A team may grant other teams access (allies). Inside a claim, a player of any
// other team may not edit terrain, build or remove parts, open chests or machines, fell plants, work farmland, use beds or
// touch animals (the server refuses each request); ops may. Graves inside their owner's team claim stay closed to others
// (combat.h gInOwnerTeamClaim). Clients get the list (NetMsg 41) for the map and a border near another team's claim.
//
// A claim is an upright cylinder on a body: its radius around the centre (body-local) across the ground, kClaimHalfHeight
// above and below along the body's up there (deep enough for basements and mines, tall enough for towers).

#include <cstdint>
#include <string>
#include <vector>
#include <cstdio>
#include <raylib.h>

constexpr float kClaimRadius = 50.0f;        // m: a command claim's (ops), as the first banner's
constexpr float kClaimHalfHeight = 128.0f;   // m up and down
constexpr int64_t kClaimDecaySeconds = 14 * 24 * 3600;

struct Claim {
    uint32_t id = 0;
    int team = 0;
    std::string body;          // the body's name
    Vector3 centre{};          // body-local
    float radius = kClaimRadius;
    uint32_t bannerPlot = 0;   // the banner that made it (0: a command claim)
    int bannerPart = -1;
    int64_t created = 0, lastSeen = 0; // wall-clock seconds; lastSeen: a member of the team was online
};

void Claims_Reset();
int Claims_Budget(int bodiesReached);           // 1 + 2 x planets and moons the team has reached
bool Claims_Contains(const Claim& c, Vector3 local);
const Claim* Claims_At(const std::string& body, Vector3 local); // the claim covering that point (nullptr: none)
std::vector<Claim> Claims_OfTeam(int team);
// A new claim for `team` centred at `local`, `radius` across the ground: refused (0, with why) over the budget or
// overlapping another team's claim. A banner's: its plot and part.
uint32_t Claims_Make(int team, int budget, const std::string& body, Vector3 local, float radius, int64_t now, std::string& why,
                     uint32_t bannerPlot = 0, int bannerPart = -1);
bool Claims_Remove(uint32_t id);
void Claims_SetBanner(uint32_t id, uint32_t plot, int part); // (the banner placed after its claim was made)
uint32_t Claims_OfBanner(const std::string& body, uint32_t plot, int part); // the claim a banner made (0: none)
void Claims_MoveTeam(int from, int to);          // a team emptied into another (#604 Teams_Move): its claims go with them
// Allies: `team` lets `other` in (build, open, ...) on its claims.
void Claims_SetAlly(int team, int other, bool on);
bool Claims_Allows(int claimTeam, int team);     // the claim's own team, or one it allies with
std::vector<int> Claims_AlliesOf(int team);
// Decay: the teams with a member online now keep their claims fresh; claims of teams away kClaimDecaySeconds go. Returns
// the claims removed.
int Claims_Update(const std::vector<int>& teamsOnline, int64_t now);
const std::vector<Claim>& Claims_All();
uint32_t Claims_Version();                       // bumps on every change (claims or allies): the server resends the list

void Claims_Save(FILE* f);
bool Claims_LoadLine(const char* tok, FILE* f);

// A client's mirror (NetMsg 41): what the map draws and the border near another team's claim.
struct ClaimView { uint32_t id = 0; int team = 0; std::string teamName; uint32_t colour = 0; std::string body; Vector3 centre{}; float radius = 0.0f; bool ours = false, allowed = false; };
void Claims_SetMirror(const std::vector<ClaimView>& list);
const std::vector<ClaimView>& Claims_Mirror();
uint32_t Claims_TeamColour(int team, uint32_t set); // the team's colour (its own, else one from its id)
class Scene; class Planetoid;
void Claims_DrawBorder(Scene& world); // (inside BeginMode3D) near another team's claim, from the camera

// (#343 tests) the pure rules: containment, budget, overlap, banners, allies, decay, and the lookup's cost. Prints
// "bench: claims: ..." and PASS/FAIL.
bool Claims_SelfTest();
