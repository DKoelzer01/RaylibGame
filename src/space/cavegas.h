#ifndef CAVEGAS_H
#define CAVEGAS_H

// Cave gas pockets (#1472, the underground plan #1465): pockets of bad air in cave passages and chambers, a pure function
// of the world seed: nothing saved, nothing sent over the network; every peer finds the same pockets.
// - CO2 (heavy) fills a pocket's lower part: the air there is not breathable and the suit's oxygen drains 3x (as by a
//   dungeon gas vent).
// - Methane (light) fills its upper part: not breathable, oxygen drains 1.5x, and a flame in it (a held torch, a placed
//   torch or campfire) sets it off: a blast that hurts by distance (no terrain damage: no cave-ins) and burns the pocket
//   off for a while.
// Every pocket warns before it harms: a hiss at its edge, a haze, and a suit warning from a few metres out (further, with a
// direction, with a gas sensor module).
// Multiplayer: the host (or single player) applies damage, in each player's survival step; clients show the warnings and
// the haze from the same field. GUESS values (the user confirms): the chances per band, sizes, the blast and the burn time.

#include <cstdint>
#include <vector>
#include <raylib.h>

class Scene;
class Planetoid;
class PlanetField;

enum class CaveGasKind : uint8_t { None = 0, CO2 = 1, Methane = 2 };

struct CaveGasPocket {
    uint32_t id = 0;          // from the seed and the cave span: the same on every peer
    Vector3 centre{};         // body-local (m)
    Vector3 up{};             // the span's floor up: CO2 lies below centre + 0.3 r along it, methane above centre - 0.3 r
    float radius = 0.0f;      // m
    CaveGasKind kind = CaveGasKind::None;
};

// The pockets whose sphere comes within `range` m of a body-local point (seed-pure; the cave graph's spans).
void CaveGas_PocketsNear(const PlanetField& f, Vector3 local, float range, std::vector<CaveGasPocket>& out);
// Whether a body-local point is in a pocket's gas (no terrain test: the caller's point is a head in the air).
bool CaveGas_InPocket(const CaveGasPocket& p, Vector3 local);
// The gas at a render-space point on `body` (None off a planet, in a burnt-off methane pocket, or with no pocket near).
CaveGasKind CaveGas_At(const Scene& world, const Planetoid* body, Vector3 renderPos);

// Per frame on the main thread, after the player moved: pockets round the player, flames in methane (the authority sets
// them off), the haze, the hiss, the suit warning.
void CaveGas_Update(Scene& world);
void CaveGas_DrawHUD(const Scene& world); // the suit warning line
void CaveGas_Reset();         // a new world: no burnt pockets, no test pockets

// A blast for a player's survival step: the damage for a player at this render-space point from blasts since `lastSeen`
// (a blast id; updated). 0 when none reaches them.
float CaveGas_BlastDamage(Vector3 renderPos, uint32_t& lastSeen);

// Tests: pockets added by hand (body-local, on the player's body; ids from 0xF0000000 up); what the warning shows.
void CaveGas_TestAdd(const CaveGasPocket& p);
struct CaveGasStats { int pockets = 0; int blasts = 0; int puffs = 0; /* haze puffs emitted (since the start) */ int warnKind = 0; float warnDist = -1.0f; float edgeDist = -1.0f; bool hissing = false; }; // edgeDist: the head to the nearest live pocket's gas (m, 0 in it)
const CaveGasStats& CaveGas_Stats();

#endif
