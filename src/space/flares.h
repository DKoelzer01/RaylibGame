#ifndef FLARES_H
#define FLARES_H

// (#2062, underground plan #1465) Thrown flares: a cheap item the player throws (RMB with one in hand) that falls on a
// ballistic arc, comes to rest on the ground and burns for a minute. While it burns at rest it is a source in the light
// grid (lightlevel.h), not a point light: it costs no point-light slot (#1739/#1576) and the creatures' and caves' light
// reads see it with no code of their own. Its reach fades to nothing over its last seconds, then it is gone.
// Local only: flares are not saved and do not replicate yet (a client cannot throw one; a host's flare is its own).

#include <raylib.h>
#include <cstdint>
#include <vector>

class Scene;
class Planetoid;

constexpr float kFlareBurnSeconds = 60.0f; // GUESS
constexpr float kFlareFadeSeconds = 10.0f; // the last seconds the reach shrinks to 0 (GUESS)
constexpr float kFlareReach = 16.0f;       // m of light-grid reach at full burn (GUESS: between a torch's 12 and a lamp's 20)

// A flare thrown from `from` (body-local) with `velocity` (body-local m/s); burnSeconds < 0: kFlareBurnSeconds.
void Flares_Throw(Planetoid& body, Vector3 from, Vector3 velocity, float burnSeconds = -1.0f);
// The local player's throw from the eye along the view (RMB with a flare in hand): takes one from the inventory.
bool Flares_ThrowFromPlayer(Scene& world, float burnSeconds = -1.0f); // burnSeconds: as Flares_Throw (tests shorten it)
void Flares_Update(Scene& world, float dt);
void Flares_Draw3D(Scene& world);
void Flares_Clear();

// A flare's light for the grid: its cell's body-local point and its reach now (m; 0 = dark).
struct FlareLight { Vector3 local; float reach; };
void Flares_Lights(const Planetoid& body, std::vector<FlareLight>& out);

// Tests: the flares in the world (any body), and where the first one is.
int Flares_Count();
bool Flares_First(Vector3& local, bool& resting, float& reach);

#endif
