#ifndef EVENTS_H
#define EVENTS_H

// World events (#177), first: meteor showers. Each 10 minutes of world time the authority (single-player or the host)
// rolls kMeteorShowerOdds on the gameplay RNG for ONE world-wide shower (#316, user 2026-11-15): two minutes of streaks
// across the sky for every player whose sky is dark (cosmetic, per viewer), and two meteors that come down 80-180 m from
// a player or an interest centre (picked by the same RNG) on dry land. Each impact flashes, throws up dust, carves a
// crater and leaves a hot, glowing meteorite (PartType::Meteorite) at its bottom: broken apart it gives nickel-iron
// and now and then platinum or iridium ore, otherwise found only in space. A client never rolls: it applies the host's
// shower (Events_ApplyShower) and receives the impacts as world edits.

#include <cstdint>
#include <cstdio>

// Tuning (user, 2026-11-15: 5% for now, to be retuned): the chance of a shower per bucket of world time.
constexpr float kMeteorShowerOdds = 0.05f;
constexpr double kMeteorShowerBucketSeconds = 600.0; // 10 game minutes

class Scene;
class Planetoid;

struct EventsShower { double start = -1.0, end = -1.0; uint32_t id = 0; }; // what the host sends its clients

void Events_Update(Scene& world, float dt);  // the authority's roll and impacts, then this viewer's streaks
void Events_Draw3D(Scene& world);          // the streaks
void Events_AddLights(Scene& world);       // impact flashes and hot meteorites into world.localLights (#139)
void Events_Reset();
void Events_Save(FILE* f);
bool Events_LoadLine(const char* tok, FILE* f);
bool Events_ShowerOn();
void Events_ForceShower();                 // tests and debugging: a shower now
EventsShower Events_GetShower();
void Events_ApplyShower(const EventsShower& s); // a client: the host's shower
uint32_t Events_ShowerCount();             // showers so far in this world
void Events_ImpactTargets(int& atPlayer, int& atCentres); // where meteorites landed so far (tests)
// Tests: a meteor comes down `dist` metres from the player on a bearing (radians): true if it landed (dry ground).
bool Events_TestImpact(Scene& world, Planetoid& body, float dist, float bearing);
int Events_Streaks();                      // streaks drawn so far (tests)

#endif
