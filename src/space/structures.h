#ifndef STRUCTURES_H
#define STRUCTURES_H

// World structures (#164, first slice): ruins, crash sites and small wonders placed by the world generator.
// Sets are authored files in resources/structures/*.st (parts in plot cells, fill/walls/clear helpers, loot chests);
// loot tables live in resources/data/loot.txt. Each 512 m tile of a planet holds at most one site, decided by a hash
// of the world seed and the tile (so every world gets the same structures in the same places). A site is built as
// its own plot when the player comes within 250 m: turned by a random 90-degree step, sunk into the ground, partly
// decayed (more near the top), its stone swapped for the biome's, the trees cleared and the ground carved to the
// parts (as building does). Its chests are filled when the player is within 60 m. Built sites and unfilled chests
// are saved; a looted or broken structure never comes back.

#include <cstdio>
#include <string>
#include <vector>
#include <raylib.h>

class Scene;
class Planetoid;

void Structures_Update(Scene& world);   // once a second is plenty; cheap when nothing is near
void Structures_Reset();                // a new or loaded world
void Structures_Save(FILE* f);
bool Structures_LoadLine(const char* tok, FILE* f);
int Structures_SetCount();              // sets loaded (tests)
// Tests (--test-structures): every set built at sites round the player; parts placed, chests filled.
bool Structures_Test(Scene& world, Planetoid& body, Vector3 nearLocal);
// Tests: the nearest site (built or not) of any set within `range` of a body-local point; its set name.
const char* Structures_LastSkip();   // tests: why the last site was left out
int Structures_SitesBuilt();          // tests: sites built by the world placement
void Structures_TestFill(Planetoid& body);          // tests: fill every waiting chest now
int Structures_TestCheck(Planetoid& body, int& items); // tests: how many of those chests hold items (and how many)
bool Structures_NearestSite(const Planetoid& body, Vector3 local, float range, Vector3& at, std::string& set);

#endif
