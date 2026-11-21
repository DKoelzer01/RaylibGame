#pragma once

// Magic machines (#1252, #1327; the magic doc sections 1 and 3.2), behind RAYLIBGAME_MAGIC:
// - the mana siphon: set on a ley node (within 8 m), it fills empty mana jars with the node's kind; anywhere else it idles;
// - the mana dynamo: an LV generator that burns jars of mana (machines.txt; nothing to steer);
// - the research table: not a machine; standing by it (6 m) lets the Tome tab spend items for their aspects.

#include <string>

class Scene;
class ManaField;
#include <raylib.h>

// The jar a siphon at this body-local point fills: "jar_<kind>" of the nearest node within 8 m (a compound outside the
// first slice gives its first base kind's jar), or "jar_none" (it idles).
std::string MagicMachines_SiphonJar(const ManaField& f, Vector3 local);
void MagicMachines_Update(Scene& world, double now); // siphons steered, the research table's reach (twice a second)
bool MagicMachines_NearResearchTable();              // the local player stood within 6 m of one at the last update
