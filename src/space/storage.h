#ifndef STORAGE_H
#define STORAGE_H

// Nearby storage: chests on build plots within reach of the player (8 m). Crafting counts and
// pulls their contents as if they were in the inventory, and the inventory can deposit into them
// (an item goes to a chest that already holds some of it).

#include <cstdint>
#include <string>
#include <vector>

class Scene;

struct NearChest { uint64_t body; uint32_t plot; int part; float dist; };

std::vector<NearChest> Storage_NearbyChests(Scene& world, float radius = 8.0f);
int Storage_CountNearby(Scene& world, const std::string& item);           // in nearby chests (last tick)
// Ask the chests for `count` of an item (it reaches the inventory next tick); false if they lack it.
bool Storage_Pull(Scene& world, const std::string& item, int count);
// Deposit: every inventory item a nearby chest already holds goes into it (tools and the items bound
// to the hotbar stay). Returns how many kinds of item moved.
int Storage_Deposit(Scene& world);

#endif
