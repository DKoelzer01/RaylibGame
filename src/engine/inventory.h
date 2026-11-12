#ifndef INVENTORY_H
#define INVENTORY_H

#include <raylib.h>
#include <string>
#include <unordered_map>
#include <cstdint>
#include <vector>

#include "voxel.h"

// Item and recipe definitions, loaded from resources/data/*.txt so content can change without
// recompiling.
// What kind of thing an item is decides how it stacks in the inventory grid.
enum class ItemKind : uint8_t {
    Resource, // raw and processed materials: stacks of 8192, amounts to 0.01
    Block,    // placeable build parts (machines, conveyors, rocket parts...) and building blocks: 256
    Tool,     // never stacks
};
struct ItemDef {
    std::string id, name, description;
    ItemKind kind = ItemKind::Resource;
    Color color = GRAY;
};

struct Recipe {
    std::string output;
    int outputCount = 1;
    std::vector<std::pair<std::string, int>> inputs;
};

class ItemDB {
public:
    bool load(const char* itemsPath, const char* recipesPath);
    const ItemDef* find(const std::string& id) const;
    const std::vector<ItemDef>& items() const { return itemList; }
    const std::vector<Recipe>& recipes() const { return recipeList; }
private:
    std::vector<ItemDef> itemList;
    std::unordered_map<std::string, size_t> byId;
    std::vector<Recipe> recipeList;
};

extern ItemDB itemDB;

// What mining a voxel material yields.
// `soil`: a living world (a planet, not a moon or an asteroid): its regolith is dirt.
const char* ItemForMaterial(VoxelMaterial m, bool soil = false);
bool HasSoil(const class Planetoid& p);
// Which ore a vein holds at a body-local point (deposits ~48 m across), from the body's seed.
const char* OreKindAt(Vector3 local, unsigned int seed);
// The ore at a body-local point: planets from their deposits and pockets (PlanetField::oreKind),
// other bodies by the old 48 m deposit hash (OreKindAt).
const char* OreKindOf(const class Planetoid& body, Vector3 local);

ItemKind KindOf(const std::string& id);
int StackLimit(const std::string& id); // whole units per slot
// Recipe wildcards: "any_planks" matches planks of any wood (oak_planks, ...), "any_log" any log.
bool ItemMatches(const std::string& want, const std::string& have);

// One grid cell. Amounts are hundredths (2 decimal places), so resources can be fractional.
struct InvSlot {
    std::string id;
    int64_t centi = 0;
    bool empty() const { return centi <= 0 || id.empty(); }
    double amount() const { return centi / 100.0; }
};

class Inventory {
public:
    static constexpr int kCols = 9, kRows = 5, kSlots = kCols * kRows;
    Inventory() : slots(kSlots) {}
    // Mining produces fractional amounts: hundredths go into the grid, smaller remainders carry over.
    void addFractional(const std::string& id, float amount);
    // Returns what did not fit (the grid is full); `overflowed` is set so the game can say so.
    double addAmount(const std::string& id, double amount);
    void add(const std::string& id, int count) { if (count > 0) addAmount(id, count); }
    bool remove(const std::string& id, int count) { return removeAmount(id, count); }
    bool removeAmount(const std::string& id, double amount);
    int count(const std::string& id) const;     // whole units
    double amount(const std::string& id) const; // with the fraction
    bool canCraft(const Recipe& r) const;
    bool craft(const Recipe& r);
    void clear() { for (auto& s : slots) s = InvSlot{}; partial.clear(); obtained.clear(); }
    // Everything ever obtained (crafted, mined, taken from machines), whatever happened to it since:
    // quest tasks count these. Saved with the game.
    std::unordered_map<std::string, double> obtained;
    double obtainedTotal(const std::string& id) const {
        if (id.rfind("any_", 0) == 0) { double t = 0.0; for (const auto& kv : obtained) if (ItemMatches(id, kv.first)) t += kv.second; return t; }
        auto it = obtained.find(id); return it == obtained.end() ? 0.0 : it->second;
    }
    // Whole units of everything held (a copy).
    std::unordered_map<std::string, int> all() const;
    std::vector<InvSlot>& grid() { return slots; }
    const std::vector<InvSlot>& grid() const { return slots; }
    void swapSlots(int a, int b);  // move/swap two cells, merging same items up to the stack limit
    void sort();                   // by kind, then name; merges partial stacks
    bool overflowed = false;       // set when something did not fit; the game clears it after warning
private:
    // Room for `id` across the grid in hundredths (existing stacks plus empty slots).
    int64_t room(const std::string& id) const;
    std::vector<InvSlot> slots;
    std::unordered_map<std::string, float> partial;
};

extern Inventory inventory;

#endif
