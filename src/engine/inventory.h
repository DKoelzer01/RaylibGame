#ifndef INVENTORY_H
#define INVENTORY_H

#include <raylib.h>
#include <string>
#include <unordered_map>
#include <vector>

#include "voxel.h"

// Item and recipe definitions, loaded from resources/data/*.txt so content can change without
// recompiling.
struct ItemDef {
    std::string id, name, description;
    int stackLimit = 999;
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
const char* ItemForMaterial(VoxelMaterial m);

class Inventory {
public:
    // Mining produces fractional amounts; whole units appear once enough has been dug.
    void addFractional(const std::string& id, float amount);
    void add(const std::string& id, int count);
    bool remove(const std::string& id, int count);
    int count(const std::string& id) const;
    bool canCraft(const Recipe& r) const;
    bool craft(const Recipe& r);
    void clear() { counts.clear(); partial.clear(); }
    const std::unordered_map<std::string, int>& all() const { return counts; }
private:
    std::unordered_map<std::string, int> counts;
    std::unordered_map<std::string, float> partial;
};

extern Inventory inventory;

#endif
