#ifndef COMPOSED_H
#define COMPOSED_H

// (#1449, #1446, #1219) Composed items: one item per class (tool_pickaxe, weapon_sword, ...); each copy carries its parts in
// its own data (itemdata.h: h head, g handle, b binding materials; q quality %; s sockets; w wear used). Classes come from
// resources/data/composed_classes.txt, materials from part_materials.txt; every stat comes from the parts here, in one place.
// Old per-material tool ids (stone_pickaxe, bronze_axe ...) convert through Composed_FromLegacy.

#include <map>
#include <string>
#include <vector>
#include "tools.h"

struct PartMaterial {
    std::string id, name;
    unsigned char r = 128, g = 128, b = 128;
    int tier = 0;
    float hardness = 0.0f, toughness = 0.0f, density = 0.0f;
    int sockets = 0;
    float factor = 1.0f; // as a handle or binding: its durability factor
};

struct ComposedClass {
    std::string id, role, toolType, weaponClass, headPart, noun;
};

struct ComposedStats {
    bool valid = false;            // a known class with a known head
    std::string cls, head, handle, binding;
    ToolSpec tool;                 // type (Hands for weapons) and tier, from the class and the head
    int durability = 0;            // strokes before it breaks: 25 x head toughness x handle and binding factors x quality
    float damage = 0.0f;           // a light hit's base damage: (3 + 3 x head hardness) x quality
    float digSpeed = 1.0f;         // against its tier's usual (1 for the default parts)
    float weight = 0.0f;           // kg
    int sockets = 0;
    std::vector<std::string> gems; // (#1219) set in its sockets ("j="), at most sockets
    float quality = 1.0f;          // the anvil's quality (q %, 90-110)
    int wear = 0;                  // strokes used (w)
    std::string name;              // "Bronze pickaxe"
};

bool Composed_Load(const char* classesPath = "resources/data/composed_classes.txt", const char* materialsPath = "resources/data/part_materials.txt");
bool Composed_Is(const std::string& id);                    // a composed class
const ComposedClass* Composed_Class(const std::string& id);
const PartMaterial* Composed_Material(const std::string& id);
ComposedStats Composed_Stats(const std::string& id, const std::string& data);
// The record for a new item: head, handle and binding materials ("" handle or binding: the head's default), quality %.
std::string Composed_MakeData(const std::string& head, const std::string& handle = {}, const std::string& binding = {}, int qualityPct = 100);
// (#1449) The assembly hook: the record for a class made from part items, as a grid or anvil recipe consumes them. A part item
// is named <material>_<part>: the class's head part with '_' for spaces (bronze_pick_head, steel_blade, flint_spearhead), or
// handle / binding (wood_handle, leather_binding); the material is a part_materials.txt row. Other ingredients are ignored.
// A missing handle or binding takes the head's default; no head: "" (the recipe is wrong: the item stays plain).
std::string Composed_FromParts(const std::string& cls, const std::vector<std::string>& ingredients, int qualityPct = 100);
// (#1219) The five gems a socket takes, and a gem's words in a tooltip ("Ruby (burns)").
bool Composed_IsGem(const std::string& id);
std::string Composed_GemText(const std::string& id);
// An old per-material tool id (stone_pickaxe, bronze_axe, titanium_hoe) as its class and record; false for anything else.
bool Composed_FromLegacy(const std::string& oldId, std::string& cls, std::string& data);
class Inventory;
// (#1449) The single conversion point for old tools: every cell holding an old per-material tool (stone_pickaxe ...) becomes its
// class with the default parts (one item per cell). wearUsed (the old per-id wear, may be null): each id's wear moves into
// the first copy converted, and its entry goes. Only classes the item table knows are converted. Returns how many.
int Composed_ConvertInventory(Inventory& inv, std::map<std::string, int>* wearUsed);
// (#1449 part 3) What the UI shows for an item: a composed copy's name from its head ("Bronze pickaxe"); "" for anything else
// (the caller shows the item table's name). Its parts in a few lines (head, handle, binding, quality, sockets, wear).
std::string Composed_DisplayName(const std::string& id, const std::string& data);
std::string Composed_PartsText(const std::string& id, const std::string& data);

#endif
