#ifndef INVENTORY_H
#define INVENTORY_H

#include <raylib.h>
#include <string>
#include <unordered_map>
#include <cstdint>
#include <vector>
#include <functional>

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
    bool formShape = false; // (#1614) made by the form table (forms.txt): a component shape, raw material to disassembly
};

// An item a save may still hold that the game no longer has, and what it becomes on load ("" = keep the id).
// #1090: lunacrete ("concrete") left the planet game; old stacks turn into stone, 1 for 1.
const char* RetiredItem(const std::string& id);

class ItemDB {
public:
    bool load(const char* itemsPath, const char* recipesPath);
    // More items in items.txt's format, added after load (#1250: magic_items.txt, only with RAYLIBGAME_MAGIC); ids already
    // known are kept as they were. False if the file can't be read.
    bool loadExtra(const char* itemsPath);
    bool loadExtraText(const std::string& text); // (#1254) the same from text (items the game makes from its data: alchemy)
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
// Wildcards (any_planks, any_log) are recipe tags matched by ItemMatches, not items: never added to an inventory or
// listed in the Creative tab (#269).
inline bool IsWildcardItem(const std::string& id) { return id.rfind("any_", 0) == 0; }

// One grid cell. Amounts are hundredths (2 decimal places), so resources can be fractional.
struct InvSlot {
    std::string id;
    int64_t centi = 0;
    float heat = 0.0f; // #580: 1 fresh out of a furnace, kiln or crucible, cooling to 0 (Inventory::cool); merges average it by amount
    std::string data;  // (#1446) this one item's own record (itemdata.h: parts, quality, sockets, wear); "" = a plain item.
                       // A cell with data holds exactly one unit and never merges with another.
    bool empty() const { return centi <= 0 || id.empty(); }
    bool hasData() const { return !data.empty(); }
    double amount() const { return centi / 100.0; }
};

// How an item glows (#580): heat (a hot ingot: dark red -> orange -> yellow-white as it gets hotter) or a constant glow
// from resources/data/item_glow.txt (radioactive materials). strength 0 = no glow. Used by the held item; the inventory,
// hotbar and machine-slot icons can draw it over the icon.
struct ItemGlow { Color color = BLACK; float strength = 0.0f; };
ItemGlow Item_GlowOf(const std::string& id, float heat);
inline ItemGlow Item_Glow(const InvSlot& s) { return s.empty() ? ItemGlow{} : Item_GlowOf(s.id, s.heat); }

class Inventory {
public:
    // The 10 x 5 grid (as wide as the hotbar, #495), then the hotbar's 10 cells (#209: the hotbar holds real stacks;
    // keys 1-0 = cells kGrid..kSlots-1).
    static constexpr int kCols = 10, kRows = 5, kGrid = kCols * kRows, kHotbar = 10, kSlots = kGrid + kHotbar;
    // Saves before #495 had a 9 x 5 grid: their hotbar and off-hand cells sit kGrid - kOldGrid lower.
    static constexpr int kOldGrid = 45;
    static int SavedCell(int i, int savedCols) { return savedCols != kCols && i >= kOldGrid ? i + kGrid - kOldGrid : i; }
    // The off-hand cell (#165) after them: what the other hand holds (a shield: RMB blocks; a torch later). Nothing is
    // filled into it or taken from it except by hand, and it is used last when paying for something.
    static constexpr int kOffHand = kSlots, kCells = kSlots + 1;
    int preferred = -1; // a cell used first when taking items out (the selected hotbar cell: the tool in hand wears)
    Inventory() : slots(kCells) {}
    const InvSlot& offHand() const { return slots[kOffHand]; }
    // Mining produces fractional amounts: hundredths go into the grid, smaller remainders carry over.
    void addFractional(const std::string& id, float amount);
    // Returns what did not fit (the grid is full); `overflowed` is set so the game can say so.
    // heat (#580): how hot what comes in is (0 cold .. 1 straight out of a furnace); stacks it joins average it by amount.
    // credit false (#1728): don't count it in `obtained` (a stack already counted when it went onto the cursor comes back in).
    double addAmount(const std::string& id, double amount, float heat = 0.0f, bool credit = true);
    // (#1728) Count what was obtained without adding it to a cell: a machine take that goes straight onto the cursor never passes
    // addAmount. Credits `obtained` exactly as addAmount / addInstance do (an item with data: one, plus the old tool id a composed
    // tool counts as); never the pickup feed or onAdded (the cursor take is not a pickup; a server credits its joined players itself).
    void noteObtained(const std::string& id, double amount, const std::string& data = std::string());
    void cool(float dt);           // hot stacks cool to nothing over kCoolSeconds (#580)
    static constexpr float kCoolSeconds = 90.0f;
    void add(const std::string& id, int count) { if (count > 0) addAmount(id, count); }
    bool remove(const std::string& id, int count) { return removeAmount(id, count); }
    // Plain items only (#1446): add, remove, count and amount never touch a cell with data, so paying for a recipe or a
    // build never takes someone's weapon. Items with data go in with addInstance and out with takeInstance (a cell).
    bool removeAmount(const std::string& id, double amount);
    int count(const std::string& id) const;     // whole units
    double amount(const std::string& id) const; // with the fraction
    // (#1446) One item with its own record into an empty cell (a tool prefers an empty hotbar cell, as plain tools do).
    // Returns the cell, or -1 when no cell is free (overflowed is set). data "" adds one plain unit instead.
    int addInstance(const std::string& id, const std::string& data, float heat = 0.0f, bool credit = true); // credit false (#1728): see addAmount
    InvSlot takeInstance(int cell);             // the cell's item, record and all; the cell is emptied (an empty slot if none)
    int countInstances(const std::string& id) const; // items of `id` that carry data
    int held(const std::string& id) const { return count(id) + countInstances(id); } // (#1449) possession: plain ones and copies with data
    bool setInstanceData(int cell, const std::string& newData); // a held item's record changed (wear, sockets); false if no such item
    int findInstance(const std::string& id, const std::string& data) const; // the cell holding exactly that item, -1 none
    // (#1446) With `journal` on (a server, for each remote player), every change to an item with data is logged here:
    // +1 added, -1 taken, 0 its record changed to newData. The server drains it into InventoryDelta messages each tick.
    struct InstanceChange { std::string id; int sign = 0; std::string data, newData; };
    std::vector<InstanceChange> instanceLog;
    bool journal = false;
    bool canCraft(const Recipe& r) const;
    bool craft(const Recipe& r);
    // (#1613) Takes n of a composed class's head part: forged heads (with a record) first, the best quality first (then the
    // lowest cell), then plain ones. Returns the quality the tool gets (the first head taken; 100 for a plain one), or -1
    // (nothing taken) when there aren't n.
    int takeHeadParts(const std::string& id, int n);
    void clear() {
        if (journal) for (const InvSlot& s : slots) if (!s.empty() && s.hasData()) instanceLog.push_back({ s.id, -1, s.data, {} }); // (#1446)
        for (auto& s : slots) s = InvSlot{};
        partial.clear(); obtained.clear();
    }
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
    void sort();                   // the main grid by kind, then name; merges partial stacks (the hotbar stays as it is)
    bool overflowed = false;       // set when something did not fit; the game clears it after warning
    // Pickup feed (#208): called with what actually went in (not what overflowed); the game sets it on the player's
    // inventory. muted > 0 (loading a save) skips it.
    std::function<void(const std::string&, double)> onAdded;
    int muted = 0;
private:
    // Room for `id` across the grid in hundredths (existing stacks plus empty slots).
    int64_t room(const std::string& id) const;
    std::vector<InvSlot> slots;
    std::unordered_map<std::string, float> partial;
};

extern Inventory& inventory;   // the local PlayerSession's (#313)

#endif
