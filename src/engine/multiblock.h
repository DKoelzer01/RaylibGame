#ifndef MULTIBLOCK_H
#define MULTIBLOCK_H

// (#1429, multiblocks PR 1) Multiblock structures: definitions from resources/data/multiblocks.txt and a validator that
// checks a built structure against one. Pure data and logic: no game code calls it yet (the part types and the controller
// come with the multiblock parts PR). The validator reads the world through a lookup function, so tests stamp a std::map
// and the game will pass its cell reader; it reports every checked cell (Ok, Missing, Wrong, Blocked) so a ghost can show
// what to build next, not just a yes or no.
//
// Definition frame: x = column (west to east), y = layer (up), z = row (north to south). The controller faces +z.

#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

struct MbInt3 {
    int x = 0, y = 0, z = 0;
    bool operator==(const MbInt3& o) const { return x == o.x && y == o.y && z == o.z; }
};

enum class MbRole : uint8_t { Air, Any, Controller, Casing, Part, Hatch };

struct MbLegend {
    char key = 0;
    MbRole role = MbRole::Any;
    std::string id;                          // Controller / Casing / Part: the part id
    std::vector<std::string> hatchKinds;     // Hatch: any hatch of one of these kinds
    std::string fallbackCasing;              // Hatch: or this casing (or a better one): a flexible hatch spot
};

struct MbDef {
    std::string id, name;
    int tier = 0;
    std::vector<MbLegend> legend;
    std::vector<std::vector<std::string>> layers; // layers[y][z] = row string (keys, no spaces), all rows equal width
    int stretchLayer = -1, stretchMin = 1, stretchMax = 1;
    std::map<std::string, std::pair<int, int>> hatchLimits; // kind -> min,max
    MbInt3 controller;                                      // its cell in the definition frame (layer y, row z, column x)
    // (#1429 PR 4) The factory: a formed structure runs this machines.txt machine's recipes ("": none, e.g. a tank), up to
    // `parallel` at once.
    std::string machine;
    int parallel = 1;
    const MbLegend* key(char c) const;
};

struct MbCatalog {
    std::vector<MbDef> defs;
    std::map<std::string, int> casingTier;                      // casing id -> tier
    std::map<std::string, std::pair<std::string, int>> hatches; // hatch id -> kind, tier
    // (#1429 PR 4) casing / hatch id -> the voltage tier name ("LV") that caps a structure's tier (an energy hatch's voltage);
    // absent: no cap (a steam-age bronze casing)
    std::map<std::string, std::string> volts;
    const MbDef* find(const std::string& id) const;
    // (#1429 PR 2) ordered part lists: the index is a placed part's Part::material, so the order is the file's (deterministic)
    std::vector<std::string> controllers, casings, hatchIds; // at most 256 each
};

// (#1429 PR 2) The multiblock parts: one PartType per role (MbController, MbCasing, MbHatch), Part::material = the index in
// the catalog's list for that role.
enum class PartType : uint8_t; // plots.h
bool Multiblock_ItemPart(const std::string& item, PartType& type, int& material); // true when the item is a controller, casing or hatch
const char* Multiblock_PartItem(PartType type, int material);                     // its item id, or "" when out of range

// Parse; errors appended as "line N: message". Bad blocks are skipped, good ones kept.
MbCatalog Multiblock_Parse(const std::string& text, std::vector<std::string>* errors);
const MbCatalog& Multiblock_All(); // resources/data/multiblocks.txt, loaded once (function-local static: thread-safe)

enum class MbCell : uint8_t { Ok, Missing, Wrong, Blocked }; // Missing: empty where a part goes; Wrong: another part; Blocked: something where air must be
struct MbCellResult { MbInt3 world; int layer = 0; MbCell status = MbCell::Ok; std::string want; char key = 0; }; // want: what goes there (an id or "air"); key: its legend key; for the ghost

struct MbResult {
    bool formed = false;
    std::vector<MbCellResult> cells;   // every checked cell (Any cells are skipped)
    int stretch = 0;                   // repeats of the stretch layer found (complete ones)
    int nextLayer = -1;                // the lowest layer index (in the stretched structure) with a cell not Ok; -1 none
    std::map<std::string, int> hatchCount;
    std::string problem;               // the first problem in words ("layer 2: 3 cells missing", "needs at least 1 fluid_in hatch"), "" when formed
};

// lookup(world cell) -> the part id there, "" = empty. controllerAt = the controller's world cell; facing = the
// controller's front as a horizontal unit vector (fx,fz) in {(0,1),(1,0),(0,-1),(-1,0)}.
MbResult Multiblock_Validate(const MbDef& def, const MbCatalog& cat, const std::function<std::string(MbInt3)>& lookup,
                             MbInt3 controllerAt, int fx, int fz);
// (#1429 PR 5) What to build where, for the ghost and the controller's parts list. A hatch spot gets a hatch while the
// structure is short of that kind's minimum (the catalog's first hatch of the kind), else its fallback casing; a shortfall
// left when every hatch spot is already filled (with casings) is still counted in `need`, as hatches to swap in.
struct MbPlan {
    std::vector<std::string> part;   // parallel to MbResult::cells: the item to place there ("" when Ok or to be cleared)
    std::map<std::string, int> need; // item -> count still to place (wrong parts count: they are replaced)
};
MbPlan Multiblock_Plan(const MbDef& def, const MbCatalog& cat, const MbResult& r);

// (#1429 PR 3) Forming. The controller is the plot's part `partIndex` (an MbController): its definition is the one whose
// controller legend id is its item, its facing is the part's local +z turned by its orientation (a controller facing up or
// down has none: not formable), and the world is the plot's cells (the part covering a cell, whatever its footprint).
// `def` (if not null) gets the definition used, nullptr when the part is not a controller with a definition.
class Plot;
struct Part;
MbResult Multiblock_CheckController(const Plot& plot, int partIndex, const MbDef** def = nullptr);
const MbDef* Multiblock_DefForController(const Part& controller); // nullptr when the part is not a controller with a definition
// Mark / unmark a controller as formed (keeps Plot::formed sorted and unique); true when the list changed.
bool Multiblock_SetFormed(Plot& plot, int partIndex, bool formed);
// Drop every formed index whose structure no longer validates (or whose part is gone or not a controller). Returns the
// dropped indices in ascending order; the caller broadcasts them.
std::vector<int> Multiblock_Recheck(Plot& plot);

// (#1429 PR 4) A formed structure as the factory runs it: its definition, its hatches (part index and kind, in the order the
// validator checked them, each once) and the voltage names capping its tier (its casings' and energy hatches', each once).
// False when the controller does not validate (the caller treats it as unformed).
struct MbMachine {
    const MbDef* def = nullptr;
    std::vector<std::pair<int, std::string>> hatches;
    std::vector<std::string> volts;
};
bool Multiblock_Members(const Plot& plot, int controllerPart, MbMachine& out);
// A hatch part's kind ("item_in", ...) from its Part::material; "" when out of range.
std::string Multiblock_HatchKind(int material);

// A definition cell in the world. Cells in layers above the stretch layer shift up by stretchExtra (repeats - 1).
MbInt3 Multiblock_ToWorld(const MbDef& def, MbInt3 defCell, MbInt3 controllerAt, int fx, int fz, int stretchExtra = 0);

#endif
