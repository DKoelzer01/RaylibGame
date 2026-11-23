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
    const MbLegend* key(char c) const;
};

struct MbCatalog {
    std::vector<MbDef> defs;
    std::map<std::string, int> casingTier;                      // casing id -> tier
    std::map<std::string, std::pair<std::string, int>> hatches; // hatch id -> kind, tier
    const MbDef* find(const std::string& id) const;
};

// Parse; errors appended as "line N: message". Bad blocks are skipped, good ones kept.
MbCatalog Multiblock_Parse(const std::string& text, std::vector<std::string>* errors);
const MbCatalog& Multiblock_All(); // resources/data/multiblocks.txt, loaded once (function-local static: thread-safe)

enum class MbCell : uint8_t { Ok, Missing, Wrong, Blocked }; // Missing: empty where a part goes; Wrong: another part; Blocked: something where air must be
struct MbCellResult { MbInt3 world; int layer = 0; MbCell status = MbCell::Ok; std::string want; }; // want: what goes there (an id or "air"); for the ghost

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
// A definition cell in the world. Cells in layers above the stretch layer shift up by stretchExtra (repeats - 1).
MbInt3 Multiblock_ToWorld(const MbDef& def, MbInt3 defCell, MbInt3 controllerAt, int fx, int fz, int stretchExtra = 0);

#endif
