#ifndef PLOTS_H
#define PLOTS_H

// Build plots (milestone 3, design doc section 4).
//
// A planet has no global grid ("up" turns as you walk), so structures live on plots: a discrete
// 0.5 m grid anchored where the first part was placed, with up = local gravity and north along the
// meridian. Parts snap to the plot's cells; the plot grows as you build. Terrain stays smooth voxels
// for digging; parts are supported by the plot, not the terrain.
//
// Parts occupy one or more cells and take any of the 24 axis-aligned orientations. Shapes are the
// unit cell cut by up to two planes (cube, slope, corner), so one convex description drives the
// mesh, collision and face culling. Structural parts merge into one mesh per 16^3-cell section
// (8 m) with faces hidden where two full faces touch, rebuilt when a part in that section changes.

#include <raylib.h>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include "dvec3.h"
#include "voxel.h" // Int3

class Planetoid;
struct RenderContext;

constexpr float kCellSize = 0.5f;   // m
constexpr int kSectionCells = 16;   // cells per section edge (8 m)

enum class PartType : uint8_t {
    Cube, Block, Slope, Corner, Beam2, Beam4, Beam8,                 // structure
    Conveyor, Cable,                                                 // item and power connectors (1 cell)
    Chest, Generator, Furnace, Macerator, Battery,                   // machines (1 m, 2x2x2 cells)
    Campfire, Kiln, CokeOven, BlastFurnace, Boiler, ForgeHammer, AlloySmelter, SteamMacerator,
    RocketEngine, FuelTank, Cockpit, Guidance, LandingLeg,           // rocket parts (milestone 7, ship_parts.txt)
    OreWasher, Centrifuge, WaterPump,                                // LV ore processing (milestone 8)
    MvGenerator, Electrolyzer, ChemicalReactor, Miner,               // MV
    HvGenerator, ArcFurnace, Refinery, PrecisionAssembler,           // HV
    Farmland,                                                        // tilled soil a crop grows on (plan 2, milestone 11)
    Fence, Gate,                                                     // 1 m wooden panels: pens (milestone 12); the player walks through a gate, animals do not
    Cabin,                                                           // rocket crew cabin: 4 seats, 30 min of air (plan 3, M19)
    StationCore,                                                     // deployed in orbit: founds a station (plan 3, M20); stores fuel
    Slab, Panel, Pane, Post,                                         // detail shapes inside one cell: half height, half thickness, a thin pane, a quarter post (M21)
    Pipe, Outlet,                                                    // water pipes from pumps, and outlets that pour it out (Plan 4, F6)
    Count
};
inline bool IsMachine(PartType t) { return (t >= PartType::Chest && t <= PartType::SteamMacerator) || (t >= PartType::OreWasher && t <= PartType::PrecisionAssembler); }
inline bool IsRocketPart(PartType t) { return (t >= PartType::RocketEngine && t <= PartType::LandingLeg) || t == PartType::Cabin || t == PartType::StationCore; }
inline bool IsShape(PartType t) { return t <= PartType::Beam8 || (t >= PartType::Slab && t <= PartType::Post); } // structure: built from a material
// The shapes in the build wheel's order (the detail shapes of M21 after the beams).
constexpr PartType kShapes[] = { PartType::Cube, PartType::Block, PartType::Slope, PartType::Corner, PartType::Beam2, PartType::Beam4, PartType::Beam8,
                                 PartType::Slab, PartType::Panel, PartType::Pane, PartType::Post };
constexpr int kShapeCount = (int)(sizeof(kShapes) / sizeof(kShapes[0]));
const char* PartName(PartType t);
int PartCostCount(PartType t);              // items one part costs (shapes: of their material)
const char* PartCostItem(PartType t);       // shapes: concrete (material 0); else the part's own item
const char* PartMachineName(PartType t);    // machines.txt / recipe name ("kiln"), or nullptr

// Building materials: any shape can be built from any of them, costing the material's item.
struct BuildMaterial { const char* item; const char* name; Color color; };
int BuildMaterialCount();
const BuildMaterial& GetBuildMaterial(int index);
int BuildMaterialForItem(const char* item);  // index, or -1

// The 24 axis-aligned rotations (index 0 = identity).
const int* Orientation(int index);          // 3x3 integer matrix, row-major
int RotateOrientation(int index, int axis); // turn a further 90 degrees about plot axis 0/1/2

struct Part {
    PartType type;
    Int3 origin;        // plot cell of the part's origin cell
    uint8_t orient;     // 0..23
    uint8_t material = 0; // shapes: BuildMaterial index
};
const char* PartItem(const Part& p);         // the item it costs and refunds

// A cell's solid shape in cell space [0,1]^3: the cube cut by up to two half-spaces
// (dot(n, p) <= d keeps the solid side).
struct CellShape {
    bool full = true;
    int planes = 0;
    Vector3 n[2];
    float d[2];
};

class Plot {
public:
    uint32_t id = 0;
    uint32_t version = 0;               // bumped on every place / remove (the factory resyncs)
    DVec3 origin;                       // body-local, the corner of cell (0,0,0)
    Vector3 east, up, north;            // body-local unit axes: plot x, y, z

    std::vector<Part> parts;            // index = part id (removed parts have type Count)
    float fuel = 0.0f;                  // rocket fuel loaded (kg), when the plot is a rocket
    std::unordered_map<int, float> damage; // rocket parts: part index -> 0 sound .. 1 out of action (saved)
    float air = 0.0f;                   // rocket life support: minutes of breathable air loaded (saved)
    struct CellRef { int part; CellShape shape; };
    std::unordered_map<Int3, CellRef> cells; // occupied cell -> part index and its shape there

    // Plot space <-> body-local.
    Vector3 toPlot(Vector3 bodyLocal) const;
    Vector3 toBody(Vector3 plotPos) const;
    Int3 cellAt(Vector3 plotPos) const;

    // The cells a part would occupy, with each cell's shape (in plot cell space).
    static void Footprint(const Part& p, std::vector<std::pair<Int3, CellShape>>& out);
    bool canPlace(const Part& p) const;
    int place(const Part& p);           // part index, or -1 if blocked
    bool remove(int partIndex);
    int partAt(Int3 cell) const;
    int partCount() const { return livePartCount; }

    // Rendering (main thread): section meshes rebuilt lazily when dirty.
    void draw(const Planetoid& body, const RenderContext& ctx, Shader shader);
    void drawDepth(const Planetoid& body, const RenderContext& ctx, Material& depth);
    ~Plot();

    // Collision: push a sphere (body-local centre, radius in m) out of parts; returns true and the
    // body-local push and normal if it touched any.
    bool collideSphere(Vector3& bodyLocalCentre, float radius, Vector3& normalOut) const;
    // Ray (body-local) against part cells: the nearest hit cell, its face normal (plot space), distance.
    bool raycast(Vector3 o, Vector3 d, float maxDist, Int3& cell, Int3& faceNormal, float& dist) const;

private:
    struct Section { Mesh mesh{}; bool dirty = true; };
    std::unordered_map<Int3, Section> sections;
    int livePartCount = 0;
    void markDirty(Int3 cell);
    void rebuildDirty();                // every dirty section in one pass over the cells
    void buildSection(Int3 key, Section& s, const std::vector<Int3>& sectionCells);
};

// All plots of one body.
class PlotSet {
public:
    std::vector<std::unique_ptr<Plot>> plots;
    // The plot whose grid should take a part at a body-local point: an existing plot within
    // `reach` of it, or a new one anchored there (up = radial).
    Plot* plotFor(Vector3 bodyLocal, float reach, bool create);
    Plot* byId(uint32_t id);
    Plot* adopt(std::unique_ptr<Plot> plot);          // move a plot in (from another body), with a fresh id
    std::unique_ptr<Plot> release(Plot* plot);        // take one out
    int partCount() const;
    void save(FILE* f) const;
    bool load(FILE* f);
private:
    uint32_t nextId = 1;
};

#endif
