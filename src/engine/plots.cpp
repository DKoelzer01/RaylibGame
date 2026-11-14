#include "plots.h"
#include "worldtypes.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>

// ---------------------------------------------------------------------------------------------
// Part definitions
// ---------------------------------------------------------------------------------------------

const char* PartName(PartType t) {
    switch (t) {
        case PartType::Cube: return "cube";
        case PartType::Block: return "block";
        case PartType::Slope: return "slope";
        case PartType::Corner: return "corner";
        case PartType::Beam2: return "beam 1 m";
        case PartType::Beam4: return "beam 2 m";
        case PartType::Beam8: return "beam 4 m";
        case PartType::Slab: return "slab";
        case PartType::Panel: return "panel";
        case PartType::Pane: return "pane";
        case PartType::Post: return "post";
        case PartType::Conveyor: return "conveyor";
        case PartType::Pipe: return "water pipe";
        case PartType::Outlet: return "water outlet";
        case PartType::Cable: return "cable";
        case PartType::Chest: return "chest";
        case PartType::Generator: return "generator";
        case PartType::Furnace: return "furnace";
        case PartType::Macerator: return "macerator";
        case PartType::Battery: return "battery";
        case PartType::Campfire: return "campfire";
        case PartType::Kiln: return "kiln";
        case PartType::CokeOven: return "coke oven";
        case PartType::BlastFurnace: return "blast furnace";
        case PartType::Boiler: return "steam boiler";
        case PartType::ForgeHammer: return "forge hammer";
        case PartType::AlloySmelter: return "alloy smelter";
        case PartType::SteamMacerator: return "steam macerator";
        case PartType::RocketEngine: return "rocket engine";
        case PartType::FuelTank: return "fuel tank";
        case PartType::Cockpit: return "cockpit";
        case PartType::Guidance: return "guidance computer";
        case PartType::LandingLeg: return "landing leg";
        case PartType::OreWasher: return "ore washer";
        case PartType::Centrifuge: return "centrifuge";
        case PartType::WaterPump: return "water pump";
        case PartType::MvGenerator: return "MV generator";
        case PartType::Electrolyzer: return "electrolyzer";
        case PartType::ChemicalReactor: return "chemical reactor";
        case PartType::Miner: return "miner";
        case PartType::HvGenerator: return "HV generator";
        case PartType::ArcFurnace: return "arc furnace";
        case PartType::Refinery: return "refinery";
        case PartType::PrecisionAssembler: return "precision assembler";
        case PartType::FusionReactor: return "fusion reactor";
        case PartType::MassDriver: return "mass driver";
        case PartType::HandPump: return "hand pump";
        case PartType::Farmland: return "farmland";
        case PartType::Fence: return "fence";
        case PartType::Gate: return "gate";
        case PartType::Cabin: return "cabin";
        case PartType::StationCore: return "station_core";
        default: return "?";
    }
}
int PartCostCount(PartType t) {
    switch (t) {
        case PartType::Block: return 4;
        case PartType::Beam2: return 1;
        case PartType::Beam4: return 2;
        case PartType::Beam8: return 4;
        default: return 1; // cube, slope, corner; connectors and machines: one of their own item
    }
}
const char* PartMachineName(PartType t) {
    switch (t) {
        case PartType::Chest: return "chest";
        case PartType::Generator: return "generator";
        case PartType::Furnace: return "furnace";
        case PartType::Macerator: return "macerator";
        case PartType::Battery: return "battery";
        case PartType::Campfire: return "campfire";
        case PartType::Kiln: return "kiln";
        case PartType::CokeOven: return "coke_oven";
        case PartType::BlastFurnace: return "blast_furnace";
        case PartType::Boiler: return "boiler";
        case PartType::ForgeHammer: return "forge_hammer";
        case PartType::AlloySmelter: return "alloy_smelter";
        case PartType::SteamMacerator: return "steam_macerator";
        case PartType::OreWasher: return "ore_washer";
        case PartType::Centrifuge: return "centrifuge";
        case PartType::WaterPump: return "water_pump";
        case PartType::MvGenerator: return "mv_generator";
        case PartType::Electrolyzer: return "electrolyzer";
        case PartType::ChemicalReactor: return "chemical_reactor";
        case PartType::Miner: return "miner";
        case PartType::HvGenerator: return "hv_generator";
        case PartType::ArcFurnace: return "arc_furnace";
        case PartType::Refinery: return "refinery";
        case PartType::PrecisionAssembler: return "precision_assembler";
        case PartType::FusionReactor: return "fusion_reactor";
        case PartType::MassDriver: return "mass_driver";
        default: return nullptr;
    }
}
const char* PartCostItem(PartType t) {
    switch (t) {

        case PartType::Conveyor: return "conveyor";
        case PartType::Pipe: return "pipe";
        case PartType::Outlet: return "water_outlet";
        case PartType::Cable: return "cable";
        case PartType::Furnace: return "electric_furnace"; // the item; "furnace" is its recipe map
        case PartType::Farmland: return "farmland";
        case PartType::Fence: return "fence";
        case PartType::Gate: return "gate";
        case PartType::Cabin: return "cabin";
        case PartType::StationCore: return "station_core";
        case PartType::RocketEngine: return "rocket_engine";
        case PartType::FuelTank: return "fuel_tank";
        case PartType::Cockpit: return "cockpit";
        case PartType::Guidance: return "guidance_computer";
        case PartType::LandingLeg: return "landing_leg";
        case PartType::HandPump: return "hand_pump"; // not a machine: no factory entry (#280)
        default: { const char* m = PartMachineName(t); return m ? m : "concrete"; }
    }
}

namespace {
    const BuildMaterial kBuildMaterials[] = {
        { "concrete", "Concrete", { 168, 166, 158, 255 } },
        { "stone", "Stone", { 122, 122, 128, 255 } },
        { "iron_plate", "Iron", { 118, 122, 130, 255 } },
        { "glass", "Glass", { 172, 212, 226, 255 } },
        { "dirt", "Packed dirt", { 120, 90, 60, 255 } },
        { "heat_shield_tile", "Heat shield", { 45, 42, 40, 255 } },  // survives re-entry on a rocket's leading face
        // Building variety (M21): planks from each tree species, stone kinds from the ground, kiln bricks.
        { "oak_planks", "Oak", { 158, 112, 66, 255 } },
        { "birch_planks", "Birch", { 214, 196, 152, 255 } },
        { "pine_planks", "Pine", { 176, 132, 82, 255 } },
        { "maple_planks", "Maple", { 182, 104, 78, 255 } },
        { "spruce_planks", "Spruce", { 112, 82, 58, 255 } },
        { "granite", "Granite", { 168, 146, 140, 255 } },
        { "basalt", "Basalt", { 62, 62, 68, 255 } },
        { "sandstone", "Sandstone", { 212, 182, 132, 255 } },
        { "limestone", "Limestone", { 206, 202, 186, 255 } },
        { "marble", "Marble", { 236, 234, 230, 255 } },
        { "brick", "Brick", { 158, 74, 52, 255 } },
        // Decoration (M22): plaster in eight dye colours, stained glass, lamps, polished gems.
        { "white_plaster", "White plaster", { 236, 234, 226, 255 } },
        { "black_plaster", "Black plaster", { 40, 40, 44, 255 } },
        { "red_plaster", "Red plaster", { 178, 52, 46, 255 } },
        { "orange_plaster", "Orange plaster", { 222, 128, 48, 255 } },
        { "yellow_plaster", "Yellow plaster", { 232, 200, 70, 255 } },
        { "green_plaster", "Green plaster", { 84, 142, 66, 255 } },
        { "blue_plaster", "Blue plaster", { 58, 92, 172, 255 } },
        { "purple_plaster", "Purple plaster", { 124, 70, 150, 255 } },
        { "red_glass", "Red glass", { 210, 110, 110, 255 } },
        { "yellow_glass", "Yellow glass", { 236, 218, 120, 255 } },
        { "green_glass", "Green glass", { 130, 200, 130, 255 } },
        { "blue_glass", "Blue glass", { 120, 150, 222, 255 } },
        { "lamp", "Lamp", { 255, 244, 200, 255 } },
        { "quartz_block", "Quartz", { 240, 236, 244, 255 } },
        { "amethyst_block", "Amethyst", { 150, 100, 200, 255 } },
        { "ruby_block", "Ruby", { 196, 40, 70, 255 } },
        { "sapphire_block", "Sapphire", { 50, 80, 196, 255 } },
        { "emerald_block", "Emerald", { 40, 176, 110, 255 } },
        // Raw resources as blocks (playtest 4): logs of each wood (bark colours), clay, ice and coal.
        { "oak_log", "Oak log", { 92, 66, 42, 255 } },
        { "birch_log", "Birch log", { 226, 222, 210, 255 } },
        { "maple_log", "Maple log", { 110, 64, 50, 255 } },
        { "pine_log", "Pine log", { 104, 74, 50, 255 } },
        { "spruce_log", "Spruce log", { 70, 52, 40, 255 } },
        { "clay", "Clay", { 170, 120, 95, 255 } },
        { "ice", "Ice", { 176, 214, 240, 255 } },
        { "coal", "Coal", { 36, 36, 40, 255 } },
    };
}
int BuildMaterialCount() { return (int)(sizeof(kBuildMaterials) / sizeof(kBuildMaterials[0])); }
const BuildMaterial& GetBuildMaterial(int i) { return kBuildMaterials[std::clamp(i, 0, BuildMaterialCount() - 1)]; }
int BuildMaterialForItem(const char* item) {
    for (int i = 0; i < BuildMaterialCount(); ++i) if (strcmp(kBuildMaterials[i].item, item) == 0) return i;
    return -1;
}
const char* PartItem(const Part& p) { return IsShape(p.type) ? GetBuildMaterial(p.material).item : PartCostItem(p.type); }
static Color PartColor(PartType t) {
    switch (t) {

        case PartType::Block: return { 150, 148, 142, 255 };
        case PartType::Conveyor: return { 60, 62, 66, 255 };
        case PartType::Pipe: return { 96, 116, 140, 255 };
        case PartType::Outlet: return { 60, 120, 190, 255 };
        case PartType::Cable: return { 184, 110, 60, 255 };       // copper
        case PartType::Chest: return { 150, 105, 60, 255 };       // wood
        case PartType::Generator: return { 180, 60, 45, 255 };
        case PartType::Furnace: return { 95, 92, 88, 255 };
        case PartType::Macerator: return { 70, 110, 150, 255 };
        case PartType::Battery: return { 210, 180, 60, 255 };
        case PartType::Campfire: return { 200, 90, 40, 255 };
        case PartType::Kiln: return { 150, 95, 70, 255 };
        case PartType::CokeOven: return { 120, 80, 60, 255 };
        case PartType::BlastFurnace: return { 90, 70, 70, 255 };
        case PartType::Boiler: return { 160, 110, 60, 255 };
        case PartType::ForgeHammer: return { 150, 120, 70, 255 };
        case PartType::AlloySmelter: return { 170, 110, 70, 255 };
        case PartType::SteamMacerator: return { 140, 120, 90, 255 };
        case PartType::RocketEngine: return { 70, 72, 80, 255 };
        case PartType::FuelTank: return { 220, 222, 226, 255 };
        case PartType::Cockpit: return { 60, 120, 170, 255 };
        case PartType::Guidance: return { 40, 140, 90, 255 };
        case PartType::LandingLeg: return { 150, 150, 150, 255 };
        case PartType::OreWasher: return { 70, 120, 170, 255 };
        case PartType::Centrifuge: return { 160, 160, 170, 255 };
        case PartType::WaterPump: return { 50, 90, 160, 255 };
        case PartType::MvGenerator: return { 200, 90, 50, 255 };
        case PartType::Electrolyzer: return { 180, 180, 80, 255 };
        case PartType::ChemicalReactor: return { 120, 160, 120, 255 };
        case PartType::Miner: return { 110, 90, 70, 255 };
        case PartType::HvGenerator: return { 210, 60, 60, 255 };
        case PartType::ArcFurnace: return { 90, 90, 110, 255 };
        case PartType::Refinery: return { 140, 110, 160, 255 };
        case PartType::PrecisionAssembler: return { 100, 140, 180, 255 };
        case PartType::FusionReactor: return { 200, 110, 240, 255 };
        case PartType::MassDriver: return { 150, 160, 185, 255 };
        case PartType::HandPump: return { 120, 90, 60, 255 };
        case PartType::Farmland: return { 86, 60, 38, 255 };
        case PartType::Fence: return { 130, 95, 60, 255 };
        case PartType::Gate: return { 165, 120, 70, 255 };
        case PartType::Cabin: return { 200, 205, 215, 255 };
        case PartType::StationCore: return { 230, 220, 120, 255 };
        default: return { 168, 166, 158, 255 };                                                          // concrete
    }
}

// Surface pattern of a part (parts.fs draws it procedurally; the id rides in the vertex colour's alpha).
enum PartPattern : unsigned char {
    PatConcrete = 0, PatStone, PatMetal, PatGlass, PatDirt, PatTiles, PatPlanks, PatPolished, PatBrick, PatPlaster,
    PatGem, PatLog, PatClay, PatIce, PatCoal, PatLamp, PatMarble,
    PatMachine = 20, PatConveyor, PatCable, PatFarmland, PatHull
};
static unsigned char PartPatternOf(const Part& p) {
    if (IsShape(p.type)) {
        const std::string item = GetBuildMaterial(p.material).item;
        auto ends = [&](const char* s) { const size_t n = strlen(s); return item.size() >= n && item.compare(item.size() - n, n, s) == 0; };
        if (ends("_planks")) return PatPlanks;
        if (ends("_log")) return PatLog;
        if (ends("_plaster")) return PatPlaster;
        if (ends("_glass") || item == "glass") return PatGlass;
        if (ends("_block")) return PatGem;
        if (item == "stone") return PatStone;
        if (item == "iron_plate") return PatMetal;
        if (item == "dirt") return PatDirt;
        if (item == "heat_shield_tile") return PatTiles;
        if (item == "brick") return PatBrick;
        if (item == "marble" || item == "limestone") return PatMarble;
        if (item == "granite" || item == "basalt" || item == "sandstone") return PatPolished;
        if (item == "clay") return PatClay;
        if (item == "ice") return PatIce;
        if (item == "coal") return PatCoal;
        if (item == "lamp") return PatLamp;
        return PatConcrete;
    }
    switch (p.type) {
        case PartType::Chest: case PartType::Fence: case PartType::Gate: return PatPlanks;
        case PartType::Kiln: case PartType::CokeOven: case PartType::BlastFurnace: case PartType::Campfire: return PatBrick;
        case PartType::Conveyor: return PatConveyor;
        case PartType::Cable: return PatCable;
        case PartType::Pipe: case PartType::Outlet: return PatMetal;
        case PartType::Farmland: return PatFarmland;
        default: return IsRocketPart(p.type) ? PatHull : PatMachine;
    }
}

// Cells of a part in part space, each with its shape.
static void PartCells(PartType t, std::vector<std::pair<Int3, CellShape>>& out) {
    out.clear();
    CellShape full;
    switch (t) {
        case PartType::Cube: case PartType::Conveyor: case PartType::Cable: case PartType::Guidance: case PartType::LandingLeg: case PartType::Pipe: case PartType::Outlet:
            out.push_back({ { 0, 0, 0 }, full }); break;
        case PartType::Block: case PartType::Chest: case PartType::Generator: case PartType::Furnace: case PartType::Macerator: case PartType::Battery:
        case PartType::Campfire: case PartType::Kiln: case PartType::CokeOven: case PartType::BlastFurnace: case PartType::Boiler:
        case PartType::ForgeHammer: case PartType::AlloySmelter: case PartType::SteamMacerator:
        case PartType::RocketEngine: case PartType::FuelTank: case PartType::Cockpit: case PartType::Cabin: case PartType::StationCore:
        case PartType::OreWasher: case PartType::Centrifuge: case PartType::WaterPump: case PartType::MvGenerator: case PartType::Electrolyzer:
        case PartType::ChemicalReactor: case PartType::Miner: case PartType::HvGenerator: case PartType::ArcFurnace: case PartType::Refinery:
        case PartType::PrecisionAssembler: case PartType::FusionReactor: case PartType::MassDriver: case PartType::HandPump:
            for (int z = 0; z < 2; ++z) for (int y = 0; y < 2; ++y) for (int x = 0; x < 2; ++x) out.push_back({ { x, y, z }, full });
            break;
        case PartType::Farmland: // 1 m square, half a metre deep
            for (int z = 0; z < 2; ++z) for (int x = 0; x < 2; ++x) out.push_back({ { x, 0, z }, full });
            break;
        case PartType::Fence: case PartType::Gate: { // a 10 cm board across the middle, 1 m long and 1 m high
            CellShape s; s.full = false; s.planes = 2; s.n[0] = { 1, 0, 0 }; s.d[0] = 0.6f; s.n[1] = { -1, 0, 0 }; s.d[1] = -0.4f;
            for (int z = 0; z < 2; ++z) for (int y = 0; y < 2; ++y) out.push_back({ { 0, y, z }, s });
            break;
        }
        case PartType::Slab: {  // the lower half of the cell
            CellShape s; s.full = false; s.planes = 1; s.n[0] = { 0, 1, 0 }; s.d[0] = 0.5f;
            out.push_back({ { 0, 0, 0 }, s });
            break;
        }
        case PartType::Panel: { // half the cell's thickness, against one face
            CellShape s; s.full = false; s.planes = 1; s.n[0] = { 1, 0, 0 }; s.d[0] = 0.5f;
            out.push_back({ { 0, 0, 0 }, s });
            break;
        }
        case PartType::Pane: {  // a thin board through the middle (glass panes, screens)
            CellShape s; s.full = false; s.planes = 2; s.n[0] = { 1, 0, 0 }; s.d[0] = 0.56f; s.n[1] = { -1, 0, 0 }; s.d[1] = -0.44f;
            out.push_back({ { 0, 0, 0 }, s });
            break;
        }
        case PartType::Post: {  // a quarter of the cell, up one edge
            CellShape s; s.full = false; s.planes = 2; s.n[0] = { 1, 0, 0 }; s.d[0] = 0.5f; s.n[1] = { 0, 0, 1 }; s.d[1] = 0.5f;
            out.push_back({ { 0, 0, 0 }, s });
            break;
        }
        case PartType::Slope: {
            CellShape s; s.full = false; s.planes = 1; s.n[0] = { 0, 1, 1 }; s.d[0] = 1.0f; // rises toward -z
            out.push_back({ { 0, 0, 0 }, s });
            break;
        }
        case PartType::Corner: {
            CellShape s; s.full = false; s.planes = 2; s.n[0] = { 0, 1, 1 }; s.d[0] = 1.0f; s.n[1] = { 1, 1, 0 }; s.d[1] = 1.0f;
            out.push_back({ { 0, 0, 0 }, s });
            break;
        }
        case PartType::Beam2: case PartType::Beam4: case PartType::Beam8: {
            const int n = t == PartType::Beam2 ? 2 : t == PartType::Beam4 ? 4 : 8;
            for (int z = 0; z < n; ++z) out.push_back({ { 0, 0, z }, full });
            break;
        }
        default: break;
    }
}

// ---------------------------------------------------------------------------------------------
// The 24 orientations
// ---------------------------------------------------------------------------------------------

namespace {
    struct OrientTable {
        int m[24][9];
        OrientTable() {
            const int perms[6][3] = { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 } };
            int n = 0;
            for (int p = 0; p < 6; ++p)
                for (int s = 0; s < 8; ++s) {
                    int r[9] = {};
                    for (int row = 0; row < 3; ++row) r[row * 3 + perms[p][row]] = (s >> row) & 1 ? -1 : 1;
                    const int det = r[0] * (r[4] * r[8] - r[5] * r[7]) - r[1] * (r[3] * r[8] - r[5] * r[6]) + r[2] * (r[3] * r[7] - r[4] * r[6]);
                    if (det != 1) continue;
                    memcpy(m[n++], r, sizeof(r));
                }
            // Identity first.
            for (int i = 0; i < 24; ++i)
                if (m[i][0] == 1 && m[i][4] == 1 && m[i][8] == 1) { int t[9]; memcpy(t, m[0], sizeof t); memcpy(m[0], m[i], sizeof t); memcpy(m[i], t, sizeof t); break; }
        }
    };
    const OrientTable& Table() { static OrientTable t; return t; }
    Int3 Mul(const int* m, Int3 v) { return { m[0] * v.x + m[1] * v.y + m[2] * v.z, m[3] * v.x + m[4] * v.y + m[5] * v.z, m[6] * v.x + m[7] * v.y + m[8] * v.z }; }
    Vector3 MulF(const int* m, Vector3 v) { return { m[0] * v.x + m[1] * v.y + m[2] * v.z, m[3] * v.x + m[4] * v.y + m[5] * v.z, m[6] * v.x + m[7] * v.y + m[8] * v.z }; }
}

const int* Orientation(int index) { return Table().m[std::clamp(index, 0, 23)]; }

int RotateOrientation(int index, int axis) {
    // 90 degrees about plot axis `axis`, applied after the current orientation.
    int R[9] = {};
    if (axis == 0) { R[0] = 1; R[5] = -1; R[7] = 1; }
    else if (axis == 1) { R[2] = 1; R[4] = 1; R[6] = -1; }
    else { R[1] = -1; R[3] = 1; R[8] = 1; }
    const int* M = Orientation(index);
    int P[9];
    for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) P[r * 3 + c] = R[r * 3] * M[c] + R[r * 3 + 1] * M[3 + c] + R[r * 3 + 2] * M[6 + c];
    for (int i = 0; i < 24; ++i) if (memcmp(Table().m[i], P, sizeof P) == 0) return i;
    return index;
}

// ---------------------------------------------------------------------------------------------
// Plot geometry
// ---------------------------------------------------------------------------------------------

Vector3 Plot::toPlot(Vector3 b) const {
    const DVec3 d = DVec3(b) - origin;
    const Vector3 v = { (float)d.x, (float)d.y, (float)d.z };
    return { Vector3DotProduct(v, east), Vector3DotProduct(v, up), Vector3DotProduct(v, north) };
}
Vector3 Plot::toBody(Vector3 p) const {
    const Vector3 v = Vector3Add(Vector3Add(Vector3Scale(east, p.x), Vector3Scale(up, p.y)), Vector3Scale(north, p.z));
    return { (float)(origin.x + v.x), (float)(origin.y + v.y), (float)(origin.z + v.z) };
}
Int3 Plot::cellAt(Vector3 p) const {
    return { (int)floorf(p.x / kCellSize), (int)floorf(p.y / kCellSize), (int)floorf(p.z / kCellSize) };
}

void Plot::Footprint(const Part& p, std::vector<std::pair<Int3, CellShape>>& out) {
    PartCells(p.type, out);
    const int* M = Orientation(p.orient);
    // Rotate about the origin cell's centre: c -> M c + off, off = M(0.5) - 0.5 (each 0 or -1).
    const Vector3 h = MulF(M, { 0.5f, 0.5f, 0.5f });
    const Int3 off = { (int)lroundf(h.x - 0.5f), (int)lroundf(h.y - 0.5f), (int)lroundf(h.z - 0.5f) };
    for (auto& c : out) {
        const Int3 r = Mul(M, c.first);
        c.first = { p.origin.x + r.x + off.x, p.origin.y + r.y + off.y, p.origin.z + r.z + off.z };
        CellShape& s = c.second;
        for (int i = 0; i < s.planes; ++i) {
            const Vector3 n2 = MulF(M, s.n[i]);
            s.d[i] = s.d[i] - 0.5f * (s.n[i].x + s.n[i].y + s.n[i].z) + 0.5f * (n2.x + n2.y + n2.z);
            s.n[i] = n2;
        }
    }
}

bool Plot::canPlace(const Part& p) const {
    std::vector<std::pair<Int3, CellShape>> fp;
    Footprint(p, fp);
    for (const auto& c : fp) if (cells.count(c.first)) return false;
    return true;
}

int Plot::place(const Part& p) {
    if (!canPlace(p)) return -1;
    std::vector<std::pair<Int3, CellShape>> fp;
    Footprint(p, fp);
    const int idx = (int)parts.size();
    parts.push_back(p);
    for (const auto& c : fp) { cells[c.first] = { idx, c.second }; markDirty(c.first); }
    ++livePartCount;
    ++version;
    return idx;
}

bool Plot::remove(int idx) {
    if (idx < 0 || idx >= (int)parts.size() || parts[idx].type == PartType::Count) return false;
    std::vector<std::pair<Int3, CellShape>> fp;
    Footprint(parts[idx], fp);
    for (const auto& c : fp) { cells.erase(c.first); markDirty(c.first); }
    parts[idx].type = PartType::Count;
    damage.erase(idx);
    --livePartCount;
    ++version;
    return true;
}

int Plot::partAt(Int3 c) const {
    auto it = cells.find(c);
    return it == cells.end() ? -1 : it->second.part;
}

static Int3 SectionOf(Int3 c) {
    auto fd = [](int a) { return a >= 0 ? a / kSectionCells : -((-a + kSectionCells - 1) / kSectionCells); };
    return { fd(c.x), fd(c.y), fd(c.z) };
}

void Plot::markDirty(Int3 cell) {
    // The cell's section, and neighbours' sections when it sits on a section face (face culling).
    for (int a = 0; a < 7; ++a) {
        Int3 c = cell;
        if (a == 1) c.x--; else if (a == 2) c.x++; else if (a == 3) c.y--; else if (a == 4) c.y++; else if (a == 5) c.z--; else if (a == 6) c.z++;
        sections[SectionOf(c)].dirty = true;
    }
}

// The shape of the part cell at c (full cube if none / for full parts).
static bool CellShapeAt(const Plot& plot, Int3 c, CellShape& shape) {
    auto it = plot.cells.find(c);
    if (it == plot.cells.end()) return false;
    shape = it->second.shape;
    return true;
}

namespace {
    // Clip a convex polygon by the half-space dot(n, p) <= d.
    void Clip(std::vector<Vector3>& poly, Vector3 n, float d) {
        if (poly.empty()) return;
        std::vector<Vector3> out;
        for (size_t i = 0; i < poly.size(); ++i) {
            const Vector3 a = poly[i], b = poly[(i + 1) % poly.size()];
            const float da = Vector3DotProduct(n, a) - d, db = Vector3DotProduct(n, b) - d;
            if (da <= 0.0f) out.push_back(a);
            if ((da < 0.0f && db > 0.0f) || (da > 0.0f && db < 0.0f)) out.push_back(Vector3Lerp(a, b, da / (da - db)));
        }
        poly.swap(out);
    }
    const Vector3 kFaceN[6] = { { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 } };
    const Int3 kFaceD[6] = { { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 } };
    std::vector<Vector3> CubeFace(int f) {
        // Corners of face f of the unit cube, counter-clockwise seen from outside.
        switch (f) {
            case 0: return { { 0, 0, 0 }, { 0, 0, 1 }, { 0, 1, 1 }, { 0, 1, 0 } };
            case 1: return { { 1, 0, 0 }, { 1, 1, 0 }, { 1, 1, 1 }, { 1, 0, 1 } };
            case 2: return { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 0, 1 }, { 0, 0, 1 } };
            case 3: return { { 0, 1, 0 }, { 0, 1, 1 }, { 1, 1, 1 }, { 1, 1, 0 } };
            case 4: return { { 0, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 }, { 1, 0, 0 } };
            default: return { { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
        }
    }
    float PolyArea(const std::vector<Vector3>& p) {
        Vector3 s = { 0, 0, 0 };
        for (size_t i = 0; i < p.size(); ++i) s = Vector3Add(s, Vector3CrossProduct(p[i], p[(i + 1) % p.size()]));
        return 0.5f * Vector3Length(s);
    }
    void EmitPoly(const std::vector<Vector3>& poly, Vector3 n, Vector3 base, Color col,
                  std::vector<float>& pos, std::vector<float>& nrm, std::vector<unsigned char>& colr) {
        if (poly.size() < 3) return;
        // Order the fan so its winding faces n in render space. The plot frame (east, up, north =
        // up x east) is left-handed, which mirrors winding, so the fan winds clockwise about n here.
        Vector3 cr = Vector3CrossProduct(Vector3Subtract(poly[1], poly[0]), Vector3Subtract(poly[2], poly[0]));
        const bool flip = Vector3DotProduct(cr, n) > 0.0f;
        for (size_t i = 1; i + 1 < poly.size(); ++i) {
            const Vector3 tri[3] = { poly[0], flip ? poly[i + 1] : poly[i], flip ? poly[i] : poly[i + 1] };
            for (const Vector3& v : tri) {
                const Vector3 p = Vector3Scale(Vector3Add(base, v), kCellSize);
                pos.insert(pos.end(), { p.x, p.y, p.z });
                nrm.insert(nrm.end(), { n.x, n.y, n.z });
                colr.insert(colr.end(), { col.r, col.g, col.b, col.a });
            }
        }
    }
}

void Plot::rebuildDirty(const PlanetField* field) {
    std::unordered_map<Int3, std::vector<Int3>> byDirty;
    bool any = false;
    for (auto& kv : sections) if (kv.second.dirty) { byDirty[kv.first]; any = true; }
    if (!any) return;
    for (const auto& kv : cells) {
        auto it = byDirty.find(SectionOf(kv.first));
        if (it != byDirty.end()) it->second.push_back(kv.first);
    }
    for (auto& kv : byDirty) buildSection(kv.first, sections[kv.first], kv.second, field);
}

void Plot::buildSection(Int3 key, Section& s, const std::vector<Int3>& sectionCells, const PlanetField* field) {
    if (s.mesh.vaoId) UnloadMesh(s.mesh);
    s.mesh = Mesh{};
    s.dirty = false;
    std::vector<float> pos, nrm, sky; // per vertex u = how buried (1 - sky exposure; parts.vs / parts.fs, #216)
    std::vector<unsigned char> col;
    const Int3 lo = { key.x * kSectionCells, key.y * kSectionCells, key.z * kSectionCells };
    for (const Int3& c : sectionCells) {
        // Sky exposure of the cell (1 open, 0 buried ~10 m+): the terrain's own rule (materialSky's depth below
        // the original surface), from the cell centre. Parts underground took the full sky ambient (#216).
        float cellSky = 1.0f;
        if (field) {
            const Vector3 centre = toBody({ (c.x + 0.5f) * kCellSize, (c.y + 0.5f) * kCellSize, (c.z + 0.5f) * kCellSize });
            const float len = Vector3Length(centre);
            if (len > 1.0f) field->materialSky(centre, 1.0f, 1.0f, Vector3Scale(centre, 1.0f / len), cellSky);
        }
        const size_t firstVertex = pos.size() / 3;
        const CellRef& ref = cells.at(c);
        const CellShape& shape = ref.shape;
        const Part& pp = parts[ref.part];
        Color color = IsShape(pp.type) ? GetBuildMaterial(pp.material).color : PartColor(pp.type);
        if (pp.type == PartType::Block && IsShape(pp.type)) color = { (unsigned char)(color.r * 0.92f), (unsigned char)(color.g * 0.92f), (unsigned char)(color.b * 0.92f), 255 };
        color.a = PartPatternOf(pp); // the surface pattern id (parts.fs)
        const Vector3 base = { (float)(c.x - lo.x), (float)(c.y - lo.y), (float)(c.z - lo.z) };
        // Cube faces, clipped by the shape's planes; full faces against a full neighbour are hidden.
        for (int f = 0; f < 6; ++f) {
            std::vector<Vector3> poly = CubeFace(f);
            for (int i = 0; i < shape.planes; ++i) Clip(poly, shape.n[i], shape.d[i]);
            if (poly.size() < 3) continue;
            if (PolyArea(poly) > 0.999f) {
                CellShape ns;
                const Int3 nc = { c.x + kFaceD[f].x, c.y + kFaceD[f].y, c.z + kFaceD[f].z };
                if (CellShapeAt(*this, nc, ns) && ns.full) continue;
            }
            const float shade = f == 3 ? 1.0f : f == 2 ? 0.7f : 0.85f + 0.05f * (f & 1);
            Color fc = { (unsigned char)(color.r * shade), (unsigned char)(color.g * shade), (unsigned char)(color.b * shade), color.a };
            EmitPoly(poly, kFaceN[f], base, fc, pos, nrm, col);
        }
        // Caps on the cutting planes.
        for (int i = 0; i < shape.planes; ++i) {
            const Vector3 n = Vector3Normalize(shape.n[i]);
            const float d = shape.d[i] / Vector3Length(shape.n[i]);
            const Vector3 t = Vector3Normalize(Vector3CrossProduct(n, fabsf(n.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
            const Vector3 b = Vector3CrossProduct(n, t);
            const Vector3 o = Vector3Scale(n, d);
            std::vector<Vector3> poly = { Vector3Add(o, Vector3Add(Vector3Scale(t, -3), Vector3Scale(b, -3))), Vector3Add(o, Vector3Add(Vector3Scale(t, 3), Vector3Scale(b, -3))),
                                          Vector3Add(o, Vector3Add(Vector3Scale(t, 3), Vector3Scale(b, 3))), Vector3Add(o, Vector3Add(Vector3Scale(t, -3), Vector3Scale(b, 3))) };
            for (int f = 0; f < 6; ++f) {
                const Vector3 fn = kFaceN[f];
                Clip(poly, fn, f & 1 ? 1.0f : 0.0f); // inside the unit cube: x <= 1 / -x <= 0 ...
            }
            for (int j = 0; j < shape.planes; ++j) if (j != i) Clip(poly, shape.n[j], shape.d[j]);
            EmitPoly(poly, n, base, color, pos, nrm, col);
        }
        for (size_t v = firstVertex; v < pos.size() / 3; ++v) { sky.push_back(1.0f - cellSky); sky.push_back(0.0f); } // as "buried": meshes without it read 0 = open sky
    }
    if (pos.empty()) return;
    s.mesh.vertexCount = (int)(pos.size() / 3);
    s.mesh.triangleCount = s.mesh.vertexCount / 3;
    s.mesh.vertices = (float*)MemAlloc((unsigned int)(pos.size() * sizeof(float)));
    s.mesh.normals = (float*)MemAlloc((unsigned int)(nrm.size() * sizeof(float)));
    s.mesh.colors = (unsigned char*)MemAlloc((unsigned int)col.size());
    memcpy(s.mesh.vertices, pos.data(), pos.size() * sizeof(float));
    memcpy(s.mesh.normals, nrm.data(), nrm.size() * sizeof(float));
    memcpy(s.mesh.colors, col.data(), col.size());
    s.mesh.texcoords = (float*)MemAlloc((unsigned int)(sky.size() * sizeof(float)));
    memcpy(s.mesh.texcoords, sky.data(), sky.size() * sizeof(float));
    UploadMesh(&s.mesh, false);
}

Plot::~Plot() {
    for (auto& kv : sections) if (kv.second.mesh.vaoId) UnloadMesh(kv.second.mesh);
}

Matrix Plot::cellModel(const Planetoid& body, const DVec3& renderOrigin) const {
    const Matrix bodyM = body.localModel(origin, renderOrigin);
    Matrix r = MatrixIdentity();
    r.m0 = east.x; r.m4 = up.x; r.m8 = north.x;
    r.m1 = east.y; r.m5 = up.y; r.m9 = north.y;
    r.m2 = east.z; r.m6 = up.z; r.m10 = north.z;
    return MatrixMultiply(r, bodyM);
}

static Matrix PlotSectionModel(const Plot& p, const Planetoid& body, const RenderContext& ctx, Int3 key) {
    const float s = kSectionCells * kCellSize;
    const Vector3 off = Vector3Add(Vector3Add(Vector3Scale(p.east, key.x * s), Vector3Scale(p.up, key.y * s)), Vector3Scale(p.north, key.z * s));
    const Matrix bodyM = body.localModel(DVec3(p.origin.x + off.x, p.origin.y + off.y, p.origin.z + off.z), ctx.renderOrigin);
    Matrix r = MatrixIdentity();
    r.m0 = p.east.x; r.m4 = p.up.x; r.m8 = p.north.x;
    r.m1 = p.east.y; r.m5 = p.up.y; r.m9 = p.north.y;
    r.m2 = p.east.z; r.m6 = p.up.z; r.m10 = p.north.z;
    return MatrixMultiply(r, bodyM);
}

void Plot::draw(const Planetoid& body, const RenderContext& ctx, Shader shader) {
    static Material material = LoadMaterialDefault();
    material.shader = shader;
    rebuildDirty(dynamic_cast<const PlanetField*>(body.field.get()));
    const float half = kSectionCells * kCellSize * 0.5f;
    for (auto& kv : sections) {
        if (!kv.second.mesh.vaoId) continue;
        const Matrix model = PlotSectionModel(*this, body, ctx, kv.first);
        // Cull sections outside the view (a big factory has hundreds of them).
        if (ctx.frustum && frustumCullingEnabled && !ctx.frustum->sphereVisible(Vector3Transform({ half, half, half }, model), half * 1.7321f)) continue;
        DrawMesh(kv.second.mesh, material, model);
        if (ctx.drawCalls) ++*ctx.drawCalls;
    }
}

void Plot::drawDepth(const Planetoid& body, const RenderContext& ctx, Material& depth) {
    rebuildDirty(dynamic_cast<const PlanetField*>(body.field.get()));
    for (auto& kv : sections) {
        if (kv.second.mesh.vaoId) DrawMesh(kv.second.mesh, depth, PlotSectionModel(*this, body, ctx, kv.first));
    }
}

// ---------------------------------------------------------------------------------------------
// Collision and picking
// ---------------------------------------------------------------------------------------------

// Closest point of a cell's solid (box [lo, lo+1] in cell units, cut by its planes) to q.
static Vector3 ClosestInCell(Vector3 q, Int3 c, const CellShape& s) {
    const Vector3 lo = { (float)c.x, (float)c.y, (float)c.z };
    auto clampBox = [&](Vector3 v) { return Vector3{ std::clamp(v.x, lo.x, lo.x + 1), std::clamp(v.y, lo.y, lo.y + 1), std::clamp(v.z, lo.z, lo.z + 1) }; };
    Vector3 p = clampBox(q);
    for (int it = 0; it < 3; ++it)
        for (int i = 0; i < s.planes; ++i) {
            const Vector3 n = s.n[i];
            const float dd = s.d[i] + Vector3DotProduct(n, lo); // plane in plot cell units
            const float over = Vector3DotProduct(n, p) - dd;
            if (over > 0.0f) p = clampBox(Vector3Subtract(p, Vector3Scale(n, over / Vector3DotProduct(n, n))));
        }
    return p;
}

bool Plot::collideSphere(Vector3& bodyCentre, float r, Vector3& normalOut) const {
    if (cells.empty()) return false;
    Vector3 p = Vector3Scale(toPlot(bodyCentre), 1.0f / kCellSize); // cell units
    const float rc = r / kCellSize;
    bool hit = false;
    Vector3 nSum = { 0, 0, 0 };
    // Deepest overlap first, then measure again: pushing out of cells in grid order let a
    // neighbour's edge shove the sphere sideways off a flat floor.
    for (int pass = 0; pass < 4; ++pass) {
        float deepest = 0.0f;
        Vector3 bestN = { 0, 0, 0 };
        for (int z = (int)floorf(p.z - rc); z <= (int)floorf(p.z + rc); ++z)
            for (int y = (int)floorf(p.y - rc); y <= (int)floorf(p.y + rc); ++y)
                for (int x = (int)floorf(p.x - rc); x <= (int)floorf(p.x + rc); ++x) {
                    auto it = cells.find({ x, y, z });
                    if (it == cells.end() || parts[it->second.part].type == PartType::Gate) continue; // people pass a gate
                    const CellShape& s = it->second.shape;
                    const Vector3 q = ClosestInCell(p, { x, y, z }, s);
                    Vector3 d = Vector3Subtract(p, q);
                    float len = Vector3Length(d);
                    if (len >= rc) continue;
                    if (len < 1e-5f) { d = { 0, 1, 0 }; len = 0.0f; } // centre inside: push up
                    const float depth = rc - len;
                    if (depth > deepest) { deepest = depth; bestN = len > 0.0f ? Vector3Scale(d, 1.0f / len) : d; }
                }
        if (deepest <= 0.0f) break;
        p = Vector3Add(p, Vector3Scale(bestN, deepest));
        nSum = Vector3Add(nSum, bestN);
        hit = true;
    }
    if (!hit) return false;
    bodyCentre = toBody(Vector3Scale(p, kCellSize));
    const Vector3 n = Vector3Normalize(nSum);
    normalOut = Vector3Normalize(Vector3Add(Vector3Add(Vector3Scale(east, n.x), Vector3Scale(up, n.y)), Vector3Scale(north, n.z)));
    return true;
}

bool Plot::raycast(Vector3 o, Vector3 d, float maxDist, Int3& cell, Int3& faceNormal, float& dist) const {
    if (cells.empty()) return false;
    // DDA through the cell grid in plot space.
    Vector3 p = Vector3Scale(toPlot(o), 1.0f / kCellSize);
    const Vector3 bd = { Vector3DotProduct(d, east), Vector3DotProduct(d, up), Vector3DotProduct(d, north) };
    Int3 c = { (int)floorf(p.x), (int)floorf(p.y), (int)floorf(p.z) };
    const int sx = bd.x > 0 ? 1 : -1, sy = bd.y > 0 ? 1 : -1, sz = bd.z > 0 ? 1 : -1;
    auto inv = [](float v) { return fabsf(v) > 1e-9f ? 1.0f / fabsf(v) : 1e30f; };
    const float tdx = inv(bd.x), tdy = inv(bd.y), tdz = inv(bd.z);
    float tx = (bd.x > 0 ? (c.x + 1 - p.x) : (p.x - c.x)) * tdx;
    float ty = (bd.y > 0 ? (c.y + 1 - p.y) : (p.y - c.y)) * tdy;
    float tz = (bd.z > 0 ? (c.z + 1 - p.z) : (p.z - c.z)) * tdz;
    const float maxT = maxDist / kCellSize;
    Int3 face = { 0, 0, 0 };
    float t = 0.0f;
    while (t <= maxT) {
        if (cells.count(c)) { cell = c; faceNormal = face; dist = t * kCellSize; return true; }
        if (tx < ty && tx < tz) { c.x += sx; t = tx; tx += tdx; face = { -sx, 0, 0 }; }
        else if (ty < tz) { c.y += sy; t = ty; ty += tdy; face = { 0, -sy, 0 }; }
        else { c.z += sz; t = tz; tz += tdz; face = { 0, 0, -sz }; }
    }
    return false;
}

// ---------------------------------------------------------------------------------------------
// PlotSet
// ---------------------------------------------------------------------------------------------

Plot* PlotSet::plotFor(Vector3 b, float reach, bool create) {
    Plot* best = nullptr;
    float bestD = reach;
    for (auto& pl : plots) {
        if (pl->cells.empty()) continue;
        const Int3 c = pl->cellAt(pl->toPlot(b));
        // Distance (cells) to the nearest occupied cell, checked coarsely over a small box.
        const int R = (int)ceilf(reach / kCellSize);
        float d = 1e30f;
        if (pl->cells.size() < 4096) {
            for (const auto& kv : pl->cells)
                d = std::min(d, (float)std::max({ abs(kv.first.x - c.x), abs(kv.first.y - c.y), abs(kv.first.z - c.z) }));
        } else {
            for (int z = -R; z <= R && d > 0; ++z) for (int y = -R; y <= R; ++y) for (int x = -R; x <= R; ++x)
                if (pl->cells.count({ c.x + x, c.y + y, c.z + z })) d = std::min(d, (float)std::max({ abs(x), abs(y), abs(z) }));
        }
        d *= kCellSize;
        if (d < bestD) { bestD = d; best = pl.get(); }
    }
    if (best || !create) return best;
    auto pl = std::make_unique<Plot>();
    pl->id = nextId++;
    pl->origin = DVec3(b);
    pl->up = Vector3Normalize(b);
    pl->east = Vector3CrossProduct({ 0, 1, 0 }, pl->up);
    if (Vector3Length(pl->east) < 0.1f) pl->east = { 1, 0, 0 };
    pl->east = Vector3Normalize(pl->east);
    pl->north = Vector3CrossProduct(pl->up, pl->east);
    plots.push_back(std::move(pl));
    return plots.back().get();
}

Plot* PlotSet::byId(uint32_t id) {
    for (auto& p : plots) if (p->id == id) return p.get();
    return nullptr;
}

int PlotSet::partCount() const {
    int n = 0;
    for (const auto& p : plots) n += p->partCount();
    return n;
}

void PlotSet::save(FILE* f) const {
    int live = 0;
    for (const auto& p : plots) if (p->partCount() > 0) ++live;
    // plots4: every part slot, removed ones as "-1" placeholders, so part indices (what the factory
    // saves machine contents by) survive a reload; parts carry their material and plots their rocket
    // fuel. (plots3: live parts only; plots2: no fuel; plots: no materials.)
    // plots5: plots4 plus a "wear" line per plot (rocket part damage by part index).
    fprintf(f, "plots5 %d\n", live);
    for (const auto& p : plots) {
        if (p->partCount() == 0) continue;
        fprintf(f, "plot %u %.4f %.4f %.4f %d %.1f\n", p->id, p->origin.x, p->origin.y, p->origin.z, (int)p->parts.size(), p->fuel);
        for (const Part& part : p->parts) {
            if (part.type == PartType::Count) fprintf(f, "-1\n");
            else fprintf(f, "%d %d %d %d %d %d\n", (int)part.type, part.origin.x, part.origin.y, part.origin.z, part.orient, part.material);
        }
        int worn = 0;
        for (const auto& kv : p->damage) worn += kv.second > 0.0f ? 1 : 0;
        fprintf(f, "wear %d %.2f", worn, p->air); // and the air its life support holds
        for (const auto& kv : p->damage) if (kv.second > 0.0f) fprintf(f, " %d %.3f", kv.first, kv.second);
        fprintf(f, "\n");
    }
}

bool PlotSet::load(FILE* f) {
    plots.clear();
    int n = 0;
    char tag[16] = {};
    if (fscanf(f, " %15s %d", tag, &n) != 2 || (strcmp(tag, "plots") != 0 && strcmp(tag, "plots2") != 0 && strcmp(tag, "plots3") != 0 && strcmp(tag, "plots4") != 0 && strcmp(tag, "plots5") != 0)) return false;
    const bool wear = strcmp(tag, "plots5") == 0;
    const bool materials = strcmp(tag, "plots") != 0, hasFuel = strcmp(tag, "plots3") == 0 || strcmp(tag, "plots4") == 0 || wear;
    const bool slots = strcmp(tag, "plots4") == 0 || wear; // removed parts kept as placeholders
    const Part dead{ PartType::Count, { 0, 0, 0 }, 0, 0 };
    for (int i = 0; i < n; ++i) {
        unsigned int id = 0;
        double ox, oy, oz;
        int count = 0;
        if (fscanf(f, " plot %u %lf %lf %lf %d", &id, &ox, &oy, &oz, &count) != 5) return false;
        float fuel = 0.0f;
        if (hasFuel && fscanf(f, " %f", &fuel) != 1) return false;
        auto p = std::make_unique<Plot>();
        p->id = id;
        p->origin = DVec3(ox, oy, oz);
        p->fuel = fuel;
        const Vector3 ob = { (float)ox, (float)oy, (float)oz };
        p->up = Vector3Normalize(ob);
        p->east = Vector3CrossProduct({ 0, 1, 0 }, p->up);
        if (Vector3Length(p->east) < 0.1f) p->east = { 1, 0, 0 };
        p->east = Vector3Normalize(p->east);
        p->north = Vector3CrossProduct(p->up, p->east);
        for (int k = 0; k < count; ++k) {
            int t, x, y, z, o, m = 0;
            if (fscanf(f, " %d", &t) != 1) return false;
            if (slots && t < 0) { p->parts.push_back(dead); continue; } // a removed part's slot
            if (fscanf(f, " %d %d %d %d", &x, &y, &z, &o) != 4) return false;
            if (materials && fscanf(f, " %d", &m) != 1) return false;
            // Old saves: beams were iron.
            if (!materials && (t == (int)PartType::Beam2 || t == (int)PartType::Beam4 || t == (int)PartType::Beam8)) m = BuildMaterialForItem("iron_plate");
            // A part that no longer fits still keeps its slot, so the indices after it don't shift.
            if (p->place({ (PartType)t, { x, y, z }, (uint8_t)o, (uint8_t)std::clamp(m, 0, BuildMaterialCount() - 1) }) < 0) p->parts.push_back(dead);
        }
        if (wear) {
            int worn = 0;
            if (fscanf(f, " wear %d %f", &worn, &p->air) != 2) return false;
            for (int k = 0; k < worn; ++k) { int part = 0; float d = 0.0f; if (fscanf(f, " %d %f", &part, &d) != 2) return false; p->damage[part] = d; }
        }
        nextId = std::max(nextId, id + 1);
        plots.push_back(std::move(p));
    }
    return true;
}

Plot* PlotSet::adopt(std::unique_ptr<Plot> plot) {
    plot->id = nextId++;
    plot->version++;
    plots.push_back(std::move(plot));
    return plots.back().get();
}

std::unique_ptr<Plot> PlotSet::release(Plot* plot) {
    for (auto it = plots.begin(); it != plots.end(); ++it)
        if (it->get() == plot) { std::unique_ptr<Plot> p = std::move(*it); plots.erase(it); return p; }
    return nullptr;
}
