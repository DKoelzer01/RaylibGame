#include "buildtool.h"
#include "survival.h"
#include "audio.h"
#include "gameui.h"
#include "inventory.h"
#include "player.h"
#include "terraintool.h"
#include "world.h"
#include "worldtypes.h"
#include "factory.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>

BuildToolState buildTool;

namespace {
    // The current aim: a placement (plot may be null = a new plot at newOrigin) and/or a part to remove.
    struct Aim {
        Planetoid* body = nullptr;
        Plot* plot = nullptr;           // existing plot to place into (null: start a new one)
        Vector3 newOrigin{};            // body-local, for a new plot
        Part part{};
        bool canPlace = false;
        Plot* removePlot = nullptr;
        int removePart = -1;
    } aim;

    // A throwaway plot frame for previewing a new plot's first part.
    Plot preview;

    Vector3 PlotDir(const Plot& p, Int3 v) {
        return Vector3Add(Vector3Add(Vector3Scale(p.east, (float)v.x), Vector3Scale(p.up, (float)v.y)), Vector3Scale(p.north, (float)v.z));
    }
}

static bool IsMachineItem(const std::string& id) {
    for (int t = 0; t < (int)PartType::Count; ++t) if (id == PartCostItem((PartType)t) && IsMachine((PartType)t)) return true;
    return false;
}

// Foundation: clear terrain inside the part's cells and fill under its bottom cells.
void BuildTool_LayFoundation(Planetoid& body, Plot& plot, const Part& part) {
    std::vector<std::pair<Int3, CellShape>> fp;
    Plot::Footprint(part, fp);
    std::vector<TerrainEdit> edits;
    for (const auto& c : fp) {
        const Vector3 centre = plot.toBody({ (c.first.x + 0.5f) * kCellSize, (c.first.y + 0.5f) * kCellSize, (c.first.z + 0.5f) * kCellSize });
        // Clear only inside the cell (a 0.55 m brush dug ~0.3 m past the part into whatever it was built against).
        if (body.densityAt(centre) > ISO_LEVEL - 0.15f) edits.push_back({ centre, 0.3f, -1.6f });
        const bool bottom = plot.partAt({ c.first.x, c.first.y - 1, c.first.z }) < 0;
        if (bottom) {
            // Only a small gap to real ground is filled (air just under the part, ground within ~1 m);
            // parts higher up are carried by the plot, not propped on mounds of terrain.
            const Vector3 below = plot.toBody({ (c.first.x + 0.5f) * kCellSize, (c.first.y - 0.35f) * kCellSize, (c.first.z + 0.5f) * kCellSize });
            const Vector3 ground = plot.toBody({ (c.first.x + 0.5f) * kCellSize, c.first.y * kCellSize - 1.0f, (c.first.z + 0.5f) * kCellSize });
            if (body.densityAt(below) < ISO_LEVEL && body.densityAt(ground) >= ISO_LEVEL) edits.push_back({ below, 0.45f, 0.8f });
        }
    }
    if (!edits.empty()) ApplyTerrainEdits(body, edits.data(), edits.size());
}

bool BuildTool_Place(Planetoid& body, Plot& plot, const Part& part, bool charge) {
    if (!plot.canPlace(part)) return false;
    const PartType t = part.type;
    const char* item = PartItem(part);
    if (charge && !terrainTool.creative && !inventory.remove(item, PartCostCount(t))) {
        GameUI_Toast(TextFormat("Need %d %s for a %s", PartCostCount(t), item, PartName(t)));
        Audio_Play(Sfx::Error, 0.5f);
        return false;
    }
    plot.place(part);
    if (IsMachine(part.type) || IsRocketPart(part.type)) QuestBook_NotePlaced(PartCostItem(part.type));
    BuildTool_LayFoundation(body, plot, part);
    return true;
}

void BuildTool_Update(Scene& world, bool inputEnabled, bool usable) {
    if (!usable) { buildTool.active = false; return; }
    if (inputEnabled && IsKeyPressed(KEY_G)) {
        buildTool.active = !buildTool.active;
        GameUI_Toast(buildTool.active ? "Build mode: wheel picks the shape, R/V rotate, RMB place, LMB remove (G to leave; the hotbar picks the material or machine)" : "Build mode off");
    }
    aim = Aim{};
    if (!buildTool.active) return;
    if (inputEnabled) {
        const float wheel = GetMouseWheelMove();
        // The wheel picks the shape (any shape from the selected material); machines come from the hotbar.
        if (IsShape((PartType)buildTool.part) && wheel > 0) buildTool.part = (buildTool.part + 1) % kShapeCount;
        if (IsShape((PartType)buildTool.part) && wheel < 0) buildTool.part = (buildTool.part + kShapeCount - 1) % kShapeCount;
        if (IsKeyPressed(KEY_R)) buildTool.orient = RotateOrientation(buildTool.orient, 1);
        if (IsKeyPressed(KEY_V)) buildTool.orient = RotateOrientation(buildTool.orient, 0);
    }
    Planetoid* body = player.nearest;
    if (!body) return;
    const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
    const Vector3 lo = body->toLocal(world.camera.position), ld = body->dirToLocal(dir);
    aim.body = body;
    aim.part = { (PartType)buildTool.part, { 0, 0, 0 }, (uint8_t)buildTool.orient, (uint8_t)buildTool.material };

    // Nearest part face along the ray, and the terrain.
    float bestPart = buildTool.reach;
    Plot* hitPlot = nullptr;
    Int3 hitCell{}, hitFace{};
    if (body->plots)
        for (auto& pl : body->plots->plots) {
            Int3 c, f; float d;
            if (pl->raycast(lo, ld, bestPart, c, f, d)) { bestPart = d; hitPlot = pl.get(); hitCell = c; hitFace = f; }
        }
    TerrainHit th;
    const bool terrain = RaycastTerrain(world, Ray{ world.camera.position, dir }, buildTool.reach, th) && th.planetoid == body;
    if (hitPlot && (!terrain || bestPart <= th.distance)) {
        aim.plot = hitPlot;
        aim.part.origin = { hitCell.x + hitFace.x, hitCell.y + hitFace.y, hitCell.z + hitFace.z };
        aim.removePlot = hitPlot;
        aim.removePart = hitPlot->partAt(hitCell);
    } else if (terrain) {
        const Vector3 ground = th.local;
        Plot* near = body->plots ? body->plots->plotFor(ground, 6.0f, false) : nullptr;
        if (near) {
            aim.plot = near;
            // The cell just above the ground at the aim.
            aim.part.origin = near->cellAt(near->toPlot(Vector3Add(ground, Vector3Scale(near->up, 0.1f))));
        } else {
            aim.newOrigin = ground;
            aim.part.origin = { 0, 0, 0 };
        }
    } else {
        return;
    }
    if (aim.plot) aim.canPlace = aim.plot->canPlace(aim.part);
    else aim.canPlace = true;

    if (!inputEnabled) return;
    if (IsMouseButtonPressed(MOUSE_BUTTON_RIGHT) && aim.canPlace) {
        Plot* target = aim.plot;
        if (!target) {
            if (!body->plots) body->plots = std::make_unique<PlotSet>();
            target = body->plots->plotFor(aim.newOrigin, 0.0f, true);
        }
        if (BuildTool_Place(*body, *target, aim.part, true)) Audio_Play(Sfx::Build, 0.5f, 1.2f);
    } else if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && aim.removePlot && aim.removePart >= 0) {
        const Part removed = aim.removePlot->parts[aim.removePart];
        if (aim.removePlot->remove(aim.removePart)) {
            if (!terrainTool.creative) inventory.add(PartItem(removed), PartCostCount(removed.type));
            Audio_Play(Sfx::Dig, 0.5f, 1.3f);
        }
    }
}

void BuildTool_Draw3D(Scene& world) {
    (void)world;
    if (!buildTool.active || !aim.body) return;
    Planetoid& body = *aim.body;
    const Plot* frame = aim.plot;
    if (!frame) {
        preview.origin = DVec3(aim.newOrigin);
        preview.up = Vector3Normalize(aim.newOrigin);
        preview.east = Vector3CrossProduct({ 0, 1, 0 }, preview.up);
        if (Vector3Length(preview.east) < 0.1f) preview.east = { 1, 0, 0 };
        preview.east = Vector3Normalize(preview.east);
        preview.north = Vector3CrossProduct(preview.up, preview.east);
        frame = &preview;
    }
    auto box = [&](const Plot& p, Int3 c, Color col) {
        Vector3 v[8];
        for (int i = 0; i < 8; ++i) {
            const Vector3 l = p.toBody({ (c.x + (i & 1)) * kCellSize, (c.y + ((i >> 1) & 1)) * kCellSize, (c.z + ((i >> 2) & 1)) * kCellSize });
            v[i] = body.toWorld(l);
        }
        const int e[12][2] = { { 0, 1 }, { 2, 3 }, { 4, 5 }, { 6, 7 }, { 0, 2 }, { 1, 3 }, { 4, 6 }, { 5, 7 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
        for (auto& ed : e) DrawLine3D(v[ed[0]], v[ed[1]], col);
    };
    std::vector<std::pair<Int3, CellShape>> fp;
    Plot::Footprint(aim.part, fp);
    for (const auto& c : fp) box(*frame, c.first, aim.canPlace ? GREEN : RED);
    if (aim.removePlot && aim.removePart >= 0) {
        Plot::Footprint(aim.removePlot->parts[aim.removePart], fp);
        for (const auto& c : fp) box(*aim.removePlot, c.first, Fade(ORANGE, 0.7f));
    }
}

void BuildTool_DrawHUD() {
    if (!buildTool.active) return;
    const PartType t = (PartType)buildTool.part;
    const int cx = GetScreenWidth() / 2, cy = GetScreenHeight() / 2;
    DrawLine(cx - 8, cy, cx + 8, cy, RAYWHITE);
    DrawLine(cx, cy - 8, cx, cy + 8, RAYWHITE);
    const Part p{ t, { 0, 0, 0 }, 0, (uint8_t)buildTool.material };
    const char* item = PartItem(p);
    const char* label = IsShape(t) ? TextFormat("%s %s", GetBuildMaterial(buildTool.material).name, PartName(t)) : PartName(t);
    const char* line = TextFormat("BUILD  %s  (%d %s, have %d)", label, PartCostCount(t), item, inventory.count(item));
    DrawText(line, cx - MeasureText(line, 22) / 2, GetScreenHeight() - 150, 22, SKYBLUE);
    DrawText(IsShape(t) ? "wheel shape   R / V rotate   RMB place   LMB remove   G leave   1-0 hotbar"
                        : "R / V rotate   RMB place   LMB remove   G leave   1-0 hotbar",
             10, GetScreenHeight() - 30, 20, Fade(RAYWHITE, 0.7f));
}

// ---------------------------------------------------------------------------------------------
// Machines: look at one to see its state; I puts matching items in, O takes the outputs.
// ---------------------------------------------------------------------------------------------

namespace {
    struct MachineAim { uint64_t body = 0; uint32_t plot = 0; int part = -1; PartType type = PartType::Count; } machineAim;
}

void Machine_Update(Scene& world, bool inputEnabled) {
    machineAim = MachineAim{};
    Planetoid* body = player.nearest;
    if (!body || !body->plots || player.flyMode) return;
    const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
    const Vector3 lo = body->toLocal(world.camera.position), ld = body->dirToLocal(dir);
    float best = 6.0f;
    for (auto& pl : body->plots->plots) {
        Int3 c, f; float d;
        if (!pl->raycast(lo, ld, best, c, f, d)) continue;
        const int part = pl->partAt(c);
        if (part < 0 || !IsMachine(pl->parts[part].type)) continue;
        best = d;
        machineAim = { body->id, pl->id, part, pl->parts[part].type };
    }
    if (machineAim.part < 0 || !inputEnabled) return;
    if (IsKeyPressed(KEY_I)) {
        // Everything in the inventory the machine might take: fuel for generators, recipe inputs,
        // anything for a chest (the simulation hands back what does not fit).
        // A chest takes the hand's raw materials only (not tools and machines).
        const auto all = inventory.all(); // copy: we remove while walking it
        for (const auto& kv : all) {
            if (kv.second <= 0 || !Factory::Accepts(machineAim.type, kv.first)) continue;
            if (machineAim.type == PartType::Chest && (kv.first == "telescope" || IsMachineItem(kv.first))) continue;
            if (inventory.remove(kv.first, kv.second)) Factory::Insert(machineAim.body, machineAim.plot, machineAim.part, kv.first, kv.second);
        }
        Audio_Play(Sfx::Pickup, 0.4f, 0.8f);
    }
    if (IsKeyPressed(KEY_O)) { Factory::Extract(machineAim.body, machineAim.plot, machineAim.part); Audio_Play(Sfx::Pickup, 0.4f, 1.2f); }
}

void Machine_DrawHUD() {
    if (machineAim.part < 0) return;
    const MachineView v = Factory::View(machineAim.body, machineAim.plot, machineAim.part);
    const int x = GetScreenWidth() / 2 + 20, y = GetScreenHeight() / 2 + 40;
    DrawRectangle(x - 8, y - 6, 420, 112, Fade(BLACK, 0.55f));
    DrawText(TextFormat("%s: %s", PartName(machineAim.type), v.valid ? MachineStatusName(v.status) : "..."), x, y, 20, SKYBLUE);
    if (v.duration > 0) {
        DrawRectangle(x, y + 26, 300, 8, Fade(RAYWHITE, 0.2f));
        DrawRectangle(x, y + 26, 300 * v.progress / v.duration, 8, ORANGE);
    }
    DrawText(TextFormat("in: %s", v.input.c_str()), x, y + 40, 16, RAYWHITE);
    DrawText(TextFormat("out: %s", v.output.c_str()), x, y + 60, 16, RAYWHITE);
    if (v.energyMax > 0) DrawText(TextFormat("charge %d / %d EU", v.energy, v.energyMax), x, y + 80, 16, YELLOW);
    else DrawText("I put items in   O take the output", x, y + 80, 16, Fade(RAYWHITE, 0.7f));
}
