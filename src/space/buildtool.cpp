#include "buildtool.h"
#include "audio.h"
#include "gameui.h"
#include "inventory.h"
#include "player.h"
#include "terraintool.h"
#include "world.h"
#include "worldtypes.h"

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

// Foundation: clear terrain inside the part's cells and fill under its bottom cells.
static void LayFoundation(Planetoid& body, Plot& plot, const Part& part) {
    std::vector<std::pair<Int3, CellShape>> fp;
    Plot::Footprint(part, fp);
    std::vector<TerrainEdit> edits;
    for (const auto& c : fp) {
        const Vector3 centre = plot.toBody({ (c.first.x + 0.5f) * kCellSize, (c.first.y + 0.5f) * kCellSize, (c.first.z + 0.5f) * kCellSize });
        if (body.densityAt(centre) > ISO_LEVEL - 0.15f) edits.push_back({ centre, 0.55f, -0.8f });
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
    if (charge && !terrainTool.creative && !inventory.remove(PartCostItem(t), PartCostCount(t))) {
        GameUI_Toast(TextFormat("Need %d %s for a %s", PartCostCount(t), PartCostItem(t), PartName(t)));
        Audio_Play(Sfx::Error, 0.5f);
        return false;
    }
    plot.place(part);
    LayFoundation(body, plot, part);
    return true;
}

void BuildTool_Update(Scene& world, bool inputEnabled, bool usable) {
    if (!usable) { buildTool.active = false; return; }
    if (inputEnabled && IsKeyPressed(KEY_G)) {
        buildTool.active = !buildTool.active;
        GameUI_Toast(buildTool.active ? "Build mode: wheel picks a part, R/V rotate, RMB place, LMB remove (G to leave)" : "Build mode off");
    }
    aim = Aim{};
    if (!buildTool.active) return;
    if (inputEnabled) {
        const float wheel = GetMouseWheelMove();
        const int n = (int)PartType::Count;
        if (wheel > 0) buildTool.part = (buildTool.part + 1) % n;
        if (wheel < 0) buildTool.part = (buildTool.part + n - 1) % n;
        if (IsKeyPressed(KEY_R)) buildTool.orient = RotateOrientation(buildTool.orient, 1);
        if (IsKeyPressed(KEY_V)) buildTool.orient = RotateOrientation(buildTool.orient, 0);
    }
    Planetoid* body = player.nearest;
    if (!body) return;
    const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
    const Vector3 lo = body->toLocal(world.camera.position), ld = body->dirToLocal(dir);
    aim.body = body;
    aim.part = { (PartType)buildTool.part, { 0, 0, 0 }, (uint8_t)buildTool.orient };

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
        const PartType t = aim.removePlot->parts[aim.removePart].type;
        if (aim.removePlot->remove(aim.removePart)) {
            if (!terrainTool.creative) inventory.add(PartCostItem(t), PartCostCount(t));
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
    DrawText(TextFormat("BUILD  %s  (%d %s, have %d)", PartName(t), PartCostCount(t), PartCostItem(t), inventory.count(PartCostItem(t))),
             cx - 200, GetScreenHeight() - 110, 22, SKYBLUE);
    DrawText("wheel part   R / V rotate   RMB place   LMB remove   G leave", 10, GetScreenHeight() - 30, 20, Fade(RAYWHITE, 0.7f));
}
