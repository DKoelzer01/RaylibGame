#include "machinemodels.h"
#include "buildtool.h"
#include "playerinput.h"
#include "waterworks.h"
#include "survival.h"
#include "farming.h"
#include "audio.h"
#include "gamerng.h"
#include "gameui.h"
#include "inventory.h"
#include "player.h"
#include "terraintool.h"
#include "world.h"
#include "worldtypes.h"
#include "factory.h"
#include "blueprint.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>

// (buildTool: the local PlayerSession's, src/space/playersession.cpp, #313)

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

bool IsMachineItem(const std::string& id) {
    for (int t = 0; t < (int)PartType::Count; ++t) if (id == PartCostItem((PartType)t) && IsMachine((PartType)t)) return true;
    return false;
}

// Foundation (#210): the terrain inside the part's cells is subtracted, exactly those boxes (a carve, as mining, but with
// no yield), and nothing else changes: no fill under it, no brush reaching past it (the plot carries parts above ground).
void BuildTool_LayFoundation(Planetoid& body, Plot& plot, const Part& part) {
    std::vector<std::pair<Int3, CellShape>> fp;
    Plot::Footprint(part, fp);
    // The plot's axes in body-local space.
    const Vector3 o = plot.toBody({ 0, 0, 0 });
    const Vector3 ax = Vector3Normalize(Vector3Subtract(plot.toBody({ 1, 0, 0 }), o));
    const Vector3 ay = Vector3Normalize(Vector3Subtract(plot.toBody({ 0, 1, 0 }), o));
    const Vector3 az = Vector3Normalize(Vector3Subtract(plot.toBody({ 0, 0, 1 }), o));
    const float h = 0.5f * kCellSize;
    std::vector<TerrainEdit> edits;
    for (const auto& c : fp) {
        const Vector3 centre = plot.toBody({ (c.first.x + 0.5f) * kCellSize, (c.first.y + 0.5f) * kCellSize, (c.first.z + 0.5f) * kCellSize });
        // Anything solid in the cell? (its centre and corners)
        bool solid = body.densityAt(centre) >= ISO_LEVEL;
        for (int k = 0; k < 8 && !solid; ++k) {
            const Vector3 corner = Vector3Add(centre, Vector3Add(Vector3Scale(ax, (k & 1) ? h * 0.95f : -h * 0.95f),
                                   Vector3Add(Vector3Scale(ay, (k & 2) ? h * 0.95f : -h * 0.95f), Vector3Scale(az, (k & 4) ? h * 0.95f : -h * 0.95f))));
            solid = body.densityAt(corner) >= ISO_LEVEL;
        }
        if (!solid) continue;
        TerrainEdit e{};
        e.center = centre;
        e.carve = true;
        e.box = true;
        e.half = { h, h, h };
        e.ax = ax; e.ay = ay; e.az = az;
        e.radius = h * 1.7320508f; // its bounding reach (dirty chunks, trees, fluid)
        edits.push_back(e);
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
    if (IsBeam(t)) { // a line of single cubes (#200): quicker to build, each removed on its own
        static std::vector<Part> cubes;
        Plot::SplitBeam(part, cubes);
        for (const Part& c : cubes) plot.place(c);
    } else plot.place(part);
    if (IsMachine(part.type) || IsRocketPart(part.type)) QuestBook_NotePlaced(PartCostItem(part.type));
    BuildTool_LayFoundation(body, plot, part);
    return true;
}

// Breaking built parts: LMB held on a part for kBreakHold seconds (in build mode or out of it; a
// part in front of the ground takes the click from the terrain tool). A short hold, so a stray click
// never deletes a wall.
namespace {
    constexpr float kBreakHold = 0.45f;
    const Plot* breakPlot = nullptr;
    int breakPart = -1;
    float breakHeld = 0.0f;
    bool breakAiming = false;
}

float BuildTool_BreakProgress() { return breakAiming && breakHeld > 0.0f ? std::min(1.0f, breakHeld / kBreakHold) : 0.0f; }
bool BuildTool_AimingPart() { return breakAiming; }
int BuildTool_AimOrient() { return aim.plot || aim.newOrigin.x != 0.0f || aim.newOrigin.y != 0.0f || aim.newOrigin.z != 0.0f ? (int)aim.part.orient : -1; }

static void UpdateBreak(Scene& world, const PlayerInput& in) {
    breakAiming = false;
    Planetoid* body = player.nearest;
    if (!body || !body->plots || player.flyMode) { breakHeld = 0.0f; breakPart = -1; return; }
    const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
    const Vector3 lo = body->toLocal(world.camera.position), ld = body->dirToLocal(dir);
    float best = buildTool.active ? buildTool.reach : terrainTool.reach;
    Plot* hitPlot = nullptr;
    Int3 hitCell{};
    for (auto& pl : body->plots->plots) {
        Int3 c, f; float d;
        if (pl->raycast(lo, ld, best, c, f, d)) { best = d; hitPlot = pl.get(); hitCell = c; }
    }
    if (!hitPlot) { breakHeld = 0.0f; breakPart = -1; return; }
    TerrainHit th;
    if (RaycastTerrain(world, Ray{ world.camera.position, dir }, best, th) && th.planetoid == body && th.distance < best * body->scale) { breakHeld = 0.0f; breakPart = -1; return; }
    const int part = hitPlot->partAt(hitCell);
    if (part < 0) { breakHeld = 0.0f; breakPart = -1; return; }
    breakAiming = true;
    if (!in.down(Act::Primary)) { breakHeld = std::min(breakHeld, 0.0f); breakPart = -1; return; }
    if (hitPlot != breakPlot || part != breakPart) { breakPlot = hitPlot; breakPart = part; breakHeld = std::max(breakHeld, 0.0f); }
    breakHeld += GameFrameTime();
    if (breakHeld < kBreakHold) return;
    const Part removed = hitPlot->parts[part];
    if (hitPlot->remove(part)) {
        if (removed.type == PartType::Meteorite && !terrainTool.creative) { // (#177) space rock: nickel-iron, now and then platinum or iridium ore
            const int n = GameRng_Int(6, 10), pt = GameRng_Int(0, 99) < 35 ? 2 : 0, ir = GameRng_Int(0, 99) < 12 ? 1 : 0; // the authority rolls (#316)
            inventory.add("nickel_iron_ore", n);
            if (pt) inventory.add("platinum_ore", pt);
            if (ir) inventory.add("iridium_ore", ir);
            GameUI_Toast(TextFormat("The meteorite breaks apart: +%d nickel-iron%s%s", n, pt ? ", +2 platinum ore" : "", ir ? ", +1 iridium ore" : ""));
        } else if (!terrainTool.creative) inventory.add(PartItem(removed), PartCostCount(removed.type));
        Audio_Play(Sfx::Dig, 0.5f, 1.3f);
    }
    breakPart = -1;
    breakHeld = -0.15f; // a breath before the next part (holding on clears a wall part by part)
}

void BuildTool_Update(Scene& world, const PlayerInput& in, bool usable) {
    if (!usable) { buildTool.active = false; breakAiming = false; breakHeld = 0.0f; return; }
    UpdateBreak(world, in);
    if (in.hit(Act::BuildMode)) {
        buildTool.active = !buildTool.active;
        buildTool.orientManual = false; // parts face you again (#196)
        GameUI_Toast(buildTool.active ? "Build mode: Ctrl+wheel shape, parts face you (R/V rotate by hand), RMB place, hold LMB remove, B save blueprint, N blueprints (G to leave)" : "Build mode off");
    }
    if (buildTool.part != buildTool.autoPart) { buildTool.autoPart = buildTool.part; buildTool.orientManual = false; } // another part: automatic facing
    aim = Aim{};
    if (!buildTool.active) return;
    if (in.enabled) {
        const float wheel = in.wheel;
        // The wheel picks the shape (any shape from the selected material); machines come from the hotbar.
        const bool modWheel = IsKeyDown(KEY_LEFT_CONTROL) || IsKeyDown(KEY_RIGHT_CONTROL) || IsKeyDown(KEY_LEFT_ALT) || IsKeyDown(KEY_RIGHT_ALT); // #227: the plain wheel scrolls the hotbar
        if (modWheel && IsShape((PartType)buildTool.part) && wheel != 0.0f) { // Ctrl/Alt + wheel: step through kShapes
            int k = 0;
            for (int i = 0; i < kShapeCount; ++i) if ((int)kShapes[i] == buildTool.part) k = i;
            k = (k + (wheel > 0 ? 1 : kShapeCount - 1)) % kShapeCount;
            buildTool.part = (int)kShapes[k];
        }
        if (in.hit(Act::Rotate) || in.hit(Act::Flip)) {
            if (!buildTool.orientManual) GameUI_Toast("Rotating by hand (parts face you again with the next part, or G)");
            buildTool.orientManual = true; // the orientation sticks from now on (#196); it starts from the automatic one shown
            buildTool.autoPart = buildTool.part;
            buildTool.orient = RotateOrientation(buildTool.orient, in.hit(Act::Rotate) ? 1 : 0);
        }
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
    // A part bigger than a cell is centred on the aimed point across the face it is placed on (#196: a 2x2 machine hung
    // off one corner of the aimed cell). `normal` is the face's axis (a plot axis); the other two axes move the
    // footprint so its middle is nearest the hit point (in cells).
    auto centreOn = [&](Plot& pl, Int3 normal, Vector3 hitBodyLocal) {
        std::vector<std::pair<Int3, CellShape>> fp;
        Plot::Footprint(aim.part, fp);
        if (fp.size() <= 1) return;
        Int3 lo{ INT32_MAX, INT32_MAX, INT32_MAX }, hi{ INT32_MIN, INT32_MIN, INT32_MIN };
        for (const auto& c : fp) {
            lo = { std::min(lo.x, c.first.x), std::min(lo.y, c.first.y), std::min(lo.z, c.first.z) };
            hi = { std::max(hi.x, c.first.x), std::max(hi.y, c.first.y), std::max(hi.z, c.first.z) };
        }
        const Vector3 h = Vector3Scale(pl.toPlot(hitBodyLocal), 1.0f / kCellSize);
        auto shift = [](float hv, int a, int b) { return (int)lroundf(hv - 0.5f * (float)(a + b + 1)); };
        if (normal.x == 0) aim.part.origin.x += shift(h.x, lo.x, hi.x);
        if (normal.y == 0) aim.part.origin.y += shift(h.y, lo.y, hi.y);
        if (normal.z == 0) aim.part.origin.z += shift(h.z, lo.z, hi.z);
    };
    // Automatic facing (#196, until R/V is pressed): the part's front (its local +z, where machines have their door)
    // turns toward the player - on a wall, out of the wall; on the ground or a floor, back along the view.
    auto faceToward = [&](Int3 target) {
        if (buildTool.orientManual) return;
        int o = 0;
        for (int k = 0; k < 4; ++k, o = RotateOrientation(o, 1)) {
            const int* m = Orientation(o);
            if (m[2] == target.x && m[5] == target.y && m[8] == target.z) break;
        }
        buildTool.orient = o;          // shown, and where a hand rotation starts from
        aim.part.orient = (uint8_t)o;
    };
    auto towardViewer = [&](const Plot& pl) {
        const float vx = Vector3DotProduct(ld, pl.east), vz = Vector3DotProduct(ld, pl.north);
        return fabsf(vx) > fabsf(vz) ? Int3{ vx > 0.0f ? -1 : 1, 0, 0 } : Int3{ 0, 0, vz > 0.0f ? -1 : 1 };
    };
    if (hitPlot && (!terrain || bestPart <= th.distance)) {
        aim.plot = hitPlot;
        faceToward(hitFace.y == 0 ? hitFace : towardViewer(*hitPlot));
        aim.part.origin = { hitCell.x + hitFace.x, hitCell.y + hitFace.y, hitCell.z + hitFace.z };
        centreOn(*hitPlot, hitFace, Vector3Add(lo, Vector3Scale(ld, bestPart)));
        {   // A part bigger than a cell grows from its origin, which could reach back into the block
            // aimed at: slide it out along the aimed face until every cell of it is in front of that face.
            std::vector<std::pair<Int3, CellShape>> fp;
            Plot::Footprint(aim.part, fp);
            int behind = 0;
            for (const auto& c : fp) {
                const int d = (c.first.x - hitCell.x) * hitFace.x + (c.first.y - hitCell.y) * hitFace.y + (c.first.z - hitCell.z) * hitFace.z;
                behind = std::max(behind, 1 - d);
            }
            aim.part.origin = { aim.part.origin.x + hitFace.x * behind, aim.part.origin.y + hitFace.y * behind, aim.part.origin.z + hitFace.z * behind };
        }
        aim.removePlot = hitPlot;
        aim.removePart = hitPlot->partAt(hitCell);
    } else if (terrain) {
        const Vector3 ground = th.local;
        Plot* near = body->plots ? body->plots->plotFor(ground, 6.0f, false) : nullptr;
        if (near) {
            aim.plot = near;
            faceToward(towardViewer(*near));
            // The cell just above the ground at the aim.
            aim.part.origin = near->cellAt(near->toPlot(Vector3Add(ground, Vector3Scale(near->up, 0.1f))));
            centreOn(*near, Int3{ 0, 1, 0 }, ground); // on the ground: the bottom face's middle at the aim
        } else {
            aim.newOrigin = ground;
            aim.part.origin = { 0, 0, 0 };
        }
    } else {
        return;
    }
    if (aim.plot) aim.canPlace = aim.plot->canPlace(aim.part);
    else aim.canPlace = true;

    if (!in.enabled) return;
    if (Blueprint_Keys(body, aim.removePlot)) return; // B saves the plot looked at, N cycles blueprints
    if (const Blueprint* bp = Blueprint_Active(); bp && in.hit(Act::Secondary)) {
        Plot* target = aim.plot;
        if (!target) {
            if (!body->plots) body->plots = std::make_unique<PlotSet>();
            target = body->plots->plotFor(aim.newOrigin, 0.0f, true);
        }
        std::string missing;
        const int n = Blueprint_Place(*body, *target, *bp, aim.part.origin, true, missing);
        GameUI_Toast(missing.empty() ? TextFormat("Placed %d parts of %s", n, bp->name.c_str()) : TextFormat("Placed %d parts; still need %s", n, missing.c_str()));
        if (n > 0) Audio_Play(Sfx::Build, 0.5f, 1.0f);
        return;
    }
    if (in.hit(Act::Secondary) && aim.canPlace) {
        Plot* target = aim.plot;
        if (!target) {
            if (!body->plots) body->plots = std::make_unique<PlotSet>();
            target = body->plots->plotFor(aim.newOrigin, 0.0f, true);
        }
        if (BuildTool_Place(*body, *target, aim.part, true)) Audio_Play(Sfx::Build, 0.5f, 1.2f);
    } // LMB on a part: UpdateBreak (a short hold)
}

void BuildTool_Draw3D(Scene& world) {
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
        // In plot coordinates under the plot's own double-precision model, as the blocks are drawn: float corners
        // 60 km out drifted off them and shook as the planet turned (#271).
        Vector3 v[8];
        for (int i = 0; i < 8; ++i) {
            // Inflated 2 cm so the lines stand off the faces of the blocks they outline (they z-fought, #199).
            constexpr float kPad = 0.02f;
            auto edge = [](int cell, int hi) { return (cell + hi) * kCellSize + (hi ? kPad : -kPad); };
            v[i] = { edge(c.x, i & 1), edge(c.y, (i >> 1) & 1), edge(c.z, (i >> 2) & 1) };
        }
        rlDrawRenderBatchActive();
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(p.cellModel(body, world.renderOrigin)));
        const int e[12][2] = { { 0, 1 }, { 2, 3 }, { 4, 5 }, { 6, 7 }, { 0, 2 }, { 1, 3 }, { 4, 6 }, { 5, 7 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
        for (auto& ed : e) DrawLine3D(v[ed[0]], v[ed[1]], col);
        rlPopMatrix();
    };
    // The part's real shape (#197): each footprint cell's faces (slopes, corners, beams as they will be), translucent
    // green (or red: can't place) with their edges, inflated 2 cm about the cell centre so they stand off the blocks
    // around them (#199), both sides drawn.
    auto ghost = [&](const Plot& p, const std::vector<std::pair<Int3, CellShape>>& cells, Color fill, Color edge) {
        std::vector<std::vector<Vector3>> faces;
        rlDrawRenderBatchActive();
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(p.cellModel(body, world.renderOrigin)));
        rlDisableBackfaceCulling();
        rlDisableDepthMask();
        constexpr float kGrow = 1.0f + 0.04f / kCellSize; // 2 cm each side
        for (const auto& c : cells) {
            CellShapeFaces(c.second, faces);
            const Vector3 lo = { (float)c.first.x, (float)c.first.y, (float)c.first.z };
            auto toPlot = [&](Vector3 v) { // cell units -> plot metres, grown about the cell centre
                const Vector3 g = Vector3Add(Vector3Scale(Vector3Subtract(v, { 0.5f, 0.5f, 0.5f }), kGrow), { 0.5f, 0.5f, 0.5f });
                return Vector3Scale(Vector3Add(lo, g), kCellSize);
            };
            for (const auto& poly : faces) {
                for (size_t i = 1; i + 1 < poly.size(); ++i) DrawTriangle3D(toPlot(poly[0]), toPlot(poly[i]), toPlot(poly[i + 1]), fill);
                for (size_t i = 0; i < poly.size(); ++i) DrawLine3D(toPlot(poly[i]), toPlot(poly[(i + 1) % poly.size()]), edge);
            }
        }
        rlDrawRenderBatchActive();
        rlEnableDepthMask();
        rlEnableBackfaceCulling();
        rlPopMatrix();
    };
    std::vector<std::pair<Int3, CellShape>> fp;
    Plot::Footprint(aim.part, fp);
    if (Blueprint_Active()) Blueprint_DrawGhost(body, *frame, aim.part.origin, true, world.renderOrigin); // the whole blueprint instead of one part
    else ghost(*frame, fp, aim.canPlace ? Fade(GREEN, 0.28f) : Fade(RED, 0.30f), aim.canPlace ? GREEN : RED);
    if (aim.removePlot && aim.removePart >= 0) {
        Plot::Footprint(aim.removePlot->parts[aim.removePart], fp);
        for (const auto& c : fp) box(*aim.removePlot, c.first, Fade(ORANGE, 0.7f));
    }
}

void BuildTool_DrawHUD() {
    {   // The bar under the crosshair: a part breaking (amber, filling) or the swing recovering (white).
        const int bx = GetScreenWidth() / 2 - 30, by = GetScreenHeight() / 2 + 20, bw = 60;
        const float brk = BuildTool_BreakProgress(), cd = buildTool.active ? 0.0f : TerrainTool_Cooldown();
        if (brk > 0.0f || cd > 0.0f) {
            DrawRectangle(bx - 1, by - 1, bw + 2, 6, Fade(BLACK, 0.55f));
            if (brk > 0.0f) DrawRectangle(bx, by, (int)(bw * brk), 4, Color{ 255, 180, 60, 255 });
            else DrawRectangle(bx, by, (int)(bw * (1.0f - cd)), 4, Fade(RAYWHITE, 0.85f));
        }
    }
    if (!buildTool.active) return;
    Blueprint_DrawHUD();
    const PartType t = (PartType)buildTool.part;
    const int cx = GetScreenWidth() / 2, cy = GetScreenHeight() / 2;
    DrawLine(cx - 8, cy, cx + 8, cy, RAYWHITE);
    DrawLine(cx, cy - 8, cx, cy + 8, RAYWHITE);
    const Part p{ t, { 0, 0, 0 }, 0, (uint8_t)buildTool.material };
    const char* item = PartItem(p);
    const char* label = IsShape(t) ? TextFormat("%s %s", GetBuildMaterial(buildTool.material).name, PartName(t)) : PartName(t);
    const char* line = TextFormat("BUILD  %s  (%d %s, have %d)", label, PartCostCount(t), item, inventory.count(item));
    DrawText(line, cx - MeasureText(line, 22) / 2, GetScreenHeight() - 150, 22, SKYBLUE);
    DrawText(IsShape(t) ? "Ctrl+wheel shape   wheel hotbar   R / V rotate   RMB place   LMB remove   G leave   1-0 hotbar"
                        : "R / V rotate   RMB place   LMB remove   G leave   1-0 hotbar",
             10, GetScreenHeight() - 30, 20, Fade(RAYWHITE, 0.7f));
}

// ---------------------------------------------------------------------------------------------
// Machines: look at one to see its state; I puts matching items in, O takes the outputs.
// ---------------------------------------------------------------------------------------------

namespace {
    struct MachineAim { uint64_t body = 0; uint32_t plot = 0; int part = -1; PartType type = PartType::Count; } machineAim;
}

void Machine_Update(Scene& world, const PlayerInput& in) {
    machineAim = MachineAim{};
    Planetoid* body = player.nearest;
    if (!body || !body->plots || player.flyMode) return;
    const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
    const Vector3 lo = body->toLocal(world.camera.position), ld = body->dirToLocal(dir);
    float best = 6.0f;
    for (auto& pl : body->plots->plots) {
        Int3 c, f; float d;
        if (pl->raycast(lo, ld, best, c, f, d)) {
            const int part = pl->partAt(c);
            if (part >= 0 && (IsMachine(pl->parts[part].type) || pl->parts[part].type == PartType::Farmland || pl->parts[part].type == PartType::HandPump)) {
                best = d;
                machineAim = { body->id, pl->id, part, pl->parts[part].type };
            }
        }
        // Chimneys stick out of the cells (#285): aimable too, when nearer.
        const Vector3 po = pl->toPlot(lo);
        const Vector3 pd = { Vector3DotProduct(ld, pl->east), Vector3DotProduct(ld, pl->up), Vector3DotProduct(ld, pl->north) };
        int dp = -1; float dd = 0.0f;
        if (Vector3Length(po) < 200.0f && Machine_RaycastDetails(*pl, po, pd, best, dp, dd)) {
            best = dd;
            machineAim = { body->id, pl->id, dp, pl->parts[dp].type };
        }
    }
    if (machineAim.part < 0 || !in.enabled) return;
    if (machineAim.type == PartType::HandPump) {
        // Hand pump (#280): hold RMB (outside build mode) to work it; the water goes down its pipes at the next pass.
        if (!buildTool.active && in.down(Act::Secondary)) {
            Waterworks_HandPump(machineAim.body, machineAim.plot, machineAim.part, GameFrameTime());
            static double nextSound = 0.0;
            if (GetTime() > nextSound) { Audio_Play(Sfx::Pickup, 0.25f, 0.6f); nextSound = GetTime() + 0.6; }
        }
        return;
    }
    if (machineAim.type == PartType::Farmland) {
        // Farmland: I plants seeds you carry, O harvests a ripe crop.
        Planetoid* fb = Planetoid::find(machineAim.body);
        Plot* fp = fb && fb->plots ? fb->plots->byId(machineAim.plot) : nullptr;
        if (fb && fp && in.hit(Act::Insert)) Farm_Plant(*fb, *fp, machineAim.part);
        if (fb && fp && in.hit(Act::Extract)) Farm_Harvest(*fb, *fp, machineAim.part);
        return;
    }
    // RMB (outside build mode, where RMB places parts) opens the machine's own screen: its slots
    // beside the inventory, like the player's.
    if (!buildTool.active && in.hit(Act::Secondary)) {
        GameUI_OpenContainer(machineAim.body, machineAim.plot, machineAim.part, (int)machineAim.type);
        Audio_Play(Sfx::Pickup, 0.3f, 1.0f);
        return;
    }
    if (in.hit(Act::Insert)) {
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
    if (in.hit(Act::Extract)) { Factory::Extract(machineAim.body, machineAim.plot, machineAim.part); Audio_Play(Sfx::Pickup, 0.4f, 1.2f); }
}

void Machine_DrawHUD() {
    if (machineAim.part < 0 || GameUI_ContainerOpen()) return; // the open screen watches its own machine
    if (machineAim.type == PartType::Farmland) {
        Planetoid* fb = Planetoid::find(machineAim.body);
        Plot* fp = fb && fb->plots ? fb->plots->byId(machineAim.plot) : nullptr;
        if (fb && fp) Farm_DrawAimHUD(*fb, *fp, machineAim.part);
        return;
    }
    if (machineAim.type == PartType::HandPump) {
        const int x = GetScreenWidth() / 2 + 20, y = GetScreenHeight() / 2 + 40;
        DrawRectangle(x - 8, y - 6, 380, 52, Fade(BLACK, 0.55f));
        DrawText("hand pump", x, y, 20, SKYBLUE);
        DrawText(TextFormat("hold RMB to pump %.0f L/s into its pipes", kHandPumpRate * 1000.0f), x, y + 24, 16, RAYWHITE);
        return;
    }
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
    else DrawText("RMB open   I put items in   O take the output", x, y + 80, 16, Fade(RAYWHITE, 0.7f));
}
