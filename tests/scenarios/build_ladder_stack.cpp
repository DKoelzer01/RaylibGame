// build/ladder-stack (#1640, the user: "Ladder placement after first ladder ... rotates around a corner, so you cannot stack
// ladders"): four one-cell brick columns, open to the east, west, north and south. Up each, through the build tool's own
// aim and RMB (creative: free): the first ladder against the wall's face, then three more each aimed at the top of the one
// below. All four must stand in one column in front of the wall, flush on it, with the first ladder's facing.
#include "scenario.h"
#include "buildtool.h"
#include "player.h"
#include "playerinput.h"
#include "plots.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
Plot* plot = nullptr;
std::string aimNote;
Int3 base{};
const Int3 kOut[4] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
const char* kName[4] = { "east", "west", "north", "south" };
Int3 WallCell(int f) { return { base.x + 6 * (f - 2), base.y, base.z + 6 }; } // the column's foot (a row of four, 3 m apart)
Vector3 At(float x, float y, float z) { return plot->toBody({ x * kCellSize, y * kCellSize, z * kCellSize }); }
}

SCENARIO("build/ladder-stack") {
    s.tags = { "world", "build" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 120.0;
    s.step("four brick columns, one facing each way", 20.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const Vector3 feet = c.body->toLocalD(player.position).toFloat();
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(feet, 0.0f, true);
        base = plot->cellAt(plot->toPlot(feet));
        base.y += 2; // (clear of the ground)
        const int mat = std::max(0, BuildMaterialForItem("stone_brick"));
        int built = 0;
        for (int f = 0; f < 4; ++f)
            for (int y = 0; y < 10; ++y) { const Int3 w = WallCell(f); built += plot->place(Part{ PartType::Cube, { w.x, w.y + y, w.z }, 0, (uint8_t)mat }) >= 0; }
        c.eq("brick cells placed", built, 40);
        terrainTool.creative = true;
        buildTool.active = true; buildTool.orientManual = false;
        buildTool.part = (int)PartType::Ladder;
        return true;
    });
    for (int f = 0; f < 4; ++f) {
        s.step(TextFormat("four ladders up the %s face", kName[f]), 20.0, [f](ScenarioCtx& c) {
            static int placed = 0, phase = 0, firstOrient = -1;
            if (c.t < 0.05) { placed = 0; phase = 0; firstOrient = -1; }
            const Int3 w = WallCell(f), o = kOut[f];
            const Int3 col{ w.x + o.x, w.y, w.z + o.z };                          // the ladders' column: a cell out from the wall
            const Vector3 up = Vector3Normalize(Vector3Subtract(At(0, 1, 0), At(0, 0, 0)));
            // Aim: the first at the wall's face (row 1), the next at the previous ladder's top (its frame, against the wall).
            Vector3 eye, target;
            if (placed == 0) {
                target = At(w.x + 0.5f + 0.5f * o.x, w.y + 1.5f, w.z + 0.5f + 0.5f * o.z);
                eye = At(w.x + 0.5f + 5.0f * o.x, w.y + 2.0f, w.z + 0.5f + 5.0f * o.z);
            } else {
                const float fx = col.x + 0.5f - 0.42f * o.x, fz = col.z + 0.5f - 0.42f * o.z; // (the frame: 15% of the cell, by the wall)
                const float top = (float)(col.y + placed + 1);                        // (rows 1..: the previous ladder's top)
                target = At(fx, top, fz);
                eye = At(fx + 1.2f * o.x, top + 3.0f, fz + 1.2f * o.z);
            }
            c.camera(eye, target, up); // (flying: the build tool aims from the camera)
            if (phase == 0) { phase = 1; return false; } // the camera in place a frame first
            if (phase == 1) { // (flying never sets the player's body, #1555: the build tool aims on it)
                player.nearest = c.body;
                buildTool.active = true; buildTool.part = (int)PartType::Ladder; // (the game's own update, flying, turns it off each frame)
                PlayerInput in; in.press(Act::Secondary); BuildTool_Update(c.world, in, true);
                Vector3 no{}; Int3 o2{};
                const Plot* ap = BuildTool_AimPlacement(&no, &o2);
                aimNote = TextFormat("aim plot %s, origin (%d, %d, %d), orient %d", ap == plot ? "the walls'" : ap ? "another" : "none", o2.x, o2.y, o2.z, BuildTool_AimOrient());
                phase = 2; return false;
            }
            // the placed ladder: in the column, the row above the last
            const int row = col.y + 1 + placed;
            const int at = plot->partAt({ col.x, row, col.z });
            const bool ladder = at >= 0 && plot->parts[at].type == PartType::Ladder;
            if (!ladder) { // (what the build tool aimed at, for the failure's story)
                int ladders = 0; for (const Part& p : plot->parts) ladders += p.type == PartType::Ladder;
                c.note(TextFormat("%s ladder %d: %s; the column (%d, %d, %d), ladders on the plot %d", kName[f], placed + 1, aimNote.c_str(), col.x, row, col.z, ladders));
            }
            c.check(TextFormat("%s: ladder %d in the column, row %d", kName[f], placed + 1, row - col.y), ladder);
            if (!ladder) return true;
            const int orient = plot->parts[at].orient;
            if (placed == 0) firstOrient = orient;
            else c.eq(TextFormat("%s: ladder %d faces as the first", kName[f], placed + 1), orient, firstOrient);
            // flush on the wall: the frame's open side (its shape's plane normal) points out of the wall
            std::vector<std::pair<Int3, CellShape>> fp;
            Plot::Footprint(plot->parts[at], fp);
            const Vector3 n = fp.empty() ? Vector3{} : fp[0].second.n[0];
            c.check(TextFormat("%s: ladder %d against the wall", kName[f], placed + 1), (int)lroundf(n.x) == o.x && (int)lroundf(n.z) == o.z,
                    TextFormat("open side (%.0f, %.0f, %.0f)", n.x, n.y, n.z));
            if (++placed == 4) { c.grab(TextFormat("ladders-%s", kName[f])); return true; }
            phase = 0;
            return false;
        });
    }
    s.step("done", 2.0, [](ScenarioCtx&) { buildTool.active = false; return true; });
}
