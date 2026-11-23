// water/shore-dig (#1590, #1091): pits dug in the beach a few metres from the sea, as a player digs them with a shovel. The
// playtests showed jagged shards of water surface hanging in the air around such holes. The nearest shore to the start (as
// water/beach finds it); three pits 0.5, 3 and 6 m inland of the waterline, each dug down past the sea level; the camera
// 3.5 m over the beach 8 m inland, looking down at them.
#include "scenario.h"
#include "fluid.h"
#include "voxel.h"
#include "waterfield.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

namespace {
struct State { bool found = false; Vector3 up{}, out{}, along{}; } S;
Vector3 At(ScenarioCtx& c, Vector3 dir, float h) { return Vector3Scale(dir, c.field->radius + h); }
// A point on the beach: u metres inland of the waterline (negative: out to sea), v along it.
Vector3 Beach(ScenarioCtx& c, float u, float v) {
    return Vector3Normalize(Vector3Add(Vector3Scale(S.up, c.field->radius), Vector3Add(Vector3Scale(S.out, -u), Vector3Scale(S.along, v))));
}
}

SCENARIO("water/shore-dig") {
    s.tags = { "water", "visual" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day" };
    s.timeout = 180.0;

    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("find the nearest shore", 1.0, [](ScenarioCtx& c) {
        const Vector3 home = c.field->findSpawnDir(c.body->seed);
        const Vector3 ax = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), ay = Vector3CrossProduct(ax, home);
        const float R = c.field->radius;
        float best = 1e30f;
        for (int k = 0; k < 16; ++k) {
            const float a = 6.2831853f * k / 16.0f;
            const Vector3 dirT = Vector3Add(Vector3Scale(ax, cosf(a)), Vector3Scale(ay, sinf(a)));
            Vector3 prev = home;
            for (float d = 2.0f; d < 4000.0f && d < best; d += 2.0f) {
                const Vector3 dd = Vector3Normalize(Vector3Add(Vector3Scale(home, R), Vector3Scale(dirT, d)));
                if (c.field->elevation(dd) < 0.0f) { // the waterline: the first ground under the sea level
                    best = d; S.found = true; S.up = dd;
                    S.out = Vector3Normalize(Vector3Subtract(Vector3Scale(dd, R), Vector3Scale(prev, R)));
                    S.along = Vector3Normalize(Vector3CrossProduct(dd, S.out));
                    break;
                }
                prev = dd;
            }
        }
        if (!c.check("a sea shore within 4 km of the start", S.found)) { c.fail("no shore"); return true; }
        c.metric("nearest shore from the start (m)", best, "m");
        for (float u : { 2.0f, 4.0f, 6.0f }) c.note(TextFormat("beach %.0f m inland: ground %.2f m", u, c.field->elevation(Beach(c, u, 0.0f))));
        return true;
    });
    auto view = [](ScenarioCtx& c) {
        const Vector3 eyeDir = Beach(c, 8.0f, -5.0f), target = Beach(c, 1.5f, 0.0f);
        c.camera(At(c, eyeDir, std::max(0.0f, c.field->elevation(eyeDir)) + 3.5f), At(c, target, 0.0f), S.up);
    };
    s.step("the shore before", 30.0, [view](ScenarioCtx& c) {
        view(c);
        c.waitingFor("the shore built (an idle world), then 2 s");
        return c.t > 2.0 && WorldIdle();
    });
    s.step("dig three pits", 1.0, [view](ScenarioCtx& c) {
        view(c);
        // A shovel's carve: 1.2 m spheres, each pit three strokes deep from the ground down.
        std::vector<TerrainEdit> cut;
        for (float u : { 0.5f, 3.0f, 6.0f }) {
            const Vector3 d = Beach(c, u, 0.0f);
            const float g = c.field->elevation(d);
            for (int k = 0; k < 3; ++k) cut.push_back({ At(c, d, g - 0.6f - 0.6f * k), 1.2f, -1.0f, true });
        }
        ApplyTerrainEdits(*c.body, cut.data(), cut.size());
        c.cue();
        return true;
    });
    s.step("watch the pits", 25.0, [view](ScenarioCtx& c) {
        view(c);
        c.waitingFor("20 s of water");
        if (c.t < 20.0) return false;
        const Fluid::Stats st = Fluid::GetStats();
        c.note(TextFormat("after 20 s: %d sites, %d bricks (%d awake), %.1f m^3", st.sites, st.bricks, st.activeBricks, st.volume));
        // Who draws the waterline (#1590): per point along it, whether a published simulation column covers it (S: the
        // simulation's rules, G: the generated water's), and the cells 0.1 m under and over the sea level.
        if (getenv("SHOREDBG")) {
            const WaterField::SnapshotPtr snap = WaterField::Current();
            for (float u : { -1.0f, 0.0f, 1.0f, 2.0f }) {
                std::string row;
                for (float v = -8.0f; v <= 8.0f; v += 1.0f) {
                    const Vector3 d = Beach(c, u, v), p = At(c, d, 0.0f);
                    WaterField::ChunkSampler cs(snap.get(), c.body->id, Vector3SubtractValue(p, 1.0f), Vector3AddValue(p, 1.0f));
                    row += cs.covered(p) ? 'S' : 'G';
                }
                c.note(TextFormat("waterline %+.0f m inland: %s", u, row.c_str()));
            }
            for (float v = -8.0f; v <= 8.0f; v += 2.0f) {
                const Vector3 d = Beach(c, 0.0f, v);
                c.note(TextFormat("at the waterline, %+.0f m along: ground %.2f m; -0.1 m: %s; +0.1 m: %s", v, c.field->elevation(d),
                                  Fluid::DebugCell(*c.body, At(c, d, -0.1f)).c_str(), Fluid::DebugCell(*c.body, At(c, d, 0.1f)).c_str()));
            }
            // The cells round each pit, column by column and top to bottom (#1590): shows which cells hold fill (and with what
            // flags) next to the dug holes, where the stray surface shards hang.
            for (float u : { 0.5f, 3.0f, 6.0f }) {
                for (float du = -1.5f; du <= 1.51f; du += 0.5f) {
                    const Vector3 col = Beach(c, u + du, 0.0f);
                    std::string row;
                    for (float h = 1.0f; h >= -2.51f; h -= 0.25f) {
                        const std::string cell = Fluid::DebugCell(*c.body, At(c, col, h));
                        int fill = 0; char flags[64] = "";
                        if (sscanf(cell.c_str(), "fill %d flags %63s", &fill, flags) == 2) row += TextFormat("%d/%s ", fill, flags);
                        else row += "- ";
                    }
                    c.note(TextFormat("pit %.1f col %+.1f: %s ground %.2f", u, du, row.c_str(), c.field->elevation(col)));
                }
            }
        }
        c.keepCapture("the look of the dug shore (#1590)");
        return true;
    });
}
