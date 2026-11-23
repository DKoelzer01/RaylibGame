// visual/ores (#1433): each ore kind in the rock, side by side. For iron, copper, tin, coal, bauxite and ilmenite the
// nearest deposit to the start is opened: a shaft down to its centre and a 4 m room there. The camera stands in the room
// with a light at the eye (RAYLIBGAME_CAMERA_LIGHT) and looks at its wall; one shot a kind, compared with its baseline.
// For a person to judge whether the kinds can be told apart (also under the colour vision simulation of #1017).
#include "scenario.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <string>
#include <vector>

namespace {
const char* kNames[] = { "iron", "copper", "tin", "coal", "bauxite", "ilmenite" };
struct Deposit { bool found = false; Vector3 centre{}, up{}, e1{}; float radius = 0.0f; };
Deposit D[6];
Vector3 Tangent(Vector3 up) { return Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 })); }
}

SCENARIO("visual/ores") {
    s.tags = { "visual", "ores" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day" };
    s.env = { { "RAYLIBGAME_CAMERA_LIGHT", "1" } };
    s.timeout = 300.0;
    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("the nearest deposit of each kind, opened", 1.0, [](ScenarioCtx& c) {
        // Deposits sit on a 96 m grid; bauxite lies below 20 m, ilmenite below 60 m: probe each grid cell at a few depths.
        const Vector3 home = c.field->findSpawnDir(c.body->seed), k1 = Tangent(home), k2 = Vector3CrossProduct(home, k1);
        const float R = c.field->radius;
        int left = 6;
        for (int ring = 0; ring <= 60 && left > 0; ++ring)
            for (int a = -ring; a <= ring && left > 0; ++a)
                for (int b = -ring; b <= ring && left > 0; ++b) {
                    if (std::max(abs(a), abs(b)) != ring) continue; // this ring's cells only
                    const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(home, R), Vector3Add(Vector3Scale(k1, 96.0f * a), Vector3Scale(k2, 96.0f * b))));
                    const float surf = R + c.field->elevation(d);
                    for (float depth : { 12.0f, 30.0f, 70.0f }) {
                        Vector3 centre{}; float rad = 0.0f; int kind = -1;
                        if (!c.field->depositAt(Vector3Scale(d, surf - depth), centre, rad, kind) || kind < 0 || kind > 5 || D[kind].found) continue;
                        // Its centre under the ground over itself (a cell probed under a hill can hold a blob centred above a
                        // valley next to it: the first bauxite sat 21 m in the air).
                        const Vector3 cd = Vector3Normalize(centre);
                        if (R + c.field->elevation(cd) - Vector3Length(centre) < 8.0f) continue;
                        D[kind] = { true, centre, Vector3Normalize(centre), Tangent(Vector3Normalize(centre)), rad };
                        --left;
                    }
                }
        std::vector<TerrainEdit> cut;
        for (int k = 0; k < 6; ++k) {
            if (!c.check(std::string("a ") + kNames[k] + " deposit near the start", D[k].found)) continue;
            const Vector3 up = D[k].up;
            const float top = R + c.field->elevation(up), at = Vector3Length(D[k].centre);
            for (float h = at; h < top + 2.0f; h += 1.5f) cut.push_back({ Vector3Scale(up, h), 1.6f, -1.0f, true }); // the shaft
            cut.push_back({ D[k].centre, 4.0f, -1.0f, true });                                                      // the room
            c.note(TextFormat("%s: %.0f m down, radius %.1f m", kNames[k], top - at, D[k].radius));
        }
        ApplyTerrainEdits(*c.body, cut.data(), cut.size());
        return true;
    });
    for (int k = 0; k < 6; ++k)
        s.step(std::string("look at the ") + kNames[k], 40.0, [k](ScenarioCtx& c) {
            if (!D[k].found) return true;
            const Vector3 eye = Vector3Add(D[k].centre, Vector3Scale(D[k].e1, -1.5f));      // off the room's centre,
            c.camera(eye, Vector3Add(D[k].centre, Vector3Scale(D[k].e1, 4.0f)), D[k].up); // facing its far wall
            c.waitingFor("the room built (an idle world)");
            if (c.t < 3.0 || (!WorldIdle() && c.t < 30.0)) return false;
            c.shot(kNames[k], 100.0); // a picture to compare by eye: the look is meant to change (#1433)
            return true;
        });
    s.step("done", 1.0, [](ScenarioCtx& c) { c.keepCapture("the ore kinds in the rock (#1433): for a person to judge"); return true; });
}
