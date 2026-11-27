#include "multiblockmodels.h"
#include "factory.h"
#include "instancebatch.h"
#include "multiblock.h"
#include "plots.h"
#include "scene.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

// Multiblock animation (#1429 PR 6). See multiblockmodels.h.

namespace {

struct PlotLooks {
    uint32_t version = 0xFFFFFFFFu;
    std::vector<int> formed;   // Plot::formed when built: forming changes no part, so the version alone misses it
    std::vector<MbLook> looks;
};
std::unordered_map<uint64_t, PlotLooks> cache; // by body id ^ plot id
int drawn = 0;

uint64_t PlotKey(uint64_t body, uint32_t plot) { return (body * 1315423911ull) ^ plot; }

Color Scale(Color c, float k) {
    k = std::clamp(k, 0.0f, 1.0f);
    return { (unsigned char)(c.r * k), (unsigned char)(c.g * k), (unsigned char)(c.b * k), 255 };
}

// A hatch's colour by its kind: what flows through it (the items' greens and oranges, the fluids' blues, power's yellow).
Color KindColor(const std::string& kind) {
    if (kind == "item_in") return { 90, 230, 110, 255 };
    if (kind == "item_out") return { 255, 150, 60, 255 };
    if (kind == "fluid_in") return { 80, 150, 255, 255 };
    if (kind == "fluid_out") return { 80, 230, 230, 255 };
    if (kind == "energy") return { 255, 220, 70, 255 };
    if (kind == "rotation") return { 220, 170, 110, 255 };
    return { 200, 200, 200, 255 };
}

// The controller's lamp: its state at a glance from across the room.
Color ControllerColor(MbLampState s, float t) {
    switch (s) {
    case MbLampState::Unformed: return { 110, 24, 20, 255 };                                         // dim red: build or wrench it
    case MbLampState::Idle: return { 40, 150, 70, 255 };                                             // steady green: formed, waiting
    case MbLampState::Working: return Scale({ 90, 255, 130, 255 }, 0.75f + 0.25f * sinf(t * 4.0f)); // bright green, pulsing
    case MbLampState::Stalled: return fmodf(t, 1.0f) < 0.5f ? Color{ 255, 170, 30, 255 } : Color{ 70, 45, 10, 255 }; // blinking amber
    }
    return { 0, 0, 0, 255 };
}

// A flat lamp on the face `out` of the cell at `centre` (plot metres, relative to the anchor), just proud of the casing.
void Lamp(Vector3 centre, MbInt3 out, float side, Color c) {
    const float half = 0.5f * kCellSize + 0.006f, thin = 0.012f;
    const Vector3 at = { centre.x + out.x * half, centre.y + out.y * half, centre.z + out.z * half };
    const Vector3 size = { out.x ? thin : side, out.y ? thin : side, out.z ? thin : side };
    IB_Cube(at, size, c, true);
}

// The seam glow: thin unlit strips along the 12 edges of the box [lo, hi] (plot metres), a hair outside its faces.
void Seams(Vector3 lo, Vector3 hi, Color c) {
    const float e = 0.012f, t = 0.035f;
    lo = { lo.x - e, lo.y - e, lo.z - e };
    hi = { hi.x + e, hi.y + e, hi.z + e };
    const Vector3 mid = Vector3Scale(Vector3Add(lo, hi), 0.5f), len = Vector3Subtract(hi, lo);
    for (float a : { lo.y, hi.y }) for (float b : { lo.z, hi.z }) IB_Cube({ mid.x, a, b }, { len.x + t, t, t }, c, true);
    for (float a : { lo.x, hi.x }) for (float b : { lo.z, hi.z }) IB_Cube({ a, mid.y, b }, { t, len.y + t, t }, c, true);
    for (float a : { lo.x, hi.x }) for (float b : { lo.y, hi.y }) IB_Cube({ a, b, mid.z }, { t, t, len.z + t }, c, true);
}

} // namespace

int MultiblockModels_DrawnCount() { return drawn; }

void MultiblockModels_Draw3D(Scene& world) {
    drawn = 0;
    const float t = (float)GameClock();
    rlDrawRenderBatchActive();
    rlDisableBackfaceCulling(); // the plot frame (east, up, north) is left-handed: mirrored winding
    for (const auto& obj : world.objects) {
        Planetoid* body = dynamic_cast<Planetoid*>(obj.get());
        if (!body || !body->isActive || !body->plots) continue;
        const Vector3 camBody = body->toLocal(world.camera.position);
        for (auto& plot : body->plots->plots) {
            if (Vector3Distance(plot->toBody({ 0, 0, 0 }), camBody) > 300.0f) continue;
            PlotLooks& pl = cache[PlotKey(body->id, plot->id)];
            if (pl.version != plot->version || pl.formed != plot->formed) {
                pl.version = plot->version;
                pl.formed = plot->formed;
                Multiblock_Looks(*plot, pl.looks);
            }
            if (pl.looks.empty()) continue;
            const Vector3 camPlot = plot->toPlot(camBody);
            const std::vector<int> stalled = Factory::StalledParts(body->id, plot->id);
            const Plot& p = *plot;
            Matrix axes = MatrixIdentity();
            axes.m0 = p.east.x; axes.m4 = p.up.x; axes.m8 = p.north.x;
            axes.m1 = p.east.y; axes.m5 = p.up.y; axes.m9 = p.north.y;
            axes.m2 = p.east.z; axes.m6 = p.up.z; axes.m10 = p.north.z;
            for (const MbLook& look : pl.looks) {
                // Anchored at the controller's cell corner: lamps and seams are drawn in plot metres relative to it.
                const MbInt3 a = look.lamps[0].cell;
                const Vector3 anchor = { a.x * kCellSize, a.y * kCellSize, a.z * kCellSize };
                if (Vector3DistanceSqr(anchor, camPlot) > 32.0f * 32.0f) continue; // the casings read on their own beyond this (as machines)
                const bool working = look.formed && Factory::Working(body->id, plot->id, look.controller);
                const bool stall = look.formed && std::find(stalled.begin(), stalled.end(), look.controller) != stalled.end();
                const MbLampState state = Multiblock_LampState(look.formed, working, stall);
                // The plot's model with its translation in double (#284), as the machine details: no shaking against the casings.
                const Vector3 off = Vector3Add(Vector3Add(Vector3Scale(p.east, anchor.x), Vector3Scale(p.up, anchor.y)), Vector3Scale(p.north, anchor.z));
                const Matrix mm = MatrixMultiply(axes, body->localModel(DVec3(p.origin.x + off.x, p.origin.y + off.y, p.origin.z + off.z), world.renderOrigin));
                rlPushMatrix();
                rlMultMatrixf(MatrixToFloat(mm));
                auto centreOf = [&](MbInt3 c) { return Vector3{ (c.x - a.x + 0.5f) * kCellSize, (c.y - a.y + 0.5f) * kCellSize, (c.z - a.z + 0.5f) * kCellSize }; };
                for (size_t k = 0; k < look.lamps.size(); ++k) {
                    const MbLamp& l = look.lamps[k];
                    if (k == 0) { Lamp(centreOf(l.cell), l.out, 0.16f, ControllerColor(state, t)); continue; }
                    // Lit while the structure works, each hatch a little out of step (a flow along the structure); dim when idle.
                    const float lit = state == MbLampState::Working ? 0.7f + 0.3f * sinf(t * 5.0f - 1.3f * (float)k) : 0.3f;
                    Lamp(centreOf(l.cell), l.out, 0.12f, Scale(KindColor(l.kind), lit));
                }
                if (look.formed) {
                    const Vector3 lo = { (look.lo.x - a.x) * kCellSize, (look.lo.y - a.y) * kCellSize, (look.lo.z - a.z) * kCellSize };
                    const Vector3 hi = { (look.hi.x + 1 - a.x) * kCellSize, (look.hi.y + 1 - a.y) * kCellSize, (look.hi.z + 1 - a.z) * kCellSize };
                    const float glow = state == MbLampState::Working ? 0.6f + 0.2f * sinf(t * 2.0f) : 0.35f; // subtle: formed at a glance
                    Seams(lo, hi, Scale({ 255, 175, 90, 255 }, glow));
                }
                rlPopMatrix();
                ++drawn;
            }
        }
    }
    rlDrawRenderBatchActive();
    rlEnableBackfaceCulling();
}
