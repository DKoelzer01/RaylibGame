#include "machinemodels.h"
#include "factory.h"
#include "instancebatch.h"
#include "plots.h"
#include "scene.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <unordered_map>

namespace {

// The machine's frame: x along the plot's east, y up, z along its north, origin at the centre of the
// footprint's floor; w, h, d its size (m); front = the part's facing in that frame (unit, horizontal).
struct Frame { float w, h, d; Vector3 front; float yawDeg; };

Color Shade(Color c, float light, float m = 1.0f) {
    return { (unsigned char)std::clamp(c.r * light * m, 0.0f, 255.0f), (unsigned char)std::clamp(c.g * light * m, 0.0f, 255.0f), (unsigned char)std::clamp(c.b * light * m, 0.0f, 255.0f), c.a };
}

struct Painter {
    float light;
    void cube(Vector3 c, Vector3 s, Color col, float m = 1.0f) const { IB_Cube(c, s, Shade(col, 1.0f, m)); } // lit by the batch shader
    void glow(Vector3 c, Vector3 s, Color col) const { IB_Cube(c, s, col, true); } // unlit (lit doors, hot coils)
    void cyl(Vector3 base, float rTop, float rBottom, float h, Color col, float m = 1.0f, int = 10) const { IB_Cylinder(base, rTop, rBottom, h, Shade(col, 1.0f, m)); }
};

const Color kSteel = { 140, 144, 150, 255 }, kDarkSteel = { 70, 72, 78, 255 }, kBrick = { 150, 80, 60, 255 }, kCopper = { 190, 110, 60, 255 };
const Color kGlowOrange = { 255, 150, 50, 255 }, kGlass = { 120, 200, 230, 255 }, kPaint = { 200, 170, 50, 255 };

void Chimney(const Painter& p, const Frame& f, float x, float z, float r, float h, Color c) {
    p.cyl({ x, f.h, z }, r, r * 1.15f, h, c);
    p.cyl({ x, f.h + h, z }, r * 1.2f, r * 1.2f, 0.06f, kDarkSteel);
}

// The door on the front face (front = +z after the yaw): dark when idle, glowing when working.
void Door(const Painter& p, const Frame& f, bool working, float wFrac = 0.45f, float hFrac = 0.45f) {
    rlPushMatrix();
    rlRotatef(f.yawDeg, 0, 1, 0);
    const float depth = (fabsf(f.front.x) > 0.5f ? f.w : f.d) * 0.5f;
    const Vector3 c = { 0, f.h * 0.35f, depth + 0.01f }, s = { f.w * wFrac, f.h * hFrac, 0.03f };
    p.cube({ c.x, c.y, c.z - 0.005f }, { s.x + 0.08f, s.y + 0.08f, 0.02f }, kDarkSteel);
    if (working) p.glow(c, s, Color{ (unsigned char)(220 + 35 * sinf((float)GetTime() * 7.0f)), 120, 40, 255 });
    else p.cube(c, s, { 30, 28, 26, 255 });
    rlPopMatrix();
}

void DrawMachine(PartType type, const Frame& f, const Painter& p, bool working) {
    const float t = (float)GetTime();
    const float spin = working ? t : 0.0f;
    const float x0 = -f.w * 0.5f, z0 = -f.d * 0.5f, top = f.h;
    switch (type) {
    case PartType::Kiln:
        Chimney(p, f, x0 + 0.2f, z0 + 0.2f, 0.1f, 0.6f, kBrick);
        Door(p, f, working);
        break;
    case PartType::CokeOven:
        Chimney(p, f, x0 + 0.2f, z0 + 0.2f, 0.09f, 0.7f, kBrick);
        Chimney(p, f, -x0 - 0.2f, z0 + 0.2f, 0.09f, 0.7f, kBrick);
        Door(p, f, working, 0.6f, 0.35f);
        break;
    case PartType::BlastFurnace:
        Chimney(p, f, 0, 0, 0.16f, 1.2f, kBrick);
        p.cyl({ 0, top + 0.5f, 0 }, 0.22f, 0.22f, 0.1f, kDarkSteel); // a charging ring
        Door(p, f, working, 0.35f, 0.3f);
        break;
    case PartType::Furnace:
        Door(p, f, working);
        for (int k = 0; k < 3; ++k) p.cube({ x0 + 0.25f + k * 0.25f, top + 0.02f, 0 }, { 0.12f, 0.04f, f.d * 0.6f }, kDarkSteel);
        break;
    case PartType::ArcFurnace:
        for (int k = 0; k < 3; ++k) {
            const float a = k * 2.094f;
            p.cyl({ 0.18f * cosf(a), top, 0.18f * sinf(a) }, 0.04f, 0.04f, 0.55f, { 60, 60, 66, 255 });
            if (working) p.glow({ 0.18f * cosf(a), top + 0.02f, 0.18f * sinf(a) }, { 0.1f, 0.04f, 0.1f }, { 200, 230, 255, 255 });
        }
        Door(p, f, working, 0.4f, 0.3f);
        break;
    case PartType::AlloySmelter:
        for (int k = 0; k < 2; ++k) {
            p.cyl({ (k ? 0.22f : -0.22f), top, 0 }, 0.17f, 0.12f, 0.2f, kBrick);
            if (working) p.glow({ (k ? 0.22f : -0.22f), top + 0.19f, 0 }, { 0.2f, 0.02f, 0.2f }, kGlowOrange);
        }
        break;
    case PartType::Boiler:
        rlPushMatrix(); rlTranslatef(0, top + 0.22f, 0); rlRotatef(90, 0, 0, 1);
        p.cyl({ 0, -f.w * 0.42f, 0 }, 0.2f, 0.2f, f.w * 0.84f, kCopper);    // a horizontal tank
        rlPopMatrix();
        Chimney(p, f, x0 + 0.15f, z0 + 0.15f, 0.07f, 0.8f, kDarkSteel);
        p.cyl({ -x0 - 0.2f, top + 0.42f, 0 }, 0.06f, 0.06f, 0.04f, { 230, 230, 230, 255 }); // gauge
        break;
    case PartType::ForgeHammer: {
        p.cube({ x0 + 0.08f, top + 0.4f, 0 }, { 0.08f, 0.8f, 0.12f }, kDarkSteel);
        p.cube({ -x0 - 0.08f, top + 0.4f, 0 }, { 0.08f, 0.8f, 0.12f }, kDarkSteel);
        p.cube({ 0, top + 0.82f, 0 }, { f.w, 0.08f, 0.14f }, kDarkSteel);
        const float lift = working ? 0.25f * (0.5f + 0.5f * sinf(t * 5.0f)) : 0.25f;
        p.cube({ 0, top + 0.2f + lift, 0 }, { 0.35f, 0.18f, 0.3f }, kSteel);   // the hammer head
        p.cube({ 0, top + 0.03f, 0 }, { 0.45f, 0.06f, 0.4f }, kDarkSteel);      // anvil
        break;
    }
    case PartType::Macerator:
    case PartType::SteamMacerator:
        p.cyl({ 0, top, 0 }, 0.38f, 0.18f, 0.35f, kSteel);                       // hopper
        p.cyl({ 0, top + 0.35f, 0 }, 0.4f, 0.4f, 0.04f, kDarkSteel);
        rlPushMatrix(); rlTranslatef(0, top + 0.12f, f.d * 0.5f + 0.02f); rlRotatef(spin * 400.0f, 0, 0, 1);
        p.cube({ 0, 0, 0 }, { 0.3f, 0.05f, 0.04f }, kDarkSteel); p.cube({ 0, 0, 0 }, { 0.05f, 0.3f, 0.04f }, kDarkSteel); // grinder wheel
        rlPopMatrix();
        if (type == PartType::SteamMacerator) Chimney(p, f, x0 + 0.15f, z0 + 0.15f, 0.06f, 0.5f, kDarkSteel);
        break;
    case PartType::Generator:
    case PartType::MvGenerator:
    case PartType::HvGenerator: {
        const Color coil = type == PartType::HvGenerator ? Color{ 200, 60, 60, 255 } : type == PartType::MvGenerator ? Color{ 200, 120, 50, 255 } : kCopper;
        rlPushMatrix(); rlTranslatef(x0 - 0.01f, top * 0.5f, 0); rlRotatef(90, 0, 0, 1);
        p.cyl({ 0, -0.02f, 0 }, 0.28f, 0.28f, 0.18f, coil);                     // the coil drum on the side
        rlPopMatrix();
        rlPushMatrix(); rlTranslatef(0, top + 0.05f, 0); rlRotatef(spin * 720.0f, 0, 1, 0);
        p.cube({ 0, 0, 0 }, { 0.6f, 0.03f, 0.1f }, kDarkSteel); p.cube({ 0, 0, 0 }, { 0.1f, 0.03f, 0.6f }, kDarkSteel); // fan
        rlPopMatrix();
        p.cyl({ 0, top, 0 }, 0.34f, 0.34f, 0.02f, kSteel);
        if (working) p.glow({ x0 - 0.12f, top * 0.5f, 0 }, { 0.02f, 0.3f, 0.3f }, { 120, 200, 255, 255 });
        break;
    }
    case PartType::Battery:
        p.cyl({ -0.2f, top, 0 }, 0.07f, 0.07f, 0.12f, { 200, 40, 40, 255 });
        p.cyl({ 0.2f, top, 0 }, 0.07f, 0.07f, 0.12f, { 30, 30, 30, 255 });
        p.cube({ 0, top * 0.5f, f.d * 0.5f + 0.01f }, { f.w * 0.7f, 0.08f, 0.02f }, kPaint);
        break;
    case PartType::OreWasher:
        rlPushMatrix(); rlTranslatef(0, top + 0.25f, 0); rlRotatef(90, 0, 0, 1); rlRotatef(spin * 180.0f, 0, 1, 0);
        p.cyl({ 0, -f.w * 0.4f, 0 }, 0.22f, 0.22f, f.w * 0.8f, kSteel);        // the washing drum, turning
        rlPopMatrix();
        p.cube({ 0, top + 0.03f, 0 }, { f.w * 0.9f, 0.06f, 0.5f }, kDarkSteel);
        break;
    case PartType::Centrifuge:
        rlPushMatrix(); rlTranslatef(0, top, 0); rlRotatef(spin * 900.0f, 0, 1, 0);
        p.cyl({ 0, 0, 0 }, 0.4f, 0.4f, 0.3f, kSteel);
        p.cube({ 0.3f, 0.15f, 0 }, { 0.14f, 0.12f, 0.14f }, kPaint);           // a marker to see it spin
        rlPopMatrix();
        break;
    case PartType::WaterPump:
        p.cyl({ x0 + 0.2f, -0.6f, z0 + 0.2f }, 0.08f, 0.08f, 0.6f + top + 0.3f, kSteel); // the intake pipe
        p.cyl({ 0.1f, top, 0.1f }, 0.28f, 0.28f, 0.45f, { 70, 120, 190, 255 }); // tank
        break;
    case PartType::Electrolyzer:
        p.cyl({ -0.22f, top, 0 }, 0.16f, 0.16f, 0.5f, { 230, 230, 230, 255 }); // oxygen
        p.cyl({ 0.22f, top, 0 }, 0.16f, 0.16f, 0.5f, { 120, 160, 230, 255 });  // hydrogen
        if (working) p.glow({ 0, top + 0.05f, 0 }, { 0.1f, 0.1f, 0.1f }, { 150, 220, 255, 255 });
        break;
    case PartType::ChemicalReactor:
        p.cyl({ 0, top, 0 }, 0.3f, 0.36f, 0.45f, { 120, 170, 120, 255 });
        p.cyl({ 0, top + 0.45f, 0 }, 0.12f, 0.3f, 0.15f, { 120, 170, 120, 255 }, 0.9f);
        p.cyl({ x0 + 0.12f, top, 0 }, 0.05f, 0.05f, 0.6f, kSteel);
        break;
    case PartType::Refinery:
        p.cyl({ 0, top, 0 }, 0.22f, 0.26f, 1.6f, kSteel);                      // the distillation column
        for (int k = 1; k <= 4; ++k) p.cyl({ 0, top + 0.35f * k, 0 }, 0.29f, 0.29f, 0.04f, kDarkSteel);
        p.cyl({ 0, top + 1.6f, 0 }, 0.05f, 0.05f, 0.25f, kDarkSteel);
        break;
    case PartType::Miner: {
        // A derrick and a drill bit under the casing's front edge, turning while it works.
        for (int k = 0; k < 4; ++k) p.cube({ (k & 1 ? 0.35f : -0.35f), top + 0.45f, (k & 2 ? 0.35f : -0.35f) }, { 0.05f, 0.9f, 0.05f }, kPaint);
        p.cube({ 0, top + 0.9f, 0 }, { 0.8f, 0.06f, 0.8f }, kPaint);
        rlPushMatrix(); rlTranslatef(0, top + 0.1f, 0); rlRotatef(spin * 300.0f, 0, 1, 0);
        p.cyl({ 0, 0, 0 }, 0.12f, 0.02f, 0.7f, kDarkSteel);
        rlPopMatrix();
        break;
    }
    case PartType::PrecisionAssembler: {
        p.cube({ 0, top + 0.02f, 0 }, { f.w * 0.8f, 0.04f, f.d * 0.8f }, kDarkSteel);
        const float a = working ? 40.0f * sinf(t * 1.6f) : 0.0f;
        rlPushMatrix(); rlTranslatef(x0 + 0.2f, top, 0); rlRotatef(a, 0, 1, 0);
        p.cube({ 0, 0.25f, 0 }, { 0.08f, 0.5f, 0.08f }, kPaint);
        rlTranslatef(0, 0.5f, 0); rlRotatef(-60.0f + (working ? 20.0f * sinf(t * 2.3f) : 0.0f), 0, 0, 1);
        p.cube({ 0.22f, 0, 0 }, { 0.44f, 0.06f, 0.06f }, kPaint);               // the arm
        rlPopMatrix();
        rlPushMatrix(); rlRotatef(f.yawDeg, 0, 1, 0);
        p.glow({ 0, top * 0.6f, (fabsf(f.front.x) > 0.5f ? f.w : f.d) * 0.5f + 0.01f }, { f.w * 0.4f, f.h * 0.25f, 0.02f }, working ? Color{ 80, 230, 140, 255 } : Color{ 20, 50, 40, 255 }); // screen
        rlPopMatrix();
        break;
    }
    case PartType::FusionReactor: {
        const float pulse = working ? 0.6f + 0.4f * sinf(t * 3.0f) : 0.25f;
        for (int k = 0; k < 16; ++k) {
            const float a = k * 0.3927f;
            const Vector3 c = { 0.42f * cosf(a), top + 0.18f, 0.42f * sinf(a) };
            p.cube(c, { 0.16f, 0.22f, 0.16f }, kSteel);
            p.glow({ c.x * 0.8f, c.y, c.z * 0.8f }, { 0.06f, 0.06f, 0.06f }, Color{ (unsigned char)(200 * pulse), (unsigned char)(120 * pulse), (unsigned char)(255 * pulse), 255 });
        }
        break;
    }
    case PartType::MassDriver:
        rlPushMatrix(); rlTranslatef(0, top + 0.1f, 0); rlRotatef(f.yawDeg, 0, 1, 0); rlRotatef(-35.0f, 1, 0, 0);
        p.cube({ 0, 0, 0.8f }, { 0.22f, 0.22f, 1.8f }, kSteel);                // the launch rail, raised
        for (int k = 0; k < 5; ++k) p.cube({ 0, 0, 0.1f + 0.35f * k }, { 0.3f, 0.3f, 0.06f }, working ? kGlass : kCopper);
        rlPopMatrix();
        break;
    case PartType::Chest:
        p.cube({ 0, top - 0.06f, 0 }, { f.w * 1.02f, 0.06f, f.d * 1.02f }, kDarkSteel);  // a banded lid
        rlPushMatrix(); rlRotatef(f.yawDeg, 0, 1, 0);
        p.cube({ 0, top * 0.7f, (fabsf(f.front.x) > 0.5f ? f.w : f.d) * 0.5f + 0.02f }, { 0.12f, 0.14f, 0.04f }, kPaint); // the latch
        rlPopMatrix();
        break;
    default: break;
    }
}

} // namespace

// Per plot, rebuilt when its version changes: each machine's floor centre (plot metres), size and facing.
namespace {
    struct MachineSpot { int part; PartType type; Vector3 floorPlot; Frame f; };
    struct PlotCache { uint32_t version = 0xFFFFFFFFu; std::vector<MachineSpot> spots; };
    std::unordered_map<uint64_t, PlotCache> plotCache; // key: body id ^ plot id
}

static int drawnLastFrame = 0;
int Machines_DrawnCount() { return drawnLastFrame; }

void Machines_Draw3D(Scene& world) {
    rlDrawRenderBatchActive();
    rlDisableBackfaceCulling(); // the plot frame (east, up, north) is left-handed: mirrored winding
    drawnLastFrame = 0;
    const Vector3 sun = Vector3Normalize(world.sunDirection);
    std::vector<std::pair<Int3, CellShape>> fp;
    struct Pending { float d2; Planetoid* body; Plot* plot; const MachineSpot* spot; Vector3 e, u, n; float light; };
    static std::vector<Pending> pending;
    pending.clear();
    for (const auto& obj : world.objects) {
        Planetoid* body = dynamic_cast<Planetoid*>(obj.get());
        if (!body || !body->isActive || !body->plots) continue;
        const Vector3 camPlotBody = body->toLocal(world.camera.position);
        for (auto& plot : body->plots->plots) {
            if (Vector3Distance(plot->toBody({ 0, 0, 0 }), camPlotBody) > 500.0f) continue;
            PlotCache& pc = plotCache[(body->id * 1315423911ull) ^ plot->id];
            if (pc.version != plot->version) {
                pc.version = plot->version;
                pc.spots.clear();
                for (size_t i = 0; i < plot->parts.size(); ++i) {
                    const Part& part = plot->parts[i];
                    if (part.type == PartType::Count || !IsMachine(part.type)) continue;
                    Plot::Footprint(part, fp);
                    if (fp.empty()) continue;
                    Int3 lo = fp[0].first, hi = fp[0].first;
                    for (auto& c : fp) { lo = { std::min(lo.x, c.first.x), std::min(lo.y, c.first.y), std::min(lo.z, c.first.z) }; hi = { std::max(hi.x, c.first.x), std::max(hi.y, c.first.y), std::max(hi.z, c.first.z) }; }
                    MachineSpot s;
                    s.part = (int)i; s.type = part.type;
                    s.floorPlot = { (lo.x + hi.x + 1) * 0.5f * kCellSize, lo.y * kCellSize, (lo.z + hi.z + 1) * 0.5f * kCellSize };
                    s.f.w = (hi.x - lo.x + 1) * kCellSize; s.f.h = (hi.y - lo.y + 1) * kCellSize; s.f.d = (hi.z - lo.z + 1) * kCellSize;
                    const int* m = Orientation(part.orient);
                    Vector3 front = { (float)m[2], 0.0f, (float)m[8] };              // the part's +z, flattened
                    if (Vector3LengthSqr(front) < 0.5f) front = { 0, 0, 1 };
                    s.f.front = Vector3Normalize(front);
                    s.f.yawDeg = atan2f(s.f.front.x, s.f.front.z) * RAD2DEG;
                    pc.spots.push_back(s);
                }
            }
            if (pc.spots.empty()) continue;
            const Vector3 camPlot = plot->toPlot(camPlotBody);                     // camera in plot metres: cheap culling
            const Vector3 e = body->dirToWorld(plot->east), u = body->dirToWorld(plot->up), n = body->dirToWorld(plot->north);
            Painter p;
            p.light = 0.45f + 0.55f * std::clamp(Vector3DotProduct(sun, u) * 2.0f + 0.5f, 0.0f, 1.0f);
            for (const MachineSpot& s : pc.spots) {
                const float d2 = Vector3DistanceSqr(s.floorPlot, camPlot);
                if (d2 > 32.0f * 32.0f) continue; // the casings read on their own beyond this
                Pending q; q.d2 = d2; q.body = body; q.plot = plot.get(); q.spot = &s; q.e = e; q.u = u; q.n = n; q.light = p.light;
                pending.push_back(q);
            }
        }
    }
    // The nearest kMaxDetailed only: a huge factory in view costs a bounded amount.
    constexpr size_t kMaxDetailed = 120;
    if (pending.size() > kMaxDetailed) {
        std::nth_element(pending.begin(), pending.begin() + kMaxDetailed, pending.end(), [](const Pending& a, const Pending& b) { return a.d2 < b.d2; });
        pending.resize(kMaxDetailed);
    }
    for (const Pending& q : pending) {
        const MachineSpot& s = *q.spot;
        const Vector3 base = q.body->toWorld(q.plot->toBody(s.floorPlot));
        const bool working = Factory::Working(q.body->id, q.plot->id, s.part);
        const Matrix mm = { q.e.x, q.u.x, q.n.x, base.x, q.e.y, q.u.y, q.n.y, base.y, q.e.z, q.u.z, q.n.z, base.z, 0, 0, 0, 1 };
        Painter pp; pp.light = q.light;
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(mm));
        DrawMachine(s.type, s.f, pp, working);
        ++drawnLastFrame;
        rlPopMatrix();
    }
    rlDrawRenderBatchActive();
    rlEnableBackfaceCulling();
}
