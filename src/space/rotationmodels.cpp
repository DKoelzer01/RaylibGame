#include "rotationmodels.h"
#include "factory.h"
#include "instancebatch.h"
#include "plots.h"
#include "rotationplot.h"
#include "scene.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

// Rotation models (#1269). See rotationmodels.h.

namespace {

const Color kWood = { 150, 112, 70, 255 }, kDarkWood = { 84, 62, 42, 255 }, kStone = { 178, 175, 168, 255 }, kLight = { 214, 182, 120, 255 }, kSteel = { 150, 154, 160, 255 };
const Color kSail = { 226, 218, 196, 255 };
constexpr int kWindmillSailsDrawn = Rotation::kWindmillSails;

struct Spot {
    int part;
    PartType type;
    Vector3 centre;   // plot metres: the cell centre (transmissions), or over the machine's top
    int axis;         // plot axis it turns about
};
struct PlotSpots { uint32_t version = 0xFFFFFFFFu; std::vector<Spot> spots; };
std::unordered_map<uint64_t, PlotSpots> cache;          // by body id ^ plot id
std::unordered_map<uint64_t, float> angles;             // by body, plot, part: degrees turned so far
int drawn = 0;

uint64_t PlotKey(uint64_t body, uint32_t plot) { return (body * 1315423911ull) ^ plot; }
uint64_t PartKey(uint64_t body, uint32_t plot, int part) { return (PlotKey(body, plot) * 2654435761ull) ^ (uint64_t)(uint32_t)part; }

void Rebuild(const Plot& plot, PlotSpots& ps) {
    ps.version = plot.version;
    ps.spots.clear();
    for (size_t i = 0; i < plot.parts.size(); ++i) {
        const Part& p = plot.parts[i];
        if (p.type == PartType::Count || !(IsSpinningPart(p.type) || (IsRotationPart(p.type) && IsMachine(p.type)))) continue;
        RotationPartMsg rm;
        if (!RotationPlot_Describe(p, (int)i, rm) || rm.geom.cells.empty()) continue;
        Vector3 lo = { 1e9f, 1e9f, 1e9f }, hi = { -1e9f, -1e9f, -1e9f };
        for (const Rotation::Cell& c : rm.geom.cells) {
            lo = { std::min(lo.x, (float)c.x), std::min(lo.y, (float)c.y), std::min(lo.z, (float)c.z) };
            hi = { std::max(hi.x, (float)c.x + 1), std::max(hi.y, (float)c.y + 1), std::max(hi.z, (float)c.z + 1) };
        }
        Spot s;
        s.part = (int)i; s.type = p.type; s.axis = rm.geom.axis;
        s.centre = { (lo.x + hi.x) * 0.5f * kCellSize, (lo.y + hi.y) * 0.5f * kCellSize, (lo.z + hi.z) * 0.5f * kCellSize };
        if (IsMachine(p.type)) s.centre.y = hi.y * kCellSize + (p.type == PartType::Millstone ? 0.0f : p.type == PartType::Saw ? 0.12f : 0.32f); // on its top
        ps.spots.push_back(s);
    }
}

// A ring of teeth round the local y axis.
void Teeth(int n, float r, float h, float size, Color c) {
    for (int k = 0; k < n; ++k) {
        rlPushMatrix();
        rlRotatef(360.0f * k / n, 0, 1, 0);
        IB_Cube({ r, 0, 0 }, { size, h, size }, c);
        rlPopMatrix();
    }
}

// The part in its own frame: y = its axis, origin at its centre, already turned by its angle.
void DrawPart(PartType t) {
    switch (t) {
    case PartType::Shaft:
        IB_Cylinder({ 0, -0.25f, 0 }, 0.06f, 0.06f, 0.5f, kWood);
        IB_Cube({ 0.055f, 0, 0 }, { 0.03f, 0.48f, 0.03f }, kLight); // a light stripe: the turning shows
        break;
    case PartType::Cogwheel:
        IB_Cylinder({ 0, -0.25f, 0 }, 0.06f, 0.06f, 0.5f, kWood);
        IB_Cylinder({ 0, -0.06f, 0 }, 0.19f, 0.19f, 0.12f, kWood);
        Teeth(8, 0.21f, 0.1f, 0.06f, kDarkWood);
        break;
    case PartType::LargeCogwheel:
        IB_Cylinder({ 0, -0.25f, 0 }, 0.06f, 0.06f, 0.5f, kWood);
        IB_Cylinder({ 0, -0.06f, 0 }, 0.44f, 0.44f, 0.12f, kWood);
        IB_Cube({ 0, 0, 0 }, { 0.8f, 0.13f, 0.08f }, kDarkWood); // spokes
        IB_Cube({ 0, 0, 0 }, { 0.08f, 0.13f, 0.8f }, kDarkWood);
        Teeth(16, 0.47f, 0.1f, 0.06f, kDarkWood);
        break;
    case PartType::HandCrank:
        IB_Cylinder({ 0, -0.25f, 0 }, 0.06f, 0.06f, 0.4f, kWood);
        IB_Cube({ 0.12f, 0.17f, 0 }, { 0.26f, 0.05f, 0.06f }, kDarkWood);  // the arm
        IB_Cylinder({ 0.23f, 0.17f, 0 }, 0.025f, 0.025f, 0.16f, kWood);    // the handle
        break;
    case PartType::Millstone: // the runner stone on top, with a notch and the grain hole
        IB_Cylinder({ 0, 0, 0 }, 0.42f, 0.44f, 0.14f, kStone);
        IB_Cube({ 0, 0.15f, 0 }, { 0.8f, 0.04f, 0.09f }, kLight); // the crossbar (rynd): the turning shows
        IB_Cylinder({ 0, 0.14f, 0 }, 0.07f, 0.07f, 0.01f, kDarkWood);
        break;
    case PartType::Saw: // the blade standing out of the table
        IB_Cylinder({ 0, -0.012f, 0 }, 0.3f, 0.3f, 0.024f, kSteel);
        Teeth(18, 0.31f, 0.02f, 0.035f, kSteel);
        IB_Cylinder({ 0, -0.04f, 0 }, 0.05f, 0.05f, 0.08f, kDarkWood);
        break;
    case PartType::Press: // the flywheel over the press
        IB_Cylinder({ 0, -0.04f, 0 }, 0.28f, 0.28f, 0.08f, kSteel);
        IB_Cube({ 0, 0, 0 }, { 0.5f, 0.1f, 0.06f }, kDarkWood);
        IB_Cylinder({ 0, -0.12f, 0 }, 0.04f, 0.04f, 0.24f, kWood);
        break;
    case PartType::WaterWheel: // (#1270) two rims on spokes, paddles between them (1.5 m across)
        IB_Cylinder({ 0, -0.25f, 0 }, 0.07f, 0.07f, 0.5f, kWood);
        for (float side : { -0.2f, 0.2f }) {
            for (int k = 0; k < 4; ++k) { rlPushMatrix(); rlRotatef(45.0f * k, 0, 1, 0); IB_Cube({ 0, side, 0 }, { 1.4f, 0.04f, 0.06f }, kDarkWood); rlPopMatrix(); }
            rlPushMatrix(); rlTranslatef(0, side, 0); Teeth(16, 0.68f, 0.04f, 0.14f, kWood); rlPopMatrix(); // the rim, as short boards round it
        }
        for (int k = 0; k < 8; ++k) { rlPushMatrix(); rlRotatef(45.0f * k, 0, 1, 0); IB_Cube({ 0.66f, 0, 0 }, { 0.16f, 0.44f, 0.05f }, kLight); rlPopMatrix(); }
        break;
    case PartType::Windmill: // (#1270) a hub and four sails, 2.5 m long
        IB_Cylinder({ 0, -0.25f, 0 }, 0.08f, 0.08f, 0.5f, kWood);
        IB_Cylinder({ 0, 0.2f, 0 }, 0.14f, 0.12f, 0.12f, kDarkWood);
        for (int k = 0; k < kWindmillSailsDrawn; ++k) {
            rlPushMatrix(); rlRotatef(90.0f * k, 0, 1, 0);
            IB_Cube({ 1.3f, 0.26f, 0 }, { 2.5f, 0.05f, 0.06f }, kDarkWood);      // the stock
            IB_Cube({ 1.5f, 0.26f, 0.2f }, { 2.0f, 0.02f, 0.36f }, kSail);       // the cloth
            rlPopMatrix();
        }
        break;
    default: break;
    }
}

} // namespace

int RotationModels_DrawnCount() { return drawn; }

void RotationModels_Draw3D(Scene& world) {
    drawn = 0;
    const float dt = std::min(0.1f, GameFrameTime());
    static std::vector<std::pair<int, float>> speeds;
    static std::unordered_map<int, float> rpmOf;
    rlDrawRenderBatchActive();
    rlDisableBackfaceCulling(); // the plot frame (east, up, north) is left-handed: mirrored winding
    for (const auto& obj : world.objects) {
        Planetoid* body = dynamic_cast<Planetoid*>(obj.get());
        if (!body || !body->isActive || !body->plots) continue;
        const Vector3 camBody = body->toLocal(world.camera.position);
        for (auto& plot : body->plots->plots) {
            if (Vector3Distance(plot->toBody({ 0, 0, 0 }), camBody) > 300.0f) continue;
            PlotSpots& ps = cache[PlotKey(body->id, plot->id)];
            if (ps.version != plot->version) Rebuild(*plot, ps);
            if (ps.spots.empty()) continue;
            const Vector3 camPlot = plot->toPlot(camBody);
            Factory::RotationSpeeds(body->id, plot->id, speeds);
            rpmOf.clear();
            for (const auto& [part, rpm] : speeds) rpmOf[part] = rpm;
            const Plot& pl = *plot;
            Matrix axes = MatrixIdentity();
            axes.m0 = pl.east.x; axes.m4 = pl.up.x; axes.m8 = pl.north.x;
            axes.m1 = pl.east.y; axes.m5 = pl.up.y; axes.m9 = pl.north.y;
            axes.m2 = pl.east.z; axes.m6 = pl.up.z; axes.m10 = pl.north.z;
            for (const Spot& s : ps.spots) {
                float& angle = angles[PartKey(body->id, plot->id, s.part)];
                auto it = rpmOf.find(s.part);
                angle = fmodf(angle + (it != rpmOf.end() ? it->second : 0.0f) * 6.0f * dt, 360.0f); // RPM x 360 / 60 degrees a second
                if (Vector3DistanceSqr(s.centre, camPlot) > 48.0f * 48.0f) continue;
                // The plot's model with its translation in double (as the machine models: no shaking against the casings, #284).
                const Vector3 off = Vector3Add(Vector3Add(Vector3Scale(pl.east, s.centre.x), Vector3Scale(pl.up, s.centre.y)), Vector3Scale(pl.north, s.centre.z));
                const Matrix mm = MatrixMultiply(axes, body->localModel(DVec3(pl.origin.x + off.x, pl.origin.y + off.y, pl.origin.z + off.z), world.renderOrigin));
                rlPushMatrix();
                rlMultMatrixf(MatrixToFloat(mm));
                if (s.axis == 0) rlRotatef(-90.0f, 0, 0, 1);  // its y onto plot x
                else if (s.axis == 2) rlRotatef(90.0f, 1, 0, 0); // its y onto plot z
                rlRotatef(angle, 0, 1, 0);
                DrawPart(s.type);
                rlPopMatrix();
                ++drawn;
            }
        }
    }
    rlDrawRenderBatchActive();
    rlEnableBackfaceCulling();
}
