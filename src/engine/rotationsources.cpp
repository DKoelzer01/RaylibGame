#include "rotationsources.h"
#include "factory.h"
#include "fluid.h"
#include "plots.h"
#include "rotationplot.h"
#include "scene.h"
#include "simclock.h"
#include "voxel.h"
#include "wind.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

// Water wheels and windmills (#1270). See rotationsources.h.

namespace RotationSources {

int WaterWheelRpm(Vector3 v, Vector3 axis, Vector3 up) {
    // The rim's bottom point sits at r = -R up from the hub; turning at w about the axis it moves at w x r = -w R (axis x up).
    // Carried by the water there: w = -(v . (axis x up)) / R, less the slip.
    const Vector3 t = Vector3CrossProduct(axis, up);
    const float omega = -Vector3DotProduct(v, t) / kWheelRadius * kWheelSlip; // rad/s
    const float rpm = omega * 60.0f / (2.0f * PI);
    const int cap = Rotation::kWaterWheel.rpm;
    return std::clamp((int)lroundf(rpm), -cap, cap);
}

int WindmillRpm(float s) {
    if (!(s >= kWindCutInMps)) return 0;
    const float f = std::min(1.0f, s / kWindFullMps);
    return std::max(1, (int)lroundf(f * Rotation::kWindmill.rpm));
}

int LowestRpm(const std::deque<int>& w) {
    if (w.empty()) return 0;
    bool pos = false, neg = false;
    int lowest = w.front();
    for (int r : w) { pos |= r > 0; neg |= r < 0; if (std::abs(r) < std::abs(lowest)) lowest = r; }
    return pos && neg ? 0 : lowest;
}

} // namespace RotationSources

namespace {
using namespace RotationSources;

struct Source {
    int part = -1;
    PartType type = PartType::Count;
    Vector3 at{};        // body-local: the bottom of a wheel's rim / a windmill's hub
    Vector3 axis{};      // body-local: the plot axis it turns about (+)
    Vector3 up{};        // body-local: the plot's up
    int sent = 0x7fffffff;
    std::deque<int> window; // recent samples (water wheels)
};
struct PlotSources { uint32_t version = 0xFFFFFFFFu; std::vector<Source> sources; };
std::unordered_map<uint64_t, PlotSources> cache; // by body id ^ plot id
double lastSample = -1e9;
bool windOverride = false;
float windOverrideSpeed = 0.0f;

void Rebuild(const Plot& plot, PlotSources& ps) {
    std::unordered_map<int, Source> old;
    for (Source& s : ps.sources) old[s.part] = std::move(s);
    ps.version = plot.version;
    ps.sources.clear();
    for (size_t i = 0; i < plot.parts.size(); ++i) {
        const Part& p = plot.parts[i];
        if (p.type != PartType::WaterWheel && p.type != PartType::Windmill) continue;
        RotationPartMsg rm;
        if (!RotationPlot_Describe(p, (int)i, rm) || rm.geom.cells.empty()) continue;
        Source s;
        if (auto it = old.find((int)i); it != old.end() && it->second.type == p.type) s = std::move(it->second);
        s.part = (int)i; s.type = p.type;
        const Rotation::Cell hub = p.type == PartType::WaterWheel && rm.geom.cells.size() > 4 ? rm.geom.cells[4] : rm.geom.cells[0];
        Vector3 c = { (hub.x + 0.5f) * kCellSize, (hub.y + 0.5f) * kCellSize, (hub.z + 0.5f) * kCellSize };
        if (p.type == PartType::WaterWheel) c.y -= kWheelRadius * 0.8f; // just inside the rim's bottom, where it dips
        s.at = plot.toBody(c);
        s.axis = rm.geom.axis == 0 ? plot.east : rm.geom.axis == 1 ? plot.up : plot.north;
        s.up = plot.up;
        ps.sources.push_back(std::move(s));
    }
}

int Sample(Planetoid& body, Source& s, bool near) {
    if (s.type == PartType::Windmill) {
        float speed = windOverrideSpeed;
        if (!windOverride) {
            float height = 10.0f;
            if (const auto* pf = dynamic_cast<const PlanetField*>(body.field.get())) {
                const float len = Vector3Length(s.at);
                if (len > 1.0f) height = std::max(1.0f, len - (pf->radius + pf->elevation(Vector3Scale(s.at, 1.0f / len))));
            }
            speed = Vector3Length(Wind_At(body, s.at, height, true));
        }
        return WindmillRpm(speed);
    }
    // A water wheel: simulated water where a simulated cell decides, else the generated river's current.
    const bool simulated = Fluid::Simulated(body, s.at);
    Vector3 v{};
    if (simulated) v = Vector3Scale(Fluid::FlowAt(body, s.at, 1.0f), 1.0f / kWheelDipArea);
    else v = body.dirToLocal(RiverCurrent(&body, body.toWorld(s.at)));
    const int rpm = WaterWheelRpm(v, s.axis, s.up);
    s.window.push_back(rpm);
    while ((int)s.window.size() > kFarWindow) s.window.pop_front();
    return near || !simulated ? rpm : LowestRpm(s.window); // (2.9) far on simulated water: the lowest of the last minute
}
} // namespace

void RotationSources_Reset() { cache.clear(); lastSample = -1e9; }
void RotationSources_SetWindOverride(bool on, float speed) { windOverride = on; windOverrideSpeed = speed; }

void RotationSources_Update(Scene& world) {
    const double now = SimTime();
    if (now < lastSample) lastSample = -1e9; // the clock went back (a load)
    if (now - lastSample < kSampleSeconds) return;
    lastSample = now;
    for (auto& obj : world.objects) {
        auto* body = dynamic_cast<Planetoid*>(obj.get());
        if (!body || !body->plots) continue;
        const Vector3 camBody = body->toLocal(world.camera.position);
        for (auto& pl : body->plots->plots) {
            PlotSources& ps = cache[(body->id * 1315423911ull) ^ pl->id];
            if (ps.version != pl->version) Rebuild(*pl, ps);
            if (ps.sources.empty()) continue;
            const bool near = Vector3Distance(pl->toBody({ 0, 0, 0 }), camBody) < kNearMetres;
            for (Source& s : ps.sources) {
                const int rpm = Sample(*body, s, near);
                if (rpm != s.sent) { Factory::SetSourceRpm(body->id, pl->id, s.part, rpm); s.sent = rpm; }
            }
        }
    }
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1270): the wheel and windmill speed rules ----
#include <gtest/gtest.h>

TEST(RotationSources, AWheelTurnsWithTheWaterAcrossItsAxis) {
    const Vector3 up = { 0, 1, 0 }, axisZ = { 0, 0, 1 };
    // 1 m/s along +x under a wheel on the z axis: w = -(v . (z x y)) / R x 0.8 = -(1 . -1) / 0.75 x 0.8 = 1.07 rad/s = 10.2 RPM.
    EXPECT_EQ(RotationSources::WaterWheelRpm({ 1, 0, 0 }, axisZ, up), 10);
    EXPECT_EQ(RotationSources::WaterWheelRpm({ -1, 0, 0 }, axisZ, up), -10) << "the other way";
    EXPECT_EQ(RotationSources::WaterWheelRpm({ 0, 0, 1 }, axisZ, up), 0) << "water along the axis turns nothing";
    EXPECT_EQ(RotationSources::WaterWheelRpm({ 1, 0, 0 }, up, up), 0) << "a wheel lying flat";
    EXPECT_EQ(RotationSources::WaterWheelRpm({ 5, 0, 0 }, axisZ, up), Rotation::kWaterWheel.rpm) << "a flood: capped";
    EXPECT_EQ(RotationSources::WaterWheelRpm({ 0, 0, 0 }, axisZ, up), 0) << "still water";
}

TEST(RotationSources, AWindmillNeedsTheCutInAndTopsOutAtFullWind) {
    EXPECT_EQ(RotationSources::WindmillRpm(0.0f), 0);
    EXPECT_EQ(RotationSources::WindmillRpm(1.4f), 0);
    EXPECT_EQ(RotationSources::WindmillRpm(4.0f), 8);
    EXPECT_EQ(RotationSources::WindmillRpm(8.0f), Rotation::kWindmill.rpm);
    EXPECT_EQ(RotationSources::WindmillRpm(30.0f), Rotation::kWindmill.rpm);
}

TEST(RotationSources, AFarWheelUsesTheLowestOfTheLastMinute) {
    EXPECT_EQ(RotationSources::LowestRpm({}), 0);
    EXPECT_EQ(RotationSources::LowestRpm({ 12, 9, 14 }), 9);
    EXPECT_EQ(RotationSources::LowestRpm({ -12, -9, -14 }), -9);
    EXPECT_EQ(RotationSources::LowestRpm({ 12, -3, 14 }), 0) << "it turned both ways: nothing counted on";
}
#endif
