#include "rocket.h"
#include "buildtool.h"
#include "gameui.h"
#include "inventory.h"
#include "audio.h"
#include "player.h"
#include "plots.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <sstream>
#include <fstream>

namespace {
    constexpr double kG0 = 9.80665;
    constexpr double kAirScaleHeight = 4000.0, kAirDensity = 1.2, kAirTop = 40000.0;
    constexpr double kMaxStep = 1.0 / 120.0;
    constexpr float kTurnAccel = 1.2f, kTurnRate = 0.7f;   // rad/s^2, rad/s (per guidance computer, capped)

    struct PartStat { float mass = 0, thrust = 0, isp = 0, fuel = 0; };
    std::map<std::string, PartStat> partStats;
    std::map<std::string, float> materialMass;  // kg per cell
    float fuelPerItem = 100.0f;
    bool loaded = false;

    // Flight state.
    RocketState state = RocketState::None;
    uint64_t padBody = 0;        // OnPad: the body and plot
    uint32_t padPlot = 0;
    uint64_t vehicleId = 0;      // Flying: the vehicle body (it carries the plot)
    DVec3 com{};                 // world, centre of mass
    Vector3 vel{};               // world, m/s
    Quaternion rot = { 0, 0, 0, 1 }; // plot-local -> world
    Vector3 angVel{};            // world, rad/s
    Vector3 comLocal{};          // plot space (m), centre of mass
    float throttle = 0.0f;
    double flightTime = 0.0;
    RocketAutopilot autopilot = RocketAutopilot::Off;
    double autoTarget = 45000.0;
    int landings = 0, damage = 0;
    float lastTouchdown = 0.0f;
    float camYaw = 0.6f, camPitch = 0.25f;
    // On-foot aim.
    uint64_t aimBody = 0;
    uint32_t aimPlot = 0;
    bool aimCockpit = false;

    void Load() {
        if (loaded) return;
        loaded = true;
        std::ifstream in("resources/data/ship_parts.txt");
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::stringstream ss(line);
            std::string first;
            ss >> first;
            if (first == "material") { std::string m; float kg; if (ss >> m >> kg) materialMass[m] = kg; continue; }
            if (first == "fuel") { std::string id; float kg; if (ss >> id >> kg) fuelPerItem = kg; continue; }
            std::vector<std::string> cols;
            size_t start = 0, bar;
            while ((bar = line.find('|', start)) != std::string::npos) { cols.push_back(line.substr(start, bar - start)); start = bar + 1; }
            cols.push_back(line.substr(start));
            if (cols.size() < 5) continue;
            for (auto& c : cols) { const size_t a = c.find_first_not_of(" \t"), b = c.find_last_not_of(" \t\r"); c = a == std::string::npos ? "" : c.substr(a, b - a + 1); }
            PartStat s;
            s.mass = (float)atof(cols[1].c_str());
            s.thrust = (float)atof(cols[2].c_str()) * 1000.0f;
            s.isp = (float)atof(cols[3].c_str());
            s.fuel = (float)atof(cols[4].c_str());
            partStats[cols[0]] = s;
        }
    }

    Planetoid* PadBody() { return Planetoid::find(padBody); }
    Plot* PadPlot() { Planetoid* b = PadBody(); return b && b->plots ? b->plots->byId(padPlot) : nullptr; }
    Planetoid* Vehicle() { return Planetoid::find(vehicleId); }
    Plot* VehiclePlot() { Planetoid* v = Vehicle(); return v && v->plots && !v->plots->plots.empty() ? v->plots->plots[0].get() : nullptr; }
    Plot* CurrentPlot() { return state == RocketState::Flying ? VehiclePlot() : state == RocketState::OnPad ? PadPlot() : nullptr; }

    // Gravity at a point from the bodies (not the rocket).
    Vector3 Gravity(Scene& world, const DVec3& p) { return GravityAt(world, (p - world.renderOrigin).toFloat()); }

    // The body whose pull dominates here (orbits are measured against it).
    Planetoid* Dominant(Scene& world, const DVec3& p) {
        Planetoid* best = nullptr;
        double bestG = -1.0;
        for (auto& obj : world.objects) {
            auto* b = dynamic_cast<Planetoid*>(obj.get());
            if (!b || b->vehicle || b->flat) continue;
            const double d = (p - b->worldPos).length();
            const double R = b->gravityRadius > 0.0f ? b->gravityRadius : b->radius * b->scale;
            const double g = b->surfaceGravity * R * R / std::max(d * d, R * R);
            if (g > bestG) { bestG = g; best = b; }
        }
        return best;
    }
    double Mu(const Planetoid& b) { const double R = b.gravityRadius > 0.0f ? b.gravityRadius : b.radius * b.scale; return b.surfaceGravity * R * R; }

    // Orbit about a body: apoapsis / periapsis radii from the state (ap < 0: hyperbolic).
    void Orbit(const Planetoid& b, const DVec3& p, Vector3 v, double& ra, double& rp) {
        const DVec3 r = p - b.worldPos;
        const DVec3 vv = DVec3(Vector3Subtract(v, b.linearVelocity));
        const double mu = Mu(b), rl = r.length(), v2 = vv.x * vv.x + vv.y * vv.y + vv.z * vv.z;
        const double energy = 0.5 * v2 - mu / rl;
        const DVec3 h(r.y * vv.z - r.z * vv.y, r.z * vv.x - r.x * vv.z, r.x * vv.y - r.y * vv.x);
        const double h2 = h.x * h.x + h.y * h.y + h.z * h.z;
        const double e = sqrt(std::max(0.0, 1.0 + 2.0 * energy * h2 / (mu * mu)));
        if (energy >= 0.0) { ra = -1.0; rp = h2 / (mu * (1.0 + e)); return; }
        const double a = -mu / (2.0 * energy);
        ra = a * (1.0 + e);
        rp = a * (1.0 - e);
    }

    // Euler angles (degrees, MatrixRotateXYZ's convention) of a rotation with local axes ex, ey, ez.
    Vector3 EulerFromBasis(Vector3 ex, Vector3 ey, Vector3 ez) {
        const float r00 = ex.x, r01 = ey.x, r02 = ez.x, r12 = ez.y, r22 = ez.z;
        const float ay = -asinf(std::clamp(r02, -1.0f, 1.0f));
        const float az = atan2f(r01, r00), ax = atan2f(r12, r22);
        return { -ax * RAD2DEG, -ay * RAD2DEG, -az * RAD2DEG };
    }

    // Plot cell centre in plot space (metres).
    Vector3 CellCentre(Int3 c) { return { (c.x + 0.5f) * kCellSize, (c.y + 0.5f) * kCellSize, (c.z + 0.5f) * kCellSize }; }

    Vector3 Up() { return Vector3RotateByQuaternion({ 0, 1, 0 }, rot); }
    // A plot-space point in world space while flying.
    DVec3 ToWorld(Vector3 plotPos) { return com + DVec3(Vector3RotateByQuaternion(Vector3Subtract(plotPos, comLocal), rot)); }

    // The vehicle body follows the rocket: its origin is plot (0,0,0), its rotation `rot`.
    void SyncVehicle(Scene& world) {
        Planetoid* v = Vehicle();
        if (!v) return;
        v->worldPos = com - DVec3(Vector3RotateByQuaternion(comLocal, rot));
        const Vector3 ex = Vector3RotateByQuaternion({ 1, 0, 0 }, rot), ey = Vector3RotateByQuaternion({ 0, 1, 0 }, rot), ez = Vector3RotateByQuaternion({ 0, 0, 1 }, rot);
        v->rotation = EulerFromBasis(ex, ey, ez);
        v->linearVelocity = vel;
        world.syncRenderPositions();
    }

    float RocketMass(const RocketStats& s, float fuel) { return s.dryMass + fuel; }

    // Pad: the plot's cockpit centre in world space.
    bool CockpitCentre(const Planetoid& body, const Plot& plot, DVec3& out) {
        for (const Part& p : plot.parts)
            if (p.type == PartType::Cockpit) {
                const Vector3 c = { (p.origin.x + 1.0f) * kCellSize, (p.origin.y + 1.0f) * kCellSize, (p.origin.z + 1.0f) * kCellSize };
                out = body.toWorldD(DVec3(plot.toBody(c)));
                return true;
            }
        return false;
    }

    // Bottom contact points (plot space): the lower corners of cells with nothing below them.
    void BottomPoints(const Plot& plot, std::vector<Vector3>& out) {
        out.clear();
        for (const auto& kv : plot.cells) {
            const Int3 c = kv.first;
            if (plot.cells.count({ c.x, c.y - 1, c.z })) continue;
            for (int k = 0; k < 4; ++k)
                out.push_back({ (c.x + 0.15f + 0.7f * (k & 1)) * kCellSize, c.y * kCellSize, (c.z + 0.15f + 0.7f * (k >> 1)) * kCellSize });
        }
    }

    // The ground body nearest a world point (not the rocket).
    Planetoid* GroundBody(Scene& world, const DVec3& p) { return NearestPlanetoid(world, (p - world.renderOrigin).toFloat()); }

    // Turn the vehicle plot back into a plot on `body`, upright on its grid, standing on the ground.
    void Land(Scene& world, Planetoid& body) {
        Planetoid* veh = Vehicle();
        Plot* plotPtr = VehiclePlot();
        if (!veh || !plotPtr) return;
        int ymin = 1 << 30;
        for (const auto& kv : plotPtr->cells) ymin = std::min(ymin, kv.first.y);
        // Plot frame on the body at the rocket's origin projected down to the ground.
        const DVec3 originW = veh->worldPos;
        const DVec3 originL = body.toLocalD(originW);
        const Vector3 up = Vector3Normalize(originL.toFloat());
        const Vector3 o = originL.toFloat();
        float groundAt = 0.0f;
        for (float s = 30.0f; s > -60.0f; s -= 0.1f)
            if (body.densityAt(Vector3Add(o, Vector3Scale(up, s))) >= ISO_LEVEL) { groundAt = s; break; }
        const Vector3 base = Vector3Add(o, Vector3Scale(up, groundAt + 0.05f - ymin * kCellSize));
        std::unique_ptr<Plot> plot = veh->plots->release(plotPtr);
        plot->origin = DVec3(base);
        plot->up = up;
        plot->east = Vector3CrossProduct({ 0, 1, 0 }, up);
        if (Vector3Length(plot->east) < 0.1f) plot->east = { 1, 0, 0 };
        plot->east = Vector3Normalize(plot->east);
        plot->north = Vector3CrossProduct(plot->up, plot->east);
        if (!body.plots) body.plots = std::make_unique<PlotSet>();
        Plot* placed = body.plots->adopt(std::move(plot));
        for (const Part& p : placed->parts)
            if (p.type != PartType::Count && !placed->cells.count({ p.origin.x, p.origin.y - 1, p.origin.z })) BuildTool_LayFoundation(body, *placed, p);
        padBody = body.id;
        padPlot = placed->id;
        // The vehicle body goes.
        for (auto it = world.objects.begin(); it != world.objects.end(); ++it)
            if (it->get() == veh) { world.objects.erase(it); break; }
        vehicleId = 0;
        state = RocketState::OnPad;
        throttle = 0.0f;
        autopilot = RocketAutopilot::Off;
        angVel = { 0, 0, 0 };
        ++landings;
    }

    // Lift off: the pad plot becomes the vehicle's.
    void Launch(Scene& world) {
        Planetoid* body = PadBody();
        Plot* plotPtr = PadPlot();
        if (!body || !plotPtr) return;
        // Frame and mass.
        const Vector3 eW = body->dirToWorld(plotPtr->east), uW = body->dirToWorld(plotPtr->up), nW = body->dirToWorld(plotPtr->north);
        const DVec3 originW = body->toWorldD(plotPtr->origin);
        // Rotation: plot local x -> east, y -> up, z -> north (the plot frame is left-handed: z = -north).
        const Matrix m = { eW.x, uW.x, -nW.x, 0, eW.y, uW.y, -nW.y, 0, eW.z, uW.z, -nW.z, 0, 0, 0, 0, 1 };
        rot = QuaternionNormalize(QuaternionFromMatrix(m));
        // Centre of mass in plot space (the vehicle plot keeps cell coordinates; its axes become x, y, -z).
        double mx = 0, my = 0, mz = 0, mt = 0;
        for (const auto& kv : plotPtr->cells) {
            const Part& p = plotPtr->parts[kv.second.part];
            const char* item = PartItem(p);
            float kg = materialMass.count(item) ? materialMass[item] : 100.0f;
            if (IsRocketPart(p.type) && partStats.count(item)) kg = partStats[item].mass / 8.0f * (p.type == PartType::Guidance || p.type == PartType::LandingLeg ? 8.0f : 1.0f);
            const Vector3 c = CellCentre(kv.first);
            mx += c.x * kg; my += c.y * kg; mz += -c.z * kg; mt += kg;
        }
        comLocal = mt > 0 ? Vector3{ (float)(mx / mt), (float)(my / mt), (float)(mz / mt) } : Vector3{ 0, 0, 0 };
        // The vehicle body takes the plot (origin 0, axes x, y, -z).
        float extent = 2.0f;
        for (const auto& kv : plotPtr->cells) extent = std::max(extent, Vector3Length(CellCentre(kv.first)) + 1.0f);
        auto veh = std::make_unique<Planetoid>("rocket", originW, std::make_shared<EmptyField>(extent), extent);
        veh->parent = &world.rootObject;
        veh->plots = std::make_unique<PlotSet>();
        std::unique_ptr<Plot> plot = body->plots->release(plotPtr);
        plot->origin = DVec3(0, 0, 0);
        plot->east = { 1, 0, 0 };
        plot->up = { 0, 1, 0 };
        plot->north = { 0, 0, -1 };
        veh->plots->adopt(std::move(plot));
        vehicleId = veh->id;
        com = originW + DVec3(Vector3RotateByQuaternion(comLocal, rot));
        vel = body->surfaceVelocity((com - world.renderOrigin).toFloat());
        angVel = { 0, 0, 0 };
        flightTime = 0.0;
        world.objects.push_back(std::move(veh));
        state = RocketState::Flying;
        SyncVehicle(world);
        Audio_Play(Sfx::Build, 0.6f, 0.6f);
    }

    // Steering toward a direction (world), as the guidance computer would.
    void SteerToward(Vector3 d, int guidance, float h) {
        if (guidance <= 0) return;
        const Vector3 u = Up();
        const Vector3 axis = Vector3CrossProduct(u, d);
        const float s = Vector3Length(axis), c = Vector3DotProduct(u, d);
        const float ang = atan2f(s, c);
        const float rate = kTurnRate * std::min(2, guidance);
        const Vector3 want = s > 1e-5f ? Vector3Scale(axis, std::min(2.0f * ang, rate) / s) : Vector3{ 0, 0, 0 };
        Vector3 dw = Vector3Subtract(want, angVel);
        const float dl = Vector3Length(dw), maxDw = kTurnAccel * std::min(2, guidance) * h;
        if (dl > maxDw) dw = Vector3Scale(dw, maxDw / dl);
        angVel = Vector3Add(angVel, dw);
    }

    struct Autopilot { Vector3 dir; float throttle; bool steer; };

    Autopilot RunAutopilot(Scene& world, const RocketStats& st, float fuel, float h) {
        (void)h;
        Autopilot a{ Up(), throttle, false };
        Planetoid* b = Dominant(world, com);
        if (!b) return a;
        const DVec3 r = com - b->worldPos;
        const double rl = r.length();
        const Vector3 up = Vector3Normalize(r.toFloat());
        const Vector3 vin = Vector3Subtract(vel, b->linearVelocity);
        double ra, rp;
        Orbit(*b, com, vel, ra, rp);
        const double R = b->radius * b->scale;
        const double alt = rl - R;
        const Vector3 g = Gravity(world, com);
        const float gm = Vector3Length(g);
        const float mass = RocketMass(st, fuel);
        const float amax = st.thrust / std::max(mass, 1.0f);
        switch (autopilot) {
        case RocketAutopilot::Ascent: {
            // Gravity turn toward the direction the surface moves (east), then along the velocity.
            Vector3 east = Vector3Subtract(b->surfaceVelocity((com - world.renderOrigin).toFloat()), b->linearVelocity);
            east = Vector3Subtract(east, Vector3Scale(up, Vector3DotProduct(east, up)));
            if (Vector3Length(east) < 1e-3f) east = Vector3CrossProduct({ 0, 1, 0 }, up);
            east = Vector3Normalize(east);
            const float pitch = (float)std::clamp(sqrt(std::max(0.0, alt - 200.0) / 30000.0), 0.0, 1.0) * 1.35f; // rad from vertical
            a.dir = Vector3Normalize(Vector3Add(Vector3Scale(up, cosf(pitch)), Vector3Scale(east, sinf(pitch))));
            a.steer = true;
            a.throttle = ra > 0 && ra - R >= autoTarget ? 0.0f : 1.0f;
            if (a.throttle == 0.0f && alt > 0.6 * autoTarget) autopilot = RocketAutopilot::Circularize;
            break;
        }
        case RocketAutopilot::Circularize: {
            Vector3 pro = Vector3Subtract(vin, Vector3Scale(up, Vector3DotProduct(vin, up)));
            pro = Vector3Normalize(pro);
            a.dir = pro;
            a.steer = true;
            const float vv = Vector3DotProduct(vin, up);
            const bool nearAp = vv < 25.0f;                       // close to apoapsis (climbing slowly or falling)
            a.throttle = (nearAp && Vector3DotProduct(Up(), pro) > 0.95f && rp - R < autoTarget - 6000.0) ? 1.0f : 0.0f;
            if (rp - R >= autoTarget - 6000.0) { a.throttle = 0.0f; autopilot = RocketAutopilot::Off; }
            break;
        }
        case RocketAutopilot::Deorbit: {
            a.dir = Vector3Normalize(Vector3Negate(vin));
            a.steer = true;
            a.throttle = Vector3DotProduct(Up(), a.dir) > 0.95f ? 1.0f : 0.0f;
            if (rp - R < -15000.0) { a.throttle = 0.0f; autopilot = RocketAutopilot::Land; }
            break;
        }
        case RocketAutopilot::Land: {
            const Vector3 vrel = Vector3Subtract(vel, b->surfaceVelocity((com - world.renderOrigin).toFloat()));
            const float vv = Vector3DotProduct(vrel, up);
            const Vector3 vh = Vector3Subtract(vrel, Vector3Scale(up, vv));
            double ground = alt;
            if (auto* pf = dynamic_cast<const PlanetField*>(b->field.get())) ground = alt - std::max(0.0f, pf->elevation(Vector3Normalize(b->dirToLocal(up))));
            else ground = alt - 0.72 * R; // asteroids: their mean surface
            ground -= Vector3Length(comLocal) + 1.0;   // the bottom of the rocket
            const float margin = std::max(0.3f, amax - gm);
            // Target descent speed: a braking curve using 30% of the spare thrust, then a slow final 30 m.
            float vt = -(float)std::min(400.0, 1.5 + sqrt(2.0 * 0.3 * margin * std::max(0.0, ground - 8.0)));
            if (ground < 30.0) vt = std::max(vt, -(float)(1.2 + ground * 0.12));
            Vector3 ades = Vector3Add(Vector3Scale(up, gm + 2.5f * (vt - vv)), Vector3Scale(vh, -0.5f));
            if (Vector3DotProduct(ades, up) < 0.0f) ades = Vector3Subtract(ades, Vector3Scale(up, Vector3DotProduct(ades, up)));
            const float need = Vector3Length(ades);
            a.dir = need > 0.05f ? Vector3Scale(ades, 1.0f / need) : up;
            if (Vector3Length(vrel) > 30.0f && Vector3DotProduct(ades, up) < gm * 0.5f) a.dir = Vector3Normalize(Vector3Negate(vrel)); // coasting: point retrograde
            a.steer = true;
            const float align = Vector3DotProduct(Up(), a.dir);
            a.throttle = align > 0.8f ? std::clamp(need / std::max(amax, 0.01f), 0.0f, 1.0f) : 0.0f;
            break;
        }
        default: break;
        }
        return a;
    }

    void Damage(Plot& plot, bool engines) {
        for (size_t i = 0; i < plot.parts.size(); ++i) {
            const PartType t = plot.parts[i].type;
            if (t == PartType::LandingLeg || (engines && t == PartType::RocketEngine)) { plot.remove((int)i); ++damage; }
        }
    }

    void FlightStep(Scene& world, float h, bool input) {
        Plot* plot = VehiclePlot();
        if (!plot) { state = RocketState::None; return; }
        const RocketStats st = Rocket_Stats(*plot, Dominant(world, com));
        // Controls.
        Autopilot ap{ Up(), throttle, false };
        if (autopilot != RocketAutopilot::Off) { ap = RunAutopilot(world, st, plot->fuel, h); throttle = ap.throttle; }
        bool turning = false;
        if (input && autopilot == RocketAutopilot::Off && st.guidance > 0) {
            const Vector3 ax = Vector3RotateByQuaternion({ 1, 0, 0 }, rot), az = Vector3RotateByQuaternion({ 0, 0, 1 }, rot), ay = Up();
            Vector3 cmd = { 0, 0, 0 };
            if (IsKeyDown(KEY_W)) cmd = Vector3Add(cmd, ax);
            if (IsKeyDown(KEY_S)) cmd = Vector3Subtract(cmd, ax);
            if (IsKeyDown(KEY_A)) cmd = Vector3Add(cmd, az);
            if (IsKeyDown(KEY_D)) cmd = Vector3Subtract(cmd, az);
            if (IsKeyDown(KEY_Z)) cmd = Vector3Add(cmd, ay);
            if (IsKeyDown(KEY_C)) cmd = Vector3Subtract(cmd, ay);
            if (Vector3Length(cmd) > 0.0f) { angVel = Vector3Add(angVel, Vector3Scale(Vector3Normalize(cmd), kTurnAccel * h)); turning = true; }
        }
        if (ap.steer) SteerToward(ap.dir, st.guidance, h);
        else if (!turning && st.guidance > 0) { // hold attitude
            const float l = Vector3Length(angVel), dec = kTurnAccel * h;
            angVel = l > dec ? Vector3Scale(angVel, 1.0f - dec / l) : Vector3{ 0, 0, 0 };
        }
        // Forces.
        const float mass = RocketMass(st, plot->fuel);
        const float thrust = plot->fuel > 0.0f ? throttle * st.thrust : 0.0f;
        if (thrust > 0.0f) plot->fuel = std::max(0.0f, plot->fuel - thrust / (float)(st.isp * kG0) * h);
        const Vector3 g = Gravity(world, com);
        Vector3 force = Vector3Scale(Up(), thrust);
        if (Planetoid* air = Dominant(world, com)) {
            const auto* pf = dynamic_cast<const PlanetField*>(air->field.get());
            const double alt = (com - air->worldPos).length() - air->radius * air->scale;
            if (pf && !pf->moon && alt < kAirTop) {
                const double rho = kAirDensity * exp(-std::max(0.0, alt) / kAirScaleHeight);
                const Vector3 vrel = Vector3Subtract(vel, air->surfaceVelocity((com - world.renderOrigin).toFloat()));
                const float area = std::max(1.0f, (float)plot->cells.size() * kCellSize * kCellSize * 0.25f);
                force = Vector3Add(force, Vector3Scale(vrel, -(float)(0.5 * rho * Vector3Length(vrel) * 0.8 * area)));
            }
        }
        vel = Vector3Add(vel, Vector3Scale(Vector3Add(g, Vector3Scale(force, 1.0f / std::max(mass, 1.0f))), h));
        com = com + DVec3(Vector3Scale(vel, h));
        const float w = Vector3Length(angVel);
        if (w > 1e-7f) rot = QuaternionNormalize(QuaternionMultiply(QuaternionFromAxisAngle(Vector3Scale(angVel, 1.0f / w), w * h), rot));
        flightTime += h;
        // Contacts: the rocket's bottom against the ground body below.
        if (flightTime < 1.0) return;
        static std::vector<Vector3> pts;
        BottomPoints(*plot, pts);
        Planetoid* ground = GroundBody(world, com);
        if (!ground) return;
        bool touching = false;
        for (const Vector3& p : pts) {
            const DVec3 wp = ToWorld({ p.x, p.y, -p.z });
            if (ground->densityAt(ground->toLocalD(wp).toFloat()) >= ISO_LEVEL) { touching = true; break; }
        }
        if (!touching) return;
        const Vector3 vrel = Vector3Subtract(vel, ground->surfaceVelocity((com - world.renderOrigin).toFloat()));
        const float speed = Vector3Length(vrel);
        const Vector3 gup = Vector3Normalize((com - ground->worldPos).toFloat());
        const float tilt = acosf(std::clamp(Vector3DotProduct(Up(), gup), -1.0f, 1.0f)) * RAD2DEG;
        const float safe = 4.0f + 1.5f * std::min(4, st.legs);
        lastTouchdown = speed;
        if (tilt < 35.0f || speed < 3.0f) {
            if (speed > 2.5f * safe) { Damage(*plot, true); GameUI_Toast("Crash landing: legs and engines broken"); }
            else if (speed > safe) { Damage(*plot, false); GameUI_Toast("Hard landing: landing legs broken"); }
            else GameUI_Toast(TextFormat("Touchdown at %.1f m/s", speed));
            Audio_Play(Sfx::Build, 0.6f, 0.5f);
            Land(world, *ground);
        } else {
            // Toppling or sideways: bounce off and keep going.
            const float vn = Vector3DotProduct(vrel, gup);
            if (vn < 0.0f) vel = Vector3Subtract(vel, Vector3Scale(gup, vn * 1.3f));
            angVel = Vector3Scale(angVel, 0.5f);
            if (speed > safe) { Damage(*plot, speed > 2.5f * safe); GameUI_Toast("Impact: parts broken"); }
        }
    }

    void ApplyCamera(Scene& world) {
        Plot* plot = CurrentPlot();
        if (!plot) return;
        DVec3 centre;
        Vector3 U;
        float size = 4.0f;
        for (const auto& kv : plot->cells) size = std::max(size, Vector3Length(CellCentre(kv.first)));
        if (state == RocketState::Flying) {
            centre = com;
            Planetoid* b = Dominant(world, com);
            U = b ? Vector3Normalize((com - b->worldPos).toFloat()) : Up();
        } else {
            Planetoid* body = PadBody();
            if (!body) return;
            Vector3 c = { 0, 0, 0 };
            int n = 0;
            for (const auto& kv : plot->cells) { c = Vector3Add(c, CellCentre(kv.first)); ++n; }
            c = Vector3Scale(c, 1.0f / std::max(1, n));
            centre = body->toWorldD(DVec3(plot->toBody(c)));
            U = body->dirToWorld(plot->up);
        }
        if (!IsKeyDown(KEY_LEFT_ALT)) {
            const Vector2 md = GetMouseDelta();
            camYaw -= md.x * 0.004f;
            camPitch = std::clamp(camPitch + md.y * 0.004f, -1.2f, 1.4f);
        }
        Vector3 F = Vector3CrossProduct(U, { 0, 0, 1 });
        if (Vector3Length(F) < 0.1f) F = Vector3CrossProduct(U, { 1, 0, 0 });
        F = Vector3Normalize(F);
        const Vector3 S = Vector3CrossProduct(U, F);
        const float dist = 3.0f * size + 10.0f;
        const Vector3 off = Vector3Scale(Vector3Add(Vector3Scale(Vector3Add(Vector3Scale(F, cosf(camYaw)), Vector3Scale(S, sinf(camYaw))), cosf(camPitch)), Vector3Scale(U, sinf(camPitch))), dist);
        world.setCameraWorld(centre + DVec3(off), centre);
        world.camera.up = U;
    }
}

Scene* gWorld = nullptr;
void Rocket_Init(Scene* world) { gWorld = world; Load(); }
bool Rocket_Boarded() { return state != RocketState::None; }
bool Rocket_AimingCockpit() { return state == RocketState::None && aimCockpit; }

RocketStats Rocket_Stats(const Plot& plot, const Planetoid* body) {
    Load();
    RocketStats s;
    float ispThrust = 0.0f;
    for (const Part& p : plot.parts) {
        if (p.type == PartType::Count) continue;
        const char* item = PartItem(p);
        if (IsRocketPart(p.type)) {
            const PartStat& ps = partStats[item];
            s.dryMass += ps.mass;
            s.thrust += ps.thrust;
            ispThrust += ps.isp * ps.thrust;
            s.fuelCapacity += ps.fuel;
            if (p.type == PartType::RocketEngine) ++s.engines;
            if (p.type == PartType::FuelTank) ++s.tanks;
            if (p.type == PartType::Cockpit) ++s.cockpits;
            if (p.type == PartType::Guidance) ++s.guidance;
            if (p.type == PartType::LandingLeg) ++s.legs;
        } else {
            std::vector<std::pair<Int3, CellShape>> fp;
            Plot::Footprint(p, fp);
            s.dryMass += (materialMass.count(item) ? materialMass[item] : 100.0f) * (float)fp.size();
        }
    }
    s.isp = s.thrust > 0.0f ? ispThrust / s.thrust : 0.0f;
    s.fuel = std::min(plot.fuel, s.fuelCapacity);
    const float m0 = s.dryMass + s.fuel;
    float g = 9.81f, needed = 700.0f;
    if (body) {
        const double R = body->radius * body->scale;
        g = body->surfaceGravity;
        const double rOrbit = R * 1.15 + 30000.0;
        needed = (float)(sqrt(Mu(*body) / rOrbit) * 1.25 + 150.0);
    }
    s.twr = m0 > 0.0f ? s.thrust / (m0 * g) : 0.0f;
    s.deltaV = s.dryMass > 0.0f && s.isp > 0.0f ? (float)(s.isp * kG0 * log(m0 / s.dryMass)) : 0.0f;
    s.deltaVNeeded = needed;
    if (s.cockpits == 0) s.missing.push_back("a cockpit");
    if (s.guidance == 0) s.missing.push_back("a guidance computer");
    if (s.engines == 0) s.missing.push_back("a rocket engine");
    if (s.tanks == 0) s.missing.push_back("a fuel tank");
    else if (s.fuel <= 0.0f) s.missing.push_back("fuel (I loads rocket fuel)");
    if (s.engines > 0 && s.twr < 1.2f) s.missing.push_back(TextFormat("thrust-to-weight %.2f (needs 1.2)", s.twr));
    if (s.engines > 0 && s.fuel > 0.0f && s.deltaV < needed) s.missing.push_back(TextFormat("delta-v %.0f of %.0f m/s", s.deltaV, needed));
    s.ready = s.missing.empty();
    return s;
}

int Rocket_LoadFuel(Plot& plot, int items) {
    const RocketStats s = Rocket_Stats(plot, nullptr);
    const int room = (int)((s.fuelCapacity - plot.fuel) / fuelPerItem);
    const int n = std::min({ items, room, inventory.count("rocket_fuel") });
    if (n <= 0 || !inventory.remove("rocket_fuel", n)) return 0;
    plot.fuel += n * fuelPerItem;
    return n;
}

bool Rocket_Board(Scene& world, Planetoid& body, Plot& plot) {
    (void)world;
    if (Rocket_Stats(plot, &body).cockpits == 0) return false;
    padBody = body.id;
    padPlot = plot.id;
    state = RocketState::OnPad;
    throttle = 0.0f;
    autopilot = RocketAutopilot::Off;
    player.flyMode = false;
    return true;
}

bool Rocket_Exit(Scene& world) {
    if (state != RocketState::OnPad) return false;
    Planetoid* body = PadBody();
    Plot* plot = PadPlot();
    state = RocketState::None;
    if (!body || !plot) return true;
    // Step out beside the rocket, on the ground.
    int xmax = 0;
    for (const auto& kv : plot->cells) xmax = std::max(xmax, kv.first.x);
    const Vector3 beside = plot->toBody({ (xmax + 3.0f) * kCellSize, 1.0f, 0.5f });
    player.spawnOn(world, *body, body->dirToWorld(Vector3Normalize(beside)));
    return true;
}

void Rocket_SetThrottle(float t) { throttle = std::clamp(t, 0.0f, 1.0f); }
void Rocket_SetAutopilot(RocketAutopilot mode, double target) { autopilot = mode; autoTarget = target; }

void Rocket_Update(Scene& world, bool input, float dt) {
    Load();
    aimCockpit = false;
    aimPlot = 0;
    if (state == RocketState::None) {
        // On foot: look at a rocket plot.
        Planetoid* body = player.nearest;
        if (!body || !body->plots || player.flyMode) return;
        const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
        const Vector3 lo = body->toLocal(world.camera.position), ld = body->dirToLocal(dir);
        float best = 7.0f;
        Plot* hit = nullptr;
        int hitPart = -1;
        for (auto& pl : body->plots->plots) {
            Int3 c, f; float d;
            if (!pl->raycast(lo, ld, best, c, f, d)) continue;
            const int part = pl->partAt(c);
            if (part < 0 || !IsRocketPart(pl->parts[part].type)) continue;
            best = d; hit = pl.get(); hitPart = part;
        }
        if (!hit) return;
        aimBody = body->id;
        aimPlot = hit->id;
        aimCockpit = hit->parts[hitPart].type == PartType::Cockpit;
        if (!input) return;
        if (IsKeyPressed(KEY_I)) {
            const int n = Rocket_LoadFuel(*hit, 1 << 20);
            GameUI_Toast(n > 0 ? TextFormat("Loaded %d rocket fuel (%.0f kg)", n, n * fuelPerItem) : "No rocket fuel to load, or the tanks are full");
        }
        if (aimCockpit && IsKeyPressed(KEY_E)) { Rocket_Board(world, *body, *hit); GameUI_Toast("Aboard: Shift throttle up to launch, Ctrl down, X cut, W/S/A/D/Z/C steer, E leave"); }
        return;
    }
    // Aboard.
    Plot* plot = CurrentPlot();
    if (!plot) { state = RocketState::None; return; }
    if (input) {
        if (IsKeyDown(KEY_LEFT_SHIFT)) throttle = std::min(1.0f, throttle + dt * 0.6f);
        if (IsKeyDown(KEY_LEFT_CONTROL)) throttle = std::max(0.0f, throttle - dt * 0.6f);
        if (IsKeyPressed(KEY_X)) throttle = 0.0f;
        if (state == RocketState::OnPad && IsKeyPressed(KEY_E)) { Rocket_Exit(world); return; }
    }
    if (state == RocketState::OnPad) {
        Planetoid* body = PadBody();
        if (!body) { state = RocketState::None; return; }
        DVec3 seat;
        if (CockpitCentre(*body, *plot, seat)) player.position = seat;
        const RocketStats st = Rocket_Stats(*plot, body);
        if (throttle > 0.0f && st.ready && throttle * st.thrust > (st.dryMass + st.fuel) * body->surfaceGravity) Launch(world);
        else if (throttle > 0.0f && !st.ready && input && IsKeyDown(KEY_LEFT_SHIFT)) GameUI_Toast(TextFormat("Not ready: needs %s", st.missing.front().c_str()));
        if (autopilot == RocketAutopilot::Ascent && st.ready) { throttle = 1.0f; Launch(world); }
    }
    if (state == RocketState::Flying) {
        int steps = std::max(1, (int)ceil(dt / kMaxStep));
        const float h = dt / steps;
        for (int i = 0; i < steps && state == RocketState::Flying; ++i) FlightStep(world, h, input);
        if (state == RocketState::Flying) {
            SyncVehicle(world);
            Plot* vp = VehiclePlot();
            DVec3 seat = com;
            if (vp) for (const Part& p : vp->parts) if (p.type == PartType::Cockpit) { seat = ToWorld({ (p.origin.x + 1.0f) * kCellSize, (p.origin.y + 1.0f) * kCellSize, -(p.origin.z + 1.0f) * kCellSize }); break; }
            player.position = seat;
            player.velocity = vel;
        }
    }
    ApplyCamera(world);
}

RocketTelemetry Rocket_Telemetry(Scene& world) {
    RocketTelemetry t;
    t.state = state;
    t.throttle = throttle;
    t.autopilot = autopilot;
    t.landings = landings;
    t.lastTouchdown = lastTouchdown;
    t.damage = damage;
    Plot* plot = CurrentPlot();
    if (!plot) return t;
    t.fuel = plot->fuel;
    DVec3 p = com;
    Vector3 v = vel;
    if (state == RocketState::OnPad) {
        Planetoid* body = PadBody();
        if (!body) return t;
        CockpitCentre(*body, *plot, p);
        v = body->surfaceVelocity((p - world.renderOrigin).toFloat());
    }
    Planetoid* b = Dominant(world, p);
    if (!b) return t;
    const RocketStats st = Rocket_Stats(*plot, b);
    const DVec3 r = p - b->worldPos;
    const double R = b->radius * b->scale;
    t.body = b->name;
    t.altitude = r.length() - R;
    double ground = t.altitude;
    if (auto* pf = dynamic_cast<const PlanetField*>(b->field.get())) ground = t.altitude - pf->elevation(Vector3Normalize(b->dirToLocal(r.toFloat())));
    t.groundAltitude = ground;
    const Vector3 up = Vector3Normalize(r.toFloat());
    const Vector3 vrel = Vector3Subtract(v, b->surfaceVelocity((p - world.renderOrigin).toFloat()));
    t.verticalSpeed = Vector3DotProduct(vrel, up);
    t.surfaceSpeed = Vector3Length(vrel);
    t.orbitalSpeed = Vector3Length(Vector3Subtract(v, b->linearVelocity));
    double ra, rp;
    Orbit(*b, p, v, ra, rp);
    t.apoapsis = ra < 0 ? -1.0 : ra - R;
    t.periapsis = rp - R;
    const float gl = Vector3Length(Gravity(world, p));
    t.twr = gl > 0.0f ? st.thrust / ((st.dryMass + plot->fuel) * gl) : 0.0f;
    return t;
}

void Rocket_DrawHUD(Scene& world) {
    const int x = 20, w = 470;
    int y = 290;
    if (state == RocketState::None) {
        if (!aimPlot) return;
        Planetoid* body = Planetoid::find(aimBody);
        Plot* plot = body && body->plots ? body->plots->byId(aimPlot) : nullptr;
        if (!plot) return;
        const RocketStats s = Rocket_Stats(*plot, body);
        DrawRectangle(x - 8, y - 8, w, 150 + 20 * (int)s.missing.size(), Fade(BLACK, 0.55f));
        DrawText("ROCKET  launch check", x, y, 20, SKYBLUE); y += 26;
        DrawText(TextFormat("mass %.1f t (dry %.1f t)   thrust %.0f kN   TWR %.2f", (s.dryMass + s.fuel) / 1000.0f, s.dryMass / 1000.0f, s.thrust / 1000.0f, s.twr), x, y, 16, RAYWHITE); y += 20;
        DrawText(TextFormat("fuel %.0f / %.0f kg   delta-v %.0f m/s (orbit needs ~%.0f)", s.fuel, s.fuelCapacity, s.deltaV, s.deltaVNeeded), x, y, 16, RAYWHITE); y += 20;
        DrawText(TextFormat("%d engine(s), %d tank(s), %d landing leg(s)", s.engines, s.tanks, s.legs), x, y, 16, RAYWHITE); y += 24;
        if (s.ready) { DrawText("READY", x, y, 18, GREEN); y += 22; }
        for (const std::string& m : s.missing) { DrawText(TextFormat("needs %s", m.c_str()), x, y, 16, ORANGE); y += 20; }
        DrawText(aimCockpit ? "E board   I load fuel" : "I load fuel (look at the cockpit to board)", x, y + 4, 16, Fade(RAYWHITE, 0.8f));
        return;
    }
    const RocketTelemetry t = Rocket_Telemetry(world);
    DrawRectangle(x - 8, y - 8, w, 190, Fade(BLACK, 0.55f));
    DrawText(state == RocketState::OnPad ? "ROCKET  on the ground" : "ROCKET  flying", x, y, 20, SKYBLUE); y += 26;
    DrawText(TextFormat("%s: altitude %.0f m (%.0f m above ground)", t.body.c_str(), t.altitude, t.groundAltitude), x, y, 16, RAYWHITE); y += 20;
    DrawText(TextFormat("vertical %+.1f m/s   surface %.1f m/s   orbital %.0f m/s", t.verticalSpeed, t.surfaceSpeed, t.orbitalSpeed), x, y, 16, RAYWHITE); y += 20;
    if (t.apoapsis < 0) DrawText(TextFormat("escaping   periapsis %.1f km", t.periapsis / 1000.0), x, y, 16, YELLOW);
    else DrawText(TextFormat("apoapsis %.1f km   periapsis %.1f km%s", t.apoapsis / 1000.0, t.periapsis / 1000.0, t.periapsis > 35000.0 ? "   IN ORBIT" : ""), x, y, 16, t.periapsis > 35000.0 ? GREEN : RAYWHITE);
    y += 20;
    DrawText(TextFormat("throttle %3.0f%%   fuel %.0f kg   TWR %.2f", t.throttle * 100.0f, t.fuel, t.twr), x, y, 16, RAYWHITE); y += 24;
    DrawText("Shift / Ctrl throttle   X cut   W/S A/D Z/C steer   mouse looks", x, y, 14, Fade(RAYWHITE, 0.75f)); y += 18;
    if (state == RocketState::OnPad) DrawText("E leave the cockpit", x, y, 14, Fade(RAYWHITE, 0.75f));
}

void Rocket_Reset() {
    state = RocketState::None;
    vehicleId = 0;
    padBody = 0;
    padPlot = 0;
    throttle = 0.0f;
    autopilot = RocketAutopilot::Off;
}

void Rocket_Save(FILE* f) {
    Plot* plot = state == RocketState::Flying ? VehiclePlot() : nullptr;
    if (!plot) { fprintf(f, "rocket none\n"); return; }
    fprintf(f, "rocket flying %.6f %.6f %.6f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %.3f\n", com.x, com.y, com.z, vel.x, vel.y, vel.z,
            rot.x, rot.y, rot.z, rot.w, angVel.x, angVel.y, angVel.z, comLocal.x, comLocal.y, comLocal.z, throttle, plot->fuel, flightTime);
    Planetoid* v = Vehicle();
    v->plots->save(f);
}

bool Rocket_LoadLine(const char* tok, FILE* f) {
    if (!gWorld) return false;
    Scene& world = *gWorld;
    if (strcmp(tok, "rocket") != 0) return false;
    char kind[16] = {};
    if (fscanf(f, " %15s", kind) != 1) return true;
    Rocket_Reset();
    if (strcmp(kind, "flying") != 0) return true;
    float fuel = 0.0f;
    double ft = 0.0;
    if (fscanf(f, " %lf %lf %lf %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %lf", &com.x, &com.y, &com.z, &vel.x, &vel.y, &vel.z, &rot.x, &rot.y, &rot.z, &rot.w,
               &angVel.x, &angVel.y, &angVel.z, &comLocal.x, &comLocal.y, &comLocal.z, &throttle, &fuel, &ft) != 19) return true;
    auto plots = std::make_unique<PlotSet>();
    if (!plots->load(f) || plots->plots.empty()) return true;
    Plot* p = plots->plots[0].get();
    p->origin = DVec3(0, 0, 0);
    p->east = { 1, 0, 0 };
    p->up = { 0, 1, 0 };
    p->north = { 0, 0, -1 };
    p->fuel = fuel;
    float extent = 2.0f;
    for (const auto& kv : p->cells) extent = std::max(extent, Vector3Length(CellCentre(kv.first)) + 1.0f);
    auto veh = std::make_unique<Planetoid>("rocket", com, std::make_shared<EmptyField>(extent), extent);
    veh->parent = &world.rootObject;
    veh->plots = std::move(plots);
    vehicleId = veh->id;
    world.objects.push_back(std::move(veh));
    state = RocketState::Flying;
    flightTime = std::max(ft, 2.0);
    SyncVehicle(world);
    return true;
}
