#include "station.h"
#include "factory.h"
#include "plots.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace {
    // A two-body orbit around the parent: elements from the state at founding (or at a load).
    struct Station {
        uint64_t body = 0, parent = 0;
        std::string parentName;
        double mu = 1.0, a = 1.0, e = 0.0, n = 0.0, M0 = 0.0, t0 = 0.0;
        DVec3 P, Q;                       // perifocal axes (world)
        Vector3 rotation{};               // the body's Euler angles (fixed: a station does not spin)
    };
    std::vector<Station> stations;
    Scene* gScene = nullptr;   // for loads (Station_Init)

    double Mu(const Planetoid& b) { const double R = b.gravityRadius > 0.0f ? b.gravityRadius : b.radius * b.scale; return b.surfaceGravity * R * R; }
    DVec3 Cross(const DVec3& a, const DVec3& b) { return DVec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }
    double Dot(const DVec3& a, const DVec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
    DVec3 Scale(const DVec3& a, double s) { return DVec3(a.x * s, a.y * s, a.z * s); }

    // Elements from a state (r, v relative to the parent); false if not a bound orbit.
    bool Elements(Station& s, const DVec3& r, const DVec3& v, double now) {
        const double rl = r.length(), v2 = Dot(v, v);
        const double energy = 0.5 * v2 - s.mu / rl;
        if (energy >= 0.0 || rl < 1.0) return false;
        const DVec3 h = Cross(r, v);
        const double hl = h.length();
        if (hl < 1e-9) return false;
        const DVec3 ev = Scale(Cross(v, h), 1.0 / s.mu) - Scale(r, 1.0 / rl);
        s.e = std::min(ev.length(), 0.99);
        s.a = -s.mu / (2.0 * energy);
        s.P = s.e > 1e-7 ? Scale(ev, 1.0 / ev.length()) : Scale(r, 1.0 / rl);
        s.Q = Cross(Scale(h, 1.0 / hl), s.P);
        const double nu = atan2(Dot(r, s.Q), Dot(r, s.P));
        const double E = 2.0 * atan(sqrt((1.0 - s.e) / (1.0 + s.e)) * tan(nu / 2.0));
        s.M0 = E - s.e * sin(E);
        s.n = sqrt(s.mu / (s.a * s.a * s.a));
        s.t0 = now;
        return true;
    }

    void State(const Station& s, double t, DVec3& r, DVec3& v) {
        const double M = s.M0 + s.n * (t - s.t0);
        double E = M;
        for (int k = 0; k < 10; ++k) E -= (E - s.e * sin(E) - M) / (1.0 - s.e * cos(E));
        const double b = s.a * sqrt(1.0 - s.e * s.e);
        const double x = s.a * (cos(E) - s.e), y = b * sin(E);
        const double f = s.n * s.a / (1.0 - s.e * cos(E));
        r = Scale(s.P, x) + Scale(s.Q, y);
        v = Scale(s.P, -f * sin(E)) + Scale(s.Q, f * sqrt(1.0 - s.e * s.e) * cos(E));
    }

    Planetoid* Body(uint64_t id) { return Planetoid::find(id); }
    Planetoid* ParentOf(Scene& world, Station& s) {
        if (Planetoid* p = Body(s.parent)) return p;
        for (auto& o : world.objects) if (auto* p = dynamic_cast<Planetoid*>(o.get()); p && p->name == s.parentName && !p->vehicle) { s.parent = p->id; return p; }
        return nullptr;
    }

    void Place(Scene& world, Station& s) {
        Planetoid* b = Body(s.body);
        Planetoid* parent = ParentOf(world, s);
        if (!b || !parent) return;
        DVec3 r, v;
        State(s, WorldTime(), r, v);
        b->worldPos = parent->worldPos + r;
        b->linearVelocity = Vector3Add(parent->linearVelocity, v.toFloat());
        b->rotation = s.rotation;
    }

    Vector3 EulerFromQuat(Quaternion q) {
        const Vector3 ex = Vector3RotateByQuaternion({ 1, 0, 0 }, q), ey = Vector3RotateByQuaternion({ 0, 1, 0 }, q), ez = Vector3RotateByQuaternion({ 0, 0, 1 }, q);
        const float r00 = ex.x, r01 = ey.x, r02 = ez.x, r12 = ez.y, r22 = ez.z;
        const float ay = -asinf(std::clamp(r02, -1.0f, 1.0f));
        const float az = atan2f(r01, r00), ax = atan2f(r12, r22);
        return { -ax * RAD2DEG, -ay * RAD2DEG, -az * RAD2DEG };
    }
}

Planetoid* Station_Create(Scene& world, Planetoid& parent, DVec3 pos, Vector3 vel, Quaternion rot, float extent) {
    Station s;
    s.parent = parent.id;
    s.parentName = parent.name;
    s.mu = Mu(parent);
    s.rotation = EulerFromQuat(rot);
    if (!Elements(s, pos - parent.worldPos, DVec3(Vector3Subtract(vel, parent.linearVelocity)), WorldTime())) return nullptr;
    auto body = std::make_unique<Planetoid>("station", pos, std::make_shared<EmptyField>(extent), extent);
    body->parent = &world.rootObject;
    body->plots = std::make_unique<PlotSet>();
    s.body = body->id;
    Planetoid* raw = body.get();
    world.objects.push_back(std::move(body));
    stations.push_back(s);
    Place(world, stations.back());
    world.syncRenderPositions();
    return raw;
}

void Station_Update(Scene& world) {
    bool moved = false;
    for (auto it = stations.begin(); it != stations.end();) {
        Planetoid* b = Body(it->body);
        if (!b) { it = stations.erase(it); continue; }
        Place(world, *it);
        moved = true;
        ++it;
    }
    if (moved) world.syncRenderPositions();
}

bool Station_Is(const Planetoid& body) { for (const Station& s : stations) if (s.body == body.id) return true; return false; }
int Station_Count() { return (int)stations.size(); }

Planetoid* Station_Nearest(Scene& world, DVec3 pos, float range) {
    (void)world;
    Planetoid* best = nullptr;
    double bestD = range;
    for (const Station& s : stations)
        if (Planetoid* b = Body(s.body)) {
            const double d = (b->worldPos - pos).length();
            if (d < bestD) { bestD = d; best = b; }
        }
    return best;
}

bool Station_Orbit(int i, double& ap, double& pe) {
    if (i < 0 || i >= (int)stations.size()) return false;
    const Station& s = stations[i];
    Planetoid* parent = Body(s.parent);
    const double R = parent ? parent->radius * parent->scale : 0.0;
    ap = s.a * (1.0 + s.e) - R;
    pe = s.a * (1.0 - s.e) - R;
    return true;
}

void Station_Reset() {
    if (!gScene) { stations.clear(); return; }
    Scene& world = *gScene;
    for (const Station& s : stations)
        for (auto it = world.objects.begin(); it != world.objects.end(); ++it)
            if (auto* p = dynamic_cast<Planetoid*>(it->get()); p && p->id == s.body) { world.objects.erase(it); break; }
    stations.clear();
}

void Station_Save(FILE* f) {
    for (Station& s : stations) {
        Planetoid* b = Body(s.body);
        if (!b || !b->plots) continue;
        DVec3 r, v;
        State(s, WorldTime(), r, v); // the state now: elements rebuilt from it at the load
        fprintf(f, "station %s %.4f %.4f %.4f %.6f %.6f %.6f %.4f %.4f %.4f %.2f\n", s.parentName.c_str(), r.x, r.y, r.z, v.x, v.y, v.z,
                s.rotation.x, s.rotation.y, s.rotation.z, b->radius);
        b->plots->save(f);
        Factory::Save(f, b->id);
    }
}

bool Station_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "station") != 0) return false;
    char parent[128] = {};
    double rx, ry, rz, vx, vy, vz;
    float ex, ey, ez, extent;
    if (fscanf(f, " %127s %lf %lf %lf %lf %lf %lf %f %f %f %f", parent, &rx, &ry, &rz, &vx, &vy, &vz, &ex, &ey, &ez, &extent) != 11) return true;
    auto plots = std::make_unique<PlotSet>();
    const bool ok = plots->load(f);
    const std::vector<std::string> machines = Factory::ReadSaved(f);
    Scene* world = gScene;
    if (!ok || !world) return true;
    Planetoid* par = nullptr;
    for (auto& o : world->objects) if (auto* p = dynamic_cast<Planetoid*>(o.get()); p && p->name == parent && !p->vehicle) par = p;
    if (!par) return true;
    Station s;
    s.parent = par->id;
    s.parentName = par->name;
    s.mu = Mu(*par);
    s.rotation = { ex, ey, ez };
    if (!Elements(s, DVec3(rx, ry, rz), DVec3(vx, vy, vz), WorldTime())) return true;
    auto body = std::make_unique<Planetoid>("station", par->worldPos + DVec3(rx, ry, rz), std::make_shared<EmptyField>(extent), extent);
    body->parent = &world->rootObject;
    body->plots = std::move(plots);
    s.body = body->id;
    Factory::LoadSaved(body->id, machines);
    world->objects.push_back(std::move(body));
    stations.push_back(s);
    Place(*world, stations.back());
    return true;
}

void Station_Init(Scene* world) { gScene = world; }
