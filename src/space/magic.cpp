#include "magic.h"

#include "envflag.h"
#include "voxel.h"
#include "weather.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>

namespace {
    std::unique_ptr<ManaField> field;
    unsigned int fieldSeed = 0;
    float fieldRadius = 0.0f;

    Planetoid* StartBody(Scene& world) { return world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get()); }
    const PlanetField* StartField(Scene& world, Planetoid** body = nullptr) {
        Planetoid* p = StartBody(world);
        if (body) *body = p;
        const auto* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
        return pf && !pf->moon && !pf->sandbox ? pf : nullptr;
    }
    Color KindColour(int kind, unsigned char alpha) {
        const uint32_t c = Mana_Kind(kind).colour;
        return { (unsigned char)(c >> 16), (unsigned char)(c >> 8), (unsigned char)c, alpha };
    }

    // The overlay's pieces near the camera, rebuilt when it moves on (heights from the planet field are not free).
    struct Overlay {
        Vector3 centre{};           // the unit direction it was built around
        unsigned int seed = 0;
        std::vector<ManaSegment> segs;
        std::vector<Vector3> segPos; // body-local endpoints, 2 a piece
        std::vector<ManaNode> nodes;
        std::vector<Vector3> nodePos;
        bool built = false;
    } overlay;
    constexpr float kOverlayRange = 300.0f, kOverlayStep = 8.0f, kOverlayRebuild = 40.0f;
}

bool Magic_Enabled() { static const bool on = EnvFlag("RAYLIBGAME_MAGIC"); return on; }

const ManaField* Magic_Field(Scene& world) {
    if (!Magic_Enabled()) return nullptr;
    const PlanetField* pf = StartField(world);
    if (!pf) return nullptr;
    if (!field || fieldSeed != pf->terrainSeed() || fieldRadius != pf->radius) {
        ManaFieldParams p;
        p.seed = pf->terrainSeed();
        p.radius = pf->radius;
        p.volcanoes = pf->gen.volcanoes;
        p.karst = pf->gen.karst;
        p.meanTempC = [pf](Vector3 dir) { return pf->climateAt(dir, std::max(pf->elevation(dir, 64.0f), 0.0f)).temp; };
        field = std::make_unique<ManaField>(p);
        fieldSeed = pf->terrainSeed();
        fieldRadius = pf->radius;
        overlay = Overlay{};
    }
    return field.get();
}

ManaEnv Magic_EnvAt(Scene& world, Vector3 dir, float depth) {
    ManaEnv env;
    Planetoid* p = nullptr;
    const PlanetField* pf = StartField(world, &p);
    if (!pf) return env;
    dir = Vector3Normalize(dir);
    const float elev = pf->elevation(dir, 4.0f);
    env.elevation = elev;
    env.depth = std::max(0.0f, depth);
    env.water = elev < 0.0f;
    env.tempC = pf->climateAt(dir, std::max(elev, 0.0f)).temp;
    env.biome = (int)pf->biomeAt(dir, elev);
    const float sunUp = Vector3DotProduct(Vector3Normalize(world.sunDirection), p->dirToWorld(dir));
    env.night = std::clamp(-4.0f * sunUp, 0.0f, 1.0f);
    const WeatherState w = Weather_At(world, dir, WorldTime());
    env.wind = Vector3Length(w.wind);
    env.storm = w.kind == WeatherKind::Storm ? 1.0f : 0.0f;
    return env;
}

namespace {
    bool DebugOn() { static const bool on = EnvFlag("RAYLIBGAME_MANA_DEBUG"); return on && Magic_Enabled(); }

    // The camera on the start planet: its body-local position and unit direction.
    bool CameraLocal(Scene& world, Planetoid*& body, const PlanetField*& pf, Vector3& local, Vector3& dir) {
        pf = StartField(world, &body);
        if (!pf) return false;
        local = body->toLocal(world.camera.position);
        const float r = Vector3Length(local);
        if (r < 1.0f) return false;
        dir = Vector3Scale(local, 1.0f / r);
        return true;
    }
}

void Magic_DrawDebug3D(Scene& world) {
    if (!DebugOn()) return;
    const ManaField* f = Magic_Field(world);
    Planetoid* body = nullptr; const PlanetField* pf = nullptr; Vector3 local, dir;
    if (!f || !CameraLocal(world, body, pf, local, dir)) return;
    const float R = pf->radius;
    if (!overlay.built || overlay.seed != f->params().seed || acosf(std::clamp(Vector3DotProduct(dir, overlay.centre), -1.0f, 1.0f)) * R > kOverlayRebuild) {
        overlay.built = true;
        overlay.centre = dir;
        overlay.seed = f->params().seed;
        overlay.segs.clear(); overlay.segPos.clear(); overlay.nodes.clear(); overlay.nodePos.clear();
        f->lineSegments(dir, kOverlayRange + kOverlayRebuild, kOverlayStep, 0.0f, overlay.segs);
        auto onGround = [&](Vector3 d, float over) { return Vector3Scale(d, R + std::max(pf->elevation(d, 2.0f), 0.0f) + over); };
        for (const ManaSegment& s : overlay.segs) { overlay.segPos.push_back(onGround(s.a, 1.5f)); overlay.segPos.push_back(onGround(s.b, 1.5f)); }
        f->nodesNear(dir, 2000.0f, overlay.nodes);
        for (const ManaNode& n : overlay.nodes) overlay.nodePos.push_back(onGround(n.dir, 0.0f));
    }
    // Draw in the body's frame around an anchor under the camera (float precision), as claims.cpp's border does.
    const Vector3 anchor = Vector3Scale(dir, R);
    const Matrix m = body->localModel(DVec3{ anchor.x, anchor.y, anchor.z }, world.renderOrigin);
    rlDrawRenderBatchActive();
    rlPushMatrix();
    rlMultMatrixf(MatrixToFloat(m));
    for (size_t i = 0; i < overlay.segs.size(); ++i) {
        const ManaSegment& s = overlay.segs[i];
        const Vector3 a = Vector3Subtract(overlay.segPos[2 * i], anchor), b = Vector3Subtract(overlay.segPos[2 * i + 1], anchor);
        if (Vector3Length(a) > kOverlayRange + kOverlayRebuild) continue;
        const Color c = KindColour(Mana_Find(s.base == kManaEarth ? "earth" : s.base == kManaWater ? "water" : s.base == kManaFire ? "fire" : "air"),
                                   (unsigned char)(60 + 195 * s.strength));
        DrawLine3D(a, b, c);
        DrawLine3D(Vector3Add(a, Vector3Scale(dir, 0.5f)), Vector3Add(b, Vector3Scale(dir, 0.5f)), c); // a little thicker
    }
    for (size_t i = 0; i < overlay.nodes.size(); ++i) {
        const ManaNode& n = overlay.nodes[i];
        const Vector3 at = Vector3Subtract(overlay.nodePos[i], anchor);
        const Color c = KindColour(n.kind >= 0 ? n.kind : Mana_Aether(), 230);
        DrawSphere(Vector3Add(at, Vector3Scale(n.dir, 3.0f)), 2.0f + 2.0f * n.strength, c);
        DrawLine3D(at, Vector3Add(at, Vector3Scale(n.dir, 80.0f)), c); // a beacon seen from afar
    }
    rlDrawRenderBatchActive();
    rlPopMatrix();
}

void Magic_DrawDebugHUD(Scene& world) {
    if (!DebugOn()) return;
    const ManaField* f = Magic_Field(world);
    Planetoid* body = nullptr; const PlanetField* pf = nullptr; Vector3 local, dir;
    if (!f || !CameraLocal(world, body, pf, local, dir)) return;
    const float ground = pf->radius + pf->elevation(dir, 4.0f);
    const ManaSample s = f->density(dir, Magic_EnvAt(world, dir, ground - Vector3Length(local)));
    char line[512];
    int n = snprintf(line, sizeof line, "mana:");
    for (int k = 0; k < Mana_Count() && n < (int)sizeof line - 32; ++k)
        if (Mana_Kind(k).slice) n += snprintf(line + n, sizeof line - n, " %s %.2f", Mana_Kind(k).name.c_str(), s.of(k));
    DrawText(line, 10, GetScreenHeight() - 70, 18, RAYWHITE);
    const ManaNode* best = nullptr; float bestD = 1e9f;
    for (const ManaNode& node : overlay.nodes) {
        const float d = acosf(std::clamp(Vector3DotProduct(dir, node.dir), -1.0f, 1.0f)) * pf->radius;
        if (d < bestD) { bestD = d; best = &node; }
    }
    if (best) {
        snprintf(line, sizeof line, "nearest node: %s (%.0f%%) %.0f m", best->kind >= 0 ? Mana_Kind(best->kind).name.c_str() : "?", 100.0f * best->strength, bestD);
        DrawText(line, 10, GetScreenHeight() - 46, 18, RAYWHITE);
    }
}
