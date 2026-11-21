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
    const void* fieldOwner = nullptr; // the PlanetField it was built from (its climate callback points into it)

    Planetoid* StartBody(Scene& world) { return world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get()); }
    const PlanetField* StartField(Scene& world, Planetoid** body = nullptr) {
        Planetoid* p = StartBody(world);
        if (body) *body = p;
        const auto* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
        return pf && !pf->moon ? pf : nullptr; // (the sandbox too: the field needs only the seed and the radius; scenarios build magic machines there)
    }
    Color KindColour(int kind, unsigned char alpha) {
        const uint32_t c = Mana_Kind(kind).colour;
        return { (unsigned char)(c >> 16), (unsigned char)(c >> 8), (unsigned char)c, alpha };
    }

    int BaseKind(int base) { return Mana_Find(base == kManaEarth ? "earth" : base == kManaWater ? "water" : base == kManaFire ? "fire" : "air"); }

    // Line pieces and nodes near the camera on the ground, rebuilt when it moves on (heights from the planet field are not
    // free). The debug overlay and magic sight (#1249) each keep one.
    struct Overlay {
        Vector3 centre{};           // the unit direction it was built around
        unsigned int seed = 0;
        float range = 0.0f;         // what it was built for
        std::vector<ManaSegment> segs;
        std::vector<Vector3> segPos; // body-local endpoints, 2 a piece
        std::vector<float> segFlow;  // the line flow at the piece (the glow's strength)
        std::vector<ManaNode> nodes;
        std::vector<Vector3> nodePos;
        bool built = false;
    } overlay, sight;
    constexpr float kOverlayRange = 300.0f, kOverlayStep = 8.0f, kOverlayRebuild = 40.0f;

    // Rebuild `o` around `dir` when it moved `rebuild` m (or the range changed): pieces within range + rebuild on a `step` m
    // grid, `over` m above the ground; nodes within nodeRange.
    void Refresh(Overlay& o, const ManaField& f, const PlanetField& pf, Vector3 dir, float range, float step, float rebuild, float over, float nodeRange, bool flow) {
        const float R = pf.radius;
        if (o.built && o.seed == f.params().seed && o.range == range && acosf(std::clamp(Vector3DotProduct(dir, o.centre), -1.0f, 1.0f)) * R <= rebuild) return;
        o.built = true;
        o.centre = dir;
        o.seed = f.params().seed;
        o.range = range;
        o.segs.clear(); o.segPos.clear(); o.segFlow.clear(); o.nodes.clear(); o.nodePos.clear();
        f.lineSegments(dir, range + rebuild, step, 0.0f, o.segs);
        auto onGround = [&](Vector3 d, float up) { return Vector3Scale(d, R + std::max(pf.elevation(d, 2.0f), 0.0f) + up); };
        for (const ManaSegment& s : o.segs) {
            o.segPos.push_back(onGround(s.a, over));
            o.segPos.push_back(onGround(s.b, over));
            o.segFlow.push_back(flow ? f.lineFlow(s.base, Vector3Normalize(Vector3Add(s.a, s.b))) : s.strength);
        }
        f.nodesNear(dir, nodeRange, o.nodes);
        for (const ManaNode& n : o.nodes) o.nodePos.push_back(onGround(n.dir, 0.0f));
    }
}

bool Magic_Enabled() { static const bool on = EnvFlag("RAYLIBGAME_MAGIC"); return on; }

const ManaField* Magic_Field(Scene& world) {
    if (!Magic_Enabled()) return nullptr;
    const PlanetField* pf = StartField(world);
    if (!pf) return nullptr;
    if (!field || fieldSeed != pf->terrainSeed() || fieldRadius != pf->radius || fieldOwner != pf) { // (a reload makes a new field: never keep a pointer into the old one)
        ManaFieldParams p;
        p.seed = pf->terrainSeed();
        p.radius = pf->radius;
        p.volcanoes = pf->gen.volcanoes;
        p.karst = pf->gen.karst;
        p.meanTempC = [pf](Vector3 dir) { return pf->climateAt(dir, std::max(pf->elevation(dir, 64.0f), 0.0f)).temp; };
        field = std::make_unique<ManaField>(p);
        fieldSeed = pf->terrainSeed();
        fieldRadius = pf->radius;
        fieldOwner = pf;
        overlay = Overlay{};
        sight = Overlay{};
    }
    return field.get();
}

const ManaField* Magic_FieldFor(const Planetoid& body) {
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    return field && pf && pf == fieldOwner ? field.get() : nullptr;
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
    Refresh(overlay, *f, *pf, dir, kOverlayRange, kOverlayStep, kOverlayRebuild, 1.5f, 2000.0f, false);
    // Draw in the body's frame around an anchor under the camera (float precision), as claims.cpp's border does.
    const Vector3 anchor = Vector3Scale(dir, R);
    const Matrix m = body->localModel(DVec3{ anchor.x, anchor.y, anchor.z }, world.renderOrigin);
    rlDrawRenderBatchActive();
    rlPushMatrix();
    rlMultMatrixf(MatrixToFloat(m));
    for (size_t i = 0; i < overlay.segs.size(); ++i) {
        const ManaSegment& s = overlay.segs[i];
        const Vector3 a = Vector3Subtract(overlay.segPos[2 * i], anchor), b = Vector3Subtract(overlay.segPos[2 * i + 1], anchor);
        // Along the ground, not from the sea-level anchor: on high ground every piece was over the range (capture 2026-11-21).
        if (Vector3Length(Vector3Subtract(a, Vector3Scale(dir, Vector3DotProduct(a, dir)))) > kOverlayRange) continue;
        const Color c = KindColour(BaseKind(s.base), (unsigned char)(60 + 195 * s.strength));
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

// ---- Magic sight's look (#1249; the magic doc section 2): faint ribbons coloured by kind, glowing with the line's flow and
// pulsing along it; nodes as pulsing glows. Additive, never writing depth; hidden by rock unless the sight sees through it.
void Magic_DrawSight3D(Scene& world, float range, bool throughRock) {
    const ManaField* f = Magic_Field(world);
    Planetoid* body = nullptr; const PlanetField* pf = nullptr; Vector3 local, dir;
    if (!f || !CameraLocal(world, body, pf, local, dir)) return;
    const float step = std::clamp(range / 25.0f, 2.0f, 16.0f);
    // 0.8 m over the height field: the meshed ground and the grass stand above it in places, and a flat ribbon under them is lost.
    Refresh(sight, *f, *pf, dir, range, step, 0.2f * range, 0.8f, range * 1.5f, true);
    const float t = (float)GetTime();
    const Vector3 anchor = Vector3Scale(dir, pf->radius);
    const Vector3 eye = Vector3Subtract(local, anchor);
    const Matrix m = body->localModel(DVec3{ anchor.x, anchor.y, anchor.z }, world.renderOrigin);
    rlDrawRenderBatchActive();
    rlPushMatrix();
    rlMultMatrixf(MatrixToFloat(m));
    BeginBlendMode(BLEND_ADDITIVE);
    rlDisableDepthMask();
    rlDisableBackfaceCulling();
    if (throughRock) rlDisableDepthTest();
    rlBegin(RL_QUADS);
    for (size_t i = 0; i < sight.segs.size(); ++i) {
        const ManaSegment& s = sight.segs[i];
        const Vector3 a = Vector3Subtract(sight.segPos[2 * i], anchor), b = Vector3Subtract(sight.segPos[2 * i + 1], anchor);
        const Vector3 mid = Vector3Scale(Vector3Add(a, b), 0.5f);
        const Vector3 flat = Vector3Subtract(Vector3Subtract(mid, eye), Vector3Scale(dir, Vector3DotProduct(Vector3Subtract(mid, eye), dir)));
        const float d = Vector3Length(flat);
        if (d > range) continue;
        const float fade = 1.0f - std::clamp((d - 0.6f * range) / (0.4f * range), 0.0f, 1.0f);
        const float phase = (mid.x + mid.y + mid.z) * 0.15f;                // the pulse travels along the line
        const float pulse = 0.65f + 0.35f * sinf(2.0f * t - phase);
        const float glow = std::clamp(0.35f + sight.segFlow[i], 0.0f, 1.2f); // mana density as glow strength
        const float alpha = std::clamp(0.45f * glow * fade * pulse, 0.0f, 1.0f);
        if (alpha < 0.01f) continue;
        Vector3 side = Vector3CrossProduct(Vector3Subtract(b, a), dir);
        const float sl = Vector3Length(side);
        if (sl < 1e-4f) continue;
        side = Vector3Scale(side, (0.3f + 0.6f * s.strength) / sl);        // a ribbon 0.6-1.8 m wide
        const Color c = KindColour(BaseKind(s.base), (unsigned char)(255 * alpha));
        rlColor4ub(c.r, c.g, c.b, c.a);
        const Vector3 p[4] = { Vector3Subtract(a, side), Vector3Add(a, side), Vector3Add(b, side), Vector3Subtract(b, side) };
        for (const Vector3& v : p) rlVertex3f(v.x, v.y, v.z);
    }
    rlEnd();
    rlDrawRenderBatchActive();
    for (size_t i = 0; i < sight.nodes.size(); ++i) { // a pulsing glow and a faint column, so a node is found from afar
        const ManaNode& n = sight.nodes[i];
        const Vector3 at = Vector3Subtract(sight.nodePos[i], anchor);
        const Vector3 flat = Vector3Subtract(Vector3Subtract(at, eye), Vector3Scale(dir, Vector3DotProduct(Vector3Subtract(at, eye), dir)));
        const float d = Vector3Length(flat);
        if (d > 1.5f * range) continue;
        const float fade = 1.0f - std::clamp((d - range) / (0.5f * range), 0.0f, 1.0f);
        const float pulse = 0.5f + 0.5f * sinf(3.0f * t + (float)(n.id & 63));
        const Color c = KindColour(n.kind >= 0 ? n.kind : Mana_Aether(), (unsigned char)(255 * fade * (0.35f + 0.25f * pulse)));
        DrawSphere(Vector3Add(at, Vector3Scale(n.dir, 1.5f)), (0.6f + 0.4f * pulse) * (0.6f + 0.6f * n.strength), c);
        DrawCylinderEx(at, Vector3Add(at, Vector3Scale(n.dir, 14.0f)), 0.5f, 0.05f, 8, Color{ c.r, c.g, c.b, (unsigned char)(c.a / 2) });
    }
    rlDrawRenderBatchActive();
    if (throughRock) rlEnableDepthTest();
    rlEnableBackfaceCulling();
    rlEnableDepthMask();
    EndBlendMode();
    rlPopMatrix();
}

bool Magic_SightNodeAhead(Scene& world, float range, ManaNode& out, float& dist) {
    Planetoid* body = nullptr; const PlanetField* pf = nullptr; Vector3 local, dir;
    if (!sight.built || !CameraLocal(world, body, pf, local, dir)) return false;
    const Vector3 ahead = Vector3Normalize(Vector3Subtract(body->toLocal(world.camera.target), local)); // (toLocal maps points)
    bool found = false;
    for (size_t i = 0; i < sight.nodes.size(); ++i) {
        const Vector3 glow = Vector3Add(sight.nodePos[i], Vector3Scale(sight.nodes[i].dir, 1.5f));
        const Vector3 to = Vector3Subtract(glow, local);
        const float along = Vector3DotProduct(to, ahead);
        if (along <= 0.0f || along > range) continue;
        const float off = Vector3Length(Vector3Subtract(to, Vector3Scale(ahead, along)));
        if (off > std::max(2.0f, 0.04f * along)) continue; // the crosshair on its glow
        if (!found || along < dist) { out = sight.nodes[i]; dist = along; found = true; }
    }
    return found;
}
