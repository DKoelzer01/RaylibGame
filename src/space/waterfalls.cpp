#include "waterfalls.h"
#include "particles.h"
#include "player.h"
#include "rivers.h"
#include "scene.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

namespace {
    constexpr float kDrawRange = 300.0f;  // m: sheets drawn this far from the camera
    constexpr float kMistRange = 150.0f;  // m: mist and spray emitted this close
    constexpr int kMaxDrawn = 24;
    constexpr int kRows = 14, kCols = 6;  // the sheet's grid (down the fall x across it)
    constexpr float kThrow = 2.2f;        // m/s: how fast the water leaves the lip (the arc's reach)
    constexpr float kGravity = 9.8f;

    const RiverNetwork* cachedNet = nullptr;
    std::vector<RiverNetwork::Fall> falls;   // the network's, cached
    std::vector<int> near;                   // falls near the camera, nearest first
    double nextPick = 0.0, nextMist = 0.0;
    int drawnLast = 0;
    HotShader shader{ "resources/waterfall.vs", "resources/waterfall.fs" };
    bool shaderReady = false;
    int locTime = -1, locSun = -1, locAmbient = -1, locFog = -1, locFogDensity = -1, locView = -1;

    Planetoid* StartBody(Scene& world) { return world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get()); }

    // A fall's frame (body-local): up at the lip, downstream along the ground, across the river.
    struct Frame { DVec3 lip; Vector3 up, down, across; float width, drop; };
    Frame FrameOf(const PlanetField& pf, const RiverNetwork::Fall& f) {
        Frame fr;
        fr.up = f.lip;
        fr.lip = DVec3(f.lip) * (double)(pf.radius + f.top);
        fr.down = Vector3Normalize(Vector3Subtract(f.down, Vector3Scale(fr.up, Vector3DotProduct(f.down, fr.up))));
        fr.across = Vector3Normalize(Vector3CrossProduct(fr.up, fr.down));
        fr.width = std::max(2.0f, 1.8f * f.halfWidth);
        fr.drop = f.drop;
        return fr;
    }
    // A point of the sheet (body-local): s from the lip (0) to the pool (1), a from one bank (-0.5) to the other (0.5).
    DVec3 SheetPoint(const Frame& fr, float s, float a) {
        const float fallen = fr.drop * s;
        const float out = -0.4f + kThrow * sqrtf(2.0f * fallen / kGravity); // starts a little behind the lip: over its water
        const float half = fr.width * (1.0f + 0.25f * s);                     // the sheet spreads a little as it falls
        return fr.lip + DVec3(Vector3Add(Vector3Add(Vector3Scale(fr.down, out), Vector3Scale(fr.up, -fallen)), Vector3Scale(fr.across, a * half)));
    }
}

void Waterfalls_Reset() { cachedNet = nullptr; falls.clear(); near.clear(); drawnLast = 0; }
int Waterfalls_Visible() { return drawnLast; }

void Waterfalls_Update(Scene& world, float dt) {
    (void)dt;
    static const bool off = getenv("RAYLIBGAME_WATERFALLS") && getenv("RAYLIBGAME_WATERFALLS")[0] == '0'; // A/B
    if (off) { near.clear(); return; }
    Planetoid* p = StartBody(world);
    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
    const RiverNetwork* net = pf ? pf->rivers() : nullptr;
    if (!net) { near.clear(); return; }
    if (net != cachedNet) { cachedNet = net; falls = net->waterfalls(); nextPick = 0.0; }
    const Vector3 cam = p->toLocal(world.camera.position);
    const float camR = Vector3Length(cam);
    if (camR < 1.0f) return;
    const double now = GetTime();
    if (now >= nextPick) {
        nextPick = now + 0.5;
        // Within the draw range of the camera (a chord on the sphere is close enough at these distances).
        const Vector3 camDir = Vector3Scale(cam, 1.0f / camR);
        const float cosMax = cosf((kDrawRange + 200.0f) / pf->radius);
        std::vector<std::pair<float, int>> cand;
        for (int i = 0; i < (int)falls.size(); ++i) {
            if (Vector3DotProduct(falls[i].lip, camDir) < cosMax) continue;
            const float d = Vector3Distance(Vector3Scale(falls[i].lip, pf->radius + falls[i].top), cam);
            if (d < kDrawRange) cand.push_back({ d, i });
        }
        std::sort(cand.begin(), cand.end());
        near.clear();
        for (size_t k = 0; k < cand.size() && (int)k < kMaxDrawn; ++k) near.push_back(cand[k].second);
    }
    // Mist rising from the pool and spray thrown off it, at the falls close by.
    if (now >= nextMist) {
        nextMist = now + 0.12;
        for (int i : near) {
            const RiverNetwork::Fall& f = falls[i];
            const Frame fr = FrameOf(*pf, f);
            const DVec3 foot = SheetPoint(fr, 1.0f, 0.0f);
            const Vector3 footL = { (float)foot.x, (float)foot.y, (float)foot.z };
            if (Vector3Distance(footL, cam) > kMistRange) continue;
            const float size = std::clamp(fr.drop * 0.1f, 0.3f, 2.5f);
            ParticleBurst mist;
            mist.kind = ParticleKind::Smoke; mist.color = { 236, 242, 248, 70 };
            mist.count = 1 + (int)(fr.width / 6.0f); mist.size = 0.8f + size; mist.life = 2.2f; mist.speed = 0.6f + 0.05f * fr.drop;
            mist.spread = 0.9f; mist.dir = fr.up;
            Particles_Emit(*p, Vector3Add(footL, Vector3Scale(fr.across, ((float)GetRandomValue(-50, 50) / 100.0f) * fr.width)), mist);
            ParticleBurst spray;
            spray.kind = ParticleKind::Splash; spray.color = { 220, 235, 245, 200 };
            spray.count = 3; spray.size = 0.06f; spray.life = 0.8f; spray.speed = 2.0f + 0.15f * fr.drop; spray.spread = 0.7f;
            spray.dir = Vector3Normalize(Vector3Add(fr.up, Vector3Scale(fr.down, 0.6f)));
            Particles_Emit(*p, footL, spray);
        }
    }
}

void Waterfalls_Draw3D(Scene& world) {
    drawnLast = 0;
    if (near.empty()) return;
    Planetoid* p = StartBody(world);
    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
    if (!pf || pf->rivers() != cachedNet) return;
    if (shader.reloadIfChanged(!shaderReady) || !shaderReady) {
        shaderReady = shader.shader.id > 0;
        locTime = GetShaderLocation(shader.shader, "time");
        locSun = GetShaderLocation(shader.shader, "sunDir");
        locAmbient = GetShaderLocation(shader.shader, "ambient");
        locFog = GetShaderLocation(shader.shader, "fogColor");
        locFogDensity = GetShaderLocation(shader.shader, "fogDensity");
        locView = GetShaderLocation(shader.shader, "viewPos");
    }
    if (!shaderReady) return;
    const float t = (float)fmod(GetTime(), 3600.0);
    SetShaderValue(shader.shader, locTime, &t, SHADER_UNIFORM_FLOAT);
    SetShaderValue(shader.shader, locSun, &world.frameSun, SHADER_UNIFORM_VEC3);
    SetShaderValue(shader.shader, locAmbient, &world.frameAmbient, SHADER_UNIFORM_VEC3);
    SetShaderValue(shader.shader, locFog, &world.frameFog, SHADER_UNIFORM_VEC3);
    SetShaderValue(shader.shader, locFogDensity, &world.frameFogDensity, SHADER_UNIFORM_FLOAT);
    SetShaderValue(shader.shader, locView, &world.camera.position, SHADER_UNIFORM_VEC3);
    rlDrawRenderBatchActive();
    rlDisableBackfaceCulling();
    rlDisableDepthMask(); // translucent: tested against the ground, writes no depth
    BeginShaderMode(shader.shader);
    auto toRender = [&](const DVec3& local) { return (p->toWorldD(local) - world.renderOrigin).toFloat(); };
    rlBegin(RL_TRIANGLES);
    for (int i : near) {
        const RiverNetwork::Fall& f = falls[i];
        const Frame fr = FrameOf(*pf, f);
        const Vector3 n = p->dirToWorld(fr.down); // the sheet faces downstream
        Vector3 grid[kRows + 1][kCols + 1];
        for (int r = 0; r <= kRows; ++r)
            for (int c = 0; c <= kCols; ++c) grid[r][c] = toRender(SheetPoint(fr, (float)r / kRows, (float)c / kCols - 0.5f));
        // The sheet: uv = (across 0..1, metres fallen), colour r = the fraction fallen, alpha 255 = sheet (the foam disc: 0).
        for (int r = 0; r < kRows; ++r)
            for (int c = 0; c < kCols; ++c) {
                const float u0 = (float)c / kCols, u1 = (float)(c + 1) / kCols;
                const float v0 = fr.drop * r / kRows, v1 = fr.drop * (r + 1) / kRows;
                const Vector3 q[4] = { grid[r][c], grid[r][c + 1], grid[r + 1][c + 1], grid[r + 1][c] };
                const float qu[4] = { u0, u1, u1, u0 }, qv[4] = { v0, v0, v1, v1 };
                const int tri[6] = { 0, 1, 2, 0, 2, 3 };
                for (int k : tri) {
                    rlColor4ub((unsigned char)(255.0f * qv[k] / std::max(fr.drop, 0.1f)), 255, 255, 255); // r: the fraction fallen
                    rlNormal3f(n.x, n.y, n.z);
                    rlTexCoord2f(qu[k], qv[k]);
                    rlVertex3f(q[k].x, q[k].y, q[k].z);
                }
            }
        // White water where the sheet meets the pool: a foam disc on the lower water.
        const DVec3 foot = SheetPoint(fr, 1.0f, 0.0f) + DVec3(Vector3Scale(fr.up, 0.06f));
        const float rad = fr.width * 0.85f;
        const Vector3 ctr = toRender(foot), upW = p->dirToWorld(fr.up);
        constexpr int kSeg = 16;
        for (int k = 0; k < kSeg; ++k) {
            const float a0 = 2.0f * PI * k / kSeg, a1 = 2.0f * PI * (k + 1) / kSeg;
            const Vector3 e0 = toRender(foot + DVec3(Vector3Add(Vector3Scale(fr.down, cosf(a0) * rad), Vector3Scale(fr.across, sinf(a0) * rad))));
            const Vector3 e1 = toRender(foot + DVec3(Vector3Add(Vector3Scale(fr.down, cosf(a1) * rad), Vector3Scale(fr.across, sinf(a1) * rad))));
            const Vector3 v[3] = { ctr, e0, e1 };
            const float uu[3] = { 0.5f, 0.5f + 0.5f * cosf(a0), 0.5f + 0.5f * cosf(a1) }, vv[3] = { 0.5f, 0.5f + 0.5f * sinf(a0), 0.5f + 0.5f * sinf(a1) };
            for (int j = 0; j < 3; ++j) {
                rlColor4ub(255, 255, 255, 0);
                rlNormal3f(upW.x, upW.y, upW.z);
                rlTexCoord2f(uu[j], vv[j]);
                rlVertex3f(v[j].x, v[j].y, v[j].z);
            }
        }
        ++drawnLast;
    }
    rlEnd();
    EndShaderMode();
    rlEnableDepthMask();
    rlEnableBackfaceCulling();
}
