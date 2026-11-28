#include "waterfalls.h"
#include "audio.h"
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
    // Near the camera a sheet is animated with its foam; beyond kNearRange it is a static sheet (a coarser grid, no foam), out
    // to kFarRange where it fades out over the last fifth (#612: falls popped in at 300 m and showed as bare rock steps
    // beyond). The animation fades out between kNearRange and kStaticRange, so nothing switches at once.
    constexpr float kNearRange = 300.0f;   // m: animated sheets and foam
    constexpr float kStaticRange = 420.0f; // m: fully static from here
    constexpr float kFarRange = 2500.0f;   // m: static sheets drawn this far, fading out from 0.8 of it
    constexpr float kMistRange = 150.0f;   // m: mist and spray emitted this close
    constexpr int kMaxNear = 24, kMaxFar = 96;
    constexpr int kRows = 14, kCols = 6;   // the near sheet's grid (down the fall x across it)
    constexpr int kFarRows = 4, kFarCols = 2;
    constexpr float kThrow = 2.2f;        // m/s: how fast the water leaves the lip (the arc's reach)
    constexpr float kGravity = 9.8f;

    const RiverNetwork* cachedNet = nullptr;
    std::vector<RiverNetwork::Fall> falls;   // the network's, cached
    std::vector<int> near;                   // falls within kFarRange of the camera, nearest first
    std::vector<float> nearDist;             // their distances at the last pick (m)
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
        fr.width = std::max(1.2f, 1.8f * f.halfWidth);
        fr.drop = f.drop;
        return fr;
    }
    // How big the fall's effects are (#614: a small fall made mist and foam as big as a large one): the water's power, by
    // its drop and its width (the flow), 0 for a 2 m drop of a 1 m stream .. 1 for a 20 m drop of a 14 m-wide river.
    float PowerOf(const Frame& fr) {
        return std::clamp(sqrtf(std::max(0.0f, fr.drop - 1.0f) / 19.0f * std::min(fr.width, 14.0f) / 14.0f), 0.0f, 1.0f);
    }
    // A point of the sheet (body-local): s from the lip (0) to the pool (1), a from one bank (-0.5) to the other (0.5).
    DVec3 SheetPoint(const Frame& fr, float s, float a) {
        const float fallen = fr.drop * s;
        const float out = -0.4f + kThrow * sqrtf(2.0f * fallen / kGravity); // starts a little behind the lip: over its water
        const float half = fr.width * (1.0f + 0.25f * s);                     // the sheet spreads a little as it falls
        return fr.lip + DVec3(Vector3Add(Vector3Add(Vector3Scale(fr.down, out), Vector3Scale(fr.up, -fallen)), Vector3Scale(fr.across, a * half)));
    }
}

namespace { constexpr int kSounding = 3; int soundEm[kSounding] = { -1, -1, -1 }; } // (#1242) the falls' sound emitters

void Waterfalls_Reset() {
    cachedNet = nullptr; falls.clear(); near.clear(); nearDist.clear(); drawnLast = 0;
    for (int& e : soundEm) { Audio_EmitterDestroy(e); e = -1; }
}
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
        const float cosMax = cosf((kFarRange + 200.0f) / pf->radius);
        std::vector<std::pair<float, int>> cand;
        for (int i = 0; i < (int)falls.size(); ++i) {
            if (Vector3DotProduct(falls[i].lip, camDir) < cosMax) continue;
            const float d = Vector3Distance(Vector3Scale(falls[i].lip, pf->radius + falls[i].top), cam);
            if (d < kFarRange) cand.push_back({ d, i });
        }
        std::sort(cand.begin(), cand.end());
        near.clear(); nearDist.clear();
        int nNear = 0;
        for (size_t k = 0; k < cand.size() && (int)near.size() < kMaxNear + kMaxFar; ++k) {
            if (cand[k].first < kStaticRange && nNear >= kMaxNear) continue; // (nearest first: the closest get the animated sheets)
            nNear += cand[k].first < kStaticRange;
            near.push_back(cand[k].second); nearDist.push_back(cand[k].first);
        }
    }
    // Mist rising from the pool and spray thrown off it, at the falls close by.
    if (now >= nextMist) {
        nextMist = now + 0.12;
        for (size_t k = 0; k < near.size(); ++k) {
            if (nearDist[k] > kMistRange + 60.0f) continue;
            const int i = near[k];
            const RiverNetwork::Fall& f = falls[i];
            const Frame fr = FrameOf(*pf, f);
            const DVec3 foot = SheetPoint(fr, 1.0f, 0.0f);
            const Vector3 footL = { (float)foot.x, (float)foot.y, (float)foot.z };
            if (Vector3Distance(footL, cam) > kMistRange) continue;
            // Scaled by the fall's power (#614): a small fall gets a little spray and a wisp of mist, a big one a cloud.
            const float power = PowerOf(fr);
            if (power > 0.15f && GetRandomValue(0, 99) < (int)(100.0f * std::min(1.0f, power * 1.6f))) {
                ParticleBurst mist;
                mist.kind = ParticleKind::Smoke; mist.color = { 236, 242, 248, (unsigned char)(30 + 40 * power) };
                mist.count = 1 + (int)(power * fr.width / 6.0f); mist.size = 0.3f + 2.2f * power; mist.life = 1.0f + 1.4f * power;
                mist.speed = 0.25f + 0.04f * fr.drop * power; mist.spread = 0.9f; mist.dir = fr.up;
                Particles_Emit(*p, Vector3Add(footL, Vector3Scale(fr.across, ((float)GetRandomValue(-50, 50) / 100.0f) * fr.width)), mist);
            }
            ParticleBurst spray;
            spray.kind = ParticleKind::Splash; spray.color = { 220, 235, 245, 200 };
            spray.count = 1 + (int)(2.0f * power + fr.width / 10.0f); spray.size = 0.035f + 0.03f * power; spray.life = 0.4f + 0.4f * power;
            spray.speed = 0.8f + 0.12f * std::min(fr.drop, 20.0f) * (0.4f + 0.6f * power); spray.spread = 0.7f;
            spray.dir = Vector3Normalize(Vector3Add(fr.up, Vector3Scale(fr.down, 0.6f)));
            Particles_Emit(*p, footL, spray);
        }
    }
    // (#1242) The roar of the nearest falls (audio v2): at the pool, as loud as the fall is powerful, deeper the bigger.
    if (Audio_V2()) {
        int* em = soundEm;
        static int snd = -2;
        if (snd == -2) snd = Audio_Sound("bed_waterfall");
        for (int k = 0; k < kSounding && snd >= 0; ++k) {
            const bool have = k < (int)near.size() && nearDist[k] < 220.0f;
            if (!have) { if (em[k] >= 0) Audio_EmitterSetPlaying(em[k], false); continue; }
            const Frame fr = FrameOf(*pf, falls[near[k]]);
            const float power = PowerOf(fr);
            if (em[k] < 0) {
                AudioAttenuation att;
                att.refDistance = 4.0f + 0.5f * fr.width;
                att.maxDistance = 220.0f;
                att.lowpassDistance = 120.0f;
                em[k] = Audio_EmitterCreate(snd, AudioBus::Ambience, AudioPriority::Low, att, true);
                if (em[k] < 0) continue;
            }
            Audio_EmitterSetPosition(em[k], p->toWorldD(SheetPoint(fr, 0.9f, 0.0f)));
            Audio_EmitterSetVolume(em[k], 0.3f + 0.7f * power);
            Audio_EmitterSetPitch(em[k], 1.05f - 0.3f * power);
            Audio_EmitterSetPlaying(em[k], true);
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
    const float t = (float)fmod(GameClock(), 3600.0); // paused with the game (#511)
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
    const Vector3 camL = p->toLocal(world.camera.position);
    for (int i : near) {
        const RiverNetwork::Fall& f = falls[i];
        const Frame fr = FrameOf(*pf, f);
        const DVec3 mid = SheetPoint(fr, 0.5f, 0.0f);
        const float dist = Vector3Distance({ (float)mid.x, (float)mid.y, (float)mid.z }, camL);
        if (dist >= kFarRange) continue;
        // (#2140) Out of view: the sheet (it spreads to 1.25 x its width and arcs out as it falls) and the foam disc at its foot.
        if (!world.inView(toRender(mid), (0.75f * fr.drop + 1.25f * fr.width + 2.0f) * p->scale)) continue;
        // Colour g: the sheet's opacity (fading in from kFarRange, #612); b: how animated it is (none beyond kStaticRange).
        const float fade = 1.0f - std::clamp((dist - 0.8f * kFarRange) / (0.2f * kFarRange), 0.0f, 1.0f);
        const float anim = 1.0f - std::clamp((dist - kNearRange) / (kStaticRange - kNearRange), 0.0f, 1.0f);
        const unsigned char g8 = (unsigned char)(255.0f * fade), b8 = (unsigned char)(255.0f * anim);
        const bool full = dist < kStaticRange;
        const int rows = full ? kRows : kFarRows, cols = full ? kCols : kFarCols;
        const Vector3 n = p->dirToWorld(fr.down); // the sheet faces downstream
        Vector3 grid[kRows + 1][kCols + 1];
        for (int r = 0; r <= rows; ++r)
            for (int c = 0; c <= cols; ++c) grid[r][c] = toRender(SheetPoint(fr, (float)r / rows, (float)c / cols - 0.5f));
        // The sheet: uv = (across 0..1, metres fallen), colour r = the fraction fallen, alpha 255 = sheet (the foam disc: 0).
        for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c) {
                const float u0 = (float)c / cols, u1 = (float)(c + 1) / cols;
                const float v0 = fr.drop * r / rows, v1 = fr.drop * (r + 1) / rows;
                const Vector3 q[4] = { grid[r][c], grid[r][c + 1], grid[r + 1][c + 1], grid[r + 1][c] };
                const float qu[4] = { u0, u1, u1, u0 }, qv[4] = { v0, v0, v1, v1 };
                const int tri[6] = { 0, 1, 2, 0, 2, 3 };
                for (int k : tri) {
                    rlColor4ub((unsigned char)(255.0f * qv[k] / std::max(fr.drop, 0.1f)), g8, b8, 255); // r: the fraction fallen
                    rlNormal3f(n.x, n.y, n.z);
                    rlTexCoord2f(qu[k], qv[k]);
                    rlVertex3f(q[k].x, q[k].y, q[k].z);
                }
            }
        ++drawnLast;
        if (!full) continue;
        // White water where the sheet meets the pool: a foam disc on the lower water, sized by the fall's power (#614).
        const DVec3 foot = SheetPoint(fr, 1.0f, 0.0f) + DVec3(Vector3Scale(fr.up, 0.06f));
        const float rad = std::min(fr.width * 0.85f, 0.5f * fr.width + 0.3f * fr.drop) * (0.45f + 0.55f * PowerOf(fr));
        const Vector3 ctr = toRender(foot), upW = p->dirToWorld(fr.up);
        constexpr int kSeg = 16;
        for (int k = 0; k < kSeg; ++k) {
            const float a0 = 2.0f * PI * k / kSeg, a1 = 2.0f * PI * (k + 1) / kSeg;
            const Vector3 e0 = toRender(foot + DVec3(Vector3Add(Vector3Scale(fr.down, cosf(a0) * rad), Vector3Scale(fr.across, sinf(a0) * rad))));
            const Vector3 e1 = toRender(foot + DVec3(Vector3Add(Vector3Scale(fr.down, cosf(a1) * rad), Vector3Scale(fr.across, sinf(a1) * rad))));
            const Vector3 v[3] = { ctr, e0, e1 };
            const float uu[3] = { 0.5f, 0.5f + 0.5f * cosf(a0), 0.5f + 0.5f * cosf(a1) }, vv[3] = { 0.5f, 0.5f + 0.5f * sinf(a0), 0.5f + 0.5f * sinf(a1) };
            for (int j = 0; j < 3; ++j) {
                rlColor4ub(255, g8, b8, 0);
                rlNormal3f(upW.x, upW.y, upW.z);
                rlTexCoord2f(uu[j], vv[j]);
                rlVertex3f(v[j].x, v[j].y, v[j].z);
            }
        }
    }
    rlEnd();
    EndShaderMode();
    rlEnableDepthMask();
    rlEnableBackfaceCulling();
}
