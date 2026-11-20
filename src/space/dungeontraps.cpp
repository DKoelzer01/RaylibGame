#include "dungeontraps.h"
#include "audio.h"
#include "combat.h"
#include "net.h"
#include "particles.h"
#include "player.h"
#include "playersession.h"
#include "plots.h"
#include "scene.h"
#include "simclock.h"
#include "survival.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <utility>
#include <vector>

namespace {
    struct Trap { Planetoid* body; uint32_t plot; int part; PartType type; uint8_t kind; Vector3 at; }; // at: plot metres, the plate's top middle
    struct Falling { Planetoid* body; uint32_t plot; int part; Vector3 at; float top; double started; bool landed = false; };
    struct Dart { Planetoid* body; uint32_t plot; Vector3 from, to; double started; float length; bool hit = false; };
    std::vector<Trap> traps;
    std::vector<Falling> falling;
    std::vector<Dart> darts;
    std::map<std::pair<uint32_t, int>, double> rearmAt;  // dart plates: (plot, part) -> when they fire again
    double nextScan = -1.0, nextHaze = 0.0;
    DungeonTrapStats stats;

    constexpr float kScanRange = 120.0f;   // m from the player to a plot's origin (a dungeon reaches ~60 m from its ruin)
    constexpr double kWarn = 0.5;          // s: the click to the block landing
    constexpr float kBlockDamage = 35.0f, kDartDamage = 12.0f;
    constexpr float kDartSpeed = 24.0f;    // m/s
    constexpr double kDartRearm = 2.5;     // s
    constexpr double kDartStuck = 3.0;     // s a spent dart stays in the far wall (it flies too fast to see)
    constexpr float kGasRadius = 3.5f;     // m round a point 1.2 m over the vent

    Plot* PlotOf(Planetoid* b, uint32_t id) { return b && b->plots ? b->plots->byId(id) : nullptr; }
    bool Rock(const Planetoid& b, const Plot& p, Vector3 plotPos) { return b.densityAt(p.toBody(plotPos)) >= ISO_LEVEL; }
    // Free metres from a plot point along a plot-space direction, up to `max` (0.1 m steps).
    float FreeAlong(const Planetoid& b, const Plot& p, Vector3 from, Vector3 dir, float max) {
        for (float d = 0.1f; d < max; d += 0.1f) if (Rock(b, p, Vector3Add(from, Vector3Scale(dir, d)))) return d;
        return max;
    }
    int StoneMaterial() {
        for (int i = 0; i < BuildMaterialCount(); ++i) if (strcmp(GetBuildMaterial(i).item, "stone_brick") == 0) return i;
        return 0;
    }
    // A plot-space box (centre, size) drawn in render space with the plot's axes.
    void DrawPlotBox(const Planetoid& b, const Plot& p, Vector3 centre, Vector3 size, Color c, bool wires) {
        const Vector3 e = b.dirToWorld(p.east), u = b.dirToWorld(p.up), n = b.dirToWorld(p.north), o = b.toWorld(p.toBody(centre));
        const Matrix m = { e.x, u.x, n.x, o.x, e.y, u.y, n.y, o.y, e.z, u.z, n.z, o.z, 0, 0, 0, 1 };
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloatV(m).v);
        DrawCubeV({ 0, 0, 0 }, size, c);
        if (wires) DrawCubeWiresV({ 0, 0, 0 }, size, Fade(BLACK, 0.35f));
        rlPopMatrix();
    }
    void Hit(Scene& world, const Planetoid& b, const Plot& p, Vector3 from, float amount, DamageKind kind, float knockback) {
        DamageEvent e;
        e.amount = amount; e.kind = kind; e.hasSource = true; e.knockback = knockback;
        e.source = world.renderOrigin + DVec3(b.toWorld(p.toBody(from)));
        Combat_Damage(LocalSession(), e);
    }
    // The local player's feet against a plot point: horizontal (plot xz) and vertical (plot y) offsets, m.
    bool FeetNear(const Scene& world, const Planetoid& b, const Plot& p, Vector3 at, float horiz, float below, float above) {
        const Vector3 me = p.toPlot(b.toLocal((player.position - world.renderOrigin).toFloat()));
        const float dy = me.y - at.y;
        return fabsf(me.x - at.x) < horiz && fabsf(me.z - at.z) < horiz && dy > -below && dy < above;
    }
}

void DungeonTraps_Reset() { traps.clear(); falling.clear(); darts.clear(); rearmAt.clear(); nextScan = -1.0; stats = {}; }
const DungeonTrapStats& DungeonTraps_Stats() { return stats; }
void DungeonTraps_Rescan() { nextScan = -1.0; }

void DungeonTraps_Update(Scene& world) {
    const double now = SimTime();
    Planetoid* body = player.nearest;
    if (now >= nextScan) { // the traps round the player
        nextScan = now + 0.5;
        traps.clear();
        stats.plates = stats.vents = 0;
        if (body && body->plots) {
            const Vector3 me = body->toLocal((player.position - world.renderOrigin).toFloat());
            for (auto& plot : body->plots->plots) {
                if (Vector3Distance(plot->toBody({ 0, 0, 0 }), me) > kScanRange) continue;
                for (int i = 0; i < (int)plot->parts.size(); ++i) {
                    const Part& part = plot->parts[i];
                    if (part.type != PartType::PressurePlate && part.type != PartType::GasVent) continue;
                    const Vector3 at = { (part.origin.x + 0.5f) * kCellSize, part.origin.y * kCellSize + 0.05f, (part.origin.z + 0.5f) * kCellSize };
                    traps.push_back({ body, plot->id, i, part.type, part.material, at });
                    (part.type == PartType::GasVent ? stats.vents : stats.plates)++;
                }
            }
        }
    }
    const bool authority = Net_IsAuthority() && !LocalSession().survival.dead;
    for (const Trap& t : traps) {
        Plot* plot = PlotOf(t.body, t.plot);
        if (!plot || t.part >= (int)plot->parts.size() || plot->parts[t.part].type != t.type) continue; // gone since the scan
        if (t.type == PartType::GasVent) {
            if (now >= nextHaze && Vector3Distance(plot->toBody(t.at), t.body->toLocal((player.position - world.renderOrigin).toFloat())) < 30.0f) {
                ParticleBurst s; // a sickly haze over the vent
                s.kind = ParticleKind::Smoke; s.color = { 132, 168, 70, 70 };
                s.count = 2; s.size = 0.45f; s.life = 3.5f; s.speed = 0.25f; s.spread = 1.0f; s.dir = plot->up;
                Particles_Emit(*t.body, plot->toBody(Vector3Add(t.at, { 0, 0.3f, 0 })), s);
            }
            continue;
        }
        if (!authority || !FeetNear(world, *t.body, *plot, t.at, 0.42f, 0.3f, 0.6f)) continue;
        if (t.kind == 0) { // a falling block: once (the plate is spent when it lands)
            bool already = false;
            for (const Falling& f : falling) already = already || (f.plot == t.plot && f.part == t.part);
            if (already) continue;
            const float top = std::min(FreeAlong(*t.body, *plot, Vector3Add(t.at, { 0, 0.3f, 0 }), { 0, 1, 0 }, 4.0f) + 0.3f, 3.5f);
            falling.push_back({ t.body, t.plot, t.part, t.at, top, now });
            ++stats.blocksDropped;
            Audio_Play(Sfx::Build, 0.7f, 0.6f); // the plate's click
            ParticleBurst d; d.kind = ParticleKind::Dust; d.color = { 140, 130, 115, 200 }; d.count = 10; d.size = 0.08f; d.life = 0.8f; d.speed = 0.6f; d.spread = 0.5f;
            d.dir = Vector3Negate(plot->up);
            Particles_Emit(*t.body, plot->toBody(Vector3Add(t.at, { 0, top - 0.1f, 0 })), d); // grit from the ceiling: look up
        } else { // darts: from the nearer wall across the plate, again after the re-arm
            auto key = std::make_pair(t.plot, t.part);
            auto it = rearmAt.find(key);
            if (it != rearmAt.end() && now < it->second) continue;
            rearmAt[key] = now + kDartRearm;
            const Vector3 mid = Vector3Add(t.at, { 0, 1.0f, 0 });
            // The corridor runs where the most room is; the dart crosses it.
            float best = -1.0f; Vector3 along{ 1, 0, 0 };
            for (int k = 0; k < 4; ++k) {
                const float a = k * PI / 4.0f;
                const Vector3 d = { cosf(a), 0, sinf(a) };
                const float room = FreeAlong(*t.body, *plot, mid, d, 6.0f) + FreeAlong(*t.body, *plot, mid, Vector3Negate(d), 6.0f);
                if (room > best) { best = room; along = d; }
            }
            Vector3 across = { -along.z, 0, along.x };
            float near = FreeAlong(*t.body, *plot, mid, across, 4.0f), far = FreeAlong(*t.body, *plot, mid, Vector3Negate(across), 4.0f);
            if (far < near) { std::swap(near, far); across = Vector3Negate(across); }
            const Vector3 from = Vector3Add(mid, Vector3Scale(across, near - 0.05f)), to = Vector3Subtract(mid, Vector3Scale(across, far - 0.05f));
            darts.push_back({ t.body, t.plot, from, to, now, near + far });
            ++stats.dartsShot;
            Audio_Play(Sfx::Jump, 0.6f, 1.8f); // the plate's click and the dart's hiss
        }
    }
    // Blocks land kWarn after the click; whoever is under it then is hit and it breaks on them, else it stays as rubble
    // (never a cube round someone's feet).
    for (Falling& f : falling) {
        if (f.landed || now - f.started < kWarn) continue;
        f.landed = true;
        Plot* plot = PlotOf(f.body, f.plot);
        if (!plot) continue;
        const bool under = FeetNear(world, *f.body, *plot, f.at, 0.6f, 0.5f, 1.9f);
        if (authority && under) { Hit(world, *f.body, *plot, Vector3Add(f.at, { 0, f.top, 0 }), kBlockDamage, DamageKind::Other, 3.0f); ++stats.blockHits; }
        if (f.part < (int)plot->parts.size() && plot->parts[f.part].type == PartType::PressurePlate) {
            const Int3 cell = plot->parts[f.part].origin;
            plot->remove(f.part); // spent
            if (!under) plot->place(Part{ PartType::Cube, cell, 0, (uint8_t)StoneMaterial() });
        }
        ParticleBurst d; d.kind = ParticleKind::Dust; d.color = { 135, 128, 118, 230 }; d.count = 22; d.size = 0.16f; d.life = 1.4f; d.speed = 1.6f; d.spread = 0.9f;
        d.dir = plot->up;
        Particles_Emit(*f.body, plot->toBody(Vector3Add(f.at, { 0, 0.2f, 0 })), d);
        Audio_Play(Sfx::Dig, 1.0f, 0.55f);
        DungeonTraps_Rescan();
    }
    falling.erase(std::remove_if(falling.begin(), falling.end(), [](const Falling& f) { return f.landed; }), falling.end());
    // Darts hit whoever stands on the plate's line when they pass it.
    for (Dart& d : darts) {
        Plot* plot = PlotOf(d.body, d.plot);
        if (!plot) { d.hit = true; continue; }
        const float flown = (float)(now - d.started) * kDartSpeed;
        if (!d.hit && authority && flown > 0.0f) {
            const Vector3 me = plot->toPlot(d.body->toLocal((player.position - world.renderOrigin).toFloat()));
            // The player's place along the dart's path; hit once it has flown that far, if they are within 0.4 m of it.
            const Vector3 path = Vector3Subtract(d.to, d.from);
            const float len = Vector3Length(path);
            const Vector3 dir = Vector3Scale(path, 1.0f / std::max(len, 0.01f));
            const Vector3 rel = Vector3Subtract({ me.x, d.from.y, me.z }, d.from);
            const float s = Vector3DotProduct(rel, dir);
            const float off = Vector3Length(Vector3Subtract(rel, Vector3Scale(dir, s)));
            const float feetToDart = d.from.y - me.y;
            if (s >= 0.0f && s <= len && flown >= s && off < 0.4f && feetToDart > 0.2f && feetToDart < 1.8f) {
                Hit(world, *d.body, *plot, d.from, kDartDamage, DamageKind::Projectile, 1.5f);
                ++stats.dartHits;
                d.hit = true;
            }
        }
        if (flown > d.length + 0.3f) d.hit = true; // stuck in the far wall
    }
    darts.erase(std::remove_if(darts.begin(), darts.end(), [&](const Dart& d) { return d.hit && now - d.started > d.length / kDartSpeed + kDartStuck; }), darts.end()); // stuck in the wall a while: seen
    if (now >= nextHaze) nextHaze = now + 0.2;
}

bool DungeonTraps_InGas(const Scene& world, Vector3 renderPos) {
    (void)world;
    for (const Trap& t : traps) {
        if (t.type != PartType::GasVent) continue;
        Plot* plot = PlotOf(t.body, t.plot);
        if (!plot) continue;
        if (Vector3Distance(renderPos, t.body->toWorld(plot->toBody(Vector3Add(t.at, { 0, 1.2f, 0 })))) < kGasRadius) return true;
    }
    return false;
}

void DungeonTraps_Draw(Scene& world) {
    (void)world;
    const double now = SimTime();
    for (const Falling& f : falling) {
        const Plot* plot = PlotOf(f.body, f.plot);
        if (!plot) continue;
        const float k = (float)std::clamp((now - f.started) / kWarn, 0.0, 1.0);
        const float y = f.top - 0.25f - (f.top - 0.25f - 0.25f) * k * k; // falls from the ceiling, gathering speed
        DrawPlotBox(*f.body, *plot, Vector3Add(f.at, { 0, y - 0.05f, 0 }), { 0.5f, 0.5f, 0.5f }, { 92, 88, 82, 255 }, true); // (unlit: dark, or it glows in a dark corridor)
    }
    for (const Dart& d : darts) {
        const Plot* plot = PlotOf(d.body, d.plot);
        if (!plot) continue;
        const float flown = std::min((float)(now - d.started) * kDartSpeed, d.length - 0.05f);
        const Vector3 dir = Vector3Normalize(Vector3Subtract(d.to, d.from));
        const Vector3 tip = Vector3Add(d.from, Vector3Scale(dir, flown));
        // A 0.35 m shaft (small boxes along it: the dart flies at any angle in the plot's plane) and red fletching at the tail.
        for (int k = 0; k < 6; ++k)
            DrawPlotBox(*d.body, *plot, Vector3Subtract(tip, Vector3Scale(dir, 0.03f + 0.055f * k)), { 0.035f, 0.035f, 0.035f }, { 150, 118, 72, 255 }, false);
        DrawPlotBox(*d.body, *plot, Vector3Subtract(tip, Vector3Scale(dir, 0.36f)), { 0.07f, 0.07f, 0.07f }, { 190, 40, 32, 255 }, false);
    }
}

// --test-dungeon-traps (#719, best with --sandbox): a 12 m trench dug 2.4 m deep at the start with a falling-block plate,
// a dart plate and a gas vent along its floor; the player is set on each in turn. Checks: the block drops, hits (health
// down), and leaves rubble where the plate was; a dart flies across and hits; by the vent the air is gas (not breathable)
// and the oxygen falls. The camera watches each trap from the trench's end.
int DungeonTraps_Test(Scene& world, Planetoid& body, int frame) {
    static Plot* plot = nullptr;
    static uint32_t plotId = 0;
    static Vector3 p0{}, plate0{}, plate1{}, vent{};
    static float fy = 0.0f, health0 = 0.0f, healthBlock = 0.0f, oxygen0 = 0.0f;
    static Int3 rubbleCell{};
    static bool gasSeen = false, breathedInGas = false, airBefore = false;
    SurvivalState& sv = LocalSession().survival;
    auto toRender = [&](Vector3 plotPos) { return body.toWorld(plot->toBody(plotPos)); };
    auto stand = [&](Vector3 plotPos) { player.placeAt(world, body, toRender({ plotPos.x, plotPos.y + 0.08f, plotPos.z })); };
    if (frame == 1) {
        const Vector3 me = body.toLocal((player.position - world.renderOrigin).toFloat());
        if (!body.plots) body.plots = std::make_unique<PlotSet>();
        plot = body.plots->plotFor(me, 0.0f, true);
        if (!plot) { printf("bench: dungeon traps: no plot -> FAIL\n"); return -1; }
        plotId = plot->id;
        p0 = plot->toPlot(me);
        // The ground under the player (it may still be dropping from the spawn): the first rock going down.
        for (float d = -2.0f; d < 40.0f; d += 0.05f)
            if (body.densityAt(plot->toBody({ p0.x, p0.y - d, p0.z })) >= ISO_LEVEL) { p0.y -= d; break; }
        fy = floorf((p0.y - 2.4f) / kCellSize) * kCellSize;
        TerrainEdit e{};
        e.carve = true; e.box = true;
        e.center = plot->toBody({ p0.x + 5.0f, fy + 1.6f, p0.z });
        e.half = { 7.0f, 1.6f, 0.75f }; e.ax = plot->east; e.ay = plot->up; e.az = plot->north;
        e.radius = Vector3Length(e.half);
        ApplyTerrainEdit(body, e);
        auto cell = [&](float x, float z) { return Int3{ (int)floorf(x / kCellSize), (int)floorf((fy + 0.1f) / kCellSize), (int)floorf(z / kCellSize) }; };
        rubbleCell = cell(p0.x + 1.0f, p0.z);
        plot->place(Part{ PartType::PressurePlate, rubbleCell, 0, 0 });
        plot->place(Part{ PartType::PressurePlate, cell(p0.x + 5.0f, p0.z), 0, 1 });
        plot->place(Part{ PartType::GasVent, cell(p0.x + 9.0f, p0.z), 0, 0 });
        auto top = [&](Int3 c) { return Vector3{ (c.x + 0.5f) * kCellSize, c.y * kCellSize + 0.05f, (c.z + 0.5f) * kCellSize }; };
        plate0 = top(rubbleCell); plate1 = top(cell(p0.x + 5.0f, p0.z)); vent = top(cell(p0.x + 9.0f, p0.z));
        sv.health = 100.0f; sv.oxygen = 100.0f;
        DungeonTraps_Reset();
        return 0;
    }
    plot = body.plots ? body.plots->byId(plotId) : nullptr;
    if (!plot) { printf("bench: dungeon traps: the plot is gone -> FAIL\n"); return -1; }
    // The schedule in sim seconds (the traps run on the sim clock; frames come at any rate): the block plate at 0.5 s,
    // the dart plate at 2 s, beside the vent at 4 s, done at 6.5 s.
    static double t0 = 0.0;
    static int step = 0;
    if (frame == 2) { t0 = SimTime(); step = 0; }
    const double t = SimTime() - t0;
    // The camera: from the trench's near end, 1.8 m over its floor, at the trap in play.
    const Vector3 subject = t < 1.9 ? plate0 : t < 3.9 ? plate1 : vent;
    world.camera.position = toRender({ subject.x - 3.6f, fy + 2.1f, p0.z - 0.4f });
    world.camera.target = toRender({ subject.x, fy + 1.1f, p0.z });
    world.camera.up = body.dirToWorld(plot->up);
    if (step == 0 && t >= 0.5) { stand(plate0); health0 = sv.health; ++step; }
    if (step == 1 && t >= 0.85) { TakeScreenshot("bench_dungeon_traps.png"); printf("bench: dungeon traps: bench_dungeon_traps.png (the block falling)\n"); ++step; }
    if (step == 2 && t >= 1.8) { healthBlock = sv.health; ++step; }
    if (step == 3 && t >= 2.0) { stand(plate1); ++step; }
    if (step == 4 && t >= 3.9) { airBefore = sv.breathable && !sv.inGas; ++step; }
    if (step == 5 && t >= 4.0) { stand({ vent.x + 0.6f, vent.y, vent.z }); oxygen0 = sv.oxygen; ++step; }
    if (step == 6 && t > 4.3) { gasSeen = gasSeen || sv.inGas; breathedInGas = breathedInGas || (sv.inGas && sv.breathable); }
    if (step == 6 && t >= 6.5) {
        const DungeonTrapStats& s = DungeonTraps_Stats();
        const int at = plot->partAt(rubbleCell);
        const bool spent = at < 0 || plot->parts[at].type != PartType::PressurePlate; // (it broke on the player: no rubble)
        const bool ok = s.blocksDropped == 1 && s.blockHits == 1 && healthBlock < health0 - 20.0f && spent && s.dartsShot >= 1 && s.dartHits >= 1 &&
                        airBefore && gasSeen && !breathedInGas && sv.oxygen < oxygen0 - 1.0f;
        printf("bench: dungeon traps (#719): block dropped %d, hit %d (health %.0f -> %.0f), plate spent %d; darts shot %d, hit %d; "
               "air breathable before %d; by the vent: gas %d, breathable there %d, oxygen %.1f -> %.1f; health now %.0f -> %s\n",
               s.blocksDropped, s.blockHits, health0, healthBlock, spent, s.dartsShot, s.dartHits, airBefore, gasSeen, breathedInGas, oxygen0, sv.oxygen, sv.health, ok ? "PASS" : "FAIL");
        return ok ? 1 : -1;
    }
    return 0;
}
