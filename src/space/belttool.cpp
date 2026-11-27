#include "belttool.h"
#include "audio.h"
#include "buildtool.h"
#include "gameui.h"
#include "inventory.h"
#include "player.h"
#include "playerinput.h"
#include "plots.h"
#include "researchgate.h"
#include "terraintool.h"
#include "uipalette.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <map>

namespace {
constexpr const char* kBeltMaterial = "basic";
constexpr int kMaxRunSlots = 32 * 8;   // GUESS: a drag reaches at most 32 cells (16 m) along each axis
constexpr float kDragReach = 24.0f;    // m: how far the drag's floor plane is followed along the view

int Get(CInt3 p, int axis) { return axis == 0 ? p.x : axis == 1 ? p.y : p.z; }
void Set(CInt3& p, int axis, int v) { (axis == 0 ? p.x : axis == 1 ? p.y : p.z) = v; }
int EvenSlot(float plotMetres) { return 2 * (int)lroundf(plotMetres / kCellSize * 4.0f); } // the nearest even slot
int FloorDiv8(int s) { return s >= 0 ? s / 8 : -((-s + 7) / 8); }

bool MachineType(PartType t) { return IsMachine(t) || t == PartType::MbController || t == PartType::MbHatch; } // as the factory's layout (factory.cpp)
int MachineAt(const Plot& plot, Int3 cell) {
    const int i = plot.partAt(cell);
    return i >= 0 && i < (int)plot.parts.size() && MachineType(plot.parts[i].type) ? i : -1;
}

// The tool's state (the local player's: a remote player's placing arrives as a finished run in a later PR).
bool beltMode = true;               // false: the conveyor item places the old cell conveyor
struct Drag {
    bool on = false, bench = false;
    uint64_t body = 0;
    uint32_t plot = 0;
    CInt3 a, b;
    bool xFirst = true, reversed = false;
} drag;
struct Hover {
    Planetoid* body = nullptr;
    Plot* plot = nullptr;
    CInt3 slot;
    bool on = false;
    int belt = 0;                   // the belt under the crosshair (LMB takes it up)
} hover;
ConduitRun planned;
ConduitPlaceResult plannedCheck;
Plot* planPlot = nullptr;
Planetoid* planBody = nullptr;

Plot* PlotById(Planetoid* body, uint32_t id) { return body && body->plots ? body->plots->byId(id) : nullptr; }
} // namespace

// ---------------------------------------------------------------------------------------------
// Pure parts
// ---------------------------------------------------------------------------------------------

std::vector<CInt3> BeltTool_Route(CInt3 a, CInt3 b, bool xFirst, int* firstAxis, int* lastAxis) {
    if (a == b) return {};
    b.y = a.y; // one floor layer
    if (a.x == b.x || a.z == b.z) {
        const int axis = a.x == b.x ? 2 : 0;
        if (firstAxis) *firstAxis = axis;
        if (lastAxis) *lastAxis = axis;
        return { a, b };
    }
    if (firstAxis) *firstAxis = xFirst ? 0 : 2;
    if (lastAxis) *lastAxis = xFirst ? 2 : 0;
    const CInt3 corner = xFirst ? CInt3{ b.x, a.y, a.z } : CInt3{ a.x, a.y, b.z };
    return { a, corner, b };
}

CInt3 BeltTool_SnapEnd(const Plot& plot, CInt3 p, int axis, int sign, bool start) {
    Int3 cell{ FloorDiv8(p.x), FloorDiv8(p.y), FloorDiv8(p.z) };
    const int m = MachineAt(plot, cell);
    if (m < 0) return p;
    // The machine's outermost cell along the travel: an end walks back against it (the face it was entered by), a start
    // forward with it (the face it leaves by).
    const int step = start ? sign : -sign;
    for (;;) {
        Int3 n = cell;
        (axis == 0 ? n.x : n.z) += step;
        if (plot.partAt(n) != m) break;
        cell = n;
    }
    const int across = axis == 0 ? 2 : 0;
    CInt3 q{ cell.x * 8, cell.y * 8, cell.z * 8 };
    // A start leaves by the face toward +sign; an end enters by the face toward -sign. Slot 0 is the -axis face's layer,
    // 6 the even slot by the +axis face's layer 7 (Conduit_BeltEndAtPort allows one slot off along the axis).
    const bool plusFace = start ? sign > 0 : sign < 0;
    Set(q, axis, Get(q, axis) + (plusFace ? 6 : 0));
    Set(q, across, Get(q, across) + 4); // the face centre across (u = 4: slot 3 or 4, one off is allowed)
    return q;
}

std::vector<CInt3> BeltTool_Plan(const Plot& plot, CInt3 a, CInt3 b, bool xFirst) {
    int fa = 0, la = 0;
    std::vector<CInt3> pts = BeltTool_Route(a, b, xFirst, &fa, &la);
    if (pts.size() < 2) return pts;
    const size_t n = pts.size();
    const int s0 = Get(pts[1], fa) > Get(pts[0], fa) ? 1 : -1;
    const int s1 = Get(pts[n - 1], la) > Get(pts[n - 2], la) ? 1 : -1;
    CInt3 a2 = BeltTool_SnapEnd(plot, pts[0], fa, s0, true);
    CInt3 b2 = BeltTool_SnapEnd(plot, pts[n - 1], la, s1, false);
    const bool aSnapped = !(a2 == pts[0]), bSnapped = !(b2 == pts[n - 1]);
    if (!aSnapped && !bSnapped) return pts;
    if (b2.y != a2.y) { // a machine on another floor layer: the start's layer for the free end, else the snapped ends as they are
        if (!aSnapped) a2.y = b2.y;
        else if (!bSnapped) b2.y = a2.y;
    }
    if (fa == la) { // a straight run: the free end follows the port end across (two snapped ends in one row share it)
        const int across = fa == 0 ? 2 : 0;
        if (aSnapped && !bSnapped) Set(b2, across, Get(a2, across));
        if (bSnapped && !aSnapped) Set(a2, across, Get(b2, across));
        if (a2 == b2) return {};
        return { a2, b2 };
    }
    std::vector<CInt3> out = BeltTool_Route(a2, b2, fa == 0);
    if (out.size() >= 2) out.back().y = b2.y; // (Route puts b on a's layer; a layer mismatch stays and the check says so)
    return out;
}

std::vector<ConduitPortAt> BeltTool_Ports(const Plot& plot) {
    std::vector<ConduitPortAt> ports;
    const Int3 side[4] = { { 0, 0, -1 }, { 1, 0, 0 }, { 0, 0, 1 }, { -1, 0, 0 } };
    std::vector<std::pair<Int3, CellShape>> fp;
    for (size_t i = 0; i < plot.parts.size(); ++i) {
        const Part& part = plot.parts[i];
        if (part.type == PartType::Count || !MachineType(part.type)) continue;
        Plot::Footprint(part, fp);
        for (const auto& [c, shape] : fp)
            for (int f = 0; f < 4; ++f) {
                const Int3 n{ c.x + side[f].x, c.y, c.z + side[f].z };
                if (plot.partAt(n) == (int)i) continue; // an inner face
                ConduitPortDef def; def.face = f; def.kind = ConduitPortKind::ItemIn;
                ConduitPortAt p; p.slot = Conduit_PortSlot({ c.x, c.y, c.z }, 0, def); p.face = f; p.kind = def.kind;
                ports.push_back(p);
            }
    }
    return ports;
}

ConduitPlaceResult BeltTool_Check(const Plot& plot, const ConduitRun& run, const std::function<bool(CInt3 cell)>& ground) {
    ConduitOccupancy occ;
    std::map<int, ConduitRun> runs;
    for (const ConduitRun& r : plot.belts) { occ.add(r); runs[r.id] = r; }
    std::map<CInt3, bool> solid; // per cell: Conduit_CheckPlace asks per slot
    auto solidCell = [&](CInt3 c) {
        const auto it = solid.find(c);
        if (it != solid.end()) return it->second;
        const bool s = plot.partAt({ c.x, c.y, c.z }) >= 0 || (ground && ground(c));
        solid[c] = s;
        return s;
    };
    ConduitPlaceResult res = Conduit_CheckPlace(run, occ, runs, solidCell, BeltTool_Ports(plot));
    if (res.ok) {
        for (int axis : { 0, 2 }) {
            int lo = INT32_MAX, hi = INT32_MIN;
            for (const CInt3& p : run.points) { lo = std::min(lo, Get(p, axis)); hi = std::max(hi, Get(p, axis)); }
            if (hi - lo > kMaxRunSlots) { res.ok = false; res.why = "too long (32 cells at most)"; }
        }
    } else if (res.why.empty() && !res.blocked.empty()) {
        res.why = "blocked: it needs open floor with something solid under it";
    }
    return res;
}

std::function<bool(CInt3 cell)> BeltTool_Ground(const Planetoid& body, const Plot& plot) {
    return [&body, &plot](CInt3 c) {
        if (!body.field) return false;
        const Vector3 local = plot.toBody({ (c.x + 0.5f) * kCellSize, (c.y + 0.5f) * kCellSize, (c.z + 0.5f) * kCellSize });
        return body.field->density(local, body.edits.get()) >= ISO_LEVEL; // (float is fine: a solid test, 4 mm steps)
    };
}

int BeltTool_RunAt(const Plot& plot, CInt3 slot) {
    for (const ConduitRun& r : plot.belts)
        for (const CInt3& s : Conduit_Slots(r))
            if (s == slot || (s.x == slot.x && s.z == slot.z && s.y == slot.y + 1)) return r.id;
    return 0;
}

int BeltTool_Cost(const ConduitRun& run) {
    const ConduitMaterial* mat = Conduit_Materials().find(run.material.empty() ? kBeltMaterial : run.material);
    return mat ? Conduit_Cost(run, *mat) : 0;
}

static const char* BeltItem(const ConduitRun& run) {
    const ConduitMaterial* mat = Conduit_Materials().find(run.material.empty() ? kBeltMaterial : run.material);
    return mat ? mat->item.c_str() : "conveyor";
}

bool BeltTool_Place(Planetoid& body, Plot& plot, ConduitRun run, bool charge, std::string* why) {
    run.kind = ConduitKind::Belt;
    run.size = kBeltWidth;
    if (run.material.empty()) run.material = kBeltMaterial;
    const ConduitPlaceResult res = BeltTool_Check(plot, run, BeltTool_Ground(body, plot));
    if (!res.ok) { if (why) *why = res.why; return false; }
    if (charge && ResearchGate_On() && !terrainTool.creative && !ResearchGate_PartAllowed(ResearchGate_Team(), PartType::Conveyor)) {
        const ResearchNode* n = ResearchGate_NodeForPart(PartType::Conveyor);
        if (why) *why = std::string("research ") + (n ? n->name : "conveyors") + " first";
        return false;
    }
    const int cost = BeltTool_Cost(run);
    const char* item = BeltItem(run);
    if (charge && !terrainTool.creative && !inventory.remove(item, cost)) {
        if (why) *why = TextFormat("needs %d %s (have %d)", cost, item, inventory.count(item));
        return false;
    }
    if (plot.addBelt(std::move(run)) == 0) { // (the check passed Conduit_ValidRun: not reached)
        if (charge && !terrainTool.creative) inventory.add(item, cost);
        if (why) *why = "not a valid belt";
        return false;
    }
    return true;
}

bool BeltTool_Remove(Planetoid& body, Plot& plot, int id, bool refund) {
    (void)body;
    const auto it = std::find_if(plot.belts.begin(), plot.belts.end(), [id](const ConduitRun& b) { return b.id == id; });
    if (it == plot.belts.end()) return false;
    const int cost = BeltTool_Cost(*it);
    const std::string item = BeltItem(*it);
    if (!plot.removeBelt(id)) return false;
    // The items riding it come back through the factory's refunds (a removed run's items go to the plot's refundTo).
    if (refund && !terrainTool.creative) inventory.add(item, cost);
    return true;
}

void BeltTool_BenchDrag(Plot* plot, CInt3 a, CInt3 b, bool xFirst) {
    if (!plot) { drag = Drag{}; return; }
    drag.on = true; drag.bench = true;
    drag.plot = plot->id; drag.a = a; drag.b = b; drag.xFirst = xFirst; drag.reversed = false;
    drag.body = 0; // the update finds the body that holds the plot
    planPlot = plot;
}
const ConduitRun& BeltTool_Planned() { return planned; }
const ConduitPlaceResult& BeltTool_PlannedCheck() { return plannedCheck; }

// ---------------------------------------------------------------------------------------------
// The tool
// ---------------------------------------------------------------------------------------------

bool BeltTool_On() { return buildTool.active && beltMode && buildTool.part == (int)PartType::Conveyor; }

void BeltTool_ForgetWorld() {
    drag = Drag{}; hover = Hover{}; planned = ConduitRun{}; plannedCheck = ConduitPlaceResult{};
    planPlot = nullptr; planBody = nullptr;
}

void BeltTool_Cancel() {
    if (drag.bench) return;
    drag = Drag{}; planned = ConduitRun{}; plannedCheck = ConduitPlaceResult{};
}

namespace {
// The floor point the crosshair is on: a part's top (the cell above), a machine's side (inside that cell, at its floor), any
// other part's side (the cell in front), or the ground by a plot (the cell whose floor is nearest the ground).
bool AimFloor(Scene& world, Planetoid& body, const DVec3& lo, Vector3 ld, Vector3 dir, Plot*& plotOut, CInt3& slot) {
    float best = buildTool.reach;
    Plot* hitPlot = nullptr;
    Int3 hitCell{}, hitFace{};
    if (body.plots)
        for (auto& pl : body.plots->plots) {
            Int3 c, f; float d;
            if (pl->raycast(lo, ld, best, c, f, d)) { best = d; hitPlot = pl.get(); hitCell = c; hitFace = f; }
        }
    TerrainHit th;
    const bool terrain = RaycastTerrain(world, Ray{ world.camera.position, dir }, buildTool.reach, th) && th.planetoid == &body;
    if (hitPlot && (!terrain || best <= th.distance)) {
        const Vector3 h = hitPlot->toPlot(lo + DVec3(Vector3Scale(ld, best)));
        Int3 cell = hitCell;
        if (hitFace.y > 0) cell.y += 1;                                         // on its top
        else if (hitFace.y == 0 && MachineAt(*hitPlot, hitCell) < 0) { cell.x += hitFace.x; cell.z += hitFace.z; } // in front of a wall
        else if (hitFace.y < 0) return false;                                   // under a part
        slot = { EvenSlot(h.x), cell.y * 8, EvenSlot(h.z) };
        if (hitFace.x != 0) slot.x = cell.x * 8 + (hitFace.x > 0 ? (MachineAt(*hitPlot, hitCell) >= 0 ? 6 : 0) : (MachineAt(*hitPlot, hitCell) >= 0 ? 0 : 6));
        if (hitFace.z != 0) slot.z = cell.z * 8 + (hitFace.z > 0 ? (MachineAt(*hitPlot, hitCell) >= 0 ? 6 : 0) : (MachineAt(*hitPlot, hitCell) >= 0 ? 0 : 6));
        plotOut = hitPlot;
        return true;
    }
    if (!terrain || !body.plots) return false;
    Plot* near = body.plots->plotFor(th.local, 6.0f, false);
    if (!near) return false;
    const Vector3 h = near->toPlot(th.local);
    slot = { EvenSlot(h.x), 8 * (int)lroundf(h.y / kCellSize), EvenSlot(h.z) };
    plotOut = near;
    return true;
}

// Where the drag's far end is: on a machine the crosshair is on, else where the view meets the start's floor layer.
bool DragEnd(Scene& world, Planetoid& body, Plot& plot, const DVec3& lo, Vector3 ld, Vector3 dir, CInt3& out) {
    Plot* p = nullptr;
    CInt3 s;
    if (AimFloor(world, body, lo, ld, dir, p, s) && p == &plot && MachineAt(plot, { FloorDiv8(s.x), FloorDiv8(s.y), FloorDiv8(s.z) }) >= 0) { out = s; return true; }
    const Vector3 o = plot.toPlot(lo);
    const Vector3 d{ Vector3DotProduct(ld, plot.east), Vector3DotProduct(ld, plot.up), Vector3DotProduct(ld, plot.north) };
    const float y = drag.a.y / 8.0f * kCellSize;
    if (fabsf(d.y) < 1e-4f) return false;
    const float t = (y - o.y) / d.y;
    if (t <= 0.0f || t > kDragReach) return false;
    out = { EvenSlot(o.x + d.x * t), drag.a.y, EvenSlot(o.z + d.z * t) };
    for (int axis : { 0, 2 }) // clamp to the run's reach
        Set(out, axis, std::clamp(Get(out, axis), Get(drag.a, axis) - kMaxRunSlots, Get(drag.a, axis) + kMaxRunSlots));
    return true;
}

void Replan() {
    planned = ConduitRun{};
    plannedCheck = ConduitPlaceResult{};
    if (!drag.on || !planPlot) return;
    planned.kind = ConduitKind::Belt; planned.size = kBeltWidth; planned.material = kBeltMaterial;
    planned.points = drag.reversed ? BeltTool_Plan(*planPlot, drag.b, drag.a, drag.xFirst) : BeltTool_Plan(*planPlot, drag.a, drag.b, drag.xFirst);
    if (planned.points.size() < 2) { plannedCheck.why = "drag to lay a belt"; return; }
    plannedCheck = planBody ? BeltTool_Check(*planPlot, planned, BeltTool_Ground(*planBody, *planPlot)) : BeltTool_Check(*planPlot, planned, nullptr);
}
} // namespace

bool BeltTool_Update(Scene& world, const PlayerInput& in) {
    // Ctrl/Alt+wheel on the conveyor: belts or the old cell conveyor (they work side by side until conveyors convert, #2052).
    if (buildTool.active && buildTool.part == (int)PartType::Conveyor && in.enabled && in.wheel != 0.0f
        && (IsKeyDown(KEY_LEFT_CONTROL) || IsKeyDown(KEY_RIGHT_CONTROL) || IsKeyDown(KEY_LEFT_ALT) || IsKeyDown(KEY_RIGHT_ALT))) {
        beltMode = !beltMode;
        drag = Drag{};
        GameUI_Toast(beltMode ? "Belts: drag RMB to lay a run" : "Cell conveyors (the old kind): RMB places one");
    }
    hover = Hover{};
    if (drag.bench) { // a scripted drag: find its plot's body, plan, nothing else
        for (auto& o : world.objects)
            if (auto* b = dynamic_cast<Planetoid*>(o.get()); b && b->plots && b->plots->byId(drag.plot) == planPlot) { planBody = b; break; }
        Replan();
        return false;
    }
    if (!BeltTool_On()) { drag = Drag{}; planned = ConduitRun{}; return false; }
    Planetoid* body = player.nearest;
    if (!body) { drag = Drag{}; return false; }
    Vector3 dir, ld; DVec3 lo;
    player.aimRay(world, *body, lo, ld, dir);
    Plot* p = nullptr;
    CInt3 s;
    if (AimFloor(world, *body, lo, ld, dir, p, s)) {
        hover.body = body; hover.plot = p; hover.slot = s; hover.on = true;
        hover.belt = BeltTool_RunAt(*p, s);
    }
    if (drag.on) {
        Plot* dp = drag.body == body->id ? PlotById(body, drag.plot) : nullptr;
        if (!dp) { drag = Drag{}; return false; }
        planPlot = dp; planBody = body;
        if (in.hit(Act::Primary)) { drag = Drag{}; planned = ConduitRun{}; GameUI_Toast("Belt cancelled"); return true; }
        if (in.hit(Act::Rotate)) drag.xFirst = !drag.xFirst;
        if (in.hit(Act::Flip)) drag.reversed = !drag.reversed;
        CInt3 e;
        if (DragEnd(world, *body, *dp, lo, ld, dir, e)) drag.b = e;
        Replan();
        if (!in.down(Act::Secondary)) { // released: lay it
            const ConduitRun run = planned;
            const ConduitPlaceResult check = plannedCheck;
            drag = Drag{};
            planned = ConduitRun{};
            if (run.points.size() < 2) return true;
            if (gBuildToolNetEdit) { GameUI_Toast("Belts can't be placed on a server yet"); Audio_Play(Sfx::Error, 0.5f); return true; }
            std::string why;
            if (check.ok && BeltTool_Place(*body, *dp, run, true, &why)) {
                Audio_Play(Sfx::Build, 0.5f, 1.2f);
                GameUI_Toast(TextFormat("Belt laid: %.1f cells", Conduit_Length(run) / 8.0f));
            } else {
                GameUI_Toast(TextFormat("Can't lay the belt: %s", check.ok ? why.c_str() : check.why.c_str()));
                Audio_Play(Sfx::Error, 0.5f);
            }
        }
        return true;
    }
    if (!in.enabled) return hover.belt != 0;
    if (hover.on && in.hit(Act::Secondary)) {
        drag = Drag{};
        drag.on = true; drag.body = body->id; drag.plot = hover.plot->id; drag.a = hover.slot; drag.b = hover.slot;
        planPlot = hover.plot; planBody = body;
        return true;
    }
    if (hover.belt && in.hit(Act::Primary)) {
        if (gBuildToolNetEdit) { GameUI_Toast("Belts can't be taken up on a server yet"); Audio_Play(Sfx::Error, 0.5f); return true; }
        if (BeltTool_Remove(*body, *hover.plot, hover.belt, true)) Audio_Play(Sfx::Dig, 0.5f, 1.3f);
        hover.belt = 0;
        return true;
    }
    return hover.belt != 0;
}

// ---------------------------------------------------------------------------------------------
// Drawing: the run being dragged (green, or red with crosses), the belts already laid (outlined, so build mode shows them
// until belts are drawn, #2053), the one under the crosshair (orange: LMB takes it up), and the flow arrows on each.
// ---------------------------------------------------------------------------------------------

namespace {
constexpr float kSlotM = kCellSize / 8.0f; // m

void Box(Vector3 lo, Vector3 hi, Color fill, Color edge, bool cross) {
    const Vector3 v[8] = { { lo.x, lo.y, lo.z }, { hi.x, lo.y, lo.z }, { lo.x, hi.y, lo.z }, { hi.x, hi.y, lo.z },
                           { lo.x, lo.y, hi.z }, { hi.x, lo.y, hi.z }, { lo.x, hi.y, hi.z }, { hi.x, hi.y, hi.z } };
    if (fill.a > 0) {
        const int q[6][4] = { { 0, 1, 3, 2 }, { 4, 6, 7, 5 }, { 0, 2, 6, 4 }, { 1, 5, 7, 3 }, { 2, 3, 7, 6 }, { 0, 4, 5, 1 } };
        for (const auto& f : q) { DrawTriangle3D(v[f[0]], v[f[1]], v[f[2]], fill); DrawTriangle3D(v[f[0]], v[f[2]], v[f[3]], fill); }
    }
    const int e[12][2] = { { 0, 1 }, { 2, 3 }, { 4, 5 }, { 6, 7 }, { 0, 2 }, { 1, 3 }, { 4, 6 }, { 5, 7 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
    for (const auto& ed : e) DrawLine3D(v[ed[0]], v[ed[1]], edge);
    if (cross) { DrawLine3D(v[2], v[7], edge); DrawLine3D(v[3], v[6], edge); } // the top's diagonals: a cue beyond the colour (#1432)
}

// One run under the plot's model: a box per segment (the belt's 4 x 2 slot cross-section, grown 1 cm), chevrons on top.
void DrawRun(const ConduitRun& run, Color fill, Color edge, bool cross) {
    constexpr float kPad = 0.01f;
    for (size_t i = 0; i + 1 < run.points.size(); ++i) {
        const CInt3 a = run.points[i], b = run.points[i + 1];
        const int axis = a.x != b.x ? 0 : 2, across = axis == 0 ? 2 : 0;
        CInt3 lo = a, hi = a;
        Set(lo, axis, std::min(Get(a, axis), Get(b, axis)));
        Set(hi, axis, std::max(Get(a, axis), Get(b, axis)) + 1);
        Set(lo, across, Get(a, across) - kBeltWidth / 2);
        Set(hi, across, Get(a, across) + kBeltWidth / 2);
        hi.y = a.y + kBeltHeight;
        const Vector3 l{ lo.x * kSlotM - kPad, lo.y * kSlotM + kPad, lo.z * kSlotM - kPad };
        const Vector3 h{ hi.x * kSlotM + kPad, (a.y + 1) * kSlotM + kPad, hi.z * kSlotM + kPad }; // the belt's own layer (items ride above)
        Box(l, h, fill, edge, cross);
        // Chevrons every half cell along the flow, on the belt's top.
        const int sign = Get(b, axis) > Get(a, axis) ? 1 : -1;
        const int len = std::abs(Get(b, axis) - Get(a, axis));
        for (int t = 2; t + 2 <= len; t += 4) {
            CInt3 c = a;
            Set(c, axis, Get(a, axis) + sign * t);
            auto at = [&](float along, float side) {
                Vector3 v{ (c.x + 0.5f) * kSlotM, h.y + 0.005f, (c.z + 0.5f) * kSlotM };
                (axis == 0 ? v.x : v.z) += along * kSlotM * sign;
                (across == 0 ? v.x : v.z) += side * kSlotM;
                return v;
            };
            DrawLine3D(at(-0.8f, -1.4f), at(0.8f, 0.0f), edge);
            DrawLine3D(at(0.8f, 0.0f), at(-0.8f, 1.4f), edge);
        }
    }
}
} // namespace

void BeltTool_Draw3D(Scene& world) {
    if (!BeltTool_On() && !drag.bench) return;
    Planetoid* body = drag.bench ? planBody : player.nearest;
    if (!body || !body->plots) return;
    auto under = [&](const Plot& p, const std::function<void()>& draw) {
        rlDrawRenderBatchActive();
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(p.cellModel(*body, world.renderOrigin))); // plot coordinates in double (#271)
        rlDisableBackfaceCulling();
        rlDisableDepthMask();
        draw();
        rlDrawRenderBatchActive();
        rlEnableDepthMask();
        rlEnableBackfaceCulling();
        rlPopMatrix();
    };
    const Plot* hovered = hover.on ? hover.plot : nullptr;
    for (const auto& pl : body->plots->plots) {
        if (pl->belts.empty()) continue;
        under(*pl, [&] {
            for (const ConduitRun& r : pl->belts) {
                const bool aimed = pl.get() == hovered && r.id == hover.belt && !drag.on;
                DrawRun(r, aimed ? Fade(ORANGE, 0.35f) : Fade(SKYBLUE, 0.18f), aimed ? ORANGE : Fade(SKYBLUE, 0.8f), false);
            }
        });
    }
    if (drag.on && planPlot && planned.points.size() >= 2) {
        const bool ok = plannedCheck.ok;
        under(*planPlot, [&] { DrawRun(planned, ok ? Fade(Ui().good, 0.30f) : Fade(Ui().bad, 0.32f), ok ? Ui().good : Ui().bad, !ok); });
    } else if (!drag.on && hover.on && !hover.belt) { // where a drag would start: one belt square
        under(*hover.plot, [&] {
            const Vector3 l{ (hover.slot.x - 2) * kSlotM, hover.slot.y * kSlotM + 0.01f, (hover.slot.z - 2) * kSlotM };
            const Vector3 h{ (hover.slot.x + 2) * kSlotM, (hover.slot.y + 1) * kSlotM + 0.01f, (hover.slot.z + 2) * kSlotM };
            Box(l, h, Fade(Ui().good, 0.25f), Ui().good, false);
        });
    }
}

void BeltTool_DrawHUD() {
    const int cx = GetScreenWidth() / 2;
    const char* item = "conveyor";
    if (const ConduitMaterial* mat = Conduit_Materials().find(kBeltMaterial)) item = mat->item.c_str();
    const char* line;
    Color col = SKYBLUE;
    if (drag.on && planned.points.size() >= 2) {
        const int cost = BeltTool_Cost(planned);
        line = TextFormat("BELT  %.1f cells  (%d %s, have %d)%s%s", Conduit_Length(planned) / 8.0f, cost, item, inventory.count(item),
                          plannedCheck.ok ? "" : "  -  ", plannedCheck.ok ? "" : plannedCheck.why.c_str());
        if (!plannedCheck.ok) col = Ui().bad;
    } else {
        line = TextFormat("BELT  (1 %s a cell, have %d)", item, inventory.count(item));
    }
    DrawText(line, cx - MeasureText(line, 22) / 2, GetScreenHeight() - 150, 22, col);
    DrawText(drag.on ? "release RMB lay   R corner   V reverse   LMB cancel"
                     : "drag RMB lay a belt   LMB on a belt take it up   Ctrl+wheel cell conveyor   G leave",
             10, GetScreenHeight() - 30, 20, Fade(RAYWHITE, 0.7f));
}

// ---------------------------------------------------------------------------------------------
#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>

TEST(BeltTool, RoutesStraightOrWithOneCorner) {
    EXPECT_TRUE(BeltTool_Route({ 0, 0, 0 }, { 0, 0, 0 }, true).empty());
    EXPECT_EQ(BeltTool_Route({ 0, 0, 4 }, { 20, 0, 4 }, true), (std::vector<CInt3>{ { 0, 0, 4 }, { 20, 0, 4 } }));
    EXPECT_EQ(BeltTool_Route({ 0, 0, 0 }, { 8, 0, 6 }, true), (std::vector<CInt3>{ { 0, 0, 0 }, { 8, 0, 0 }, { 8, 0, 6 } }));
    EXPECT_EQ(BeltTool_Route({ 0, 0, 0 }, { 8, 0, 6 }, false), (std::vector<CInt3>{ { 0, 0, 0 }, { 0, 0, 6 }, { 8, 0, 6 } }));
    EXPECT_EQ(BeltTool_Route({ 0, 8, 0 }, { 8, 0, 0 }, true), (std::vector<CInt3>{ { 0, 8, 0 }, { 8, 8, 0 } })) << "one floor layer: the start's";
}

namespace {
// Two chests (2 x 2 x 2 cells) on a plot floor of cubes (y -1, z -3 .. 5), at cells x 0 and x 6 (z 0).
struct TwoChests {
    Plot plot;
    int a = -1, b = -1, xa = 0, xb = 6;
    TwoChests() {
        for (int x = -2; x <= 10; ++x) for (int z = -3; z <= 5; ++z) plot.place({ PartType::Cube, { x, -1, z }, 0, 0 });
        a = plot.place({ PartType::Chest, { 0, 0, 0 }, 0 });
        b = plot.place({ PartType::Chest, { 6, 0, 0 }, 0 });
        while (plot.partAt({ xa + 1, 0, 0 }) == a) ++xa;
        while (xb > 0 && plot.partAt({ xb - 1, 0, 0 }) == b) --xb;
    }
};
}

TEST(BeltTool, EndsOnChestsSnapToTheirFacePorts) {
    TwoChests t;
    ASSERT_GE(t.a, 0); ASSERT_GE(t.b, 0);
    // Dragged from somewhere in A to somewhere in B, off the face centres: both ends move to the faces they cross.
    const std::vector<CInt3> pts = BeltTool_Plan(t.plot, { 2, 0, 2 }, { t.xb * 8 + 4, 0, 2 }, true);
    ASSERT_EQ(pts.size(), 2u);
    EXPECT_EQ(pts[0], (CInt3{ t.xa * 8 + 6, 0, 4 })) << "A's east face, the face centre across";
    EXPECT_EQ(pts[1], (CInt3{ t.xb * 8, 0, 4 })) << "B's west face";
    ConduitRun run; run.kind = ConduitKind::Belt; run.size = kBeltWidth; run.material = "basic"; run.points = pts;
    const ConduitPlaceResult res = BeltTool_Check(t.plot, run, nullptr);
    EXPECT_TRUE(res.ok) << res.why;
    EXPECT_EQ(res.portIdx.size(), 2u) << "both ends on a port";
    // Reversed (B to A): the same faces, the other way.
    const std::vector<CInt3> back = BeltTool_Plan(t.plot, { t.xb * 8 + 4, 0, 2 }, { 2, 0, 4 }, true);
    ASSERT_EQ(back.size(), 2u);
    EXPECT_EQ(back[0], (CInt3{ t.xb * 8, 0, 4 }));
    EXPECT_EQ(back[1], (CInt3{ t.xa * 8 + 6, 0, 4 }));
}

TEST(BeltTool, CheckRefusesFloatingCrossingAndPartCells) {
    TwoChests t;
    ConduitRun run; run.kind = ConduitKind::Belt; run.size = kBeltWidth; run.material = "basic";
    run.points = { { 0, 0, 24 }, { 40, 0, 24 } }; // z cells 2-3: on the cube floor, clear of the chests (z 0-1)
    EXPECT_TRUE(BeltTool_Check(t.plot, run, nullptr).ok);
    run.points = { { 0, 0, 56 }, { 40, 0, 56 } }; // z cell 7: past the floor's edge (z 5), nothing under it
    EXPECT_FALSE(BeltTool_Check(t.plot, run, nullptr).ok);
    EXPECT_TRUE(BeltTool_Check(t.plot, run, [](CInt3 c) { return c.y < 0; }).ok) << "the ground under it holds it up";
    run.points = { { 0, 8, 24 }, { 40, 8, 24 } }; // one layer up: nothing under it but air
    EXPECT_FALSE(BeltTool_Check(t.plot, run, nullptr).ok);
    run.points = { { -8, 0, 4 }, { 40, 0, 4 } }; // straight through chest A
    EXPECT_FALSE(BeltTool_Check(t.plot, run, nullptr).ok);
    ConduitRun laid; laid.material = "basic"; laid.points = { { 0, 0, 24 }, { 40, 0, 24 } };
    const int id = t.plot.addBelt(laid);
    ASSERT_GT(id, 0);
    run.points = { { 20, 0, 8 }, { 20, 0, 40 } }; // across the laid belt's middle
    EXPECT_FALSE(BeltTool_Check(t.plot, run, nullptr).ok);
    run.points = { { 20, 0, 8 }, { 20, 0, 24 } }; // ending on its centreline: merges into it
    const ConduitPlaceResult merge = BeltTool_Check(t.plot, run, nullptr);
    EXPECT_TRUE(merge.ok) << merge.why;
    EXPECT_EQ(merge.joinsRuns, std::vector<int>{ id });
    EXPECT_EQ(BeltTool_RunAt(t.plot, { 20, 0, 24 }), id);
    EXPECT_EQ(BeltTool_RunAt(t.plot, { 20, 0, 32 }), 0);
    run.points = { { 0, 0, -16 }, { 400, 0, -16 } };
    EXPECT_FALSE(BeltTool_Check(t.plot, run, [](CInt3 c) { return c.y < 0; }).ok) << "longer than 32 cells";
}

#endif
