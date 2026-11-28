#include "beltmodels.h"
#include "beltsim.h"
#include "conduit.h"
#include "factory.h"
#include "instancebatch.h"
#include "itemicons.h"
#include "plots.h"
#include "scene.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

// Plot coordinates: x east, y up, z north, in metres from the plot origin (Plot::cellModel); conduit slots are 1/8 cell.
// A belt run's points are its centreline on a floor layer: the belt is kBeltWidth slots across, centred on it.
constexpr float kSlot = kCellSize / 8.0f;          // m
constexpr float kHalfW = kBeltWidth * 0.5f;        // slots
constexpr float kTop = 0.85f;                      // the belt surface over its floor layer (slots): just under 1/8 cell
constexpr float kRailW = 0.3f, kRailH = 1.1f;      // the side rails (slots): a lip the items sit between
constexpr float kChevronEvery = 4.0f;              // slots: a chevron every half cell
constexpr float kItemSize = 0.09f;                 // m: items ride 1/4 cell (0.125 m) apart at the closest
constexpr float kDrawRange = 96.0f, kItemRange = 48.0f, kIconRange = 20.0f; // m from the camera
const Color kSlab{ 44, 44, 48, 255 }, kRail{ 118, 122, 128, 255 }, kChevron{ 214, 176, 48, 255 };

BeltDrawStats stats;

struct Decal { Vector3 c; std::string id; float light; };
struct PlotDecals { Matrix model; std::vector<Decal> decals; };
std::vector<PlotDecals> decalPlots; // this frame's, for Belts_DrawIcons

struct PlotState { uint64_t version = ~0ull; float buried = 0.0f; };
std::unordered_map<uint64_t, PlotState> plotState;

int Get(CInt3 p, int axis) { return axis == 0 ? p.x : axis == 1 ? p.y : p.z; }
float& At(Vector3& v, int axis) { return axis == 0 ? v.x : axis == 1 ? v.y : v.z; }

// The icon's mean colour (its opaque pixels), from a copy of the atlas read back once.
Color ItemColour(const std::string& id) {
    static std::unordered_map<std::string, Color> cache;
    static Image atlas{};
    static unsigned atlasId = 0;
    const Texture2D tex = ItemIcons_Texture();
    if (tex.id != atlasId) { // first use, or the atlas was rebuilt
        if (atlas.data) UnloadImage(atlas);
        atlas = tex.id ? LoadImageFromTexture(tex) : Image{};
        atlasId = tex.id;
        cache.clear();
    }
    if (auto it = cache.find(id); it != cache.end()) return it->second;
    Color c{ 150, 150, 150, 255 };
    if (atlas.data && ItemIcons_Has(id)) {
        const Rectangle r = ItemIcons_Source(id);
        double sr = 0, sg = 0, sb = 0; int n = 0;
        for (int y = (int)r.y; y < (int)(r.y + r.height) && y < atlas.height; ++y)
            for (int x = (int)r.x; x < (int)(r.x + r.width) && x < atlas.width; ++x) {
                const Color p = GetImageColor(atlas, x, y);
                if (p.a < 128) continue;
                sr += p.r; sg += p.g; sb += p.b; ++n;
            }
        if (n > 0) c = { (unsigned char)(sr / n), (unsigned char)(sg / n), (unsigned char)(sb / n), 255 };
    }
    cache[id] = c;
    return c;
}

// The point s slots along the run's centreline (clamped to its ends), in plot metres on its floor layer.
Vector3 PointAt(const ConduitRun& run, float s) {
    for (size_t i = 0; i + 1 < run.points.size(); ++i) {
        const CInt3 a = run.points[i], b = run.points[i + 1];
        const int axis = a.x != b.x ? 0 : 2;
        const float len = (float)std::abs(Get(b, axis) - Get(a, axis));
        if (s <= len || i + 2 == run.points.size()) {
            Vector3 v{ (float)a.x, (float)a.y, (float)a.z };
            At(v, axis) += (Get(b, axis) > Get(a, axis) ? 1.0f : -1.0f) * std::clamp(s, 0.0f, len);
            return Vector3Scale(v, kSlot);
        }
        s -= len;
    }
    return run.points.empty() ? Vector3{} : Vector3Scale({ (float)run.points[0].x, (float)run.points[0].y, (float)run.points[0].z }, kSlot);
}

// A box in slot units (lo/hi corners) as one instance.
void SlotBox(Vector3 lo, Vector3 hi, Color c) {
    IB_Cube(Vector3Scale(Vector3Add(lo, hi), 0.5f * kSlot), Vector3Scale(Vector3Subtract(hi, lo), kSlot), c);
}

// One run: per segment the slab and its two rails, then the chevrons (scrolling `phase` slots along the flow). An inner
// corner point gets the slab's square; the outer rail runs on around it, the inner one stops short.
void DrawRun(const ConduitRun& run, float phase, Vector3 camPlot) {
    const size_t n = run.points.size();
    float segStart = 0.0f;
    for (size_t i = 0; i + 1 < n; ++i) {
        const CInt3 a = run.points[i], b = run.points[i + 1];
        const int axis = a.x != b.x ? 0 : 2, across = axis == 0 ? 2 : 0;
        const int dir = Get(b, axis) > Get(a, axis) ? 1 : -1;
        const float len = (float)std::abs(Get(b, axis) - Get(a, axis));
        const float c = (float)Get(a, across), y0 = (float)a.y;
        const float lo = (float)std::min(Get(a, axis), Get(b, axis)), hi = (float)std::max(Get(a, axis), Get(b, axis));
        Vector3 mid{ (float)a.x, (float)a.y, (float)a.z };
        At(mid, axis) = 0.5f * (lo + hi);
        if (Vector3Distance(Vector3Scale(mid, kSlot), camPlot) > kDrawRange + len * 0.5f * kSlot) { segStart += len; continue; }
        // The ends at an inner corner: the other segment's direction away from the corner, across this one (+1 / -1), else 0.
        auto turn = [&](size_t corner, size_t other) {
            if (corner == 0 || corner + 1 == n) return 0;
            const int d = Get(run.points[other], across) - Get(run.points[corner], across);
            return d > 0 ? 1 : d < 0 ? -1 : 0;
        };
        const int turnLo = dir > 0 ? turn(i, i - 1) : turn(i + 1, i + 2), turnHi = dir > 0 ? turn(i + 1, i + 2) : turn(i, i - 1);
        auto box = [&](float l, float h, float c0, float c1, float y1, Color col) {
            Vector3 vlo{}, vhi{};
            At(vlo, axis) = l; At(vhi, axis) = h; At(vlo, across) = c0; At(vhi, across) = c1; vlo.y = y0; vhi.y = y0 + y1;
            SlotBox(vlo, vhi, col);
        };
        box(lo - (turnLo ? kHalfW : 0.0f), hi + (turnHi ? kHalfW : 0.0f), c - kHalfW, c + kHalfW, kTop, kSlab);
        for (int side = -1; side <= 1; side += 2) {
            auto ext = [&](int t) { return t == 0 ? 0.0f : t == side ? -kHalfW : kHalfW; }; // inner: short; outer: around
            const float r0 = side < 0 ? c - kHalfW : c + kHalfW - kRailW;
            box(lo - ext(turnLo), hi + ext(turnHi), r0, r0 + kRailW, kRailH, kRail);
        }
        stats.segments++;
        // Chevrons: at run distances phase + k * kChevronEvery, clear of the segment's ends.
        const float yaw = axis == 0 ? (dir > 0 ? 0.0f : 180.0f) : (dir > 0 ? -90.0f : 90.0f); // local +x along the flow
        const float first = std::ceil((segStart + 1.5f - phase) / kChevronEvery) * kChevronEvery + phase;
        for (float s = first; s <= segStart + len - 1.5f; s += kChevronEvery) {
            Vector3 p{ (float)a.x, y0 + kTop + 0.03f, (float)a.z };
            At(p, axis) += dir * (s - segStart);
            rlPushMatrix();
            rlTranslatef(p.x * kSlot, p.y * kSlot, p.z * kSlot);
            rlRotatef(yaw, 0, 1, 0);
            for (int side = -1; side <= 1; side += 2) { // a ">" of two bars pointing along +x
                rlPushMatrix();
                rlTranslatef(-0.25f * kSlot, 0, side * 0.5f * kSlot);
                rlRotatef(side * 45.0f, 0, 1, 0);
                IB_Cube({ 0, 0, 0 }, { 1.5f * kSlot, 0.06f * kSlot, 0.32f * kSlot }, kChevron);
                rlPopMatrix();
            }
            rlPopMatrix();
            stats.chevrons++;
        }
        segStart += len;
    }
}

} // namespace

const BeltDrawStats& Belts_LastStats() { return stats; }

void Belts_Draw3D(Scene& world) {
    stats = BeltDrawStats{};
    decalPlots.clear();
    const Vector3 sun = Vector3Normalize(world.sunDirection);
    const double clock = Factory::TickClock();
    static Factory::BeltPlotView view;
    rlDrawRenderBatchActive();
    for (const auto& obj : world.objects) {
        Planetoid* body = dynamic_cast<Planetoid*>(obj.get());
        if (!body || !body->isActive || !body->plots) continue;
        const Vector3 camBody = body->toLocal(world.camera.position);
        for (auto& plot : body->plots->plots) {
            if (plot->belts.empty()) continue;
            if (Vector3Distance(plot->toBody({ 0, 0, 0 }), camBody) > 500.0f) continue;
            const Vector3 camPlot = plot->toPlot(camBody);
            PlotState& st = plotState[(body->id * 1315423911ull) ^ plot->id];
            if (st.version != (uint64_t)plot->version) { // sky exposure 1 m over the plot's floor, as the machines (#284)
                st.version = (uint64_t)plot->version;
                st.buried = 0.0f;
                if (body->field) {
                    const Vector3 at = plot->toBody({ 0.0f, 1.0f, 0.0f });
                    const float len = Vector3Length(at);
                    float sky = 1.0f;
                    if (len > 1.0f) body->field->materialSky(at, 1.0f, 1.0f, Vector3Scale(at, 1.0f / len), sky);
                    st.buried = 1.0f - sky;
                }
            }
            const Matrix model = plot->cellModel(*body, world.renderOrigin); // plot coordinates in double (#271)
            const float light = 0.45f + 0.55f * std::clamp(Vector3DotProduct(sun, body->dirToWorld(plot->up)) * 2.0f + 0.5f, 0.0f, 1.0f);
            rlPushMatrix();
            rlMultMatrixf(MatrixToFloat(model));
            IB_SetBuried(st.buried);
            for (const ConduitRun& run : plot->belts) {
                const ConduitMaterial* mat = Conduit_Materials().find(run.material);
                // The factory's clock (pauses with the game): 1/64 cell a tick -> slots, wrapped in double before the float
                const float phase = (float)std::fmod(clock * (mat ? mat->speed : 0) / kBeltSub, (double)kChevronEvery);
                DrawRun(run, phase, camPlot);
                stats.runs++;
            }
            // The items, from the factory's last published tick, carried on by the sim time since (up to two ticks, while the
            // factory thread is one behind), held back behind the item ahead and at the end as the sim moves them.
            if (Factory::BeltItems(body->id, plot->id, view)) {
                const float frac = (float)std::clamp(clock - (double)view.tick, 0.0, 2.0);
                PlotDecals pd; pd.model = model;
                for (const Factory::BeltLineView& line : view.lines) {
                    const ConduitRun* run = nullptr;
                    for (const ConduitRun& r : plot->belts) if (r.id == line.run) { run = &r; break; }
                    if (!run || run->points.size() < 2) continue;
                    float ahead = (float)line.length + kBeltGap;
                    for (size_t k = 0; k < line.pos.size() && k < line.items.size(); ++k) {
                        const float p = (float)line.pos[k];
                        const float d = std::max(p, std::min(p + (line.asleep ? 0.0f : line.speed * frac), std::min((float)line.length, ahead - kBeltGap)));
                        ahead = d;
                        Vector3 c = PointAt(*run, d / kBeltSub);
                        if (Vector3Distance(c, camPlot) > kItemRange) continue;
                        c.y += (kTop * kSlot) + kItemSize * 0.5f;
                        const std::string* name = nullptr;
                        for (const auto& nm : view.names) if (nm.first == line.items[k]) { name = &nm.second; break; }
                        IB_Cube(c, { kItemSize, kItemSize, kItemSize }, name ? ItemColour(*name) : Color{ 150, 150, 150, 255 });
                        stats.items++;
                        if (name && Vector3Distance(c, camPlot) < kIconRange) pd.decals.push_back({ c, *name, light });
                    }
                }
                if (!pd.decals.empty()) decalPlots.push_back(std::move(pd));
            }
            IB_SetBuried(0.0f);
            rlPopMatrix();
        }
    }
}

void Belts_DrawIcons(Scene&) {
    if (decalPlots.empty()) return;
    const Texture2D tex = ItemIcons_Texture();
    if (!tex.id) return;
    rlDrawRenderBatchActive();
    rlDisableBackfaceCulling(); // the plot frame is left-handed
    const float h = kItemSize * 0.45f, up = kItemSize * 0.5f + 0.002f;
    for (const PlotDecals& pd : decalPlots) {
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(pd.model));
        rlSetTexture(tex.id);
        rlBegin(RL_QUADS);
        for (const Decal& d : pd.decals) {
            float uv0[2], uv1[2];
            ItemIcons_UV(d.id, uv0, uv1);
            const unsigned char l = (unsigned char)std::clamp(d.light * 255.0f, 0.0f, 255.0f);
            rlColor4ub(l, l, l, 255);
            rlNormal3f(0, 1, 0);
            const float y = d.c.y + up;
            rlTexCoord2f(uv0[0], uv0[1]); rlVertex3f(d.c.x - h, y, d.c.z - h);
            rlTexCoord2f(uv0[0], uv1[1]); rlVertex3f(d.c.x - h, y, d.c.z + h);
            rlTexCoord2f(uv1[0], uv1[1]); rlVertex3f(d.c.x + h, y, d.c.z + h);
            rlTexCoord2f(uv1[0], uv0[1]); rlVertex3f(d.c.x + h, y, d.c.z - h);
            stats.icons++;
        }
        rlEnd();
        rlSetTexture(0);
        rlPopMatrix();
    }
    rlDrawRenderBatchActive();
    rlEnableBackfaceCulling();
    decalPlots.clear();
}
