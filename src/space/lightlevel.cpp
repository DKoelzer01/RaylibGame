#include "lightlevel.h"
#include "envflag.h"
#include "locallights.h"
#include "player.h"
#include "plots.h"
#include "scene.h"
#include "terraindeltas.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <vector>

// Torch and lamp reach and colour (user, #1739: candle 6 / torch 12 / lamp 20 / big 30 m: GUESS) and the shader gain:
// a torch's own cell then lights about as its point light did at its centre (1.1). Campfires and fire machines stay
// point lights: they flicker and go on and off with the fire.
namespace {
constexpr float kTorchReach = 12.0f, kLampReach = 20.0f; // m (GUESS)
constexpr float kGain = 2.0f;                            // GUESS
constexpr int kRelightBudget = 20000;                    // cells a frame
constexpr float kPlotRange = 160.0f;                     // m: plots whose lights and walls the grid takes

uint64_t BrickKey(int bx, int by, int bz) {
    const uint64_t m = (1u << 21) - 1;
    return ((uint64_t)(bx & m)) | ((uint64_t)(by & m) << 21) | ((uint64_t)(bz & m) << 42);
}
int FloorI(float v) { return (int)std::floor(v); }

struct Source { int x, y, z; LightCell l; };
bool operator<(const Source& a, const Source& b) {
    return std::tie(a.x, a.y, a.z, a.l.level, a.l.r, a.l.g, a.l.b) < std::tie(b.x, b.y, b.z, b.l.level, b.l.r, b.l.g, b.l.b);
}

// One plot's contribution, rebuilt when its version changes: its light sources and its opaque 1 m cells by brick.
struct PlotCache {
    uint32_t version = 0xFFFFFFFFu;
    std::vector<Source> sources;
    std::unordered_map<uint64_t, std::vector<uint16_t>> opaque; // brick key -> cell indices in it
    bool inRange = false;
};

struct State {
    uint64_t bodyId = 0;
    const Planetoid* body = nullptr;
    std::unique_ptr<LightGrid> grid;
    std::shared_ptr<const TerrainDeltas> edits;   // the snapshot the grid's opacity was read from
    std::unordered_map<uint32_t, PlotCache> plots;
    std::map<Source, int> current;                // sources in the grid, with their counts
    LightLevelStats stats;
};
State S;

void FillOpacity(int bx, int by, int bz, uint8_t* out) {
    const Planetoid* body = S.body;
    float dens[LightGrid::kBrickCells];
    const Vector3 o = { bx * 16 + 0.5f, by * 16 + 0.5f, bz * 16 + 0.5f };   // cell centres
    body->field->baseDensityBlock(o, 1.0f, 16, 0.0f, dens);
    const TerrainDeltas* ed = S.edits.get();
    bool edited = false;
    if (ed && !ed->empty()) {
        static std::vector<std::pair<TerrainDeltas::Key, const DeltaChunk*>> in;
        in.clear();
        ed->chunksIn({ o.x - 1.0f, o.y - 1.0f, o.z - 1.0f }, { o.x + 16.0f, o.y + 16.0f, o.z + 16.0f }, in);
        edited = !in.empty();
    }
    for (int z = 0, i = 0; z < 16; ++z)
        for (int y = 0; y < 16; ++y)
            for (int x = 0; x < 16; ++x, ++i) {
                float d = dens[i];
                if (edited) d += ed->sample({ o.x + x, o.y + y, o.z + z });
                out[i] = d >= ISO_LEVEL ? 1 : 0;
            }
    const uint64_t key = BrickKey(bx, by, bz);
    for (auto& kv : S.plots) {
        auto it = kv.second.opaque.find(key);
        if (it != kv.second.opaque.end()) for (uint16_t i : it->second) out[i] = 1;
    }
}

void Reset(const Planetoid* body) {
    S.plots.clear();
    S.current.clear();
    S.grid.reset();
    S.body = body;
    S.bodyId = body ? body->id : 0;
    S.edits = body ? body->edits : nullptr;
    if (body && body->field) S.grid = std::make_unique<LightGrid>(FillOpacity);
}

void RebuildPlot(const Planetoid& body, Plot& plot, PlotCache& c, std::vector<uint64_t>& changedBricks) {
    static int lampMat = -2, glassMat = -2;
    if (lampMat == -2) lampMat = BuildMaterialForItem("lamp");
    if (glassMat == -2) glassMat = BuildMaterialForItem("glass");
    for (auto& kv : c.opaque) changedBricks.push_back(kv.first);
    c.version = plot.version;
    c.sources.clear();
    c.opaque.clear();
    static std::vector<std::pair<Int3, CellShape>> fp;
    const LightCell torch = LightGrid::Emission(kTorchReach, 1.0f, 0.62f, 0.3f), lamp = LightGrid::Emission(kLampReach, 1.0f, 0.92f, 0.76f);
    for (const Part& p : plot.parts) {
        const bool isTorch = p.type == PartType::Torch, isLamp = lampMat >= 0 && IsShape(p.type) && p.material == lampMat;
        const bool solid = IsShape(p.type) && p.material != glassMat;
        if (!isTorch && !solid) continue;
        Plot::Footprint(p, fp);
        if (fp.empty()) continue;
        Vector3 sum{};
        for (auto& cell : fp) {
            const Vector3 centre = { (cell.first.x + 0.5f) * kCellSize, (cell.first.y + 0.5f) * kCellSize, (cell.first.z + 0.5f) * kCellSize };
            sum = Vector3Add(sum, centre);
            if (!solid || !cell.second.full) continue;
            // A 1 m cell holding a full 0.5 m building cell's centre stops light (walls of half-metre blocks still block).
            const Vector3 at = plot.toBody(centre);
            const int x = FloorI(at.x), y = FloorI(at.y), z = FloorI(at.z);
            c.opaque[BrickKey(x >> 4, y >> 4, z >> 4)].push_back((uint16_t)((x & 15) + (y & 15) * 16 + (z & 15) * 256));
        }
        if (!isTorch && !isLamp) continue;
        Vector3 at = plot.toBody(Vector3Scale(sum, 1.0f / fp.size()));
        if (isTorch) at = Vector3Add(at, Vector3Scale(plot.up, 0.4f)); // the flame, at the top of the post
        c.sources.push_back({ FloorI(at.x), FloorI(at.y), FloorI(at.z), isTorch ? torch : lamp });
    }
    for (auto& kv : c.opaque) changedBricks.push_back(kv.first);
    (void)body;
}
}

bool LightLevel_On() {
    static const bool on = !getenv("RAYLIBGAME_LIGHT_GRID") || EnvFlag("RAYLIBGAME_LIGHT_GRID");
    return on;
}

void LightLevel_Update(Scene& world, bool inGame) {
    const double t0 = GetTime();
    const int off[3] = { 0, 0, 0 };
    Planetoid* body = player.nearest;
    if (!inGame || !LightLevel_On() || !body || !body->field) { LocalLights_UploadGrid(nullptr, off, off, MatrixIdentity(), 0.0f); return; }
    if (S.bodyId != body->id || S.body != body) Reset(body);
    LightGrid& g = *S.grid;
    const Vector3 cam = body->toLocal(world.camera.position);

    // Plots: sources and walls, from the plots near the camera; a plot whose version moved relights the bricks it touched.
    std::map<Source, int> want;
    static std::vector<uint64_t> changed;
    changed.clear();
    if (body->plots) {
        for (auto& kv : S.plots) kv.second.inRange = false;
        for (auto& plot : body->plots->plots) {
            if (Vector3Distance(plot->toBody({ 0, 0, 0 }), cam) > kPlotRange) continue;
            PlotCache& c = S.plots[plot->id];
            c.inRange = true;
            if (c.version != plot->version) RebuildPlot(*body, *plot, c, changed);
            for (const Source& s : c.sources) ++want[s];
        }
        for (auto it = S.plots.begin(); it != S.plots.end();) {
            if (it->second.inRange) { ++it; continue; }
            for (auto& kv : it->second.opaque) changed.push_back(kv.first);
            it = S.plots.erase(it);
        }
    } else S.plots.clear();
    std::sort(changed.begin(), changed.end());
    changed.erase(std::unique(changed.begin(), changed.end()), changed.end());
    for (uint64_t k : changed) {
        const int bx = (int)((int64_t)(k << 43) >> 43), by = (int)((int64_t)(k >> 21 << 43) >> 43), bz = (int)((int64_t)(k >> 42 << 43) >> 43);
        g.RefreshOpacity(bx * 16, by * 16, bz * 16, bx * 16 + 15, by * 16 + 15, bz * 16 + 15);
    }
    // Sources: the difference between what the plots hold and what the grid has.
    for (auto it = S.current.begin(); it != S.current.end();) {
        const int have = it->second, need = want.count(it->first) ? want[it->first] : 0;
        for (int k = need; k < have; ++k) g.RemoveSource(it->first.x, it->first.y, it->first.z, it->first.l);
        if (need == 0) it = S.current.erase(it); else { it->second = std::min(have, need); ++it; }
    }
    for (auto& kv : want) {
        int& have = S.current[kv.first];
        for (int k = have; k < kv.second; ++k) g.AddSource(kv.first.x, kv.first.y, kv.first.z, kv.first.l);
        have = kv.second;
    }

    // Terrain edits: each delta chunk that changed (a new snapshot shares the untouched ones) re-reads its cells.
    if (body->edits != S.edits) {
        const std::shared_ptr<const TerrainDeltas> was = S.edits, now = body->edits;
        S.edits = now;
        auto refresh = [&](TerrainDeltas::Key k) {
            int cx, cy, cz;
            TerrainDeltas::Unpack(k, cx, cy, cz);
            const float span = kDeltaChunk * kVoxelSize;
            g.RefreshOpacity(FloorI(cx * span) - 1, FloorI(cy * span) - 1, FloorI(cz * span) - 1,
                             FloorI((cx + 1) * span) + 1, FloorI((cy + 1) * span) + 1, FloorI((cz + 1) * span) + 1);
        };
        if (now) for (auto& kv : now->chunks) {
            if (!was) { refresh(kv.first); continue; }
            auto it = was->chunks.find(kv.first);
            if (it == was->chunks.end() || it->second != kv.second) refresh(kv.first);
        }
        if (was) for (auto& kv : was->chunks) if (!now || !now->chunks.count(kv.first)) refresh(kv.first);
    }

    g.Update(kRelightBudget);

    // The GPU window: anchored at a multiple of 128 m near the camera (the texture wraps the same way), so the matrix
    // carries small numbers; localModel's translation is computed in double.
    const int camCell[3] = { FloorI(cam.x), FloorI(cam.y), FloorI(cam.z) };
    const int anchor[3] = { (camCell[0] >> 7) << 7, (camCell[1] >> 7) << 7, (camCell[2] >> 7) << 7 };
    const Matrix model = body->localModel(DVec3(anchor[0], anchor[1], anchor[2]), world.renderOrigin);
    LocalLights_UploadGrid(&g, camCell, anchor, MatrixInvert(model), kGain);

    S.stats.sources = (int)want.size();
    S.stats.bricks = (int)g.BrickCount();
    S.stats.busy = g.Busy();
    S.stats.uploaded = LocalLights_GridBricksUploaded();
    S.stats.ms = (float)((GetTime() - t0) * 1000.0);
}

LightCell LightLevel_At(const Planetoid& body, Vector3 local) {
    if (!S.grid || S.body != &body) return {};
    return S.grid->At(FloorI(local.x), FloorI(local.y), FloorI(local.z));
}

LightLevelStats LightLevel_Stats() { return S.stats; }
