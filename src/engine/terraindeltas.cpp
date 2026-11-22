#include "terraindeltas.h"
#include "voxel.h"
#include "tools.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>

namespace {
    constexpr int kBias = 1 << 20; // chunk coordinates are packed as 21-bit unsigned values
    inline int FloorDiv(int a, int b) { return a >= 0 ? a / b : -((-a + b - 1) / b); }
    inline int Idx(int lx, int ly, int lz) { return (lz * kDeltaChunk + ly) * kDeltaChunk + lx; }
}

TerrainDeltas::Key TerrainDeltas::KeyOf(int cx, int cy, int cz) {
    return ((Key)(uint32_t)(cx + kBias) << 42) | ((Key)(uint32_t)(cy + kBias) << 21) | (Key)(uint32_t)(cz + kBias);
}

void TerrainDeltas::Unpack(Key k, int& cx, int& cy, int& cz) {
    cx = (int)((k >> 42) & 0x1FFFFF) - kBias;
    cy = (int)((k >> 21) & 0x1FFFFF) - kBias;
    cz = (int)(k & 0x1FFFFF) - kBias;
}

TerrainDeltas::Key TerrainDeltas::RegionOf(Key chunk) {
    int x, y, z;
    Unpack(chunk, x, y, z);
    return KeyOf(FloorDiv(x, kRegionChunks), FloorDiv(y, kRegionChunks), FloorDiv(z, kRegionChunks));
}

float TerrainDeltas::at(int x, int y, int z) const {
    if (chunks.empty()) return 0.0f;
    const int cx = FloorDiv(x, kDeltaChunk), cy = FloorDiv(y, kDeltaChunk), cz = FloorDiv(z, kDeltaChunk);
    auto it = chunks.find(KeyOf(cx, cy, cz));
    if (it == chunks.end()) return 0.0f;
    return it->second->get(Idx(x - cx * kDeltaChunk, y - cy * kDeltaChunk, z - cz * kDeltaChunk)) * kDeltaStep;
}

float TerrainDeltas::sample(Vector3 pos) const {
    if (chunks.empty()) return 0.0f;
    const Vector3 p = Vector3Scale(pos, 1.0f / kVoxelSize); // lattice coordinates
    const int x0 = (int)floorf(p.x), y0 = (int)floorf(p.y), z0 = (int)floorf(p.z);
    const float fx = p.x - x0, fy = p.y - y0, fz = p.z - z0;
    const float c00 = at(x0, y0, z0) * (1 - fx) + at(x0 + 1, y0, z0) * fx;
    const float c10 = at(x0, y0 + 1, z0) * (1 - fx) + at(x0 + 1, y0 + 1, z0) * fx;
    const float c01 = at(x0, y0, z0 + 1) * (1 - fx) + at(x0 + 1, y0, z0 + 1) * fx;
    const float c11 = at(x0, y0 + 1, z0 + 1) * (1 - fx) + at(x0 + 1, y0 + 1, z0 + 1) * fx;
    return (c00 * (1 - fy) + c10 * fy) * (1 - fz) + (c01 * (1 - fy) + c11 * fy) * fz;
}

int TerrainDeltas::materialAt(Vector3 pos) const {
    if (chunks.empty()) return -1;
    const Vector3 p = Vector3Scale(pos, 1.0f / kVoxelSize);
    auto mat = [&](int x, int y, int z) -> int {
        const int cx = FloorDiv(x, kDeltaChunk), cy = FloorDiv(y, kDeltaChunk), cz = FloorDiv(z, kDeltaChunk);
        auto it = chunks.find(KeyOf(cx, cy, cz));
        if (it == chunks.end()) return -1;
        return MatOf(it->second->mat(Idx(x - cx * kDeltaChunk, y - cy * kDeltaChunk, z - cz * kDeltaChunk)));
    };
    const int m = mat((int)lroundf(p.x), (int)lroundf(p.y), (int)lroundf(p.z));
    if (m >= 0) return m;
    const int x0 = (int)floorf(p.x), y0 = (int)floorf(p.y), z0 = (int)floorf(p.z);
    for (int k = 0; k < 8; ++k) {
        const int c = mat(x0 + (k & 1), y0 + ((k >> 1) & 1), z0 + ((k >> 2) & 1));
        if (c >= 0) return c;
    }
    return -1;
}

int TerrainDeltas::placedAt(Vector3 pos) const {
    if (chunks.empty()) return -1;
    const int x = (int)lroundf(pos.x / kVoxelSize), y = (int)lroundf(pos.y / kVoxelSize), z = (int)lroundf(pos.z / kVoxelSize);
    const int cx = FloorDiv(x, kDeltaChunk), cy = FloorDiv(y, kDeltaChunk), cz = FloorDiv(z, kDeltaChunk);
    auto it = chunks.find(KeyOf(cx, cy, cz));
    if (it == chunks.end()) return -1;
    const uint8_t v = it->second->mat(Idx(x - cx * kDeltaChunk, y - cy * kDeltaChunk, z - cz * kDeltaChunk));
    return MatPlaced(v) ? MatOf(v) : -1;
}

void TerrainDeltas::chunksIn(Vector3 loPos, Vector3 hiPos, std::vector<std::pair<Key, const DeltaChunk*>>& out) const {
    const Vector3 lo = Vector3Scale(loPos, 1.0f / kVoxelSize), hi = Vector3Scale(hiPos, 1.0f / kVoxelSize); // lattice coordinates
    out.clear();
    if (chunks.empty()) return;
    const int x0 = FloorDiv((int)floorf(lo.x), kDeltaChunk), x1 = FloorDiv((int)floorf(hi.x), kDeltaChunk);
    const int y0 = FloorDiv((int)floorf(lo.y), kDeltaChunk), y1 = FloorDiv((int)floorf(hi.y), kDeltaChunk);
    const int z0 = FloorDiv((int)floorf(lo.z), kDeltaChunk), z1 = FloorDiv((int)floorf(hi.z), kDeltaChunk);
    const long long boxChunks = (long long)(x1 - x0 + 1) * (y1 - y0 + 1) * (z1 - z0 + 1);
    if (boxChunks < (long long)chunks.size()) {
        for (int z = z0; z <= z1; ++z)
            for (int y = y0; y <= y1; ++y)
                for (int x = x0; x <= x1; ++x) {
                    auto it = chunks.find(KeyOf(x, y, z));
                    if (it != chunks.end()) out.push_back({ it->first, it->second.get() });
                }
    } else {
        for (const auto& kv : chunks) {
            int x, y, z;
            Unpack(kv.first, x, y, z);
            if (x >= x0 && x <= x1 && y >= y0 && y <= y1 && z >= z0 && z <= z1) out.push_back({ kv.first, kv.second.get() });
        }
    }
}

// #693 A/B: RAYLIBGAME_CARVE_ROUND_OLD=1 = carved points round to the nearest step even when that lands on iso.
static const bool gCarveRoundOld = [] { const char* s = getenv("RAYLIBGAME_CARVE_ROUND_OLD"); return s && s[0] == '1'; }();

void TerrainDeltasBuilder::apply(const TerrainEdit& e, const DensityField& field, std::unordered_set<TerrainDeltas::Key>* touched) {
    next->strokes++;
    // A carve also lowers solid points just outside the sphere, so the meshed surface sits on it.
    const float reach = e.box ? sqrtf(e.half.x * e.half.x + e.half.y * e.half.y + e.half.z * e.half.z) : e.radius;
    const bool csg = e.carve || e.fill; // exact shapes (carve, #228 fill) also move points just outside them
    const float r = csg ? reach + 1.0f : reach; // 1 m: enough on the 0.25 m lattice (2 m touched 8 delta chunks a stroke)
    // Signed distance to the carved shape (the sphere, or the oriented box, #210): < 0 inside.
    auto shapeDist = [&](float dx, float dy, float dz, float dist) {
        if (!e.box) return dist - e.radius;
        const float lx = fabsf(dx * e.ax.x + dy * e.ax.y + dz * e.ax.z) - e.half.x;
        const float ly = fabsf(dx * e.ay.x + dy * e.ay.y + dz * e.ay.z) - e.half.y;
        const float lz = fabsf(dx * e.az.x + dy * e.az.y + dz * e.az.z) - e.half.z;
        const float ox = std::max(lx, 0.0f), oy = std::max(ly, 0.0f), oz = std::max(lz, 0.0f);
        return sqrtf(ox * ox + oy * oy + oz * oz) + std::min(std::max(lx, std::max(ly, lz)), 0.0f);
    };
    const int x0 = (int)floorf((e.center.x - r) / kVoxelSize), x1 = (int)ceilf((e.center.x + r) / kVoxelSize); // lattice indices
    const int y0 = (int)floorf((e.center.y - r) / kVoxelSize), y1 = (int)ceilf((e.center.y + r) / kVoxelSize); // lattice indices
    const int z0 = (int)floorf((e.center.z - r) / kVoxelSize), z1 = (int)ceilf((e.center.z + r) / kVoxelSize); // lattice indices
    DeltaChunk* chunk = nullptr;
    TerrainDeltas::Key chunkKey = ~0ull;
    for (int z = z0; z <= z1; ++z)
        for (int y = y0; y <= y1; ++y)
            for (int x = x0; x <= x1; ++x) {
                const Vector3 p = { x * kVoxelSize, y * kVoxelSize, z * kVoxelSize }; // the lattice point, body-local
                const float dx = p.x - e.center.x, dy = p.y - e.center.y, dz = p.z - e.center.z;
                const float d2 = dx * dx + dy * dy + dz * dz;
                if (d2 >= r * r) continue;
                if (!field.editable(p)) continue; // e.g. a planet's mantle
                const float dist = sqrtf(d2);
                const float t = 1.0f - dist / r;
                float add = e.strength * t * t * (3.0f - 2.0f * t);
                float carveTo = 0.0f; // carve: the density this point may have at most
                const float sd = csg ? shapeDist(dx, dy, dz, dist) : 0.0f;
                // #693: an eighth of a storage step (8 mm) under the shape's own surface. Ground exactly at iso (a point within float
                // precision of the ground, a few mm at a planet's radius, reads exactly 0.5 = solid) just outside a face
                // flush with the ground stayed solid over the cut: specks of ground in a lattice pattern.
                if (e.carve) carveTo = ISO_LEVEL + DENSITY_GRADIENT * sd - (gCarveRoundOld ? 0.0f : 0.125f * kDeltaStep); // 8 mm: float precision at a 60 km radius is ~4 mm
                const int cx = (int)floorf((float)x / kDeltaChunk), cy = (int)floorf((float)y / kDeltaChunk), cz = (int)floorf((float)z / kDeltaChunk);
                const TerrainDeltas::Key k = TerrainDeltas::KeyOf(cx, cy, cz);
                if (k != chunkKey) {
                    chunkKey = k;
                    auto own = owned.find(k);
                    if (own != owned.end()) chunk = own->second;
                    else {
                        // First touch in this build: copy (or create) the chunk so the published
                        // snapshot workers may hold is never written.
                        auto it = next->chunks.find(k);
                        auto copy = it != next->chunks.end() ? std::make_shared<DeltaChunk>(*it->second) : std::make_shared<DeltaChunk>();
                        chunk = copy.get();
                        next->chunks[k] = copy;
                        owned[k] = chunk;
                    }
                    if (touched) touched->insert(k);
                }
                const int i = Idx(x - cx * kDeltaChunk, y - cy * kDeltaChunk, z - cz * kDeltaChunk);
                const int16_t v = chunk->get(i);
                // Grass is a state of the generated surface, not a layer (#1085): ground an edit exposes or makes is bare dirt.
                // The generator calls everything in the ~1 m band of its surface grass, so a point the edit leaves on the new
                // surface takes Regolith in the material channel where the field would say Grass. Placed ground keeps its own.
                auto bareGrass = [&](float base) {
                    // (#1097) natural ground in another state (kMatState): it digs as natural ground, extras and perks included
                    if (!chunk->mat(i) && field.material(p, base) == VoxelMaterial::Grass) chunk->setMat(i, (uint8_t)(((uint8_t)VoxelMaterial::Regolith + 1) | kMatState));
                };
                if (e.till) {
                    // The surface band of soil: a lattice point at most one step over the ground (the mesh's vertices take
                    // the nearest point's material) down to half a metre under it. Grass goes: one channel, one state.
                    if (dist >= e.radius) continue;
                    const float base = field.baseDensity(p), d = base + v * kDeltaStep;
                    if (d < ISO_LEVEL - DENSITY_GRADIENT * kVoxelSize || d > ISO_LEVEL + DENSITY_GRADIENT * 0.5f) continue;
                    const uint8_t cur = chunk->mat(i);
                    const int placed = MatOf(cur); // (#1097: the material without its state bit)
                    const VoxelMaterial m = placed >= 0 ? (VoxelMaterial)placed : field.material(p, base);
                    // (#1097) tilled natural ground is a state of it (kMatState); tilled placed dirt stays placed (it digs back at
                    // most what it cost: tilling can't turn placed ground into natural)
                    if (IsTillable(m)) chunk->setMat(i, (uint8_t)(((uint8_t)VoxelMaterial::Farmland + 1) | (MatPlaced(cur) ? 0 : kMatState)));
                    continue;
                }
                if (e.carve) {
                    const float base = field.baseDensity(p);
                    // Placed ground digs like what it is (placed gravel with a shovel), not the rock below.
                    const int placed = MatOf(chunk->mat(i));
                    const VoxelMaterial was = placed >= 0 ? (VoxelMaterial)placed : base >= ISO_LEVEL ? field.material(p, base) : VoxelMaterial::Regolith;
                    if (e.softOnly && (base + v * kDeltaStep) >= ISO_LEVEL && !IsSoft(was)) continue;
                    if (base + v * kDeltaStep > carveTo) {
                        int q = (int)lroundf((carveTo - base) / kDeltaStep);
                        // #693: inside the shape the point must end below iso. A face passing within half a step (~3 cm) of
                        // the point (a box top flush with the ground, two boxes sharing a face) rounded it back up to iso:
                        // solid, so films and specks of ground stayed in a lattice pattern over the cut.
                        if (carveTo < ISO_LEVEL && base + q * kDeltaStep >= ISO_LEVEL && !gCarveRoundOld) --q;
                        chunk->set(i, (int16_t)std::clamp(q, -32767, 32767));
                        if (sd < 0.0f) chunk->setMat(i, 0); // dug out: what is exposed below shows its own material
                    }
                    if (fabsf(sd) < kVoxelSize) bareGrass(base); // the points either side of the cut's face: bare dirt, not turf
                } else if (e.fill) {
                    const float base = field.baseDensity(p);
                    const float fillTo = ISO_LEVEL - DENSITY_GRADIENT * sd;
                    if (base + v * kDeltaStep < fillTo) {
                        // Placed ground is what was air before (#1097: it is what the fill paid for, 2 items a m^3); ground
                        // that was there already keeps its own material and digs as itself.
                        const bool wasAir = base + v * kDeltaStep < ISO_LEVEL;
                        const int q = std::clamp((int)ceilf((fillTo - base) / kDeltaStep), -32767, 32767); // up: never just under iso inside
                        chunk->set(i, (int16_t)q);
                        // (#1157) Only the ground the fill made takes its material: a shell point left in the air beside it handed
                        // that material to the natural surface next to it (a surface vertex takes a marked corner's). Grass the
                        // fill covers (inside the shape) goes bare; the grass beside it stays.
                        const bool made = wasAir && base + q * kDeltaStep >= ISO_LEVEL;
                        if (e.paint >= 0 && made) chunk->setMat(i, (uint8_t)(e.paint + 1)); // placed (#1097)
                        else if (sd < 0.0f || made) bareGrass(base); // (#1085) grass under it, or ground made with no material: bare dirt
                    }
                } else {
                    const float base = field.baseDensity(p);
                    const bool wasAir = base + v * kDeltaStep < ISO_LEVEL; // (#1097, as a fill)
                    const int nv = std::clamp((int)lroundf(v + add / kDeltaStep), -32767, 32767);
                    chunk->set(i, (int16_t)nv);
                    // (#1157) As the fill: what the blob made solid is placed ground; a point it raised that is still air isn't
                    // (marked, it re-textured the natural ground beside the blob). Grass goes bare only under its core.
                    if (e.paint >= 0 && wasAir && add > 0.0f && base + nv * kDeltaStep >= ISO_LEVEL) chunk->setMat(i, (uint8_t)(e.paint + 1)); // (#1097)
                    else if (!wasAir && t > 0.5f) bareGrass(base); // (#1085) grass buried under the blob: bare dirt, natural (#1097)
                }
            }
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#869): carve and fill rounding (#693), keys, the sparse chunk store ----
#include "unit_helpers.h"
#include <random>

namespace {
// Flat ground: solid below y = 0, the lattice row at y = 0 exactly at iso (the worst case for #693's rounding).
class FlatField : public DensityField {
public:
    float baseDensity(Vector3 p) const override { return ISO_LEVEL - DENSITY_GRADIENT * p.y; }
    VoxelMaterial material(Vector3, float) const override { return VoxelMaterial::Rock; }
    Vector3 fallbackUp(Vector3) const override { return { 0, 1, 0 }; }
    float maxRadius() const override { return 1e6f; }
};
float LatticeDensity(const TerrainDeltas& d, const DensityField& f, int x, int y, int z) {
    return f.baseDensity({ x * kVoxelSize, y * kVoxelSize, z * kVoxelSize }) + d.at(x, y, z);
}
}

// #693: a box carved flush with the ground leaves no lattice point inside it (or on its faces) solid: no specks.
TEST(TerrainDeltas, FlushBoxCarveLeavesNoSolidPointInside) {
    const FlatField field;
    const float offsets[] = { 0.0f, 0.001f, -0.001f, 0.0625f, 0.124f }; // the box top at, just over and just under the ground
    for (float off : offsets) {
        TerrainDeltas empty;
        TerrainDeltasBuilder b(empty);
        TerrainEdit e{};
        e.carve = true; e.box = true;
        e.half = { 1.0f, 0.5f, 1.0f };
        e.center = { 0.1f, -0.5f + off, -0.2f };
        e.radius = Vector3Length(e.half);
        b.apply(e, field);
        const auto d = b.finish();
        int solidInside = 0, inside = 0;
        for (int z = -8; z <= 8; ++z) for (int y = -8; y <= 4; ++y) for (int x = -8; x <= 8; ++x) {
            const Vector3 p = { x * kVoxelSize, y * kVoxelSize, z * kVoxelSize };
            const bool in = fabsf(p.x - e.center.x) <= e.half.x && fabsf(p.y - e.center.y) <= e.half.y && fabsf(p.z - e.center.z) <= e.half.z;
            if (!in) continue;
            ++inside;
            if (LatticeDensity(*d, field, x, y, z) >= ISO_LEVEL) ++solidInside;
        }
        EXPECT_GT(inside, 100);
        EXPECT_EQ(solidInside, 0) << "box top " << off << " m over the ground";
    }
}

TEST(TerrainDeltas, SphereCarveAndFillAreExact) {
    const FlatField field;
    std::mt19937 rng(693);
    std::uniform_real_distribution<float> c(-0.3f, 0.3f), r(0.4f, 1.6f);
    for (int k = 0; k < 40; ++k) {
        TerrainDeltas empty;
        TerrainDeltasBuilder b(empty);
        TerrainEdit carve{}; carve.carve = true; carve.center = { c(rng), c(rng), c(rng) }; carve.radius = r(rng);
        b.apply(carve, field);
        TerrainEdit fill{}; fill.fill = true; fill.center = { carve.center.x + 4.0f, 1.0f + c(rng), carve.center.z }; fill.radius = r(rng); fill.paint = (int)VoxelMaterial::Gravel;
        b.apply(fill, field);
        const auto d = b.finish();
        EXPECT_EQ(d->strokes, 2u);
        for (int z = -10; z <= 10; ++z) for (int y = -10; y <= 12; ++y) for (int x = -10; x <= 26; ++x) {
            const Vector3 p = { x * kVoxelSize, y * kVoxelSize, z * kVoxelSize };
            const float dens = LatticeDensity(*d, field, x, y, z);
            if (Vector3Distance(p, carve.center) < carve.radius) ASSERT_LT(dens, ISO_LEVEL) << "carved point solid at " << x << " " << y << " " << z;
            if (Vector3Distance(p, fill.center) < fill.radius) ASSERT_GE(dens, ISO_LEVEL) << "filled point empty at " << x << " " << y << " " << z;
        }
        EXPECT_EQ(d->materialAt(fill.center), (int)VoxelMaterial::Gravel) << "filled ground keeps its material";
    }
}

// #1085: grass is a state of the generated surface. A field that, like the generator, calls the ~1 m band under its surface
// grass: a dig there, or a fill without a material, leaves bare dirt on every face point it makes; rock stays rock.
namespace {
struct TurfField : FlatField {
    VoxelMaterial material(Vector3, float base) const override {
        return base < ISO_LEVEL + DENSITY_GRADIENT * 1.0f ? VoxelMaterial::Grass : VoxelMaterial::Rock;
    }
};
// The material a lattice point shows: the placed one, else the field's.
VoxelMaterial Shown(const TerrainDeltas& d, const DensityField& f, Vector3 p) {
    const int m = d.materialAt(p);
    return m >= 0 ? (VoxelMaterial)m : f.material(p, f.baseDensity(p));
}
// Solid lattice points with an empty neighbour: the points the surface's vertices take their material from.
template <class Fn> void ForSurfacePoints(const TerrainDeltas& d, const DensityField& f, Fn fn) {
    for (int z = -12; z <= 12; ++z) for (int y = -12; y <= 12; ++y) for (int x = -12; x <= 12; ++x) {
        if (LatticeDensity(d, f, x, y, z) < ISO_LEVEL) continue;
        bool face = false;
        const int n[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
        for (auto& o : n) face |= LatticeDensity(d, f, x + o[0], y + o[1], z + o[2]) < ISO_LEVEL;
        if (face) fn(Vector3{ x * kVoxelSize, y * kVoxelSize, z * kVoxelSize });
    }
}
}

TEST(TerrainDeltas, DiggingGrassLeavesBareDirt) {
    const TurfField field;
    TerrainDeltas empty;
    TerrainDeltasBuilder b(empty);
    TerrainEdit dig{}; dig.carve = true; dig.center = { 0.1f, 0.0f, -0.1f }; dig.radius = 0.7f; // the hole's floor 0.7 m down: still in the turf band
    b.apply(dig, field);
    const auto d = b.finish();
    int cut = 0, grassInCut = 0, untouchedGrass = 0;
    ForSurfacePoints(*d, field, [&](Vector3 p) {
        const float r = Vector3Distance(p, dig.center);
        if (r < dig.radius + kVoxelSize) { ++cut; if (Shown(*d, field, p) == VoxelMaterial::Grass) ++grassInCut; }
        else if (r > dig.radius + 1.5f && Shown(*d, field, p) == VoxelMaterial::Grass) ++untouchedGrass;
    });
    EXPECT_GT(cut, 30);
    EXPECT_EQ(grassInCut, 0) << "the hole's faces show turf";
    EXPECT_GT(untouchedGrass, 100) << "the untouched surface keeps its grass";
}

TEST(TerrainDeltas, FillsOnGrassAreBareDirtAndRockDigsStayRock) {
    { // a fill without a material (no filler in hand): bare dirt, not turf
        const TurfField field;
        TerrainDeltas empty;
        TerrainDeltasBuilder b(empty);
        TerrainEdit fill{}; fill.fill = true; fill.center = { 0.0f, 0.5f, 0.0f }; fill.radius = 0.6f;
        b.apply(fill, field);
        const auto d = b.finish();
        int mound = 0, grass = 0;
        ForSurfacePoints(*d, field, [&](Vector3 p) {
            if (p.y > 0.2f && Vector3Distance(p, fill.center) < fill.radius + kVoxelSize) { ++mound; if (Shown(*d, field, p) == VoxelMaterial::Grass) ++grass; }
        });
        EXPECT_GT(mound, 20);
        EXPECT_EQ(grass, 0);
    }
    { // a dig in rock writes no material: the rock shows its own
        const FlatField field;
        TerrainDeltas empty;
        TerrainDeltasBuilder b(empty);
        TerrainEdit dig{}; dig.carve = true; dig.center = { 0.0f, 0.0f, 0.0f }; dig.radius = 0.7f;
        b.apply(dig, field);
        const auto d = b.finish();
        int placed = 0;
        ForSurfacePoints(*d, field, [&](Vector3 p) { placed += d->materialAt(p) >= 0; });
        EXPECT_EQ(placed, 0);
    }
}

// #1085: a till paints the surface band of soil Farmland inside its radius, leaves the ground's shape alone, and leaves
// rock (and anything outside the radius) as it was.
TEST(TerrainDeltas, TillingPaintsFarmlandWithoutChangingTheGround) {
    {
        const TurfField field;
        TerrainDeltas empty;
        TerrainDeltasBuilder b(empty);
        TerrainEdit till{}; till.till = true; till.center = { 0.05f, 0.0f, 0.1f }; till.radius = 1.0f;
        b.apply(till, field);
        const auto d = b.finish();
        int inside = 0, farmland = 0, outsideFarmland = 0;
        ForSurfacePoints(*d, field, [&](Vector3 p) {
            const float r = Vector3Distance(p, till.center);
            const VoxelMaterial m = Shown(*d, field, p);
            if (r < till.radius - kVoxelSize) { ++inside; farmland += m == VoxelMaterial::Farmland; }
            // materialAt falls back to a cell corner's placed material: the rim reads one lattice step wider (as drawn)
            else if (r > till.radius + 2.0f * kVoxelSize) outsideFarmland += m == VoxelMaterial::Farmland;
        });
        EXPECT_GT(inside, 10);
        EXPECT_EQ(farmland, inside) << "every surface point in the brush is farmland";
        EXPECT_EQ(outsideFarmland, 0);
        for (int z = -8; z <= 8; ++z) for (int y = -8; y <= 8; ++y) for (int x = -8; x <= 8; ++x)
            ASSERT_EQ(d->at(x, y, z), 0.0f) << "a till changed the ground's density";
    }
    { // rock does not till
        const FlatField field;
        TerrainDeltas empty;
        TerrainDeltasBuilder b(empty);
        TerrainEdit till{}; till.till = true; till.center = { 0.0f, 0.0f, 0.0f }; till.radius = 0.6f;
        b.apply(till, field);
        const auto d = b.finish();
        int placed = 0;
        ForSurfacePoints(*d, field, [&](Vector3 p) { placed += d->materialAt(p) >= 0; });
        EXPECT_EQ(placed, 0);
    }
}

// #1157 (the playtest: a limestone blob re-textured the grass beside it): only the ground a placement makes takes its
// material. The natural surface beside it shows its own (grass) beyond the one lattice step the mesh's corner lookup reaches.
TEST(TerrainDeltas, PlacedGroundDoesNotRetextureTheGroundBesideIt) {
    const TurfField field;
    for (int kind = 0; kind < 2; ++kind) { // the soft blob (RMB), then an exact sphere fill (#228)
        TerrainDeltas empty;
        TerrainDeltasBuilder b(empty);
        TerrainEdit e{};
        e.paint = (int)VoxelMaterial::Limestone;
        if (kind == 0) { e.center = { 0.05f, 0.3f, -0.05f }; e.radius = 1.3f; e.strength = 0.6f; }
        else { e.fill = true; e.center = { 0.05f, 0.4f, -0.05f }; e.radius = 0.7f; }
        b.apply(e, field);
        const auto d = b.finish();
        std::vector<Vector3> placed;
        for (int z = -12; z <= 12; ++z) for (int y = -12; y <= 12; ++y) for (int x = -12; x <= 12; ++x) {
            const Vector3 p = { x * kVoxelSize, y * kVoxelSize, z * kVoxelSize };
            if (d->placedAt(p) >= 0) {
                placed.push_back(p);
                ASSERT_GE(LatticeDensity(*d, field, x, y, z), ISO_LEVEL) << "a placed point is solid ground (kind " << kind << ")";
            }
        }
        ASSERT_GT(placed.size(), 10u) << "the placement made ground (kind " << kind << ")";
        int beside = 0, bleed = 0, grassFar = 0, farPoints = 0;
        ForSurfacePoints(*d, field, [&](Vector3 p) {
            if (d->placedAt(p) >= 0) return;
            float near = 1e9f;
            for (const Vector3& q : placed) near = std::min(near, Vector3Distance(p, q));
            if (near <= 1.01f * kVoxelSize * 1.7321f) return; // a corner of a placed point's cell: drawn as the rim
            ++beside; bleed += Shown(*d, field, p) == VoxelMaterial::Limestone;
            if (Vector3Distance(p, e.center) > 2.0f && p.y > -1.0f) { ++farPoints; grassFar += Shown(*d, field, p) == VoxelMaterial::Grass; }
        });
        EXPECT_GT(beside, 50);
        EXPECT_EQ(bleed, 0) << "the natural surface beside the placement shows the placed material (kind " << kind << ")";
        EXPECT_EQ(grassFar, farPoints) << "the grass beside it stays grass (kind " << kind << ")";
    }
}

TEST(TerrainDeltas, EditsMakeNewSnapshotsAndShareUntouchedChunks) {
    const FlatField field;
    TerrainDeltas empty;
    TerrainDeltasBuilder b1(empty);
    TerrainEdit a{}; a.carve = true; a.center = { 0, 0, 0 }; a.radius = 1.0f;
    TerrainEdit far{}; far.carve = true; far.center = { 100, 0, 0 }; far.radius = 1.0f;
    b1.apply(a, field); b1.apply(far, field);
    const auto s1 = b1.finish();
    TerrainDeltasBuilder b2(*s1);
    std::unordered_set<TerrainDeltas::Key> touched;
    TerrainEdit again = a; again.center = { 0.5f, 0, 0 };
    b2.apply(again, field, &touched);
    const auto s2 = b2.finish();
    EXPECT_FALSE(touched.empty());
    int shared = 0, copied = 0;
    for (const auto& kv : s2->chunks) {
        auto it = s1->chunks.find(kv.first);
        if (it == s1->chunks.end()) continue;
        (it->second == kv.second ? shared : copied)++;
        if (it->second != kv.second) EXPECT_TRUE(touched.count(kv.first)) << "a changed chunk is reported";
    }
    EXPECT_GT(shared, 0) << "the far stroke's chunks are shared, not copied";
    EXPECT_GT(copied, 0);
    EXPECT_LT(LatticeDensity(*s1, field, 0, 0, 0), ISO_LEVEL);
}

TEST(TerrainDeltas, KeysPackAndRegionsFloor) {
    std::mt19937 rng(5);
    std::uniform_int_distribution<int> d(-1000000, 1000000);
    for (int i = 0; i < 10000; ++i) {
        const int x = d(rng), y = d(rng), z = d(rng);
        int a, b, c;
        TerrainDeltas::Unpack(TerrainDeltas::KeyOf(x, y, z), a, b, c);
        ASSERT_EQ(a, x); ASSERT_EQ(b, y); ASSERT_EQ(c, z);
    }
    int rx, ry, rz;
    TerrainDeltas::Unpack(TerrainDeltas::RegionOf(TerrainDeltas::KeyOf(-1, 0, 15)), rx, ry, rz);
    EXPECT_EQ(rx, -1); EXPECT_EQ(ry, 0); EXPECT_EQ(rz, 0);
    TerrainDeltas::Unpack(TerrainDeltas::RegionOf(TerrainDeltas::KeyOf(-16, -17, 16)), rx, ry, rz);
    EXPECT_EQ(rx, -1); EXPECT_EQ(ry, -2); EXPECT_EQ(rz, 1);
}

TEST(TerrainDeltas, SparseChunkDenseRoundTrip) {
    std::mt19937 rng(11);
    std::vector<int16_t> dense(DeltaChunk::kN, 0);
    for (int k = 0; k < 3000; ++k) dense[rng() % DeltaChunk::kN] = (int16_t)((int)(rng() % 2001) - 1000);
    DeltaChunk c;
    c.fromDense(dense.data());
    for (int i = 0; i < DeltaChunk::kN; ++i) ASSERT_EQ(c.get(i), dense[i]) << i;
    std::vector<int16_t> out(DeltaChunk::kN, 7);
    c.toDense(out.data());
    EXPECT_EQ(out, dense);
    DeltaChunk e;
    for (int i = 0; i < DeltaChunk::kN; i += 97) e.set(i, 0);
    EXPECT_TRUE(e.empty()) << "zeros allocate no brick";
}
#endif
