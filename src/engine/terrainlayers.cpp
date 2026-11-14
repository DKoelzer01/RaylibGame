#include "terrainlayers.h"
#include "glcaps.h"
#include "utils/logger.h"

#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

// The GL entry points we need, from the glad loader inside raylib (initialised by InitWindow).
extern "C" {
    typedef void (*PfnGenTextures)(int, unsigned int*);
    typedef void (*PfnDeleteTextures)(int, const unsigned int*);
    typedef void (*PfnBindTexture)(unsigned int, unsigned int);
    typedef void (*PfnTexParameteri)(unsigned int, unsigned int, int);
    typedef void (*PfnTexParameterf)(unsigned int, unsigned int, float);
    typedef void (*PfnTexImage3D)(unsigned int, int, int, int, int, int, int, unsigned int, unsigned int, const void*);
    typedef void (*PfnTexSubImage3D)(unsigned int, int, int, int, int, int, int, int, unsigned int, unsigned int, const void*);
    typedef void (*PfnGenerateMipmap)(unsigned int);
    typedef void (*PfnActiveTexture)(unsigned int);
    extern PfnGenTextures glad_glGenTextures;
    extern PfnDeleteTextures glad_glDeleteTextures;
    extern PfnBindTexture glad_glBindTexture;
    extern PfnTexParameteri glad_glTexParameteri;
    extern PfnTexParameterf glad_glTexParameterf;
    extern PfnTexImage3D glad_glTexImage3D;
    extern PfnTexSubImage3D glad_glTexSubImage3D;
    extern PfnGenerateMipmap glad_glGenerateMipmap;
    extern PfnActiveTexture glad_glActiveTexture;
}

namespace {
    constexpr unsigned int GL_TEXTURE_2D_ARRAY_ = 0x8C1A, GL_RGBA8_ = 0x8058, GL_RGBA_ = 0x1908, GL_UNSIGNED_BYTE_ = 0x1401;
    constexpr unsigned int GL_TEXTURE_MIN_FILTER_ = 0x2801, GL_TEXTURE_MAG_FILTER_ = 0x2800, GL_LINEAR_ = 0x2601, GL_LINEAR_MIPMAP_LINEAR_ = 0x2703;
    constexpr unsigned int GL_TEXTURE_WRAP_S_ = 0x2802, GL_TEXTURE_WRAP_T_ = 0x2803, GL_REPEAT_ = 0x2901, GL_TEXTURE_MAX_ANISOTROPY_ = 0x84FE;
    constexpr unsigned int GL_TEXTURE0_ = 0x84C0;

    // --- Seamless procedural patterns (lattice noise with an integer period, so the image tiles) ---
    uint32_t Hash(int x, int y, uint32_t seed) {
        uint32_t h = (uint32_t)x * 374761393u + (uint32_t)y * 668265263u + seed * 2246822519u;
        h = (h ^ (h >> 13)) * 1274126177u;
        return h ^ (h >> 16);
    }
    float Rand01(int x, int y, uint32_t seed) { return (Hash(x, y, seed) & 0xFFFFFF) / 16777215.0f; }
    float Fade(float t) { return t * t * (3.0f - 2.0f * t); }

    // Value noise in [0,1] at (u,v) in [0,1), `cells` lattice cells across, wrapping.
    float ValueNoise(float u, float v, int cells, uint32_t seed) {
        const float x = u * cells, y = v * cells;
        const int x0 = (int)floorf(x), y0 = (int)floorf(y);
        const float fx = Fade(x - x0), fy = Fade(y - y0);
        auto at = [&](int i, int j) { return Rand01(((i % cells) + cells) % cells, ((j % cells) + cells) % cells, seed); };
        const float a = at(x0, y0), b = at(x0 + 1, y0), c = at(x0, y0 + 1), d = at(x0 + 1, y0 + 1);
        return (a * (1 - fx) + b * fx) * (1 - fy) + (c * (1 - fx) + d * fx) * fy;
    }
    float Fbm(float u, float v, int cells, int octaves, uint32_t seed) {
        float sum = 0.0f, amp = 0.5f, norm = 0.0f;
        for (int o = 0; o < octaves; ++o) { sum += amp * ValueNoise(u, v, cells << o, seed + o * 101u); norm += amp; amp *= 0.5f; }
        return sum / norm;
    }
    // Worley (cell) noise: distance to the nearest and second-nearest feature point, and the nearest
    // cell's id, wrapping with `cells` cells across.
    void Worley(float u, float v, int cells, uint32_t seed, float& f1, float& f2, uint32_t& id) {
        const float x = u * cells, y = v * cells;
        const int cx = (int)floorf(x), cy = (int)floorf(y);
        f1 = f2 = 1e9f; id = 0;
        for (int j = -1; j <= 1; ++j)
            for (int i = -1; i <= 1; ++i) {
                const int gx = cx + i, gy = cy + j;
                const int wx = ((gx % cells) + cells) % cells, wy = ((gy % cells) + cells) % cells;
                const float px = gx + Rand01(wx, wy, seed), py = gy + Rand01(wx, wy, seed + 7u);
                const float d = sqrtf((px - x) * (px - x) + (py - y) * (py - y));
                if (d < f1) { f2 = f1; f1 = d; id = Hash(wx, wy, seed + 13u); } else if (d < f2) f2 = d;
            }
    }
    Color Mix(Color a, Color b, float t) {
        t = std::clamp(t, 0.0f, 1.0f);
        return { (unsigned char)(a.r + (b.r - a.r) * t), (unsigned char)(a.g + (b.g - a.g) * t), (unsigned char)(a.b + (b.b - a.b) * t), 255 };
    }
    Color Scale(Color c, float s) {
        return { (unsigned char)std::clamp(c.r * s, 0.0f, 255.0f), (unsigned char)std::clamp(c.g * s, 0.0f, 255.0f), (unsigned char)std::clamp(c.b * s, 0.0f, 255.0f), 255 };
    }
}

const char* TerrainLayerName(int layer) {
    static const char* names[LayerCount] = { "dirt", "sand", "gravel", "grass", "snow", "ice", "stone" };
    return layer >= 0 && layer < LayerCount ? names[layer] : "?";
}

Image GenerateTerrainLayerImage(int layer, int size) {
    Image img = GenImageColor(size, size, BLACK);
    Color* px = (Color*)img.data;
    for (int y = 0; y < size; ++y)
        for (int x = 0; x < size; ++x) {
            const float u = (x + 0.5f) / size, v = (y + 0.5f) / size;
            Color c = BLACK;
            switch (layer) {
                case LayerDirt: { // brown soil: lumpy tone, dark clods, small stones and light grains
                    const float n = Fbm(u, v, 8, 5, 11u);
                    c = Mix({ 74, 52, 34, 255 }, { 128, 94, 62, 255 }, n * 1.4f - 0.2f);
                    if (ValueNoise(u, v, 96, 12u) > 0.82f) c = Scale(c, 0.7f);
                    if (Rand01(x, y, 13u) > 0.985f) c = Scale(c, 1.35f);
                    float f1, f2; uint32_t id;
                    Worley(u, v, 36, 14u, f1, f2, id);
                    if ((id % 9u) == 0 && f1 < 0.3f) { // a pebble in every ninth cell, shaded round, a dark rim
                        const Color stone = Mix({ 96, 90, 84, 255 }, { 150, 142, 130, 255 }, (id >> 8 & 0xFF) / 255.0f);
                        c = Mix(Scale(stone, 1.12f - f1 * 1.2f), Scale(c, 0.6f), std::clamp((f1 - 0.22f) * 12.0f, 0.0f, 1.0f));
                    }
                    break;
                }
                case LayerSand: { // pale tan: fine grain and soft wind ripples
                    const float ripple = 0.5f + 0.5f * sinf((v + 0.08f * Fbm(u, v, 4, 3, 21u)) * 6.2831853f * 14.0f);
                    c = Mix({ 196, 170, 122, 255 }, { 226, 206, 160, 255 }, 0.55f * Fbm(u, v, 16, 3, 22u) + 0.25f * ripple);
                    c = Scale(c, 0.92f + 0.16f * Rand01(x, y, 23u));
                    break;
                }
                case LayerGravel: { // small rounded pebbles of varied tone, shaded, with soil between them
                    float f1, f2; uint32_t id;
                    Worley(u, v, 56, 31u, f1, f2, id);
                    const float tone = 0.8f + 0.4f * ((id & 0xFF) / 255.0f);
                    const float edge = std::clamp((f2 - f1) * 7.0f, 0.0f, 1.0f);
                    Color pebble = Scale(Mix({ 98, 94, 88, 255 }, { 158, 152, 142, 255 }, (id >> 8 & 0xFF) / 255.0f), tone);
                    pebble = Scale(pebble, 1.15f - 0.55f * f1); // rounded: lighter in the middle
                    c = Mix({ 62, 54, 46, 255 }, pebble, edge);
                    c = Scale(c, 0.92f + 0.16f * ValueNoise(u, v, 128, 32u));
                    break;
                }
                case LayerGrass: { // green turf: patchy tone, streaky blades, dry flecks
                    const float patch = Fbm(u, v, 6, 4, 41u);
                    const float blades = ValueNoise(u * 1.0f, v, 160, 42u) * 0.6f + ValueNoise(u, v, 48, 43u) * 0.4f;
                    c = Mix({ 46, 86, 30, 255 }, { 96, 140, 52, 255 }, patch * 0.9f + blades * 0.35f - 0.15f);
                    if (Rand01(x, y, 44u) > 0.992f) c = { 170, 160, 80, 255 };
                    break;
                }
                case LayerSnow: { // near-white, faint blue shading, rare sparkles
                    const float n = Fbm(u, v, 8, 4, 51u);
                    c = Mix({ 214, 222, 236, 255 }, { 248, 250, 255, 255 }, n);
                    if (Rand01(x, y, 52u) > 0.995f) c = { 255, 255, 255, 255 };
                    break;
                }
                case LayerIce: { // pale blue, smooth depth variation, thin white cracks
                    float f1, f2; uint32_t id;
                    Worley(u, v, 10, 61u, f1, f2, id);
                    const float crack = 1.0f - std::clamp((f2 - f1) * 18.0f, 0.0f, 1.0f);
                    c = Mix({ 120, 170, 210, 255 }, { 180, 215, 240, 255 }, Fbm(u, v, 4, 4, 62u));
                    c = Mix(c, { 235, 245, 255, 255 }, crack * 0.8f);
                    break;
                }
                case LayerStone: { // grey rock: warped strata bands, fine grain, thin cracks, lichen
                    const float warp = Fbm(u, v, 3, 4, 73u);
                    const float band = v * 9.0f + warp * 1.6f;                       // 9 wavy layers per tile (seamless)
                    const float strata = 0.5f + 0.5f * sinf(band * 6.2831853f);
                    const float layerTone = Rand01((int)floorf(band) & 15, 0, 75u);  // each layer a little lighter or darker
                    const float grain = Fbm(u, v, 24, 4, 72u);
                    c = Mix({ 92, 91, 95, 255 }, { 158, 154, 148, 255 }, 0.5f * grain + 0.22f * strata + 0.25f * layerTone);
                    float f1, f2; uint32_t id;
                    Worley(u, v, 5, 71u, f1, f2, id);
                    const float crack = 1.0f - std::clamp((f2 - f1) * 70.0f, 0.0f, 1.0f);
                    c = Mix(c, { 52, 51, 55, 255 }, crack * 0.8f);
                    const float lichen = Fbm(u, v, 10, 3, 74u);
                    if (lichen > 0.68f) c = Mix(c, { 118, 128, 86, 255 }, std::min(1.0f, (lichen - 0.68f) * 6.0f) * 0.45f);
                    break;
                }
                default: c = MAGENTA;
            }
            px[y * size + x] = c;
        }
    return img;
}


// --- Detail maps (Plan 5 D: stylized PBR) ---
// Each layer's height is its albedo brightness, softened (a 3x3 blur drops the per-pixel speckle) and
// signed per layer (ice cracks are grooves although they are bright). From it: a tangent-space normal
// (RG), roughness (B: per layer, rougher in the hollows) and ambient occlusion (A: hollows below their
// surroundings darken).
namespace {
    struct LayerSurface { float bump, rough, roughVar, aoStrength, invert; };
    const LayerSurface kSurface[LayerCount] = {
        //        bump  rough  var   ao    invert
        /*dirt*/ { 5.0f, 0.92f, 0.06f, 2.5f, 0.0f },
        /*sand*/ { 3.0f, 0.86f, 0.06f, 1.5f, 0.0f },
        /*grav*/ { 7.0f, 0.78f, 0.15f, 3.5f, 0.0f },
        /*gras*/ { 4.0f, 0.88f, 0.08f, 2.0f, 0.0f },
        /*snow*/ { 2.0f, 0.62f, 0.10f, 1.0f, 0.0f },
        /*ice */ { 2.5f, 0.10f, 0.12f, 0.8f, 1.0f },
        /*ston*/ { 6.0f, 0.68f, 0.20f, 3.0f, 0.0f },
    };

    // Periodic box blur (separable), radius r.
    std::vector<float> BoxBlur(const std::vector<float>& h, int size, int r) {
        std::vector<float> tmp(h.size()), out(h.size());
        const float inv = 1.0f / (2 * r + 1);
        for (int y = 0; y < size; ++y)
            for (int x = 0; x < size; ++x) {
                float s = 0.0f;
                for (int k = -r; k <= r; ++k) s += h[y * size + ((x + k + size) % size)];
                tmp[y * size + x] = s * inv;
            }
        for (int y = 0; y < size; ++y)
            for (int x = 0; x < size; ++x) {
                float s = 0.0f;
                for (int k = -r; k <= r; ++k) s += tmp[((y + k + size) % size) * size + x];
                out[y * size + x] = s * inv;
            }
        return out;
    }
}

Image GenerateTerrainLayerDetail(int layer, const Image& albedo) {
    const int size = albedo.width;
    const LayerSurface& S = kSurface[std::clamp(layer, 0, LayerCount - 1)];
    const Color* px = (const Color*)albedo.data;
    std::vector<float> h(size * size);
    for (int i = 0; i < size * size; ++i) {
        const float l = (0.299f * px[i].r + 0.587f * px[i].g + 0.114f * px[i].b) / 255.0f;
        h[i] = S.invert > 0.5f ? 1.0f - l : l;
    }
    h = BoxBlur(h, size, 1);
    const std::vector<float> wide = BoxBlur(h, size, 6);
    // Normalise the height to its own range so the bump strength means the same for every layer.
    float lo = 1e9f, hi = -1e9f;
    for (float v : h) { lo = std::min(lo, v); hi = std::max(hi, v); }
    const float span = std::max(hi - lo, 1e-3f);
    Image img = GenImageColor(size, size, BLACK);
    Color* out = (Color*)img.data;
    auto H = [&](int x, int y) { return (h[((y + size) % size) * size + ((x + size) % size)] - lo) / span; };
    for (int y = 0; y < size; ++y)
        for (int x = 0; x < size; ++x) {
            // Slope per texel times the bump; the texture covers ~1 texel per 1/size of a tile.
            const float dx = (H(x + 1, y) - H(x - 1, y)) * 0.5f * S.bump;
            const float dy = (H(x, y + 1) - H(x, y - 1)) * 0.5f * S.bump;
            const float inv = 1.0f / sqrtf(dx * dx + dy * dy + 1.0f);
            const float nx = -dx * inv, ny = -dy * inv;
            const float hv = H(x, y);
            const float hollow = std::max(0.0f, (wide[y * size + x] - lo) / span - hv); // below the neighbourhood
            const float rough = std::clamp(S.rough + S.roughVar * (0.5f - hv) * 2.0f, 0.05f, 1.0f);
            const float ao = std::clamp(1.0f - hollow * S.aoStrength, 0.35f, 1.0f);
            out[y * size + x] = { (unsigned char)std::lround((nx * 0.5f + 0.5f) * 255.0f), (unsigned char)std::lround((ny * 0.5f + 0.5f) * 255.0f),
                                  (unsigned char)std::lround(rough * 255.0f), (unsigned char)std::lround(ao * 255.0f) };
        }
    return img;
}

TextureArray LoadTerrainLayers(int size) {
    TextureArray t;
    if (!glad_glTexImage3D || !glad_glGenTextures) { logger.logf("[Terrain] texture arrays unavailable\n"); return t; }
    t.size = size;
    t.layers = LayerCount;
    glad_glGenTextures(1, &t.id);
    glad_glBindTexture(GL_TEXTURE_2D_ARRAY_, t.id);
    glad_glTexImage3D(GL_TEXTURE_2D_ARRAY_, 0, (int)GL_RGBA8_, size, size, LayerCount, 0, GL_RGBA_, GL_UNSIGNED_BYTE_, nullptr);
    std::vector<Image> albedos;
    for (int layer = 0; layer < LayerCount; ++layer) {
        const char* path = TextFormat("assets/textures/terrain/%s.png", TerrainLayerName(layer));
        Image img = FileExists(path) ? LoadImage(path) : GenerateTerrainLayerImage(layer, size);
        ImageFormat(&img, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
        if (img.width != size || img.height != size) ImageResize(&img, size, size);
        glad_glTexSubImage3D(GL_TEXTURE_2D_ARRAY_, 0, 0, 0, layer, size, size, 1, GL_RGBA_, GL_UNSIGNED_BYTE_, img.data);
        albedos.push_back(img);
    }
    glad_glGenerateMipmap(GL_TEXTURE_2D_ARRAY_);
    glad_glTexParameteri(GL_TEXTURE_2D_ARRAY_, GL_TEXTURE_MIN_FILTER_, (int)GL_LINEAR_MIPMAP_LINEAR_);
    glad_glTexParameteri(GL_TEXTURE_2D_ARRAY_, GL_TEXTURE_MAG_FILTER_, (int)GL_LINEAR_);
    glad_glTexParameteri(GL_TEXTURE_2D_ARRAY_, GL_TEXTURE_WRAP_S_, (int)GL_REPEAT_);
    glad_glTexParameteri(GL_TEXTURE_2D_ARRAY_, GL_TEXTURE_WRAP_T_, (int)GL_REPEAT_);
    if (glad_glTexParameterf && GL_Caps().anisotropic) glad_glTexParameterf(GL_TEXTURE_2D_ARRAY_, GL_TEXTURE_MAX_ANISOTROPY_, std::min(8.0f, GL_Caps().maxAnisotropy)); // only where supported (#286)
    // Detail array: normal XY, roughness, AO per layer (assets/textures/terrain/<name>_detail.png overrides).
    const double t0 = GetTime();
    glad_glGenTextures(1, &t.detailId);
    glad_glBindTexture(GL_TEXTURE_2D_ARRAY_, t.detailId);
    glad_glTexImage3D(GL_TEXTURE_2D_ARRAY_, 0, (int)GL_RGBA8_, size, size, LayerCount, 0, GL_RGBA_, GL_UNSIGNED_BYTE_, nullptr);
    for (int layer = 0; layer < LayerCount; ++layer) {
        const char* path = TextFormat("assets/textures/terrain/%s_detail.png", TerrainLayerName(layer));
        Image d = FileExists(path) ? LoadImage(path) : GenerateTerrainLayerDetail(layer, albedos[layer]);
        ImageFormat(&d, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
        if (d.width != size || d.height != size) ImageResize(&d, size, size);
        glad_glTexSubImage3D(GL_TEXTURE_2D_ARRAY_, 0, 0, 0, layer, size, size, 1, GL_RGBA_, GL_UNSIGNED_BYTE_, d.data);
        UnloadImage(d);
        UnloadImage(albedos[layer]);
    }
    glad_glGenerateMipmap(GL_TEXTURE_2D_ARRAY_);
    glad_glTexParameteri(GL_TEXTURE_2D_ARRAY_, GL_TEXTURE_MIN_FILTER_, (int)GL_LINEAR_MIPMAP_LINEAR_);
    glad_glTexParameteri(GL_TEXTURE_2D_ARRAY_, GL_TEXTURE_MAG_FILTER_, (int)GL_LINEAR_);
    glad_glTexParameteri(GL_TEXTURE_2D_ARRAY_, GL_TEXTURE_WRAP_S_, (int)GL_REPEAT_);
    glad_glTexParameteri(GL_TEXTURE_2D_ARRAY_, GL_TEXTURE_WRAP_T_, (int)GL_REPEAT_);
    if (glad_glTexParameterf && GL_Caps().anisotropic) glad_glTexParameterf(GL_TEXTURE_2D_ARRAY_, GL_TEXTURE_MAX_ANISOTROPY_, std::min(8.0f, GL_Caps().maxAnisotropy)); // only where supported (#286)
    glad_glBindTexture(GL_TEXTURE_2D_ARRAY_, 0);
    logger.logf("[Terrain] %d terrain layers (%dx%d) in texture array %u, detail maps %u (%.0f ms)\n", LayerCount, size, size, t.id, t.detailId,
                (GetTime() - t0) * 1000.0);
    return t;
}

void UnloadTextureArray(TextureArray& t) {
    if (t.id && glad_glDeleteTextures) glad_glDeleteTextures(1, &t.id);
    if (t.detailId && glad_glDeleteTextures) glad_glDeleteTextures(1, &t.detailId);
    t = TextureArray{};
}

void BindTextureArray(int unit, const TextureArray& t) {
    glad_glActiveTexture(GL_TEXTURE0_ + unit);
    glad_glBindTexture(GL_TEXTURE_2D_ARRAY_, t.id);
    if (t.detailId) {
        glad_glActiveTexture(GL_TEXTURE0_ + kTerrainDetailUnit);
        glad_glBindTexture(GL_TEXTURE_2D_ARRAY_, t.detailId);
    }
}

void UnbindTextureArray(int unit) {
    glad_glActiveTexture(GL_TEXTURE0_ + unit);
    glad_glBindTexture(GL_TEXTURE_2D_ARRAY_, 0);
    glad_glActiveTexture(GL_TEXTURE0_ + kTerrainDetailUnit);
    glad_glBindTexture(GL_TEXTURE_2D_ARRAY_, 0);
}
