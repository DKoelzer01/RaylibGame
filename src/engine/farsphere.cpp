#include "farsphere.h"
#include "jobsystem.h"
#include "trees.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>

namespace {
    constexpr int kRowsPerJob = 16;

    Vector3 TexelDir(int x, int y, int W, int H) {
        const float lon = ((x + 0.5f) / W - 0.5f) * 2.0f * PI;
        const float lat = ((y + 0.5f) / H - 0.5f) * PI;
        return { cosf(lat) * cosf(lon), sinf(lat), cosf(lat) * sinf(lon) };
    }
    unsigned char Enc(float v) { return (unsigned char)std::clamp((int)lroundf(v * 255.0f), 0, 255); }
    // Linear -> sRGB bytes (dark colours keep their precision); the shader decodes.
    void PutColor(unsigned char* px, Vector3 c, bool water) {
        px[0] = Enc(powf(std::max(c.x, 0.0f), 1.0f / 2.2f));
        px[1] = Enc(powf(std::max(c.y, 0.0f), 1.0f / 2.2f));
        px[2] = Enc(powf(std::max(c.z, 0.0f), 1.0f / 2.2f));
        px[3] = water ? 0 : 255;
    }
    Vector3 Mix(Vector3 a, Vector3 b, float t) { return Vector3Lerp(a, b, std::clamp(t, 0.0f, 1.0f)); }

    // Average linear colours of the terrain layers (planet palette) and the moon's dust and rock.
    const Vector3 kGrass = { 0.13f, 0.27f, 0.06f }, kDirt = { 0.28f, 0.19f, 0.11f }, kSand = { 0.60f, 0.52f, 0.34f };
    const Vector3 kGravel = { 0.34f, 0.32f, 0.29f }, kSnow = { 0.84f, 0.87f, 0.92f }, kIce = { 0.58f, 0.72f, 0.86f };
    const Vector3 kStone = { 0.38f, 0.37f, 0.36f }, kMoonDust = { 0.42f, 0.40f, 0.38f }, kMoonRock = { 0.28f, 0.27f, 0.26f };
}

struct FarSphere::Bake {
    std::shared_ptr<const PlanetField> field;
    bool foliage = false;
    unsigned int foliageSeed = 0;
    int W = 0, H = 0;
    std::vector<float> elev;
    std::vector<unsigned char> albedo, normals;
    std::atomic<int> elevRowsDone{ 0 }, shadeRowsDone{ 0 };
};

static void ShadeRows(FarSphere::Bake& b, int y0, int y1);

FarSphere::FarSphere(std::shared_ptr<const PlanetField> field, bool foliage, unsigned int foliageSeed) : bake(std::make_shared<Bake>()) {
    Bake& b = *bake;
    b.field = std::move(field);
    b.foliage = foliage;
    b.foliageSeed = foliageSeed;
    b.W = b.field->radius > 30000.0f ? 1024 : 512;
    b.H = b.W / 2;
    b.elev.assign((size_t)b.W * b.H, 0.0f);
    b.albedo.assign((size_t)b.W * b.H * 4, 0);
    b.normals.assign((size_t)b.W * b.H * 4, 0);
    // Pass 1: heights, a block of rows per job. The job finishing the last block queues pass 2
    // (colours and normals need the neighbouring rows).
    auto bk = bake;
    for (int y0 = 0; y0 < b.H; y0 += kRowsPerJob) {
        WorldJobs().submit(2000.0f, [bk, y0] {
            Bake& b = *bk;
            const float footprint = 2.0f * PI * b.field->radius / b.W;
            const int y1 = std::min(b.H, y0 + kRowsPerJob);
            for (int y = y0; y < y1; ++y)
                for (int x = 0; x < b.W; ++x) b.elev[(size_t)y * b.W + x] = b.field->elevation(TexelDir(x, y, b.W, b.H), footprint);
            if (b.elevRowsDone.fetch_add(y1 - y0) + (y1 - y0) == b.H)
                for (int s0 = 0; s0 < b.H; s0 += kRowsPerJob)
                    WorldJobs().submit(2000.0f, [bk, s0] { ShadeRows(*bk, s0, std::min(bk->H, s0 + kRowsPerJob)); });
        });
    }
}

static void ShadeRows(FarSphere::Bake& b, int y0, int y1) {
    const PlanetField& f = *b.field;
    const float R = f.radius;
    const int W = b.W, H = b.H;
    auto E = [&](int x, int y) { x = (x + W) % W; y = std::clamp(y, 0, H - 1); return b.elev[(size_t)y * W + x]; };
    for (int y = y0; y < y1; ++y)
        for (int x = 0; x < W; ++x) {
            const Vector3 up = TexelDir(x, y, W, H);
            const float e = E(x, y);
            const float lat = fabsf(up.y);
            // Normal from the height grid: east-west and north-south slopes.
            const float cl = std::max(cosf(((y + 0.5f) / H - 0.5f) * PI), 0.02f);
            const float dx = 2.0f * PI * R * cl / W, dy = PI * R / H;
            const float dedx = (E(x + 1, y) - E(x - 1, y)) / (2.0f * dx), dedy = (E(x, y + 1) - E(x, y - 1)) / (2.0f * dy);
            Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, up));
            if (Vector3Length(east) < 0.5f) east = { 1, 0, 0 };
            const Vector3 north = Vector3CrossProduct(up, east);
            Vector3 n = Vector3Normalize(Vector3Subtract(up, Vector3Add(Vector3Scale(east, dedx), Vector3Scale(north, dedy))));
            const float slope = 1.0f - Vector3DotProduct(n, up);
            Vector3 c;
            bool water = false;
            if (f.moon) {
                c = Mix(kMoonDust, kMoonRock, slope * 4.0f);
            } else if (e < 0.0f && lat < 0.9f) {
                // Sea: clear in the shallows, deep blue offshore (as the water shader).
                const float clarity = expf(e / 5.0f);
                c = Mix({ 0.01f, 0.05f, 0.13f }, { 0.03f, 0.22f, 0.24f }, clarity * 0.6f);
                water = true;
                n = up;
            } else {
                const Vector3 p = Vector3Scale(up, R + e);
                const float w = f.climate.fractal(2, p.x, p.y, p.z);
                if (lat > 0.9f && e < 0.0f) c = kIce;                                     // frozen polar sea
                else if (lat > 0.93f + 0.03f * w || e > 1350.0f + 100.0f * w) c = kIce;
                else if (lat > 0.82f + 0.08f * w || e > 1050.0f + 150.0f * w) c = kSnow;
                else if (slope > 0.45f) c = kStone;
                else if (e < 6.0f + 4.0f * w) c = kSand;
                else if (slope > 0.3f || e > 800.0f + 100.0f * w) c = kGravel;
                else {
                    c = lat > 0.62f + 0.06f * w ? kDirt : kGrass;
                    if (b.foliage) {
                        const float canopy = TreeCanopy(up, R, e, b.foliageSeed);
                        const Vector3 crowns = lat < 0.35f ? Mix({ 0.22f, 0.27f, 0.07f }, { 0.09f, 0.27f, 0.05f }, (lat - 0.15f) / 0.2f)
                                                           : Mix({ 0.09f, 0.27f, 0.05f }, { 0.05f, 0.17f, 0.08f }, (lat - 0.45f) / 0.25f);
                        c = Mix(c, Vector3Scale(crowns, 0.5f), canopy);
                    }
                }
            }
            unsigned char* px = &b.albedo[((size_t)y * W + x) * 4];
            PutColor(px, c, water);
            unsigned char* pn = &b.normals[((size_t)y * W + x) * 4];
            pn[0] = Enc(n.x * 0.5f + 0.5f); pn[1] = Enc(n.y * 0.5f + 0.5f); pn[2] = Enc(n.z * 0.5f + 0.5f); pn[3] = 255;
        }
    b.shadeRowsDone.fetch_add(y1 - y0);
}

FarSphere::~FarSphere() {
    if (uploaded) {
        UnloadTexture(albedo);
        UnloadTexture(normals);
        UnloadMesh(mesh);
    }
}

bool FarSphere::ready() {
    if (uploaded) return true;
    if (bake->shadeRowsDone.load() < bake->H) return false;
    auto upload = [&](std::vector<unsigned char>& px) {
        Image img{ px.data(), bake->W, bake->H, 1, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8 };
        Texture2D t = LoadTextureFromImage(img);
        GenTextureMipmaps(&t);
        SetTextureFilter(t, TEXTURE_FILTER_TRILINEAR);
        SetTextureWrap(t, TEXTURE_WRAP_REPEAT);
        return t;
    };
    albedo = upload(bake->albedo);
    normals = upload(bake->normals);
    mesh = GenMeshSphere(1.0f, 64, 128);
    uploaded = true;
    std::vector<unsigned char>().swap(bake->albedo);
    std::vector<unsigned char>().swap(bake->normals);
    std::vector<float>().swap(bake->elev);
    return true;
}

void FarSphere::draw(const Planetoid& body, const RenderContext& ctx, Shader shader) {
    if (!uploaded) return;
    static Material material = LoadMaterialDefault();
    material.shader = shader;
    material.maps[MATERIAL_MAP_ALBEDO].texture = albedo;
    material.maps[MATERIAL_MAP_NORMAL].texture = normals;
    const float r = body.radius * body.scale;
    const Matrix model = MatrixMultiply(MatrixScale(r, r, r), body.localModel(DVec3(0.0, 0.0, 0.0), ctx.renderOrigin));
    DrawMesh(mesh, material, model);
    if (ctx.drawCalls) ++*ctx.drawCalls;
}
