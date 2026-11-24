#include "locallights.h"
#include "lightgrid.h"
#include "gpumem.h"
#include "external/glad.h"
#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <unordered_set>

namespace {
    GLuint lightTex = 0, tileTex = 0;
    int tilesX = 0, tilesY = 0, uploaded = 0, maxPerTile = 0;
    std::vector<float> lightData;
    std::vector<int> tileData;
    constexpr int kLightUnit = 13, kTileUnit = 14; // past the material maps (raylib uses 0-10) and the shadow slots (6-7)

    // (#1739) The light grid window: 8^3 bricks of 16^3 cells, brick b in slot b mod 8 on each axis, so the texture wraps
    // (GL_REPEAT) the same way the body-local cells do and moving the camera re-sends only the bricks that came into view.
    constexpr int kGridUnit = 18;               // 13-14 local lights, 15 the haze LUT, 16-17 scene.cpp
    constexpr int kGridBricks = 8, kGridCells = kGridBricks * 16;
    constexpr int kGridUploadsPerFrame = 96;    // 1.5 MB a frame at most; the rest follow on the next frames
    GLuint gridTex = 0;
    bool gridOn = false, gridUsable = true;
    float gridMat[16] = {};
    float gridInfo[4] = {};
    int gridSent = 0;
    struct GridSlot { int bx = 0, by = 0, bz = 0; bool valid = false; };
    GridSlot gridSlots[kGridBricks * kGridBricks * kGridBricks];
    std::unordered_set<uint64_t> gridPending;   // bricks whose light changed and are not sent yet
    std::vector<int> gridDirty;
    uint64_t BrickKey(int bx, int by, int bz) {
        const uint64_t m = (1u << 21) - 1;
        return ((uint64_t)(bx & m)) | ((uint64_t)(by & m) << 21) | ((uint64_t)(bz & m) << 42);
    }
    int PosMod(int v, int m) { const int r = v % m; return r < 0 ? r + m : r; }
}

static_assert(sizeof(LightCell) == 4, "the grid uploads LightCell arrays as GL_RGBA bytes");

void LocalLights_UploadGrid(LightGrid* grid, const int cam[3], const int anchor[3], const Matrix& gridFromRender, float gain) {
    gridSent = 0;
    gridOn = false;
    if (!grid || gain <= 0.0f || !gridUsable) return;
    if (!gridTex) {
        GLint units = 0;
        glGetIntegerv(GL_MAX_TEXTURE_IMAGE_UNITS, &units);
        if (units <= kGridUnit) { gridUsable = false; TraceLog(LOG_WARNING, "light grid: %d fragment texture units, needs %d: off", (int)units, kGridUnit + 1); return; }
        glGenTextures(1, &gridTex);
        glBindTexture(GL_TEXTURE_3D, gridTex);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        std::vector<uint8_t> zero((size_t)kGridCells * kGridCells * kGridCells * 4, 0);
        GpuMemTag gm_("light grid"); // VRAM accounting (#302)
        glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA8, kGridCells, kGridCells, kGridCells, 0, GL_RGBA, GL_UNSIGNED_BYTE, zero.data());
        glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
        // Every slot starts not valid: the first frames send the whole window (lit bricks lit before the texture existed too).
    }
    grid->TakeDirty(gridDirty);
    for (size_t i = 0; i + 2 < gridDirty.size(); i += 3) gridPending.insert(BrickKey(gridDirty[i], gridDirty[i + 1], gridDirty[i + 2]));
    static const LightCell zero[LightGrid::kBrickCells] = {};
    const int cb[3] = { cam[0] >> 4, cam[1] >> 4, cam[2] >> 4 };
    glBindTexture(GL_TEXTURE_3D, gridTex);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    // Nearest bricks first, so a budget-limited frame fills in from the camera outward.
    for (int ring = 0; ring <= kGridBricks / 2; ++ring)
        for (int dz = -ring; dz <= ring; ++dz)
            for (int dy = -ring; dy <= ring; ++dy)
                for (int dx = -ring; dx <= ring; ++dx) {
                    if (std::max(std::abs(dx), std::max(std::abs(dy), std::abs(dz))) != ring) continue;
                    if (dx >= kGridBricks / 2 || dy >= kGridBricks / 2 || dz >= kGridBricks / 2) continue; // the window: -4..+3
                    if (gridSent >= kGridUploadsPerFrame) continue;
                    const int bx = cb[0] + dx, by = cb[1] + dy, bz = cb[2] + dz;
                    const int sx = PosMod(bx, kGridBricks), sy = PosMod(by, kGridBricks), sz = PosMod(bz, kGridBricks);
                    GridSlot& s = gridSlots[sx + sy * kGridBricks + sz * kGridBricks * kGridBricks];
                    const uint64_t key = BrickKey(bx, by, bz);
                    const bool here = s.valid && s.bx == bx && s.by == by && s.bz == bz;
                    if (here && !gridPending.count(key)) continue;
                    const LightCell* data = grid->BrickLight(bx, by, bz);
                    if (!data && here) { gridPending.erase(key); continue; }
                    glTexSubImage3D(GL_TEXTURE_3D, 0, sx * 16, sy * 16, sz * 16, 16, 16, 16, GL_RGBA, GL_UNSIGNED_BYTE, data ? data : zero);
                    s = { bx, by, bz, true };
                    gridPending.erase(key);
                    ++gridSent;
                }
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
    // A changed brick outside the window is sent when it comes into view (its slot holds another brick by then).
    for (auto it = gridPending.begin(); it != gridPending.end();) {
        const int bx = (int)((int64_t)(*it << 43) >> 43), by = (int)((int64_t)(*it >> 21 << 43) >> 43), bz = (int)((int64_t)(*it >> 42 << 43) >> 43);
        const bool in = bx - cb[0] >= -kGridBricks / 2 && bx - cb[0] < kGridBricks / 2 && by - cb[1] >= -kGridBricks / 2 && by - cb[1] < kGridBricks / 2 &&
                        bz - cb[2] >= -kGridBricks / 2 && bz - cb[2] < kGridBricks / 2;
        it = in ? std::next(it) : gridPending.erase(it);
    }
    glActiveTexture(GL_TEXTURE0 + kGridUnit); glBindTexture(GL_TEXTURE_3D, gridTex);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, 0);
    // Render space -> texture coordinates (body-local cells minus the anchor, / 128); the window's centre and the gain.
    const Matrix m = MatrixMultiply(gridFromRender, MatrixScale(1.0f / kGridCells, 1.0f / kGridCells, 1.0f / kGridCells));
    const float16 f = MatrixToFloatV(m);
    std::copy(f.v, f.v + 16, gridMat);
    for (int a = 0; a < 3; ++a) gridInfo[a] = (float)(cb[a] * 16 - anchor[a]) / kGridCells;
    gridInfo[3] = gain;
    gridOn = true;
}

int LocalLights_GridBricksUploaded() { return gridSent; }

void LocalLights_Upload(const std::vector<LocalLight>& lights, Vector3 cameraPos) {
    GLint vp[4] = { 0, 0, 0, 0 };
    glGetIntegerv(GL_VIEWPORT, vp);
    const int w = std::max(1, (int)vp[2]), h = std::max(1, (int)vp[3]);
    const int tx = (w + kLightTileSize - 1) / kLightTileSize, ty = (h + kLightTileSize - 1) / kLightTileSize;
    // The nearest lights that can reach the view.
    std::vector<std::pair<float, int>> order;
    order.reserve(lights.size());
    for (int i = 0; i < (int)lights.size(); ++i) {
        const float d = Vector3Distance(lights[i].position, cameraPos) - lights[i].range;
        if (d < 120.0f && lights[i].intensity > 0.0f) order.push_back({ d, i });
    }
    static const bool on = !(getenv("RAYLIBGAME_LOCAL_LIGHTS") && getenv("RAYLIBGAME_LOCAL_LIGHTS")[0] == '0'); // A/B
    if (!on) order.clear();
    std::sort(order.begin(), order.end());
    if ((int)order.size() > kMaxLocalLights) order.resize(kMaxLocalLights);
    uploaded = (int)order.size();
    static bool wasEmpty = false;
    if (uploaded == 0 && wasEmpty && lightTex && tileTex && tx == tilesX && ty == tilesY) { maxPerTile = 0; return; } // nothing near, again: the shaders skip on count 0
    wasEmpty = uploaded == 0;
    lightData.assign(kMaxLocalLights * 8, 0.0f);
    tileData.assign((size_t)tx * (kLightsPerTile + 1) * ty, 0);
    maxPerTile = 0;
    const Matrix viewProj = MatrixMultiply(rlGetMatrixModelview(), rlGetMatrixProjection());
    for (int k = 0; k < uploaded; ++k) {
        const LocalLight& L = lights[order[k].second];
        float* o = &lightData[k * 8];
        o[0] = L.position.x; o[1] = L.position.y; o[2] = L.position.z; o[3] = L.range;
        o[4] = L.color.x; o[5] = L.color.y; o[6] = L.color.z; o[7] = L.intensity;
        // The sphere's screen rectangle from its bounding box corners; a corner behind the eye: the whole screen.
        float x0 = 1e9f, y0 = 1e9f, x1 = -1e9f, y1 = -1e9f;
        bool behind = false;
        for (int c = 0; c < 8 && !behind; ++c) {
            const Vector3 p = { L.position.x + ((c & 1) ? L.range : -L.range), L.position.y + ((c & 2) ? L.range : -L.range), L.position.z + ((c & 4) ? L.range : -L.range) };
            const Quaternion q = QuaternionTransform({ p.x, p.y, p.z, 1.0f }, viewProj);
            if (q.w <= 1e-3f) { behind = true; break; }
            const float nx = q.x / q.w, ny = q.y / q.w;
            x0 = std::min(x0, nx); x1 = std::max(x1, nx); y0 = std::min(y0, ny); y1 = std::max(y1, ny);
        }
        int ax = 0, bx = tx - 1, ay = 0, by = ty - 1;
        if (!behind) {
            if (x1 < -1.0f || x0 > 1.0f || y1 < -1.0f || y0 > 1.0f) continue; // off screen
            ax = std::clamp((int)((x0 * 0.5f + 0.5f) * w) / kLightTileSize, 0, tx - 1);
            bx = std::clamp((int)((x1 * 0.5f + 0.5f) * w) / kLightTileSize, 0, tx - 1);
            ay = std::clamp((int)((y0 * 0.5f + 0.5f) * h) / kLightTileSize, 0, ty - 1); // gl_FragCoord.y counts up from the bottom, as NDC y
            by = std::clamp((int)((y1 * 0.5f + 0.5f) * h) / kLightTileSize, 0, ty - 1);
        }
        for (int y = ay; y <= by; ++y)
            for (int x = ax; x <= bx; ++x) {
                int* t = &tileData[((size_t)y * tx + x) * (kLightsPerTile + 1)];
                if (t[0] < kLightsPerTile) t[1 + t[0]++] = k; // nearest first: a crowded tile keeps the nearest
                maxPerTile = std::max(maxPerTile, t[0]);
            }
    }
    if (!lightTex) {
        glGenTextures(1, &lightTex);
        glBindTexture(GL_TEXTURE_2D, lightTex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        GpuMemTag gm73_("local lights"); // VRAM accounting (#302)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, kMaxLocalLights * 2, 1, 0, GL_RGBA, GL_FLOAT, nullptr);
    }
    glBindTexture(GL_TEXTURE_2D, lightTex);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, kMaxLocalLights * 2, 1, GL_RGBA, GL_FLOAT, lightData.data());
    if (!tileTex || tx != tilesX || ty != tilesY) {
        if (tileTex) glDeleteTextures(1, &tileTex);
        glGenTextures(1, &tileTex);
        glBindTexture(GL_TEXTURE_2D, tileTex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        GpuMemTag gm83_("local lights"); // VRAM accounting (#302)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32I, tx * (kLightsPerTile + 1), ty, 0, GL_RED_INTEGER, GL_INT, nullptr);
        tilesX = tx; tilesY = ty;
    }
    glBindTexture(GL_TEXTURE_2D, tileTex);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, tx * (kLightsPerTile + 1), ty, GL_RED_INTEGER, GL_INT, tileData.data());
    glBindTexture(GL_TEXTURE_2D, 0);
    // Bound once for the frame: nothing else uses these units.
    glActiveTexture(GL_TEXTURE0 + kLightUnit); glBindTexture(GL_TEXTURE_2D, lightTex);
    glActiveTexture(GL_TEXTURE0 + kTileUnit); glBindTexture(GL_TEXTURE_2D, tileTex);
    glActiveTexture(GL_TEXTURE0);
}

void LocalLights_Bind(Shader& sh) {
    const int lu = kLightUnit, tu = kTileUnit;
    const int info[4] = { tilesX, kLightTileSize, uploaded, kLightsPerTile + 1 };
    const int locL = GetShaderLocation(sh, "localLightTex"), locT = GetShaderLocation(sh, "localTileTex"), locI = GetShaderLocation(sh, "localTileInfo");
    if (locL >= 0) SetShaderValue(sh, locL, &lu, SHADER_UNIFORM_INT);
    if (locT >= 0) SetShaderValue(sh, locT, &tu, SHADER_UNIFORM_INT);
    if (locI >= 0) SetShaderValue(sh, locI, info, SHADER_UNIFORM_IVEC4);
    // (#1739) The light grid. Its sampler gets unit 18 even while off: two sampler types on unit 0 fail on strict drivers.
    const int gu = kGridUnit;
    const float off[4] = { 0, 0, 0, 0 };
    const int locG = GetShaderLocation(sh, "lightGridTex"), locM = GetShaderLocation(sh, "lightGridFromRender"), locN = GetShaderLocation(sh, "lightGridInfo");
    if (locG >= 0) SetShaderValue(sh, locG, &gu, SHADER_UNIFORM_INT);
    if (locM >= 0) glProgramUniformMatrix4fv(sh.id, locM, 1, GL_FALSE, gridMat);
    if (locN >= 0) SetShaderValue(sh, locN, gridOn ? gridInfo : off, SHADER_UNIFORM_VEC4);
}

int LocalLights_Count() { return uploaded; }
int LocalLights_MaxPerTile() { return maxPerTile; }
