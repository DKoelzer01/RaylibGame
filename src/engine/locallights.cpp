#include "locallights.h"
#include "external/glad.h"
#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace {
    GLuint lightTex = 0, tileTex = 0;
    int tilesX = 0, tilesY = 0, uploaded = 0, maxPerTile = 0;
    std::vector<float> lightData;
    std::vector<int> tileData;
    constexpr int kLightUnit = 13, kTileUnit = 14; // past the material maps (raylib uses 0-10) and the shadow slots (6-7)
}

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
}

int LocalLights_Count() { return uploaded; }
int LocalLights_MaxPerTile() { return maxPerTile; }
