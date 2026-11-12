#include "viewmodel.h"
#include "buildtool.h"
#include "gameui.h"
#include "hotbar.h"
#include "inventory.h"
#include "menus.h"
#include "navigation.h"
#include "player.h"
#include "rocket.h"
#include "survival.h"
#include "telescope.h"
#include "terraintool.h"
#include "tools.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <string>

extern "C" { typedef void (*PfnClearVm)(unsigned int); extern PfnClearVm glad_glClear; }

namespace {
    // Bob: phase advances half a turn every 1.8 m walked (it was 0.75 m: two to three times too
    // busy in play); the hand traces a figure eight.
    constexpr float kStep = 1.8f;
    constexpr float kBobUp = 0.015f, kBobSide = 0.010f;
    float phase = 0.0f, amp = 0.0f;
    Vector3 bob{ 0, 0, 0 };
    // Swing (a mining stroke or a chop) and dip (switching what is held), 0..1 progress.
    int seenSwings = 0;
    float swing = 1.0f, dip = 1.0f, lowered = 0.0f;
    std::string held;            // what the hand holds: "" hands, a tool id, or an item id (a cube)
    bool visible = false;

    std::string Held() {
        if (buildTool.active) {
            const std::string& s = Hotbar_Slot(Hotbar_Selected());
            return s == kHotbarTerrainTool ? std::string() : s;
        }
        return inventory.count(terrainTool.heldTool) > 0 ? terrainTool.heldTool : std::string();
    }

    Color Scale(Color c, float k) {
        return { (unsigned char)std::clamp((int)(c.r * k), 0, 255), (unsigned char)std::clamp((int)(c.g * k), 0, 255), (unsigned char)std::clamp((int)(c.b * k), 0, 255), 255 };
    }

    // A box of `size` centred at the origin of `m`, its faces shaded by a light from above left.
    void Box(const Matrix& m, Vector3 size, Color col, float light) {
        const Vector3 h = Vector3Scale(size, 0.5f);
        Vector3 v[8];
        for (int i = 0; i < 8; ++i) v[i] = Vector3Transform({ (i & 1) ? h.x : -h.x, (i & 2) ? h.y : -h.y, (i & 4) ? h.z : -h.z }, m);
        static const int faces[6][4] = { { 0, 2, 6, 4 }, { 1, 5, 7, 3 }, { 0, 4, 5, 1 }, { 2, 3, 7, 6 }, { 0, 1, 3, 2 }, { 4, 6, 7, 5 } };
        static const Vector3 normals[6] = { { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 } };
        const Vector3 lightDir = Vector3Normalize({ -0.45f, 0.8f, 0.4f });
        Matrix rot = m; rot.m12 = rot.m13 = rot.m14 = 0.0f;
        rlBegin(RL_QUADS);
        for (int f = 0; f < 6; ++f) {
            const Vector3 n = Vector3Normalize(Vector3Transform(normals[f], rot));
            const float k = light * (0.55f + 0.45f * std::max(0.0f, Vector3DotProduct(n, lightDir)));
            const Color c = Scale(col, k);
            rlColor4ub(c.r, c.g, c.b, 255);
            for (int j = 0; j < 4; ++j) rlVertex3f(v[faces[f][j]].x, v[faces[f][j]].y, v[faces[f][j]].z);
        }
        rlEnd();
    }
    Matrix At(Vector3 p, Vector3 rotDeg = { 0, 0, 0 }) {
        return MatrixMultiply(MatrixRotateXYZ(Vector3Scale(rotDeg, DEG2RAD)), MatrixTranslate(p.x, p.y, p.z));
    }

    Color TierColour(int tier) {
        switch (tier) {
            case 1: return { 135, 128, 120, 255 };
            case 2: return { 180, 180, 186, 255 };
            case 3: return { 120, 126, 140, 255 };
            case 4: return { 205, 180, 215, 255 };
            default: return GRAY;
        }
    }
}

Vector3 Viewmodel_BobOffset() { return bob; }
bool Viewmodel_Visible() { return visible; }

void Viewmodel_Update(Scene& world, float dt) {
    (void)world;
    // Speed along the ground (grounded velocity is relative to the surface).
    const Vector3 flat = Vector3Subtract(player.velocity, Vector3Scale(player.up, Vector3DotProduct(player.velocity, player.up)));
    const float speed = player.grounded && !player.flyMode ? Vector3Length(flat) : 0.0f;
    phase += PI * speed * dt / kStep;
    if (phase > 2.0f * PI * 1000.0f) phase = fmodf(phase, 2.0f * PI);
    const float setting = settings.handBob == 2 ? 0.0f : settings.handBob == 1 ? 0.5f : 1.0f;
    const float target = std::clamp(speed / player.walkSpeed, 0.0f, 1.0f) * setting * (settings.reduceMotion ? 0.5f : 1.0f);
    amp += (target - amp) * std::min(1.0f, 6.0f * dt); // eases in and out: no snap when stopping
    if (amp < 1e-4f && target == 0.0f) amp = 0.0f;
    bob = { kBobSide * sinf(phase) * amp, -kBobUp * sinf(phase) * sinf(phase) * amp, 0.0f };
    // Swing on each stroke or chop; dip when what is held changes.
    if (terrainTool.swings != seenSwings) { seenSwings = terrainTool.swings; swing = 0.0f; }
    swing = std::min(1.0f, swing + dt / 0.28f);
    const std::string now = Held();
    if (now != held) { held = now; dip = 0.0f; }
    dip = std::min(1.0f, dip + dt / 0.25f);
    const float lowerTarget = (GameUI_InventoryOpen() || Nav_MapOpen()) ? 1.0f : 0.0f;
    lowered += (lowerTarget - lowered) * std::min(1.0f, 8.0f * dt);
}

void Viewmodel_Draw(Scene& world) {
    visible = settings.showHands && !player.flyMode && !Rocket_Boarded() && !telescope.active && !survival.dead;
    if (!visible) return;
    // Light: daylight from the sun's height over the player's horizon, or the headlamp.
    const float sunUp = Vector3DotProduct(Vector3Normalize(world.sunDirection), player.up);
    float light = 0.28f + 0.72f * std::clamp((sunUp + 0.15f) / 0.4f, 0.0f, 1.0f);
    if (survival.headlamp) light = std::max(light, 0.8f);

    rlDrawRenderBatchActive();
    if (glad_glClear) glad_glClear(0x00000100); // depth only: the hands draw over the world, never into it
    const double nearDefault = rlGetCullDistanceNear(), farDefault = rlGetCullDistanceFar();
    rlSetClipPlanes(0.01, 10.0);
    Camera3D cam{};
    cam.position = { 0, 0, 0 };
    cam.target = { 0, 0, -1 };
    cam.up = { 0, 1, 0 };
    cam.fovy = 55.0f; // fixed: the hands don't stretch with the FOV setting
    cam.projection = CAMERA_PERSPECTIVE;
    BeginMode3D(cam);
    rlDisableBackfaceCulling();

    // The hand's pose: low right, bobbing, swinging, dipping, lowered while a screen is open.
    const float s = sinf(PI * swing), d = sinf(PI * dip);
    const Vector3 base = { 0.27f + bob.x, -0.235f + bob.y - 0.12f * d - 0.3f * lowered, -0.52f - 0.06f * s };
    const Matrix hand = MatrixMultiply(MatrixRotateXYZ({ DEG2RAD * (8.0f - 55.0f * s), DEG2RAD * -12.0f, DEG2RAD * 4.0f }), MatrixTranslate(base.x, base.y, base.z));
    auto local = [&](Vector3 p, Vector3 rotDeg = { 0, 0, 0 }) { return MatrixMultiply(At(p, rotDeg), hand); };

    const Color glove = { 176, 186, 198, 255 }, cuff = { 120, 130, 145, 255 };
    Box(local({ 0, 0, 0.11f }), { 0.10f, 0.09f, 0.07f }, cuff, light);          // suit cuff
    Box(local({ 0, 0, 0 }), { 0.085f, 0.07f, 0.15f }, glove, light);             // glove
    const ToolSpec tool = ToolOf(held);
    if (held.empty() || tool.type == ToolType::Hands) {
        if (held.empty()) {
            Box(local({ 0, 0.005f, -0.1f }, { -12, 0, 0 }), { 0.08f, 0.03f, 0.07f }, glove, light); // open fingers
            Box(local({ -0.05f, 0.01f, -0.03f }, { 0, 0, 30 }), { 0.025f, 0.025f, 0.06f }, glove, light); // thumb
        } else {
            // An item: a small cube of its colour in the hand.
            const ItemDef* def = itemDB.find(held);
            Box(local({ 0, 0.07f, -0.08f }, { 10, 20, 0 }), { 0.11f, 0.11f, 0.11f }, def ? def->color : GRAY, light);
        }
    } else {
        // A tool: a wooden handle rising nearly upright out of the fist (leaning a little forward and
        // in), its head across the top. It leaned 52 degrees forward and read as lying along the view.
        const Color wood = { 120, 84, 50, 255 }, head = TierColour(tool.tier);
        // Scaled to 0.65: full size, the head filled the upper right of the view.
        const Matrix handle = MatrixMultiply(MatrixScale(0.65f, 0.65f, 0.65f), local({ 0, 0.09f, -0.07f }, { -14, 0, 10 }));
        Box(handle, { 0.028f, 0.42f, 0.028f }, wood, light);
        auto onHandle = [&](Vector3 p, Vector3 rotDeg = { 0, 0, 0 }) { return MatrixMultiply(At(p, rotDeg), handle); };
        switch (tool.type) {
            case ToolType::Pickaxe:
                Box(onHandle({ 0, 0.2f, 0 }), { 0.32f, 0.04f, 0.045f }, head, light);
                Box(onHandle({ 0.17f, 0.18f, 0 }, { 0, 0, -25 }), { 0.06f, 0.03f, 0.035f }, head, light);
                Box(onHandle({ -0.17f, 0.18f, 0 }, { 0, 0, 25 }), { 0.06f, 0.03f, 0.035f }, head, light);
                break;
            case ToolType::Axe:
                Box(onHandle({ 0.06f, 0.17f, 0 }), { 0.12f, 0.1f, 0.02f }, head, light);
                break;
            case ToolType::Shovel:
                Box(onHandle({ 0, 0.27f, 0 }), { 0.15f, 0.17f, 0.018f }, head, light);
                break;
            case ToolType::Hoe: // a flat blade at right angles to the top of the handle
                Box(onHandle({ 0, 0.2f, 0.06f }), { 0.14f, 0.025f, 0.12f }, head, light);
                break;
            default: break;
        }
    }

    rlEnableBackfaceCulling();
    EndMode3D();
    rlSetClipPlanes(nearDefault, farDefault);
}
