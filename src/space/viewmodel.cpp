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
#include <vector>

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
            return s == kHotbarTerrainTool || inventory.count(s) <= 0 ? std::string() : s; // only what is carried (#272)
        }
        if (Hotbar_Selected() >= 0) { // a torch in hand (#558: it showed an empty glove, selecting it clears heldTool)
            const std::string& s = Hotbar_Slot(Hotbar_Selected());
            if (Hotbar_IsLight(s) && inventory.count(s) > 0) return s;
        }
        return inventory.count(terrainTool.heldTool) > 0 ? terrainTool.heldTool : std::string();
    }

    // The held torch's flame (#558): a few additive sprites in the hands' own view space, rising from the head.
    struct Flame { Vector3 p, v; float age, life, size; };
    std::vector<Flame> flames;
    float flameSpawn = 0.0f;
    bool IsTorch(const std::string& id) { return id == "torch" || id == "fibre_torch"; }

    Color Scale(Color c, float k) {
        return { (unsigned char)std::clamp((int)(c.r * k), 0, 255), (unsigned char)std::clamp((int)(c.g * k), 0, 255), (unsigned char)std::clamp((int)(c.b * k), 0, 255), 255 };
    }

    // A box of `size` centred at the origin of `m`, its faces shaded by a light from above left.
    void Box(const Matrix& m, Vector3 size, Color col, float light) {
        const Vector3 h = Vector3Scale(size, 0.5f);
        Vector3 v[8];
        for (int i = 0; i < 8; ++i) v[i] = Vector3Transform({ (i & 1) ? h.x : -h.x, (i & 2) ? h.y : -h.y, (i & 4) ? h.z : -h.z }, m);
        static const int faces[6][4] = { { 0, 4, 6, 2 }, { 1, 3, 7, 5 }, { 0, 1, 5, 4 }, { 2, 6, 7, 3 }, { 0, 2, 3, 1 }, { 4, 5, 7, 6 } }; // counter-clockwise from outside (#273: they were inside out)
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
    swing = std::min(1.0f, swing + dt / std::max(terrainTool.swingTime, 0.1f)); // one swing per cooldown
    const std::string now = Held();
    if (now != held) { held = now; dip = 0.0f; flames.clear(); }
    dip = std::min(1.0f, dip + dt / 0.25f);
    const float lowerTarget = (GameUI_InventoryOpen() || Nav_MapOpen()) ? 1.0f : 0.0f;
    lowered += (lowerTarget - lowered) * std::min(1.0f, 8.0f * dt);
}

void Viewmodel_Draw(Scene& world) {
    visible = settings.showHands && !player.flyMode && !Rocket_Boarded() && !telescope.active && !survival.dead;
    if (!visible) return;
    // Light (#203): it came only from the sun's height (0.28 at night at least), so the hands glowed in caves, in the
    // shade and at night. Now: sky light scaled by how buried the eye is (the terrain's own rule), sunlight only when
    // no terrain stands between the eye and the sun (one ray), the headlamp; eased so walking into shade doesn't pop.
    const Vector3 sunDir = Vector3Normalize(world.sunDirection);
    const float sunUp = Vector3DotProduct(sunDir, player.up);
    const float day = std::clamp((sunUp + 0.15f) / 0.4f, 0.0f, 1.0f);
    float sky = 1.0f, sunVis = day > 0.0f ? 1.0f : 0.0f;
    if (Planetoid* b = player.nearest) {
        const Vector3 eye = world.camera.position, l = b->toLocal(eye);
        const float len = Vector3Length(l);
        if (b->field && len > 1.0f) b->field->materialSky(l, 1.0f, 1.0f, Vector3Scale(l, 1.0f / len), sky);
        TerrainHit th;
        if (sunVis > 0.0f && RaycastTerrain(world, Ray{ Vector3Add(eye, Vector3Scale(sunDir, 0.3f)), sunDir }, 120.0f, th)) sunVis = 0.0f;
    }
    float target = (0.08f + 0.27f * day) * sky + 0.65f * day * sunVis * sky;
    if (survival.headlamp) target = std::max(target, 0.8f);
    static float light = -1.0f;
    const float dtl = GameFrameTime();
    light = light < 0.0f || dtl <= 0.0f ? target : light + (target - light) * (1.0f - expf(-dtl / 0.15f));
    const float handLight = IsTorch(held) ? std::max(light, 0.55f) : light; // a held torch lights its own hand (#558)

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
    Box(local({ 0, 0, 0.11f }), { 0.10f, 0.09f, 0.07f }, cuff, handLight);          // suit cuff
    Box(local({ 0, 0, 0 }), { 0.085f, 0.07f, 0.15f }, glove, handLight);             // glove
    const ToolSpec tool = ToolOf(held);
    if (IsTorch(held)) {
        // A thin stick held upright (it was an item cube, #558), a wrapped head, and the particle flame over it.
        const bool fibre = held == "fibre_torch";
        const float lit = std::max(light, 0.55f); // its own flame lights the stick and the glove (they read black at night)
        const Matrix stick = local({ 0, 0.1f, -0.05f }, { -10, 0, 8 });
        Box(stick, { 0.016f, 0.3f, 0.016f }, fibre ? Color{ 150, 135, 80, 255 } : Color{ 120, 84, 50, 255 }, lit);
        const Matrix head = MatrixMultiply(At({ 0, 0.15f, 0 }), stick);
        Box(head, { 0.03f, 0.05f, 0.03f }, fibre ? Color{ 170, 150, 90, 255 } : Color{ 92, 70, 52, 255 }, lit);
        const Vector3 foot = Vector3Transform({ 0, 0.18f, 0 }, stick); // the flame's foot in view space
        const float dt = GameFrameTime();
        flameSpawn += dt * (fibre ? 45.0f : 60.0f);
        while (flameSpawn >= 1.0f && flames.size() < 64) {
            flameSpawn -= 1.0f;
            const float a = GetRandomValue(0, 1000) * 0.001f * 2.0f * PI, r = GetRandomValue(0, 1000) * 0.001f * 0.008f;
            flames.push_back({ { foot.x + cosf(a) * r, foot.y, foot.z + sinf(a) * r },
                               { cosf(a) * 0.01f, 0.1f + GetRandomValue(0, 100) * 0.0006f, sinf(a) * 0.01f },
                               0.0f, 0.22f + GetRandomValue(0, 100) * 0.0015f, (fibre ? 0.026f : 0.032f) * (0.8f + GetRandomValue(0, 100) * 0.004f) });
        }
        if (flameSpawn > 1.0f) flameSpawn = 1.0f;
        for (Flame& f : flames) { f.age += dt; f.p = Vector3Add(f.p, Vector3Scale(f.v, dt)); f.v.y += 0.15f * dt; }
        flames.erase(std::remove_if(flames.begin(), flames.end(), [](const Flame& f) { return f.age >= f.life; }), flames.end());
        rlDrawRenderBatchActive();
        rlSetBlendMode(BLEND_ADDITIVE);
        rlDisableDepthMask();
        rlBegin(RL_TRIANGLES); // each sprite a soft disc: its colour at the centre, black (nothing, additive) at the rim
        for (const Flame& f : flames) {
            const float k = f.age / f.life;                                    // 0 new .. 1 gone
            const float sz = f.size * (1.0f - 0.7f * k);                       // shrinks as it rises
            const float fade = (1.0f - k) * std::min(1.0f, f.age * 25.0f);
            const unsigned char r = (unsigned char)(255 * fade), g = (unsigned char)((220 - 150 * k) * fade), b = (unsigned char)((120 - 110 * k) * fade); // yellow -> deep orange
            constexpr int kSeg = 8;
            for (int j = 0; j < kSeg; ++j) {
                const float a0 = 2.0f * PI * j / kSeg, a1 = 2.0f * PI * (j + 1) / kSeg;
                rlColor4ub(r, g, b, 255); rlVertex3f(f.p.x, f.p.y, f.p.z);
                rlColor4ub(0, 0, 0, 255);
                rlVertex3f(f.p.x + cosf(a0) * sz, f.p.y + sinf(a0) * sz * 1.4f, f.p.z); // a little taller than wide
                rlVertex3f(f.p.x + cosf(a1) * sz, f.p.y + sinf(a1) * sz * 1.4f, f.p.z);
            }
        }
        rlEnd();
        rlDrawRenderBatchActive();
        rlEnableDepthMask();
        rlSetBlendMode(BLEND_ALPHA);
    } else if (held.empty() || tool.type == ToolType::Hands) {
        if (held.empty()) {
            Box(local({ 0, 0.005f, -0.1f }, { -12, 0, 0 }), { 0.08f, 0.03f, 0.07f }, glove, light); // open fingers
            Box(local({ -0.05f, 0.01f, -0.03f }, { 0, 0, 30 }), { 0.025f, 0.025f, 0.06f }, glove, light); // thumb
        } else {
            // An item: a small cube of its colour in the hand.
            const ItemDef* def = itemDB.find(held);
            Color col = def ? def->color : GRAY;
            float lit = light;
            // A hot (or radioactive) stack glows (#580): its colour pulled toward the glow, lit by itself.
            const int sel = Hotbar_Selected();
            if (sel >= 0) {
                const InvSlot& cell = inventory.grid()[Inventory::kGrid + sel];
                if (cell.id == held) {
                    const ItemGlow g = Item_Glow(cell);
                    if (g.strength > 0.0f) {
                        const float k = g.strength;
                        col = { (unsigned char)(col.r + (g.color.r - col.r) * k), (unsigned char)(col.g + (g.color.g - col.g) * k), (unsigned char)(col.b + (g.color.b - col.b) * k), 255 };
                        lit = std::max(light, 0.35f + 0.65f * k);
                    }
                }
            }
            Box(local({ 0, 0.07f, -0.08f }, { 10, 20, 0 }), { 0.11f, 0.11f, 0.11f }, col, lit);
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
                Box(onHandle({ 0, 0.2f, 0 }), { 0.045f, 0.04f, 0.32f }, head, light); // the head runs along the view (#273)
                Box(onHandle({ 0, 0.18f, 0.17f }, { 25, 0, 0 }), { 0.035f, 0.03f, 0.06f }, head, light);
                Box(onHandle({ 0, 0.18f, -0.17f }, { -25, 0, 0 }), { 0.035f, 0.03f, 0.06f }, head, light);
                break;
            case ToolType::Axe:
                Box(onHandle({ 0, 0.17f, -0.06f }), { 0.02f, 0.1f, 0.12f }, head, light); // the blade points away (#273)
                break;
            case ToolType::Shovel:
                Box(onHandle({ 0, 0.27f, 0 }), { 0.15f, 0.17f, 0.018f }, head, light);
                break;
            case ToolType::Hoe: // a flat blade at right angles to the top of the handle
                Box(onHandle({ 0, 0.2f, -0.06f }), { 0.14f, 0.025f, 0.12f }, head, light); // the blade reaches away from you (#273)
                break;
            default: break;
        }
    }

    rlDrawRenderBatchActive(); // the quads are batched: draw them before culling comes back on (#273)
    rlEnableBackfaceCulling();
    EndMode3D();
    rlSetClipPlanes(nearDefault, farDefault);
}
