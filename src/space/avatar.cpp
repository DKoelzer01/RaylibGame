#include "avatar.h"
#include "instancebatch.h"
#include "inventory.h"
#include "tools.h"
#include "world.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <unordered_map>

// Remote player avatars (#328): see avatar.h.
namespace {
    // Proportions (m) for the 1.7 m player capsule: feet at 0, hips 0.88, shoulders 1.40, top of the helmet 1.72.
    constexpr float kHip = 0.88f, kThigh = 0.44f, kShin = 0.42f, kTorso = 0.52f, kUpperArm = 0.30f, kForearm = 0.28f;
    constexpr float kHipHalf = 0.11f, kShoulderHalf = 0.21f;
    const Color kSuit = { 188, 194, 204, 255 }, kSuitDark = { 120, 128, 142, 255 }, kVisor = { 40, 70, 110, 255 };
    const Color kBoot = { 70, 72, 78, 255 }, kPack = { 96, 104, 118, 255 };

    struct Anim { float gait = 0.0f; double last = -1.0; float speed = 0.0f; };
    std::unordered_map<uint32_t, Anim> anims;

    Color TierAccent(int tier) {
        switch (tier) {
            case 1: return { 200, 150, 60, 255 };   // bronze-ish
            case 2: return { 90, 150, 210, 255 };
            case 3: return { 120, 200, 120, 255 };
            case 4: return { 200, 110, 200, 255 };
            default: return { 210, 90, 60, 255 };   // the starting suit's orange stripes
        }
    }
    // The suit colour a player picked when joining (#341): the suit takes most of it (a little of the default grey keeps it lit
    // like the rest of the scene); trims are a darker shade of it. Alpha 0: the default suit.
    Color SuitOf(const AvatarState& a) {
        if (a.suitColor.a == 0) return kSuit;
        auto mix = [](unsigned char d, unsigned char c) { return (unsigned char)std::clamp((int)lroundf(d * 0.25f + c * 0.75f), 0, 255); };
        return { mix(kSuit.r, a.suitColor.r), mix(kSuit.g, a.suitColor.g), mix(kSuit.b, a.suitColor.b), 255 };
    }
    Color SuitDarkOf(const AvatarState& a) {
        if (a.suitColor.a == 0) return kSuitDark;
        const Color s = SuitOf(a);
        return { (unsigned char)(s.r * 0.62f), (unsigned char)(s.g * 0.62f), (unsigned char)(s.b * 0.62f), 255 };
    }
    Color ToolTier(int tier) {
        switch (tier) { case 1: return { 135, 128, 120, 255 }; case 2: return { 180, 180, 186, 255 }; case 3: return { 120, 126, 140, 255 }; case 4: return { 205, 180, 215, 255 }; default: return GRAY; }
    }
    // A box hanging from the current joint (its top at the origin), as the animal models do.
    void Limb(float w, float len, float d, Color c) { IB_Cube({ 0, -len * 0.5f, 0 }, { w, len, d }, c); }
    void Box(Vector3 at, Vector3 size, Color c, bool glow = false) { IB_Cube(at, size, c, glow); }

    // What the right hand holds, drawn from the grip (the hand's origin, pointing down the forearm).
    void HeldItem(const std::string& held) {
        if (held.empty()) return;
        if (held == "torch" || held == "fibre_torch") {
            rlPushMatrix(); rlRotatef(80.0f, 1, 0, 0); // upright out of the fist (the hand frame: y back up the forearm)
            Box({ 0, 0.18f, 0 }, { 0.03f, 0.40f, 0.03f }, { 120, 84, 50, 255 });
            Box({ 0, 0.40f, 0 }, { 0.06f, 0.07f, 0.06f }, { 60, 50, 42, 255 });
            Box({ 0, 0.47f, 0 }, { 0.07f, 0.09f, 0.07f }, { 255, 170, 70, 255 }, true);
            rlPopMatrix();
            return;
        }
        const ToolSpec tool = ToolOf(held);
        if (tool.type == ToolType::Hands) { // an item: a small cube of its colour
            const ItemDef* def = itemDB.find(held);
            Box({ 0, -0.06f, 0.04f }, { 0.12f, 0.12f, 0.12f }, def ? def->color : GRAY);
            return;
        }
        rlPushMatrix(); rlRotatef(80.0f, 1, 0, 0); // the handle upright out of the fist (tools are held upright)
        Box({ 0, 0.20f, 0 }, { 0.035f, 0.56f, 0.035f }, { 120, 84, 50, 255 });
        const Color head = ToolTier(tool.tier);
        switch (tool.type) {
            case ToolType::Pickaxe: Box({ 0, 0.46f, 0 }, { 0.05f, 0.05f, 0.40f }, head); break;
            case ToolType::Axe:     Box({ 0, 0.42f, -0.07f }, { 0.025f, 0.14f, 0.15f }, head); break;
            case ToolType::Shovel:  Box({ 0, 0.52f, 0 }, { 0.17f, 0.2f, 0.025f }, head); break;
            case ToolType::Hoe:     Box({ 0, 0.46f, -0.07f }, { 0.16f, 0.03f, 0.14f }, head); break;
            default: break;
        }
        rlPopMatrix();
    }

    // The humanoid in its own frame (x right, y up, z forward; feet at the origin). gait: radians; stride 0..1.
    void DrawBody(const AvatarState& a, float gait, float stride, bool lowDetail) {
        const Color kSuit = SuitOf(a), kSuitDark = SuitDarkOf(a); // (#341) the picked colour, else the defaults
        const Color accent = TierAccent(a.suitTier);
        const bool crouch = a.posture == Player::Posture::Crouch, prone = a.posture == Player::Posture::Prone || a.swimming;
        const bool airborne = !a.grounded && !a.flying && !a.swimming && !prone;
        rlPushMatrix();
        // Prone and swimming: the whole figure lies forward along the ground, its feet where they stood.
        if (prone) { rlTranslatef(0, 0.18f, -0.75f); rlRotatef(90.0f, 1, 0, 0); }
        // Flying: leaning into the motion.
        if (a.flying) rlRotatef(std::min(35.0f, stride * 35.0f), 1, 0, 0);
        const float hipDrop = crouch ? 0.32f : 0.0f;
        const float hipY = kHip - hipDrop;
        const float lean = crouch ? 22.0f : stride * 10.0f;           // forward lean (deg): running, crouching
        if (lowDetail) { // far away: a body, a head and two legs (a few boxes)
            Box({ 0, hipY * 0.5f, 0 }, { 0.3f, hipY, 0.2f }, kSuitDark);
            Box({ 0, hipY + kTorso * 0.5f, 0 }, { 0.44f, kTorso, 0.26f }, kSuit);
            Box({ 0, hipY + kTorso + 0.15f, 0 }, { 0.26f, 0.28f, 0.28f }, kSuit);
            rlPopMatrix();
            return;
        }
        // Legs: thigh and shin from each hip; walking swings them (opposite phase), crouching folds them.
        for (int s = 0; s < 2; ++s) {
            const float side = s == 0 ? 1.0f : -1.0f;
            const float phase = gait + (s == 0 ? 0.0f : PI);
            float hipAng = stride * 38.0f * sinf(phase);                 // + forward
            float knee = std::max(0.0f, -sinf(phase)) * stride * 55.0f;   // the back leg bends as it comes through
            if (crouch) { hipAng = 55.0f + hipAng * 0.4f; knee = 105.0f; }
            if (airborne) { hipAng = 35.0f + 10.0f * side; knee = 60.0f; }
            if (a.flying) { hipAng = -8.0f + 6.0f * side; knee = 12.0f; }
            if (a.swimming) { hipAng = 14.0f * sinf(gait * 2.0f + (s ? PI : 0.0f)); knee = 10.0f; } // a flutter kick
            rlPushMatrix();
            rlTranslatef(side * kHipHalf, hipY, 0);
            rlRotatef(-hipAng, 1, 0, 0);
            Limb(0.15f, kThigh, 0.17f, kSuit);
            rlTranslatef(0, -kThigh, 0);
            rlRotatef(knee, 1, 0, 0);
            Limb(0.13f, kShin, 0.15f, kSuit);
            Box({ 0, -kShin - 0.03f, 0.04f }, { 0.14f, 0.08f, 0.26f }, kBoot);   // the boot
            rlPopMatrix();
        }
        // Torso from the hips, leaning; the arms and head ride on it.
        rlPushMatrix();
        rlTranslatef(0, hipY, 0);
        rlRotatef(lean, 1, 0, 0);
        Box({ 0, 0.06f, 0 }, { 0.34f, 0.16f, 0.22f }, kSuitDark);                   // belt and hips
        Box({ 0, kTorso * 0.55f, 0 }, { 0.42f, kTorso * 0.8f, 0.26f }, kSuit);      // chest
        Box({ 0, kTorso * 0.55f, 0.135f }, { 0.30f, 0.05f, 0.01f }, accent);       // a stripe across the chest
        Box({ 0, kTorso * 0.5f, -0.19f }, { 0.34f, kTorso * 0.7f, 0.14f }, kPack);  // the life-support pack
        if (a.jetpack) { // jet nozzles under the pack, glowing
            Box({ 0.09f, 0.04f, -0.2f }, { 0.07f, 0.12f, 0.07f }, { 255, 160, 60, 255 }, true);
            Box({ -0.09f, 0.04f, -0.2f }, { 0.07f, 0.12f, 0.07f }, { 255, 160, 60, 255 }, true);
        }
        // Head: the helmet with its visor, turned up or down with the look.
        rlPushMatrix();
        rlTranslatef(0, kTorso + 0.02f, 0);
        rlRotatef(-std::clamp(a.lookPitch * RAD2DEG, -60.0f, 60.0f) - (prone ? 60.0f : 0.0f), 1, 0, 0);
        Box({ 0, 0.15f, 0 }, { 0.27f, 0.29f, 0.29f }, kSuit);
        Box({ 0, 0.16f, 0.13f }, { 0.21f, 0.13f, 0.04f }, kVisor);
        Box({ 0, 0.30f, -0.02f }, { 0.12f, 0.03f, 0.18f }, accent);
        if (a.headlamp) Box({ 0, 0.27f, 0.14f }, { 0.07f, 0.05f, 0.03f }, { 255, 245, 210, 255 }, true);
        rlPopMatrix();
        // Arms: swinging against the legs when walking; the right one swings the held tool when acting.
        for (int s = 0; s < 2; ++s) {
            const float side = s == 0 ? 1.0f : -1.0f;
            const bool right = s == 0;
            float shoulder = -stride * 32.0f * sinf(gait + (s == 0 ? 0.0f : PI));  // opposite the leg on its side
            float elbow = 15.0f + stride * 20.0f;
            float spread = 6.0f;
            if (crouch) { shoulder = 25.0f; elbow = 40.0f; }
            if (airborne) { shoulder = 40.0f; spread = 30.0f; elbow = 20.0f; }
            if (a.flying) { shoulder = -10.0f; spread = 12.0f; elbow = 5.0f; }
            if (a.swimming) { shoulder = 90.0f + 60.0f * sinf(gait + (s ? PI : 0.0f)); elbow = 15.0f; } // strokes
            if (prone && !a.swimming) { shoulder = 150.0f; elbow = 20.0f; }              // arms ahead, on the ground
            if (right && !a.held.empty() && a.action == 0 && !a.swimming && !prone) { shoulder = 35.0f; elbow = 55.0f; } // holding it ready
            if (right && a.action > 0) { // a swing: raised back, then down through the stroke
                const float p = std::clamp(a.actionPhase, 0.0f, 1.0f);
                const float up = p < 0.35f ? p / 0.35f : 1.0f - (p - 0.35f) / 0.65f;
                shoulder = (a.action == 3 ? 60.0f : 30.0f) + up * (a.action == 2 ? 110.0f : 130.0f);
                elbow = 20.0f + up * 30.0f;
            }
            rlPushMatrix();
            rlTranslatef(side * kShoulderHalf, kTorso * 0.92f, 0);
            rlRotatef(side * spread, 0, 0, 1);
            rlRotatef(-shoulder, 1, 0, 0);
            Limb(0.11f, kUpperArm, 0.12f, kSuit);
            rlTranslatef(0, -kUpperArm, 0);
            rlRotatef(-elbow, 1, 0, 0);
            Limb(0.10f, kForearm, 0.11f, kSuit);
            Box({ 0, -kForearm - 0.04f, 0 }, { 0.11f, 0.09f, 0.12f }, kSuitDark);    // the glove
            if (right) { rlTranslatef(0, -kForearm - 0.05f, 0); HeldItem(a.held); }
            rlPopMatrix();
        }
        rlPopMatrix(); // torso
        rlPopMatrix(); // whole figure
    }

    // RAYLIBGAME_AVATAR_GLTF: a glTF / glb body instead of the boxes (static for now: no animation), scaled to 1.7 m.
    Model gltfModel{}; bool gltfTried = false; float gltfScale = 1.0f;
    bool GltfReady() {
        if (gltfTried) return gltfModel.meshCount > 0;
        gltfTried = true;
        const char* path = getenv("RAYLIBGAME_AVATAR_GLTF");
        if (!path || !*path || !FileExists(path)) return false;
        gltfModel = LoadModel(path);
        if (gltfModel.meshCount <= 0) return false;
        const BoundingBox bb = GetModelBoundingBox(gltfModel);
        gltfScale = 1.7f / std::max(bb.max.y - bb.min.y, 0.01f);
        return true;
    }
}

AvatarState Avatar_FromPlayer(const Player& p, uint32_t session) {
    AvatarState a;
    a.session = session;
    a.position = p.position;
    a.up = p.up;
    Vector3 f = Vector3Subtract(p.forward, Vector3Scale(p.up, Vector3DotProduct(p.forward, p.up)));
    a.facing = Vector3LengthSqr(f) > 1e-6f ? Vector3Normalize(f) : Vector3{ 0, 0, 1 };
    a.lookPitch = asinf(std::clamp(Vector3DotProduct(Vector3Normalize(p.forward), p.up), -1.0f, 1.0f));
    a.velocity = p.velocity;
    a.posture = p.posture;
    a.grounded = p.grounded;
    a.flying = p.flyMode || p.hovering;
    a.swimming = p.inWater && !p.grounded;
    return a;
}

void Avatar_Draw3D(Scene& world, const std::vector<AvatarState>& avatars) {
    if (avatars.empty()) return;
    const double now = GameClock();
    const bool gltf = GltfReady();
    for (const AvatarState& a : avatars) {
        const Vector3 base = (a.position - world.renderOrigin).toFloat();
        const float dist = Vector3Distance(base, world.camera.position);
        if (dist > 160.0f) continue; // like the animals
        // The gait: phase advances with the speed along the ground (a stride ~1.4 m walking, longer running).
        Anim& an = anims[a.session];
        const float dt = an.last < 0.0 ? 0.0f : (float)std::clamp(now - an.last, 0.0, 0.1);
        an.last = now;
        const Vector3 flat = Vector3Subtract(a.velocity, Vector3Scale(a.up, Vector3DotProduct(a.velocity, a.up)));
        const float speed = Vector3Length(a.flying ? a.velocity : flat);
        an.speed += (speed - an.speed) * std::min(1.0f, dt * 8.0f);  // eased: no snap starting or stopping
        an.gait += dt * an.speed / 1.4f * PI;
        const float stride = std::clamp(an.speed / 6.0f, 0.0f, 1.0f); // 0 standing .. 1 running
        const Vector3 u = Vector3Normalize(a.up), f = Vector3Normalize(a.facing);
        // (#2140) Out of view: a 2 m figure with a held tool (after the gait update, so the stride stays continuous).
        if (!world.inView(Vector3Add(base, u), 1.8f)) continue;
        const Vector3 r = Vector3CrossProduct(u, f);
        const Matrix m = { r.x, u.x, f.x, base.x, r.y, u.y, f.y, base.y, r.z, u.z, f.z, base.z, 0, 0, 0, 1 };
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(m));
        if (gltf) {
            rlDrawRenderBatchActive();
            DrawModel(gltfModel, { 0, 0, 0 }, gltfScale, WHITE);
        } else {
            DrawBody(a, an.gait, stride, dist > 60.0f);
        }
        rlPopMatrix();
    }
}

void Avatar_DrawNameplates(Scene& world, const std::vector<AvatarState>& avatars) {
    for (const AvatarState& a : avatars) {
        if (a.name.empty()) continue;
        const Vector3 base = (a.position - world.renderOrigin).toFloat();
        const Vector3 head = Vector3Add(base, Vector3Scale(Vector3Normalize(a.up), a.posture == Player::Posture::Prone ? 0.8f : a.posture == Player::Posture::Crouch ? 1.55f : 2.0f));
        const float dist = Vector3Distance(head, world.camera.position);
        const float fade = 1.0f - std::clamp((dist - 20.0f) / 40.0f, 0.0f, 1.0f); // full to 20 m, gone at 60 m
        if (fade <= 0.0f) continue;
        const Vector3 toHead = Vector3Subtract(head, world.camera.position);
        if (Vector3DotProduct(toHead, Vector3Subtract(world.camera.target, world.camera.position)) <= 0.0f) continue; // behind
        const Vector2 sp = GetWorldToScreen(head, world.camera);
        const int fs = 20, w = MeasureText(a.name.c_str(), fs);
        DrawRectangle((int)sp.x - w / 2 - 6, (int)sp.y - fs / 2 - 3, w + 12, fs + 6, Fade(BLACK, 0.45f * fade));
        DrawText(a.name.c_str(), (int)sp.x - w / 2, (int)sp.y - fs / 2, fs, Fade(RAYWHITE, fade));
        if (a.suitColor.a > 0) // (#341) a thin bar of the player's colour under the name
            DrawRectangle((int)sp.x - w / 2 - 6, (int)sp.y + fs / 2 + 3, w + 12, 3, Fade(Color{ a.suitColor.r, a.suitColor.g, a.suitColor.b, 255 }, fade));
    }
}

// (#822) The inventory screen's player preview: one avatar alone in `target`, lit like the world where it stands (sun,
// shadows, ambient), seen from `eye` towards `at` (render space). Call outside the world's 3D pass and outside ImGui.
void Avatar_DrawPreview(Scene& world, const AvatarState& a, RenderTexture2D& target, Vector3 eye, Vector3 at, Vector3 up) {
    BeginTextureMode(target);
    ClearBackground(Color{ 18, 25, 36, 255 });
    Camera3D cam{};
    cam.position = eye; cam.target = at; cam.up = up; cam.fovy = 32.0f; cam.projection = CAMERA_PERSPECTIVE;
    const Camera3D saved = world.camera;
    world.camera = cam; // (the batch lights from the camera's place: the preview's)
    // (#913 review) A studio light, not the world's: at night the world-lit figure was near-black on the dark panel. A key
    // light from the camera's side and above, bright ambient, no fog; the world's lighting comes back after.
    const Vector3 sun0 = world.frameSun, amb0 = world.frameAmbient, ambUp0 = world.frameAmbientUp, fade0 = world.sunFade, light0 = world.lightDirection;
    const float fog0 = world.frameFogDensity, moon0 = world.moonLight;
    const Vector3 key = Vector3Normalize(Vector3Add(Vector3Normalize(Vector3Subtract(eye, at)), Vector3Scale(up, 0.8f)));
    world.frameSun = key; world.lightDirection = key; world.sunFade = { 0, 0, 0 }; world.moonLight = 0.0f;
    world.frameAmbient = { 0.32f, 0.33f, 0.36f }; world.frameAmbientUp = { 0.40f, 0.42f, 0.46f }; world.frameFogDensity = 0.0f;
    BeginMode3D(cam);
    const std::vector<AvatarState> one{ a };
    Avatar_Draw3D(world, one);
    InstanceBatch_Flush(world);
    EndMode3D();
    world.frameSun = sun0; world.frameAmbient = amb0; world.frameAmbientUp = ambUp0; world.sunFade = fade0; world.lightDirection = light0;
    world.frameFogDensity = fog0; world.moonLight = moon0;
    world.camera = saved;
    EndTextureMode();
}
