#include "flares.h"
#include "audio.h"
#include "fluid.h"
#include "gameui.h"
#include "inventory.h"
#include "net.h"
#include "player.h"
#include "scene.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"
#include "instancebatch.h"
#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>

namespace {
constexpr float kThrowSpeed = 12.0f; // m/s off the hand (GUESS: a few metres' lob, ~15 m on the flat)
constexpr float kLift = 0.25f;       // m: the light's cell is taken this far above the flare, off the floor it lies on

struct Flare {
    uint64_t bodyId = 0;
    Vector3 pos{}, vel{};      // body-local (the frame of the body it lies on: it turns with the ground)
    float age = 0.0f, burn = kFlareBurnSeconds;
    bool landed = false;       // touched the ground: from then on it lights the grid
};
std::vector<Flare> flares;

Planetoid* BodyOf(uint64_t id) { return Planetoid::find(id); }

// The reach now: full, then down to 0 over the fade (whole metres, so the fade is ~16 relights, not one a frame).
float ReachOf(const Flare& f) {
    const float left = f.burn - f.age;
    if (left <= 0.0f) return 0.0f;
    const float fade = std::min(kFlareFadeSeconds, 0.5f * f.burn);
    return std::ceil(kFlareReach * std::min(1.0f, left / fade));
}
}

void Flares_Throw(Planetoid& body, Vector3 from, Vector3 velocity, float burnSeconds) {
    Flare f;
    f.bodyId = body.id;
    f.pos = from; f.vel = velocity;
    if (burnSeconds > 0.0f) f.burn = burnSeconds;
    flares.push_back(f);
}

bool Flares_ThrowFromPlayer(Scene& world, float burnSeconds) {
    Planetoid* body = player.nearest;
    if (!body || !body->field || inventory.count("flare") <= 0) return false;
    if (gNetRole == NetRole::Client) { // (#2062) not replicated yet: a client's inventory is the server's
        GameUI_Toast("Flares work in single player and on a host for now");
        Audio_Play(Sfx::Error, 0.3f);
        return false;
    }
    const Vector3 eye = body->toLocal(world.camera.position);
    const Vector3 look = body->dirToLocal(Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position)));
    const Vector3 up = Vector3Normalize(eye);
    // A little upward lob on the view, from just in front of the eye.
    const Vector3 dir = Vector3Normalize(Vector3Add(look, Vector3Scale(up, 0.2f)));
    if (!inventory.remove("flare", 1)) return false;
    Flares_Throw(*body, Vector3Add(eye, Vector3Scale(look, 0.4f)), Vector3Scale(dir, kThrowSpeed), burnSeconds);
    Audio_Play(Sfx::Build, 0.4f, 1.6f);
    return true;
}

void Flares_Update(Scene&, float dt) {
    if (dt <= 0.0f) return;
    for (size_t k = 0; k < flares.size();) {
        Flare& f = flares[k];
        Planetoid* body = BodyOf(f.bodyId);
        f.age += dt;
        bool gone = !body || !body->field || f.age >= f.burn;
        if (!gone) {
            // Ballistic under the body's pull; a step into the ground is not taken and the flare stops there. A resting
            // flare keeps trying the next step down, so it settles onto the floor and falls again if the floor is dug.
            const Vector3 up = Vector3Normalize(f.pos);
            f.vel = Vector3Add(f.vel, Vector3Scale(up, -9.81f * Fluid::GravityScale(*body) * dt));
            const Vector3 next = Vector3Add(f.pos, Vector3Scale(f.vel, dt));
            if (body->densityAt(next) >= ISO_LEVEL) { f.vel = {}; f.landed = true; }
            else f.pos = next;
        }
        if (gone) { flares[k] = flares.back(); flares.pop_back(); }
        else ++k;
    }
}

void Flares_Draw3D(Scene& world) {
    for (const Flare& f : flares) {
        Planetoid* body = BodyOf(f.bodyId);
        if (!body) continue;
        // Under the body's localModel (double translation): a float toWorld 60 km out steps ~4 mm and jitters (#271/#284).
        const Matrix m = body->localModel(DVec3(f.pos), world.renderOrigin);
        if (Vector3Distance({ m.m12, m.m13, m.m14 }, world.camera.position) > 150.0f) continue;
        const Vector3 up = Vector3Normalize(f.pos);
        // Lying on its side once down; end over end in flight.
        Vector3 side = Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 });
        side = Vector3Normalize(side);
        if (!f.landed) side = Vector3RotateByAxisAngle(side, up, f.age * 9.0f);
        const Quaternion q = QuaternionFromVector3ToVector3({ 0, 1, 0 }, side);
        const float burning = f.landed ? std::min(1.0f, ReachOf(f) / kFlareReach) : 1.0f;
        const float flick = 0.85f + 0.15f * sinf(f.age * 31.0f) * sinf(f.age * 17.0f);
        const unsigned char glow = (unsigned char)(255.0f * std::clamp(burning * flick, 0.15f, 1.0f));
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(m));
        rlMultMatrixf(MatrixToFloat(QuaternionToMatrix(q)));
        IB_Cylinder({ 0, -0.11f, 0 }, 0.025f, 0.025f, 0.2f, Color{ 170, 40, 35, 255 });                        // the red tube
        IB_Cube({ 0, 0.11f, 0 }, { 0.06f, 0.05f, 0.06f }, Color{ glow, (unsigned char)(glow * 0.45f), (unsigned char)(glow * 0.3f), 255 }, true); // the burning end
        rlPopMatrix();
    }
}

void Flares_Clear() { flares.clear(); }

void Flares_Lights(const Planetoid& body, std::vector<FlareLight>& out) {
    for (const Flare& f : flares) {
        if (f.bodyId != body.id || !f.landed) continue; // in flight it only glows: a moving source would relight every frame
        const float reach = ReachOf(f);
        if (reach <= 0.0f) continue;
        out.push_back({ Vector3Add(f.pos, Vector3Scale(Vector3Normalize(f.pos), kLift)), reach });
    }
}

int Flares_Count() { return (int)flares.size(); }

bool Flares_First(Vector3& local, bool& resting, float& reach) {
    if (flares.empty()) return false;
    local = flares[0].pos; resting = flares[0].landed; reach = ReachOf(flares[0]);
    return true;
}
