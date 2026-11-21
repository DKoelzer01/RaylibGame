#include "fishing.h"
#include "audio.h"
#include "gameui.h"
#include "inventory.h"
#include "particles.h"
#include "player.h"
#include "scene.h"
#include "terraintool.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <rlgl.h>
#include <cmath>

namespace {
    enum class St { Idle, Waiting, Bite };
    St st = St::Idle;
    Planetoid* bodyP = nullptr;
    Vector3 bobber{}, castFrom{};   // body-local
    float wait = 0.0f, biteLeft = 0.0f, fishClock = 0.0f;
    uint32_t rng = 0x2545F491u;
    float Rand() { rng = rng * 1664525u + 1013904223u; return (rng >> 8) * (1.0f / 16777216.0f); }

    Vector3 PlayerLocal(Scene& world, const Planetoid& b) { return b.toLocalOf(player.position, world.renderOrigin); }

    // The first water along the view within reach (render space -> body-local), as the bucket finds it.
    bool WaterAlongView(Scene& world, Planetoid& b, Vector3& outLocal) {
        const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
        for (float s = 0.8f; s <= 14.0f; s += 0.25f) {
            const Vector3 q = Vector3Add(world.camera.position, Vector3Scale(dir, s));
            if (!IsUnderwater(world, q)) continue;
            outLocal = b.toLocal(q);
            return true;
        }
        return false;
    }

    void Splash(Planetoid& b, Vector3 at, int n) {
        ParticleBurst s;
        s.kind = ParticleKind::Dust; s.color = { 200, 225, 245, 200 }; s.count = n;
        s.size = 0.12f; s.life = 0.7f; s.speed = 1.6f; s.spread = 0.6f;
        Particles_Emit(b, at, s);
    }

    void Cast(Planetoid& b, Vector3 water, Vector3 from) {
        bodyP = &b; bobber = water; castFrom = from;
        st = St::Waiting;
        wait = 4.0f + 8.0f * Rand();
        Splash(b, water, 4);
        Audio_Play(Sfx::Pickup, 0.4f, 0.6f);
        GameUI_Toast("Cast: wait for a bite, then RMB");
    }

    // RMB while fishing: land a biting fish, or reel in empty.
    std::string Reel() {
        std::string got;
        if (st == St::Bite && bodyP) {
            got = Fishing_CatchAt(*bodyP, bobber);
            inventory.add(got, 1);
            const ItemDef* d = itemDB.find(got);
            GameUI_Toast(TextFormat("Caught a %s!", d ? d->name.c_str() : got.c_str()));
            Audio_Play(Sfx::Pickup, 0.7f, 1.2f);
            Splash(*bodyP, bobber, 8);
        } else if (st == St::Waiting) GameUI_Toast("Nothing on the line");
        st = St::Idle;
        return got;
    }

    void Tick(float dt) {
        fishClock += dt;
        if (st == St::Waiting && (wait -= dt) <= 0.0f) {
            st = St::Bite; biteLeft = 1.5f;
            if (bodyP) Splash(*bodyP, bobber, 6);
            Audio_Play(Sfx::Error, 0.35f, 1.6f);
            GameUI_Toast("Bite! RMB now");
        } else if (st == St::Bite && (biteLeft -= dt) <= 0.0f) {
            st = St::Waiting; wait = 3.0f + 5.0f * Rand();
            GameUI_Toast("It got away");
        }
    }
}

std::string Fishing_CatchAt(const Planetoid& body, Vector3 local) {
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    if (!pf) return "raw_trout";
    const float r = Vector3Length(local);
    const Vector3 d = Vector3Scale(local, 1.0f / std::max(r, 1.0f));
    const float e = pf->elevation(d);
    if (e < -0.5f) return "raw_cod";                                     // the sea's floor under it: sea fish
    const Biome b = pf->biomeAt(d, e);
    if (b == Biome::Jungle || b == Biome::Swamp || b == Biome::Savanna) return "raw_catfish";
    return "raw_trout";
}

void Fishing_Reset() { st = St::Idle; bodyP = nullptr; }

void Fishing_Update(Scene& world, float dt, bool inputEnabled) {
    if (terrainTool.heldTool != "fishing_rod") { if (st != St::Idle) st = St::Idle; return; }
    Planetoid* b = player.nearest;
    if (!b) return;
    if (st != St::Idle && (bodyP != b || Vector3Distance(PlayerLocal(world, *b), castFrom) > 4.0f)) { st = St::Idle; GameUI_Toast("The line is lost"); }
    Tick(dt);
    if (!inputEnabled || !IsMouseButtonPressed(MOUSE_BUTTON_RIGHT)) return;
    if (st != St::Idle) { Reel(); return; }
    Vector3 water;
    if (WaterAlongView(world, *b, water)) Cast(*b, water, PlayerLocal(world, *b));
    else GameUI_Toast("Cast into water: the sea, a river or a lake");
}

void Fishing_Draw3D(Scene& world) {
    if (st == St::Idle || !bodyP || terrainTool.heldTool != "fishing_rod") return;
    const Vector3 fwd = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
    const Vector3 right = Vector3Normalize(Vector3CrossProduct(fwd, world.camera.up));
    const Vector3 tip = Vector3Add(world.camera.position, Vector3Add(Vector3Scale(fwd, 0.9f), Vector3Add(Vector3Scale(right, 0.35f), Vector3Scale(world.camera.up, 0.25f))));
    const Vector3 up = Vector3Normalize(bobber);
    const float bob = st == St::Bite ? -0.12f + 0.05f * sinf(fishClock * 30.0f) : 0.03f * sinf(fishClock * 3.0f);
    const Vector3 at = bodyP->toWorld(Vector3Add(bobber, Vector3Scale(up, 0.05f + bob)));
    rlDrawRenderBatchActive();
    DrawLine3D(tip, at, Fade(RAYWHITE, 0.7f));
    DrawSphere(at, 0.07f, st == St::Bite ? Color{ 255, 60, 40, 255 } : Color{ 230, 60, 50, 255 });
    DrawSphere(Vector3Add(at, Vector3Scale(bodyP->dirToWorld(up), 0.06f)), 0.045f, RAYWHITE);
    rlDrawRenderBatchActive();
}

void Fishing_DrawHUD() {
    if (st == St::Idle || terrainTool.heldTool != "fishing_rod") return;
    const int cx = GetScreenWidth() / 2, cy = GetScreenHeight() / 2;
    if (st == St::Bite) DrawText("BITE!  RMB", cx - 40, cy - 60, 28, Color{ 255, 220, 80, 255 });
    else DrawText("fishing...  (RMB reels in)", cx + 14, cy - 34, 18, Fade(RAYWHITE, 0.7f));
}

std::string Fishing_TestCast(Scene& world, Planetoid& body, Vector3 waterLocal, bool early) {
    Cast(body, waterLocal, PlayerLocal(world, body));
    for (float t = 0.0f; t < 15.0f; t += 0.1f) {
        Tick(0.1f);
        if (early && t >= 1.0f) return Reel();
        if (st == St::Bite) return Reel();
    }
    st = St::Idle;
    return "";
}
