#include "footsteps.h"

#include "audio.h"
#include "avatar.h"
#include <raymath.h>
#include "scene.h"
#include "world.h"
#include "worldtypes.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <random>
#include <string>
#include <unordered_map>

const char* StepSurface_Sound(StepSurface s) {
    static const char* kNames[] = { "step_stone", "step_dirt", "step_grass", "step_sand", "step_gravel", "step_snow",
                                    "step_wood", "step_metal", "step_water" };
    return (int)s < (int)StepSurface::Count ? kNames[(int)s] : "step_stone";
}

StepSurface StepSurfaceOfVoxel(VoxelMaterial m) {
    switch (m) {
        case VoxelMaterial::Grass: return StepSurface::Grass;
        case VoxelMaterial::Sand: return StepSurface::Sand;
        case VoxelMaterial::Gravel: case VoxelMaterial::Pebbles: case VoxelMaterial::Regolith: return StepSurface::Gravel;
        case VoxelMaterial::Snow: return StepSurface::Snow;
        case VoxelMaterial::Clay: case VoxelMaterial::Silt: case VoxelMaterial::Mud: case VoxelMaterial::Farmland: return StepSurface::Dirt;
        default: return StepSurface::Stone; // rock and its strata, ore, ice
    }
}

StepSurface StepSurfaceOfPart(const Part& p) {
    if (p.type == PartType::Farmland) return StepSurface::Dirt;
    if (p.type == PartType::Fence || p.type == PartType::Gate || p.type == PartType::Bed) return StepSurface::Wood;
    if (!IsShape(p.type)) return StepSurface::Metal; // machines, rocket parts, conveyors, pipes
    const std::string item = GetBuildMaterial(p.material).item;
    auto has = [&](const char* s) { return item.find(s) != std::string::npos; };
    if (has("planks") || has("_log")) return StepSurface::Wood;
    if (has("iron") || has("plate") || has("heat_shield")) return StepSurface::Metal;
    return StepSurface::Stone; // brick, concrete, plaster, glass, gem blocks, ice
}

float StepLoudness(StepSurface s, Player::Posture posture, bool running) {
    static const float kSurface[] = { 1.0f, 0.85f, 0.75f, 0.7f, 1.2f, 0.9f, 1.1f, 1.3f, 1.2f }; // as StepSurface
    const float pace = posture == Player::Posture::Prone ? 0.25f : posture == Player::Posture::Crouch ? 0.45f : (running ? 1.6f : 1.0f);
    return pace * kSurface[(int)s < (int)StepSurface::Count ? (int)s : 0];
}

StepTracker::Event StepTracker::update(bool grounded, float moved, float downSpeed, Player::Posture posture, float& landing) {
    landing = 0.0f;
    if (!grounded) {
        wasGrounded = false;
        lastDownSpeed = std::max(0.0f, downSpeed);
        return None;
    }
    if (!wasGrounded) { // touchdown
        wasGrounded = true;
        const float s = std::max(downSpeed, lastDownSpeed);
        lastDownSpeed = 0.0f;
        distance = 0.0f;
        if (s > kLandingSpeed) { landing = s; return Landing; }
        return Step; // a hop: one soft step
    }
    const float stride = kStride * (posture == Player::Posture::Prone ? 0.4f : posture == Player::Posture::Crouch ? 0.6f : 1.0f);
    distance += std::max(0.0f, moved);
    if (distance < stride) return None;
    distance = std::min(distance - stride, stride * 0.5f); // (a teleport-sized move is one step, not many)
    return Step;
}

namespace {
    std::vector<StepNoise> noises;
    std::mt19937 rng(77);
    std::uniform_real_distribution<float> jitter(-1.0f, 1.0f);

    struct Walker { StepTracker t; DVec3 lastLocal{}; const Planetoid* body = nullptr; bool known = false; };
    Walker self;
    std::unordered_map<uint32_t, Walker> remote;

    // What is underfoot at `feet` (world): the part stood on, else the ground's material.
    StepSurface SurfaceAt(Scene& world, const DVec3& feet, Vector3 up, bool wading) {
        if (wading) return StepSurface::Water;
        const DVec3 probe = feet - DVec3(Vector3Scale(up, 0.2f));
        Planetoid* body = NearestPlanetoid(world, (probe - world.renderOrigin).toFloat());
        if (!body || !body->field) return StepSurface::Stone;
        const Vector3 local = body->toLocalD(probe).toFloat();
        if (body->plots)
            if (Plot* pl = body->plots->plotFor(local, 1.0f, false)) {
                const int idx = pl->partAt(pl->cellAt(pl->toPlot(local)));
                if (idx >= 0 && idx < (int)pl->parts.size()) return StepSurfaceOfPart(pl->parts[idx]);
            }
        return StepSurfaceOfVoxel(body->materialAt(local));
    }

    void Play(Scene& world, const DVec3& feet, Vector3 up, bool wading, StepTracker::Event ev, float landing,
              Player::Posture posture, bool running, uint32_t who, bool isSelf, bool audible = true) {
        const StepSurface s = SurfaceAt(world, feet, up, wading);
        float loud = StepLoudness(s, posture, running);
        if (ev == StepTracker::Landing) loud = std::max(loud, 0.8f + 0.12f * landing);
        noises.push_back({ feet, loud, GetTime(), who });
        if (!audible) return; // (the dedicated server: the noise for creature hearing, no sound)
        static int ids[(int)StepSurface::Count];
        static bool looked[(int)StepSurface::Count] = {};
        if (!looked[(int)s]) { ids[(int)s] = Audio_Sound(StepSurface_Sound(s)); looked[(int)s] = true; }
        if (ids[(int)s] < 0) return; // the old engine (or no such sound)
        AudioAttenuation att;
        att.refDistance = 1.5f;
        att.maxDistance = 35.0f;
        att.lowpassDistance = 25.0f;
        const float volume = std::min(1.2f, 0.4f * loud * (isSelf ? 0.7f : 1.0f));
        const float pitch = (ev == StepTracker::Landing ? 0.8f : 1.0f) * (1.0f + 0.07f * jitter(rng));
        Audio_PlaySoundAt(ids[(int)s], feet, AudioBus::Effects, AudioPriority::Low, volume, pitch, att);
    }

    void Prune() {
        const double now = GetTime();
        noises.erase(std::remove_if(noises.begin(), noises.end(), [&](const StepNoise& n) { return now - n.time > 2.0; }), noises.end());
    }

    // Walkers known only by their avatar state (remote players on a client; every player on the server).
    void UpdateWalkers(Scene& world, const std::vector<AvatarState>& walkers, float dt) {
        for (const AvatarState& a : walkers) {
            Walker& w = remote[a.session];
            const Vector3 horiz = Vector3Subtract(a.velocity, Vector3Scale(a.up, Vector3DotProduct(a.velocity, a.up)));
            const float speed = Vector3Length(horiz);
            const float down = std::max(0.0f, -Vector3DotProduct(a.velocity, a.up));
            float landing = 0.0f;
            const bool grounded = a.grounded && !a.flying && !a.swimming;
            const StepTracker::Event ev = w.t.update(grounded, speed * dt, down, a.posture, landing);
            if (ev != StepTracker::None) Play(world, a.position, a.up, false, ev, landing, a.posture, speed > 7.5f, a.session, false);
        }
        if (remote.size() > walkers.size() + 8) { // players who left
            std::unordered_map<uint32_t, Walker> kept;
            for (const AvatarState& a : walkers) kept[a.session] = remote[a.session];
            remote.swap(kept);
        }
    }

    // A full Player (the local one, or a session's on the server): the distance walked measured in its body's frame (the
    // planet's spin is not walking; Player::velocity includes it), its touchdown speed from landingSpeed.
    void UpdatePlayer(Scene& world, Walker& w, const Player& p, float dt, uint32_t who, bool isSelf, bool audible) {
        const bool grounded = p.grounded && !p.flyMode && !p.hovering;
        float moved = 0.0f;
        const Planetoid* body = p.nearest;
        if (body) {
            const DVec3 local = body->toLocalD(p.position);
            if (w.known && w.body == body) {
                const double d = (local - w.lastLocal).length();
                moved = d < 5.0 ? (float)d : 0.0f; // (a teleport or respawn is not a step)
            }
            w.lastLocal = local;
            w.body = body;
            w.known = true;
        } else w.known = false;
        float landing = 0.0f;
        const bool swimming = p.inWater && !grounded;
        const StepTracker::Event ev = w.t.update(grounded && !swimming, moved, p.landingSpeed, p.posture, landing);
        if (ev != StepTracker::None) {
            const bool running = moved / dt > p.walkSpeed * 1.25f;
            Play(world, p.position, p.up, p.inWater, ev, landing, p.posture, running, who, isSelf, audible);
        }
    }
    std::unordered_map<uint32_t, Walker> serverWalkers;
}

const std::vector<StepNoise>& Footsteps_RecentNoises() { return noises; }

void Footsteps_UpdateServer(Scene& world, const std::vector<std::pair<uint32_t, const Player*>>& players, float dt) {
    Prune();
    if (dt <= 0.0f) return;
    for (const auto& [session, p] : players)
        if (p) UpdatePlayer(world, serverWalkers[session], *p, dt, session, false, false);
    if (serverWalkers.size() > players.size() + 8) { // players who left
        std::unordered_map<uint32_t, Walker> kept;
        for (const auto& kv : players) kept[kv.first] = serverWalkers[kv.first];
        serverWalkers.swap(kept);
    }
}

void Footsteps_Update(Scene& world, const Player& p, const std::vector<AvatarState>& others, float dt) {
    Prune();
    if (dt <= 0.0f) return;
    UpdatePlayer(world, self, p, dt, 0, true, true);
    // Remote players (multiplayer): their replicated state, already interpolated (#327).
    UpdateWalkers(world, others, dt);
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1240): strides, landings, loudness, surfaces ----
#include "unit_helpers.h"

TEST(Footsteps, OneStepPerStrideAndShorterCrouching) {
    StepTracker t;
    float landing = 0.0f;
    int steps = 0;
    for (int i = 0; i < 1800; ++i) steps += t.update(true, 0.01f, 0.0f, Player::Posture::Stand, landing) == StepTracker::Step; // 18 m
    EXPECT_EQ(steps, 10); // the first after 0.72 m (it starts 60% into a stride), then every 1.8 m
    StepTracker c;
    int crouched = 0;
    for (int i = 0; i < 1800; ++i) crouched += c.update(true, 0.01f, 0.0f, Player::Posture::Crouch, landing) == StepTracker::Step;
    EXPECT_GT(crouched, steps);
    // Standing still: nothing.
    StepTracker s;
    int still = 0;
    for (int i = 0; i < 600; ++i) still += s.update(true, 0.0f, 0.0f, Player::Posture::Stand, landing) != StepTracker::None;
    EXPECT_EQ(still, 0);
}

TEST(Footsteps, LandingsByFallSpeedAndHopsAreSteps) {
    StepTracker t;
    float landing = 0.0f;
    EXPECT_EQ(t.update(false, 0.0f, 6.0f, Player::Posture::Stand, landing), StepTracker::None);    // falling at 6 m/s
    EXPECT_EQ(t.update(true, 0.0f, 0.0f, Player::Posture::Stand, landing), StepTracker::Landing);  // (remote: the last airborne speed)
    EXPECT_FLOAT_EQ(landing, 6.0f);
    EXPECT_EQ(t.update(false, 0.0f, 1.0f, Player::Posture::Stand, landing), StepTracker::None);    // a little hop
    EXPECT_EQ(t.update(true, 0.0f, 1.2f, Player::Posture::Stand, landing), StepTracker::Step);
    EXPECT_EQ(t.update(true, 0.0f, 0.0f, Player::Posture::Stand, landing), StepTracker::None);
    // The local player's touchdown speed comes with the touchdown frame.
    StepTracker s;
    s.update(false, 0.0f, 0.0f, Player::Posture::Stand, landing);
    EXPECT_EQ(s.update(true, 0.0f, 9.0f, Player::Posture::Stand, landing), StepTracker::Landing);
    EXPECT_FLOAT_EQ(landing, 9.0f);
}

TEST(Footsteps, LoudnessAndSurfaces) {
    const float walk = StepLoudness(StepSurface::Stone, Player::Posture::Stand, false);
    EXPECT_GT(StepLoudness(StepSurface::Stone, Player::Posture::Stand, true), walk);
    EXPECT_LT(StepLoudness(StepSurface::Stone, Player::Posture::Crouch, false), walk);
    EXPECT_LT(StepLoudness(StepSurface::Stone, Player::Posture::Prone, false), StepLoudness(StepSurface::Stone, Player::Posture::Crouch, false));
    EXPECT_LT(StepLoudness(StepSurface::Grass, Player::Posture::Stand, false), StepLoudness(StepSurface::Metal, Player::Posture::Stand, false));
    EXPECT_EQ(StepSurfaceOfVoxel(VoxelMaterial::Grass), StepSurface::Grass);
    EXPECT_EQ(StepSurfaceOfVoxel(VoxelMaterial::Granite), StepSurface::Stone);
    EXPECT_EQ(StepSurfaceOfVoxel(VoxelMaterial::Pebbles), StepSurface::Gravel);
    EXPECT_EQ(StepSurfaceOfVoxel(VoxelMaterial::Mud), StepSurface::Dirt);
    Part plank{ PartType::Cube, { 0, 0, 0 }, 0, (uint8_t)BuildMaterialForItem("oak_planks") };
    Part brick{ PartType::Cube, { 0, 0, 0 }, 0, (uint8_t)BuildMaterialForItem("stone_brick") };
    Part plate{ PartType::Slab, { 0, 0, 0 }, 0, (uint8_t)BuildMaterialForItem("iron_plate") };
    Part machine{ PartType::Furnace, { 0, 0, 0 }, 0, 0 };
    EXPECT_EQ(StepSurfaceOfPart(plank), StepSurface::Wood);
    EXPECT_EQ(StepSurfaceOfPart(brick), StepSurface::Stone);
    EXPECT_EQ(StepSurfaceOfPart(plate), StepSurface::Metal);
    EXPECT_EQ(StepSurfaceOfPart(machine), StepSurface::Metal);
    for (int s = 0; s < (int)StepSurface::Count; ++s) EXPECT_EQ(std::string(StepSurface_Sound((StepSurface)s)).rfind("step_", 0), 0u);
}
#endif
