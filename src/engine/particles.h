#ifndef PARTICLES_H
#define PARTICLES_H

// Particles: camera-facing soft sprites for rocket exhaust, smoke, dust and debris from digging,
// wood chips, embers, rain splashes and fireflies. Each particle lives in the local frame of a body
// (a planet or moon), so it stays put on the spinning ground and across floating-origin shifts;
// gravity pulls toward that body's centre (scaled per kind; smoke rises). One pooled array (oldest
// replaced when full), updated on the main thread; drawn in two batches (alpha-blended, sorted far to
// near, lit by the day/night light; then additive glow for fire, sparks and fireflies), depth-tested
// without depth writes.

#include <raylib.h>
#include "dvec3.h"

class Scene;
class Planetoid;

enum class ParticleKind : unsigned char {
    Smoke,      // grows, fades, rises and drifts
    Fire,       // additive, shrinks, yellow -> red
    Dust,       // a puff of the ground's colour, settles slowly
    Debris,     // small chunks thrown out, fall and bounce once
    Spark,      // additive, fast, falls
    Splash,     // water droplets: fall, gone on landing
    Firefly,    // additive, slow wander, blinks
    Count
};

struct ParticleBurst {
    ParticleKind kind = ParticleKind::Dust;
    Color color = WHITE;
    int count = 8;
    float size = 0.15f;            // m (start; smoke grows, fire shrinks)
    float life = 1.0f;             // s (each particle +-30%)
    float speed = 1.0f;            // m/s, along `dir` plus a random cone of `spread`
    float spread = 0.6f;           // 0 = straight along dir, 1 = any direction
    Vector3 dir = { 0, 1, 0 };     // body-local; zero = away from the body centre ("up")
    Vector3 baseVel = { 0, 0, 0 }; // body-local velocity added to every particle (a moving emitter)
};

void Particles_Emit(const Planetoid& body, Vector3 local, const ParticleBurst& b);
// The same from a world-space point (absolute, double) and a world-space direction / velocity.
void Particles_EmitWorld(const Planetoid& body, const DVec3& world, const ParticleBurst& b, Vector3 worldDir, Vector3 worldBaseVel);
void Particles_Update(Scene& world, float dt);
void Particles_Draw(Scene& world);  // inside BeginMode3D (render space)
void Particles_Clear();
int Particles_Count();
void Particles_Unload();

#endif
