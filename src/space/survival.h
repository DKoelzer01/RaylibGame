#ifndef SURVIVAL_H
#define SURVIVAL_H

class Scene;

// Suit survival loop: oxygen drains outside the ship's life-support field, power drains with
// jetpack use, health drops when suffocating or from hard landings. Death respawns you at the ship.
struct SurvivalState {
    float oxygen = 100.0f, power = 100.0f, health = 100.0f;   // 0..100
    float oxygenDrainPerSec = 0.35f;                          // ~5 minutes per tank
    float powerDrainPerSec = 0.05f;                           // life support baseline
    float jetpackPowerPerSec = 4.0f;
    float suffocationDamagePerSec = 6.0f;
    float lifeSupportRefillPerSec = 12.0f;                    // inside / near the ship
    bool inLifeSupport = false;
    bool dead = false;
    double deathTime = 0.0;
    int deaths = 0;
};
extern SurvivalState survival;

void Survival_Update(Scene& world, float dt, bool jetpackFiring, bool inputEnabled);
void Survival_OnLanding(float impactSpeed);   // fall damage
void Survival_Respawn(Scene& world);
void Survival_DrawHUD();

// Lightweight objective tracker for the vertical slice.
void Objectives_Update(Scene& world);
void Objectives_DrawHUD();
void Objectives_MarkVisitedBody(int bodyIndex);

#endif
