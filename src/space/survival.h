#ifndef SURVIVAL_H
#define SURVIVAL_H

class Scene;

// Suit survival loop: oxygen drains outside the ship's life-support field, power drains with
// recharging the jetpack and the headlamp, health drops when suffocating or from hard landings. Death
// respawns you at the ship.
struct SurvivalState {
    float oxygen = 100.0f, power = 100.0f, health = 100.0f;   // 0..100
    float oxygenDrainPerSec = 0.35f;                          // ~5 minutes per tank
    float powerDrainPerSec = 0.05f;                           // life support baseline
    float jetpackPowerPerTank = 25.0f;                        // suit power to recharge an empty jetpack tank
    float lastJetpackFuel = 1.0f;                             // fuel at the end of the last update (recharge tracking)
    float suffocationDamagePerSec = 6.0f;
    float lifeSupportRefillPerSec = 12.0f;                    // inside / near the ship
    float headlampPowerPerSec = 0.15f;                        // ~11 min of lamp on a full suit
    bool headlamp = false;                                    // L toggles; goes out at 0 power
    bool inLifeSupport = false;
    float breathableAltitude = 15000.0f;                      // a planet's air is breathable below this (m above sea level)
    float airRefillPerSec = 4.0f;                             // suit oxygen top-up from breathable air
    bool breathable = false;
    bool dead = false;
    double deathTime = 0.0;
    int deaths = 0;
};
extern SurvivalState survival;

void Survival_Update(Scene& world, float dt, bool inputEnabled);
void Survival_OnLanding(float impactSpeed);   // fall damage
void Survival_Respawn(Scene& world);
void Survival_DrawHUD();

// Lightweight objective tracker for the vertical slice.
void Objectives_Update(Scene& world);
void Objectives_DrawHUD();
void Objectives_MarkVisitedBody(int bodyIndex);

#endif
