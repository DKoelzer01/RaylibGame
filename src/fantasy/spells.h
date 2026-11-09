#ifndef SPELLS_H
#define SPELLS_H

#include <string>
#include <vector>
#include "sealrecognizer.h"

class Scene;

// Drawn-seal magic. Press G to open the grimoire canvas, draw a ring with a core sign in the
// middle and modifier signs toward the ring (dots add power), then press G / Enter to cast.
// The seal's drawing quality scales the spell's power; a sloppy ring makes it fizzle.
//
// Core signs:  flame (triangle)  earth (square)  water (wave)  wind (spiral)  light (asterisk)
// Modifiers:   arrow (chevron) = launch / forward   orb (small circle) = burst / sphere
//              line = wall / span                   dots = +25% power each
struct SpellResult {
    bool success = false;
    std::string name;       // e.g. "Fireball"
    std::string message;    // what happened / why it failed
    float power = 0.0f;     // 0..~2
    float manaCost = 0.0f;
    SealAnalysis analysis;
};

struct SpellState {
    float mana = 100.0f, maxMana = 100.0f, manaRegenPerSec = 6.0f;
    float health = 100.0f;
    bool grimoireOpen = false;
    std::vector<Stroke> strokes;   // current drawing (screen space)
    SpellResult last;
    double lastTime = -10.0;
    int casts = 0, craters = 0, fizzles = 0;
};
extern SpellState spells;

void Spells_Update(Scene& world, float dt, bool inputEnabled);  // grimoire input + active effects
void Spells_Draw3D(const Scene& world);                         // projectiles, orbs (inside BeginMode3D)
void Spells_DrawOverlay(Scene& world);                          // grimoire canvas + HUD (2D)
SpellResult Spells_Cast(Scene& world, const std::vector<Stroke>& strokes); // analyse + cast
int Spells_ActiveProjectiles();
void Spells_ShowInGrimoire(const std::vector<Stroke>& strokes); // tests/demo: open the canvas with a drawing

#endif
