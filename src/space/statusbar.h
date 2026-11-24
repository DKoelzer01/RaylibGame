#pragma once
// Status effects on the HUD (#1777): a stack of entries near the middle of the screen's left edge, each a label (and an
// optional item icon) over a bar that starts full and drains over the effect's duration. The first user is the parry window
// (#1776); buffs, debuffs, wetness and poison come later. Timed on the sim clock (SimTime): bars freeze under the pause menu.
// An entry is removed when it runs out. The label is always drawn (no colour-only cue); the widget follows the HUD's size.
#include <raylib.h>

// Starts the effect `id` full, draining over `seconds`. Re-adding an id restarts it in place (it keeps its slot in the stack).
// icon: an item id drawn from the item icon atlas, or nullptr.
void StatusBar_Add(const char* id, float seconds, const char* label, Color color, const char* icon = nullptr);
void StatusBar_Remove(const char* id);      // an early end (the parry landed)
bool StatusBar_Active(const char* id);
float StatusBar_Remaining(const char* id);  // 1 full -> 0 empty; 0 when not active
int StatusBar_Count();                      // entries showing
void StatusBar_Clear();                     // all gone (a new world: its sim clock starts again)
void StatusBar_Draw();                      // from GameUI_DrawHUD
