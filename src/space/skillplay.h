#pragma once

// Skills in play (#883): the local player's SkillState (PlayerSession::skills) fed by the game's events and read for its
// multipliers. skills.h is the model (no game state); this connects it: the terrain tool's hooks (gSkillEvent,
// gSkillStat), the world's XP rate, the anti-farming cell share for strokes, level-up toasts and pick wear.
// A joined client's strokes are the server's (no local XP yet): server-side XP is slice 2.

#include <string>

void SkillPlay_Init();                                              // once at start-up: sets the terrain tool's hooks
// Tool wear for `strokes` strokes with `toolId` in hand, after the perks (pick_wear for picks): the strokes that wear it.
int SkillPlay_WearStrokes(const std::string& toolId, int strokes);
int SkillPlay_LevelUps();                                           // (tests) level-ups the local player has had since start-up
