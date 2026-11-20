#pragma once

// Skills in play (#883): the local player's SkillState (PlayerSession::skills) fed by the game's events and read for its
// multipliers. skills.h is the model (no game state); this connects it: the terrain tool's hooks (gSkillEvent,
// gSkillStat), the world's XP rate, the anti-farming cell share for strokes, level-up toasts and pick wear.
// A joined client's strokes are the server's (no local XP yet): server-side XP is slice 2.

#include <string>

void SkillPlay_Init();                                              // once at start-up: sets the terrain tool's hooks
// Tool wear for `strokes` strokes with `toolId` in hand, after the perks (pick_wear for picks; #919: shovel_wear for shovels and hoes,
// axe_wear for axes): the strokes that wear it.
int SkillPlay_WearStrokes(const std::string& toolId, int strokes);
int SkillPlay_LevelUps();                                           // (tests) level-ups the local player has had since start-up
// (#919) Once a sim step on foot: Athletics' walk_10m (every 10 m walked on the ground) and Survival's night (a night
// lived through: the sun was down a minute or more and came up again with the player alive).
class Scene;
void SkillPlay_Step(Scene& world, float dt);
