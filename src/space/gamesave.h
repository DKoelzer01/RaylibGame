#ifndef GAMESAVE_H
#define GAMESAVE_H

#include <functional>

// Registers save/load hooks for gameplay state with the world save file: the world's own state in the save, the player's in
// their record beside it (#346, playerrecord.h).
void GameSave_Register();
// Whether the last load found this game's player (their record, or a format 12 save's player lines).
bool GameSave_PlayerLoaded();
// Called by a load that found no player of ours (a world from a dedicated server, say): put a new player at the world's
// start with the starting kit. The world's own state (quests, difficulty) is the save's.
void GameSave_SetNewPlayerHook(std::function<void()> hook);
class Scene;
bool GameSave_TestSplit(Scene& world); // --test-save-split (#346)

#endif
