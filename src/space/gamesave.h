#ifndef GAMESAVE_H
#define GAMESAVE_H

#include <functional>
#include "dvec3.h"
#include "frames.h"

// Registers save/load hooks for gameplay state with the world save file: the world's own state in the save, the player's in
// their record beside it (#346, playerrecord.h).
void GameSave_Register();
// Whether the last load found this game's player (their record, or a format 12 save's player lines).
bool GameSave_PlayerLoaded();
double GameSave_PlayTime();   // (#827) seconds this player has played in this world (their record's, plus this session)
// (#901) The local player's place from the world just loaded, in its body's frame (the body by NetId), once: the loading
// screen holds the player there exactly (a frame later the planet has turned ~1.5 m).
// (Plan 11 F7a, #2074) framePos: the saved frame state too (frame.body 0: the record had none), for the hold to restore.
bool GameSave_TakeLoadedFrame(unsigned long long& bodyNetId, DVec3& local, FramePos* framePos = nullptr);
// Called by a load that found no player of ours (a world from a dedicated server, say): put a new player at the world's
// start with the starting kit. The world's own state (quests, difficulty) is the save's.
void GameSave_SetNewPlayerHook(std::function<void()> hook);
class Scene;
bool GameSave_TestSplit(Scene& world); // --test-save-split (#346)

// A player's record on a server (#346 step 2): read into and written from a remote player's state (position in their
// body's frame, inventory, survival); the lines a server doesn't hold (hotbar, map pins, suit, ...) ride along as text.
#include <string>
#include <vector>
struct Player; class Inventory; struct SurvivalState;
struct PlayerRecordMeta { std::string uuid, name; double firstSeen = 0.0, playTime = 0.0; std::string kept; };
struct SkillState;
bool PlayerRecord_Read(const std::string& path, Player& p, Inventory& inv, SurvivalState& sv, PlayerRecordMeta& meta, SkillState* skills = nullptr);
bool PlayerRecord_Write(const std::string& path, Player& p, Inventory& inv, SurvivalState& sv, PlayerRecordMeta& meta, double sessionSeconds, SkillState* skills = nullptr); // skills (#920): a remote player's
std::vector<PlayerRecordMeta> PlayerRecord_List(); // the uuid and name of every record beside the save (a server's name claims)

#endif
