#ifndef PLAYERSESSION_H
#define PLAYERSESSION_H

// Per-player state (#313, multiplayer phase 1): everything that belongs to one player rather than to the world, in one
// struct. Single-player has one session, the local one; a server will hold one per connected player.
// The old global names (player, inventory, survival, buildTool, terrainTool) are references to the local session's
// members, so existing code keeps working; the hotbar, suit, tool wear and recipe book read their state from it
// (LocalSession()). Shared, not per-player (decided on #309): the quest book and tiers, the world, machines.
// Not moved yet (see #313): navigation (explored map, waypoints), the telescope, rocket crew state, weather "here",
// the terrain tool's aim/swing statics (terraintool.cpp) and the death-hold statics in main.cpp.

#include "player.h"
#include "inventory.h"
#include "playerinput.h"
#include "terraintool.h"
#include "survival.h"
#include "buildtool.h"
#include "hotbar.h"
#include "suit.h"
#include "toolwear.h"
#include "recipebook.h"
#include "skills.h"
#include "netid.h"

#include <cstdint>

struct PlayerSession {
    NetId netId = NetId_Make(NetKind::Player, 0); // 0 = the local player; a server issues ids (or SteamIDs) later (#317)
    Player player;
    Inventory inventory;
    SurvivalState survival;
    TerrainToolSettings terrainTool; // the tool in hand, creative/terraform modes, reach
    BuildToolState buildTool;
    HotbarState hotbar;
    SuitState suit;
    ToolWearState toolWear;
    RecipeBookState recipes;
    SkillState skills;               // (#883) XP and perks (skills.h)
    PlayerInput input;               // the last step's input (filled by the server from the network later)
};

PlayerSession& LocalSession();

#endif
