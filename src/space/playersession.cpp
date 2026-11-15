#include "playersession.h"

namespace {
PlayerSession gLocal; // single-player: the one session
}

PlayerSession& LocalSession() { return gLocal; }

// The old globals, now the local session's members (#313).
Player& player = gLocal.player;
Inventory& inventory = gLocal.inventory;
SurvivalState& survival = gLocal.survival;
TerrainToolSettings& terrainTool = gLocal.terrainTool;
BuildToolState& buildTool = gLocal.buildTool;
