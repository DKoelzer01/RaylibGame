#pragma once

// Discovery (#1250, magic step 3; the magic doc section 1): magic starts the first time the player stumbles on it, by any
// route, and that first find gives the First Sight page (reading it teaches magic sight, #1249) and opens the Arcane quest
// line (quests.txt, hidden until then), which guides loosely: any route satisfies its steps.
// Routes in this step:
//   a location - a ley node: faint motes rise over it (no sight needed); standing on it wakes it and a page rises;
//   a formation - a resonant crystal: rock dug on a strong ley line now and then holds one;
//   loot - a First Sight page in a ruin's, camp's or dungeon's chest, only in chests near a ley line (so the first cast of
//          the sight shows a line right there; #1249's first-range question);
//   a creature drop (a Glimmer Wisp's mote) waits for the magic creatures.
// Everything behind RAYLIBGAME_MAGIC: the magic items (magic_items.txt) are only loaded with it, so loot and drops can't
// hold them otherwise.

#include <string>

class Scene;
class Planetoid;
class ManaField;
struct MagicPlayerState;
#include <raylib.h>

enum class MagicRoute { Node, Crystal, Loot, Creature };

// The rules (pure). A find: the first one marks the player as having found magic (true) and gives the page unless they
// already know the sight or carry a page.
struct MagicFind { bool first = false, givePage = false; };
MagicFind MagicDiscovery_OnFind(MagicPlayerState& s, bool carriesPage);
// Near a ley line: some base kind's line within `range` m of the direction, at a strength of at least 0.3.
bool MagicDiscovery_NearLine(const ManaField& f, Vector3 dir, float range);
// Rock dug at a place with this much line flow (0..2) holds a resonant crystal one stroke in N (0: never).
int MagicDiscovery_CrystalOdds(float flow);
// The distinct node kinds studied (the Arcane line's last step).
int MagicDiscovery_NodeKindsStudied(const MagicPlayerState& s);

// The game.
void MagicDiscovery_Install();                    // the hooks (the crystal in rock); once, with RAYLIBGAME_MAGIC
void MagicDiscovery_Update(Scene& world, float dt); // node motes and touch, finds by item, quest events
bool MagicDiscovery_IsPage(const std::string& id);
bool MagicDiscovery_ReadPage();                    // the local player reads a First Sight page: true if it taught the sight
// May a chest at this body-local point hold this item? Magic items only near a ley line; everything else yes.
bool MagicDiscovery_LootAllowed(const Planetoid& body, Vector3 local, const std::string& item);
