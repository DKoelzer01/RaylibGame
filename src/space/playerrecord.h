#ifndef PLAYERRECORD_H
#define PLAYERRECORD_H

// Multiplayer save format (#346, save format 13): the world file holds the world; each player who ever played in it has a
// record of their own, <save>.players/<uuid>.player (position, inventory, vitals, suit, hotbar, map pins, help). Players are
// known by a UUID (this install's: settings.txt player_uuid, made once); their display name is a separate field (Steam's
// later). Quests and tiers belong to teams (a solo player is a team of one): team lines in the world file, each team's
// quest book in a block of its own. Files are written whole or not at all (a .tmp renamed over the old one).

#include <cstdint>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

std::string PlayerId_New();                         // a random UUID (version 4)
std::string PlayerId_Local();                       // this install's player UUID (made and kept in settings on first use)
bool PlayerId_Valid(const std::string& id);         // the 8-4-4-4-12 hex form
std::string PlayerRecord_Dir();                     // <save>.players (beside <save>.regions)
std::string PlayerRecord_Path(const std::string& uuid);

// Writes `path` through path.tmp, renamed over it when complete (a crash mid-write leaves the old file whole).
bool Save_WriteAtomic(const std::string& path, const std::function<bool(FILE*)>& write);

// Teams (world file lines): the owner (who made the world), each team (id, suit-style colour, name), its members (UUIDs)
// and its quest book. This game's own team's quest book loads into the quest book; other teams' blocks are kept as they
// were and written back unchanged (a server keeps everyone's; quests for joined players run on the server later, #337).
struct TeamRecord {
    int id = 0;
    uint32_t colour = 0;               // 0xRRGGBB
    std::string name;
    std::vector<std::string> members;  // UUIDs
    std::string quests;                // a quest block not loaded into this game's quest book, as its lines
};
void Teams_Reset();
// The world file's lines: owner, team, team_member, team_quests ... team_quests_end. `localId` is the UUID whose team's
// quests go into the quest book ("" none: a server keeps all as text).
bool Teams_LoadLine(const char* tok, FILE* f, const std::string& localId);
// Writes the teams; `writeLocalQuests` writes the local player's team's block from the quest book (else kept text).
void Teams_Save(FILE* f, const std::string& localId, const std::function<void(FILE*)>& writeLocalQuests);
int Teams_Of(const std::string& uuid);              // 0: none
int Teams_Join(const std::string& uuid, const std::string& name, uint32_t colour); // a team of one if it has none; its id
const std::string& Teams_Owner();
void Teams_SetOwner(const std::string& uuid);
void Teams_ReplaceMember(const std::string& from, const std::string& to);
const std::vector<TeamRecord>& Teams_All();

#endif
