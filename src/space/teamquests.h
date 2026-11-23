#pragma once

// Team quest books on a server (#338; the user's rules on #303: quests and tiers belong to teams, rewards go to the
// completing player). Each team's book is a QuestBookState, loaded on first use from the team's kept block in the world
// file (#346: TeamRecord::quests, its save lines) and written back before a save. Events the server decides are counted
// into the player's team's book; evaluation runs here, its rewards and news through the given callbacks (the server pays
// the completing player with Host_Give, so the client gets inventory deltas: #816). A client mirrors its team's book.

#include <functional>
#include <string>
#include <raylib.h>

struct QuestBookState;

QuestBookState& TeamQuests_Book(int team);
void TeamQuests_Count(int team, const std::string& event, double n = 1.0);
void TeamQuests_AtLeast(int team, const std::string& event, double value);
void TeamQuests_NotePlaced(int team, const std::string& item);
void TeamQuests_Obtained(int team, const std::string& item, double n);
void TeamQuests_Made(int team, const std::string& item, double n); // (#1728) a team's machine finished n of item (counts for "get" as the larger of made and obtained)
// (#338 slice 3) a player of `team` on `body` (body-local unit `up`, radius R) reveals map tiles as the local map does
// (Nav_RevealTiles); a team's tiles are one set per body (teammates add up, a tile counts once), saved with its book and
// never sent to clients; "explored" rises to their total.
int TeamQuests_Explore(int team, const std::string& body, Vector3 up, float radius); // the team's total when it grew, else 0
double TeamQuests_Explored(int team);
void TeamQuests_Visit(int team, bool moon);           // a moon visited (else the ring)
int TeamQuests_BodiesReached(int team);              // (#343) planets and moons the team has explored ground on (claim budget)
// Teams forming and splitting (#604): a team's book and explored tiles copied to another (a join, a split), or dropped (an emptied
// team; its record goes too). Use these rather than copying TeamQuests_Book: the tiles live beside the book.
void TeamQuests_Copy(int from, int to);
// (#861) a joiner's own progress added to the team's (the user, on #604): tier the higher, tasks / nodes / placed items /
// visits the union, counts and obtained totals the larger of the two (not the sum), explored tiles the union.
void TeamQuests_Merge(int from, int into);
void TeamQuests_Forget(int team);
// Evaluates the book (tiers, side lines); reward(item, count) for each reward, news(text) for each tier-up / node done.
void TeamQuests_Update(int team, const std::function<void(const std::string&, int)>& reward, const std::function<void(const std::string&)>& news);
std::string TeamQuests_Text(int team);   // the book as its save lines (what a client gets)
void TeamQuests_StoreAll();              // every loaded book back into its team's record (before the world is saved)
void TeamQuests_Reset();                 // a new or reloaded world: none loaded
// A client: the server's book for its team replaces the local one (the Quests tab and the HUD show it).
void TeamQuests_ApplyText(const std::string& text);
bool TeamQuests_SelfTest();             // --test-team-quests
