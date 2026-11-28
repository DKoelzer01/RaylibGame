#pragma once

// The team progression store (#1886): one place per team for the progression state a late joiner needs - research,
// the team's explored map; the blueprint library next. Each module keeps its own data and its own text form; the store
// names them as sections, joins them for the net (JoinInfo, a resync) and writes them into the world save as
// `team_<section> <id>` ... `team_<section>_end` beside the team record (so research keeps its #1030 save lines). The
// map's tiles stay in the team's quest block in the save (team_tiles, #338), so the store doesn't write them again.

#include <cstdio>
#include <string>

// Every section of a team, as text:  "section <name>\n" <the module's lines> "section_end\n"  (empty sections left out).
std::string TeamStore_Text(int team);
// A client applies a TeamStore_Text: every known section takes its block, and a known section the text leaves out gets
// "" (research is cleared: the client mirrors the host). The map is added to the local map, never cleared: a player keeps
// what they explored (the user, on #1886). Unknown sections are skipped (a newer host). False: no "section" line.
bool TeamStore_SetText(int team, const std::string& text);
// A client: the host's TeamStore message for our team. Applies it (an empty text, a team with nothing yet, clears
// research and the library) and remembers the team: the host's id of our team, which the research tab reads (a client's
// own team registry doesn't know the host's ids). TeamStore_LocalTeam: 0 = none (offline, a host, or not joined yet).
void TeamStore_ApplyFromHost(int team, const std::string& text);
int TeamStore_LocalTeam();
void TeamStore_ResetLocalTeam(); // (left the server)

// World save: writes the team's non-empty sections (called from Teams_Save for each team).
void TeamStore_SaveTeam(FILE* f, int team);
// World load: handles a `team_<section> <id>` line by reading its block up to `team_<section>_end`. False: not ours.
bool TeamStore_LoadLine(const char* tok, FILE* f);
