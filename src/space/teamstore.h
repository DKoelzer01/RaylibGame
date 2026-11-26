#pragma once

// The team progression store (#1886): one place per team for the progression state a late joiner needs - research
// today; the team's explored map and blueprint library next. Each module keeps its own data and its own text form; the
// store names them as sections, joins them for the net (JoinInfo, a resync) and writes them into the world save as
// `team_<section> <id>` ... `team_<section>_end` beside the team record (so research keeps its #1030 save lines).

#include <cstdio>
#include <string>

// Every section of a team, as text:  "section <name>\n" <the module's lines> "section_end\n"  (empty sections left out).
std::string TeamStore_Text(int team);
// Applies a TeamStore_Text: every known section takes its block, and a known section the text leaves out is cleared, so
// a client mirrors the host's state exactly. Unknown sections are skipped (a newer host). False: no "section" line.
bool TeamStore_SetText(int team, const std::string& text);

// World save: writes the team's non-empty sections (called from Teams_Save for each team).
void TeamStore_SaveTeam(FILE* f, int team);
// World load: handles a `team_<section> <id>` line by reading its block up to `team_<section>_end`. False: not ours.
bool TeamStore_LoadLine(const char* tok, FILE* f);
