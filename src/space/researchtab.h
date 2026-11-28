#pragma once

// The Research tab (#1192, #1030; the factory doc's UI, decided): the research nodes per tier as a pannable, zoomable graph
// (core nodes on each tier's spine, branches beside it; lines to their prerequisites), what the team may see (the current
// tier in full, the next tier names only, later tiers hidden), the team's queue and the milestone checklist. Inert in play
// for now (user decision 2026-11-21): no lab researches yet and nothing is gated; a click only queues.

#include "research.h"

#include <raylib.h>
#include <functional>
#include <string>

// (#1886) Set by the game while joined to a server: a queue change made in the tab goes to the host (MsgResearchQueue),
// which keeps the team's research; the tab applies it to the mirrored state at once and the host's next store replaces it.
extern std::function<void(ResearchQueueOp op, const std::string& node)> gResearchQueueRequest;

void ResearchTab_Draw();            // inside the game window's screen (ImGui)
// Tests: the screen-space centre of a node's box as last drawn (false: not drawn), and the team whose state the tab shows.
bool ResearchTab_TestNodeCenter(const std::string& id, Vector2& center);
int ResearchTab_Team();
