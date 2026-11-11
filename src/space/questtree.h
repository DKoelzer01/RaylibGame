#ifndef QUESTTREE_H
#define QUESTTREE_H

// The quest tree (deferred plan-2 work): side lines (farming, husbandry, the suit, exploring) that
// branch off the main line (tiers.txt), from resources/data/quests.txt. A node opens once its
// requirements (other nodes, or a main-line tier) are done; it completes when its task is met (the
// same task kinds as the main line) and pays its reward into the inventory. The Quests tab draws the
// main line and the side lines as one tree.

#include <cstdio>
#include <string>

void QuestTree_Update();                 // after the main line (QuestBook_Update calls it)
void QuestTree_DrawTab();                // ImGui: the tree
void QuestTree_Reset();
void QuestTree_Save(FILE* f);
bool QuestTree_LoadLine(const char* tok, FILE* f);

// Tests.
int QuestTree_NodeCount();
bool QuestTree_IsDone(const std::string& id);
bool QuestTree_IsOpen(const std::string& id);   // requirements done, task not yet

#endif
