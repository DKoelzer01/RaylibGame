#ifndef HELPMANUAL_H
#define HELPMANUAL_H

// In-game database (#354): the inventory's Help tab, like Factorio's Factoriopedia. Hand-written guides from
// resources/data/help.txt ("# Title" lines, then text; "##" comments) and a page per item built from the data files
// (hand and machine recipes, uses, fuels), with links between pages, Back, search, and "Discovered only".

#include <cstdio>
#include <string>

void HelpManual_DrawTab();               // inside the inventory's ImGui tab bar
int HelpManual_TopicCount();             // tests: guides loaded, plus the item database
void HelpManual_OpenItem(const std::string& id); // opens an item's page (the next time the tab draws)
const std::string& HelpManual_PageItem(); // the item page shown ("" a guide) (tests)
void HelpManual_Save(FILE* f);           // "Discovered only", with the world
bool HelpManual_LoadLine(const char* tok, FILE* f);

#endif
