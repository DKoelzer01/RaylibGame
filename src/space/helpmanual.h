#ifndef HELPMANUAL_H
#define HELPMANUAL_H

// In-game help manual (#354): the inventory's Help tab. Hand-written topics from resources/data/help.txt ("# Title"
// lines, then text; "##" comments) and an Items topic generated from items.txt and recipes.txt, so it stays current.
// A search box filters topics and items.

void HelpManual_DrawTab();               // inside the inventory's ImGui tab bar
int HelpManual_TopicCount();             // tests: topics loaded (the Items topic included)

#endif
