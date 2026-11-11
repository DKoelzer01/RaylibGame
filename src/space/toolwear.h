#ifndef TOOLWEAR_H
#define TOOLWEAR_H

// Tool durability (plan 2, milestone 14): every stroke, chop or hit with a tool wears it (tools.h
// ToolDurability: stone 150 strokes ... titanium 2000). A worn-out tool breaks (one is removed).
// Repairing restores it for a little of its head material (Crafting tab). Wear is kept per tool
// item id: tools never stack, and a second copy only starts wearing once the first has broken.

#include <cstdio>
#include <string>

void ToolWear_OnStroke(const std::string& toolId, int strokes = 1); // true use only (not creative)
float ToolWear_Left(const std::string& toolId);   // 1 new .. 0 about to break (1 for non-tools)
int ToolWear_Used(const std::string& toolId);
bool ToolWear_RepairCost(const std::string& toolId, std::string& item, int& count);
bool ToolWear_Repair(const std::string& toolId);  // pays the cost; false if it can't
void ToolWear_DrawRepairUI();                     // ImGui, in the Crafting tab
void ToolWear_Reset();
void ToolWear_Save(FILE* f);
bool ToolWear_LoadLine(const char* tok, FILE* f);

#endif
