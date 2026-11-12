#ifndef GAMEUI_H
#define GAMEUI_H

#include <cstdint>

class Scene;

// In-game panels (Dear ImGui for now): Tab = inventory + crafting.
bool GameUI_InventoryOpen();
void GameUI_ToggleInventory();
void GameUI_ShowTab(const char* tab);     // open the inventory window on a tab ("Quest tree", "Suit")
void GameUI_Draw(Scene& world);  // inside rlImGuiBegin()/rlImGuiEnd()
void GameUI_DrawHUD(Scene& world); // plain raylib 2D HUD (status, toasts)
void GameUI_Toast(const char* text);
// A machine's or chest's screen (RMB on it): its slots in the inventory window's right column.
// Clicking an inventory cell puts that stack in; clicking one of its slots takes it out.
void GameUI_OpenContainer(uint64_t body, uint32_t plot, int part, int partType);
bool GameUI_ContainerOpen();
void GameUI_CloseContainer();

#endif
