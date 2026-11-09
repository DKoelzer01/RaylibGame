#ifndef GAMEUI_H
#define GAMEUI_H

class Scene;

// In-game panels (Dear ImGui for now): Tab = inventory + crafting.
bool GameUI_InventoryOpen();
void GameUI_ToggleInventory();
void GameUI_Draw(Scene& world);  // inside rlImGuiBegin()/rlImGuiEnd()
void GameUI_DrawHUD(Scene& world); // plain raylib 2D HUD (status, toasts)
void GameUI_Toast(const char* text);

#endif
