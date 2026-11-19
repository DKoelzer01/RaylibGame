#ifndef GAMEUI_H
#define GAMEUI_H

#include <cstdint>

class Scene;

// The game window (#821, Dear ImGui for now): one window filling the screen, its screens switched by large tabs across
// the top (Inventory, Crafting, Quests, Skills, Team). Tab or E opens it at the screen last used.
enum class GameScreen { Inventory, Crafting, Quests, Skills, Team, Catalog, Count }; // Catalog: only in creative (#823)
bool GameUI_InventoryOpen();              // the game window is open (any screen)
void GameUI_ToggleInventory();
void GameUI_ShowScreen(GameScreen s);     // open the window on a screen
GameScreen GameUI_Screen();               // the screen shown (or shown last)
void GameUI_ShowTab(const char* tab);     // older names: "Crafting", "Quest tree" (Quests, tree view), "Suit" (Inventory), "Help"
bool GameUI_TextInputActive();            // a text field in the window has the keyboard (the window's keys don't apply)
float GameUI_Scale();                     // the interface size (Settings; auto by the screen's height)
void GameUI_PrepareFrame(Scene& world);  // (#822) before ImGui, after the world: renders the inventory screen's player preview
bool GameUI_PreviewDrawn();              // (#822) the preview was rendered this frame (tests)
void GameUI_QueueCraft(int recipe, int count);   // (#824) into the crafting queue (one every 0.3 s, through Crafting_Request)
int GameUI_CraftQueueLength();                    // crafts still queued
void GameUI_SelectRecipe(int recipe);             // the Crafting screen's selected recipe
#include <string>
#include <functional>
// (#866) Set while joined (main.cpp): creative items from the server (the Catalog, a creative middle-click; "" clears the inventory).
inline std::function<void(const std::string& item, int count)> gCreativeGive;
void GameUI_OpenDatabase(const std::string& itemId); // (#826) the database pop-up, at an item's page ("" where it was)
bool GameUI_DatabaseOpen();
void GameUI_Draw(Scene& world);  // inside rlImGuiBegin()/rlImGuiEnd()
void GameUI_DrawHUD(Scene& world); // plain raylib 2D HUD (status, toasts)
void GameUI_Toast(const char* text);
void GameUI_InitPickupFeed(); // #208: "+N item" lines in the bottom right for everything that enters the inventory
void GameUI_InstallKeyHints(); // (#827) the key hints line (the terrain tool's hook): bindings, size, fade
bool GameUI_TestCursorStack();  // --test-inv-cursor (#192): pick up, put down, halves, ones, merge, swap, quick-move, return
bool GameUI_TestContainerCursor(Scene& world, uint64_t body, uint32_t plot, int part, int partType); // --test-container-cursor (#192)
// A machine's or chest's screen (RMB on it): its slots in the inventory window's right column.
// Clicking an inventory cell puts that stack in; clicking one of its slots takes it out.
void GameUI_OpenContainer(uint64_t body, uint32_t plot, int part, int partType);
bool GameUI_ContainerOpen();
void GameUI_CloseContainer();

#endif
