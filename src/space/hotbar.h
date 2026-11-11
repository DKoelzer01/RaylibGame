#ifndef HOTBAR_H
#define HOTBAR_H

// The hotbar: keys 1-9 and 0 select one of ten slots, each bound to a tool or an item.
//   tool:terrain     the terrain tool (LMB dig, RMB fill) - build mode off
//   a material item  build mode with that material (concrete, stone, wood, iron plate, glass,
//                    dirt): the wheel picks any shape (cube, block, slope, corner, beams)
//   a machine item   build mode placing that machine (or a conveyor / cable)
//   a consumable     used at once (oxygen canister, power cell)
// In the inventory (Tab), hover an item and press 1-0 to bind it. Bindings are saved.

#include <cstdio>
#include <string>

class Scene;

constexpr int kHotbarSlots = 10;
extern const char* const kHotbarTerrainTool; // "tool:terrain"

void Hotbar_Reset();                                  // the default bindings
void Hotbar_Update(Scene& world, bool inputEnabled, bool piloting);
void Hotbar_DrawHUD(bool piloting);
void Hotbar_Bind(int slot, const std::string& id);    // slot 0..9 (key 1 .. key 0)
const std::string& Hotbar_Slot(int slot);
int Hotbar_Selected();
std::string Hotbar_HeldFood();                        // the food or animal feed in hand (a selected slot), or ""
// A number key pressed this frame: its slot (0..9), or -1.
int Hotbar_KeySlot();
void Hotbar_Press(int slot, bool piloting);           // as if its key were pressed (tests)
void Hotbar_Save(FILE* f);                            // "hotbar <slot> <id>" lines
bool Hotbar_LoadLine(const char* tok, FILE* f);       // true if it consumed a hotbar line

#endif
