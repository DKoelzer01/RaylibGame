#ifndef ITEMICONS_H
#define ITEMICONS_H

// Item icons: one 64 px picture per item in items.txt, painted at startup from anti-aliased shapes
// (no image files): ores as rock with crystals of the ore's colour, ingots and plates, gems, logs,
// planks, bottles of liquid, foods, circuits, tools with their material's head, a machine front with
// a symbol of what it does, rocket parts, suits and modules. The item's colour drives each picture;
// unknown items fall back to a shape by kind (a cube for blocks, a lump for resources). All icons
// share one mipmapped atlas texture.

#include <raylib.h>
#include <string>

void ItemIcons_Init();                                    // after InitWindow and itemDB.load
void ItemIcons_Unload();
bool ItemIcons_Has(const std::string& id);
Texture2D ItemIcons_Texture();
Rectangle ItemIcons_Source(const std::string& id);        // the icon's cell in the atlas (pixels)
void ItemIcons_Draw(const std::string& id, Rectangle dst, Color tint = WHITE); // raylib 2D
// ImGui: uv corners of the icon (ImGui::Image / ImDrawList::AddImage with ItemIcons_Texture().id).
void ItemIcons_UV(const std::string& id, float uv0[2], float uv1[2]);
bool ItemIcons_SaveSheet(const char* path);              // tests: every icon in one PNG

#endif
