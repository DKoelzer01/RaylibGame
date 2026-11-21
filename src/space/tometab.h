#pragma once

// The Tome tab (#1251, magic step 4; the magic doc section 2): the game window's page for magic, once the player has found
// it (RAYLIBGAME_MAGIC). The reader's aspects, the team's pages in four branches (owned, ready to write, locked and why),
// writing a page with a click, and spending carried items for their aspects once the team has the Research Table page.

#include <raylib.h>
#include <string>

bool TomeTab_Shown();  // the tab is in the window: magic on and found
void TomeTab_Draw();   // inside the game window's screen (ImGui)
// Tests: where a page's box was drawn this frame (screen pixels); false if it wasn't.
bool TomeTab_TestPageCenter(const std::string& id, Vector2& center);
