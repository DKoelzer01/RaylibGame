#ifndef UNIT_HELPERS_H
#define UNIT_HELPERS_H

// Shared helpers for the Tier 1 unit tests (#869). Included only inside `#ifdef RAYLIBGAME_UNIT_TESTS` blocks: the game
// never sees it. Tests run from the repo root (CTest sets the working directory), so resources/ is where the game has it.

#include <gtest/gtest.h>
#include <string>

class Scene;
class ItemDB;

namespace UnitTest {

// The item and recipe tables (resources/data), loaded once per process.
void EnsureItemDB();

// A fresh, empty directory for this test under the system temp dir (unique per process and test); removed with its
// contents at the end of the process. Never inside the repo or saves/.
std::string ScratchDir();

// A Scene for code that only reads its plain fields (renderOrigin, objects). Scene's constructor loads shaders and
// meshes and needs a GL context, which Tier 1 never opens; this is zero-filled storage reinterpreted as a Scene, whose
// constructor never ran. Use it only for functions that read renderOrigin and the like, never for anything that touches
// GL resources, strings or containers of the scene.
Scene& FieldOnlyScene();

} // namespace UnitTest

#endif
