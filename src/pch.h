#ifndef RAYLIBGAME_PCH_H
#define RAYLIBGAME_PCH_H

// Precompiled header (#1758): the heavy, stable headers that most files include. CMake force-includes it first in every
// file of engine, game_objs, game and game_tests (option RAYLIBGAME_PCH). Keep only third-party and standard headers
// here: a project header in this list would rebuild everything on every edit. Files that need <windows.h> before
// raylib's names appear are excluded in CMakeLists.txt (SKIP_PRECOMPILE_HEADERS).

#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#endif
