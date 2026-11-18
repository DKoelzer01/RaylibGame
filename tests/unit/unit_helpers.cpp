// Shared helpers for the Tier 1 unit tests (#869); see unit_helpers.h.
#include "unit_helpers.h"
#include "inventory.h"
#include "scene.h"

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <mutex>
#include <new>

namespace UnitTest {

void EnsureItemDB() {
    static std::once_flag once;
    std::call_once(once, [] {
        const bool ok = itemDB.load("resources/data/items.txt", "resources/data/recipes.txt");
        ASSERT_TRUE(ok) << "resources/data not found: run game_tests from the repo root (ctest does)";
    });
}

namespace {
    // The root and its removal at exit in one function-local static (a separate static's destructor ran after the path's).
    struct ScratchRoot { std::filesystem::path path; ~ScratchRoot() { std::error_code ec; if (!path.empty()) std::filesystem::remove_all(path, ec); } };
    std::filesystem::path& Root() { static ScratchRoot root; return root.path; }
}

std::string ScratchDir() {
    namespace fs = std::filesystem;
    static std::atomic<int> counter{ 0 };
    std::error_code ec;
    if (Root().empty()) {
        // Unique per process (several game_tests may run at once): the clock's ticks at the first call.
        Root() = fs::temp_directory_path(ec) / ("raylibgame_unit_" + std::to_string((long long)std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::remove_all(Root(), ec);
        fs::create_directories(Root(), ec);
    }
    const ::testing::TestInfo* info = ::testing::UnitTest::GetInstance()->current_test_info();
    std::string name = info ? std::string(info->test_suite_name()) + "." + info->name() : "test";
    for (char& c : name) if (c == '/' || c == '\\' || c == ':') c = '_';
    const fs::path dir = Root() / (name + "_" + std::to_string(counter++));
    fs::create_directories(dir, ec);
    return dir.generic_string();
}

Scene& FieldOnlyScene() {
    alignas(Scene) static unsigned char storage[sizeof(Scene)] = {};
    return *std::launder(reinterpret_cast<Scene*>(storage));
}

} // namespace UnitTest
