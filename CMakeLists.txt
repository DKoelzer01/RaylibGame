cmake_minimum_required(VERSION 3.24)
project(RaylibGame LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # for clangd / VS Code IntelliSense

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Precompiled headers (#1758): src/pch.h (raylib + the standard headers) is force-included in engine, game_objs, game and
# game_tests. One PCH per target (their flags and defines differ, so REUSE_FROM does not apply).
option(RAYLIBGAME_PCH "Precompiled headers (src/pch.h) for faster builds" ON)
# ccache only caches a PCH build with these sloppiness settings (and GCC's -fpch-preprocess, below): add them when the
# launcher is ccache (the sandbox passes CMAKE_CXX_COMPILER_LAUNCHER=ccache).
if(RAYLIBGAME_PCH AND CMAKE_CXX_COMPILER_LAUNCHER MATCHES "ccache" AND NOT CMAKE_CXX_COMPILER_LAUNCHER MATCHES "CCACHE_SLOPPINESS")
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CMAKE_COMMAND} -E env CCACHE_SLOPPINESS=pch_defines,time_macros,include_file_mtime,include_file_ctime ${CMAKE_CXX_COMPILER_LAUNCHER})
endif()
# Files that include <windows.h> / <winsock2.h> before raylib's names (Rectangle, CloseWindow...) cannot take the force-included raylib.h.
set(PCH_SKIP_FILES src/engine/dynlib.cpp src/engine/profiler.cpp src/engine/utils/logger.cpp src/space/hostprocess.cpp src/space/scenarioproc.cpp)
function(raylibgame_pch target)
    if(NOT RAYLIBGAME_PCH)
        return()
    endif()
    target_precompile_headers(${target} PRIVATE "${CMAKE_SOURCE_DIR}/src/pch.h")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-fpch-preprocess>)
    endif()
endfunction()
foreach(f ${PCH_SKIP_FILES})
    set_source_files_properties("${CMAKE_SOURCE_DIR}/${f}" PROPERTIES SKIP_PRECOMPILE_HEADERS ON)
endforeach()

option(USE_SYSTEM_RAYLIB "Link a prebuilt raylib at RAYLIB_PATH instead of fetching raylib 5.5" OFF)
set(RAYLIB_PATH "C:/raylib/raylib" CACHE PATH "raylib source tree containing src/libraylib.a (USE_SYSTEM_RAYLIB only)")
option(ENABLE_TRACY "Build with the Tracy profiler client (connect with tracy-profiler v0.11.1)" OFF)
# Multiplayer transport (#310): Valve's GameNetworkingSockets, built with this compiler (static, BCrypt crypto, protobuf
# 3.21) by tools/build_gns.sh into the gitignored deps/gns/install (versions pinned in the script): GNS_ROOT holds
# include/steam/*.h, lib/libGameNetworkingSockets_s.a and lib/libprotobuf.a. Off: src/engine/net.* builds as stubs
# (Net_Available() false).
option(ENABLE_GNS "Build with GameNetworkingSockets (multiplayer transport) from GNS_ROOT" OFF)
set(GNS_ROOT "${CMAKE_SOURCE_DIR}/deps/gns/install" CACHE PATH "GameNetworkingSockets + protobuf built by tools/build_gns.sh (ENABLE_GNS only)")
set(GAME_OUTPUT_DIR "${CMAKE_SOURCE_DIR}" CACHE PATH "Where game.exe is written (default: repo root, which is also the working directory they expect)")
# Logging (#514): OFF compiles every LOG_* macro to nothing (arguments not evaluated, no strings in the binary, no
# log.txt) for shipping builds. RAYLIBGAME_LOG_MIN strips the levels below it at compile time (the runtime level,
# --log-level or settings.txt's log_level, filters what is left).
option(RAYLIBGAME_LOGGING "Build with logging (log.txt); OFF for shipping builds" ON)
set(RAYLIBGAME_LOG_MIN "TRACE" CACHE STRING "Lowest log level compiled in: TRACE, DEBUG, INFO, WARN, ERROR or FATAL")
set_property(CACHE RAYLIBGAME_LOG_MIN PROPERTY STRINGS TRACE DEBUG INFO WARN ERROR FATAL)
# Linux dedicated server (#792, #326): raylib built with src/platform/rcore_headless.c instead of GLFW: no X11 or
# Wayland, an off-screen OpenGL 4.3 context from EGL (a GPU in the container) or a null GL (no GPU: CPU paths only).
# The Linux build is the dedicated server; there is no Linux desktop client yet. tools/build_linux.sh builds it.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(HEADLESS_DEFAULT ON)
else()
    set(HEADLESS_DEFAULT OFF)
endif()
option(RAYLIBGAME_HEADLESS "Build raylib with the headless EGL / null-GL backend (Linux dedicated servers)" ${HEADLESS_DEFAULT})

include(FetchContent)

# ---------------------------------------------------------------------------------------------
# raylib
# ---------------------------------------------------------------------------------------------
if(USE_SYSTEM_RAYLIB)
    add_library(raylib STATIC IMPORTED)
    set_target_properties(raylib PROPERTIES
        IMPORTED_LOCATION "${RAYLIB_PATH}/src/libraylib.a"
        INTERFACE_INCLUDE_DIRECTORIES "${RAYLIB_PATH}/src")
    if(WIN32)
        set_property(TARGET raylib PROPERTY INTERFACE_LINK_LIBRARIES opengl32 gdi32 winmm)
    endif()
elseif(RAYLIBGAME_HEADLESS)
    # raylib's sources only (SOURCE_SUBDIR names no CMakeLists, so raylib's own build with GLFW is not added); rcore.c
    # is compiled inside src/platform/rcore_headless.c, which adds the platform backend.
    FetchContent_Declare(raylib
        GIT_REPOSITORY https://github.com/raysan5/raylib.git
        GIT_TAG 5.5
        GIT_SHALLOW TRUE
        SOURCE_SUBDIR no-cmake-build)
    FetchContent_MakeAvailable(raylib)
    set(RL_SRC "${raylib_SOURCE_DIR}/src")
    add_library(raylib STATIC src/platform/rcore_headless.c
        ${RL_SRC}/rshapes.c ${RL_SRC}/rtextures.c ${RL_SRC}/rtext.c ${RL_SRC}/rmodels.c ${RL_SRC}/utils.c ${RL_SRC}/raudio.c)
    set_target_properties(raylib PROPERTIES C_STANDARD 99 C_EXTENSIONS ON)
    target_include_directories(raylib PUBLIC "${RL_SRC}")
    target_compile_definitions(raylib PUBLIC GRAPHICS_API_OPENGL_43)
    target_compile_options(raylib PRIVATE -w) # raylib's own warnings are not ours to fix
    target_link_libraries(raylib PUBLIC m dl pthread)
else()
    set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(BUILD_GAMES OFF CACHE BOOL "" FORCE)
    # OpenGL 4.3 core (Plan: GPU compute for the fluid sim and terrain): compute shaders and SSBOs through
    # rlgl. The game checks rlGetVersion() and keeps CPU paths for 3.3-only GPUs.
    set(OPENGL_VERSION "4.3" CACHE STRING "" FORCE)
    FetchContent_Declare(raylib
        GIT_REPOSITORY https://github.com/raysan5/raylib.git
        GIT_TAG 5.5
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(raylib)
endif()

# ---------------------------------------------------------------------------------------------
# Dear ImGui + rlImGui (vendored, see third_party/README.md)
# ---------------------------------------------------------------------------------------------
add_library(imgui STATIC
    third_party/imgui/imgui.cpp
    third_party/imgui/imgui_draw.cpp
    third_party/imgui/imgui_tables.cpp
    third_party/imgui/imgui_widgets.cpp
    third_party/imgui/imgui_demo.cpp
    third_party/rlImGui/rlImGui.cpp)
target_include_directories(imgui PUBLIC third_party/imgui third_party/rlImGui)
target_compile_definitions(imgui PUBLIC NO_FONT_AWESOME)
target_link_libraries(imgui PUBLIC raylib)

# ---------------------------------------------------------------------------------------------
# Optional Tracy profiler client
# ---------------------------------------------------------------------------------------------
if(ENABLE_TRACY)
    set(TRACY_ENABLE ON CACHE BOOL "" FORCE)
    set(TRACY_ON_DEMAND ON CACHE BOOL "" FORCE)
    FetchContent_Declare(tracy
        GIT_REPOSITORY https://github.com/wolfpld/tracy.git
        GIT_TAG v0.11.1
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(tracy)
endif()

# ---------------------------------------------------------------------------------------------
# Engine (shared by both games) + the two game executables
# ---------------------------------------------------------------------------------------------
file(GLOB_RECURSE ENGINE_SOURCES CONFIGURE_DEPENDS src/engine/*.cpp)
add_library(engine STATIC ${ENGINE_SOURCES})
raylibgame_pch(engine)
target_include_directories(engine PUBLIC src/engine)
target_link_libraries(engine PUBLIC raylib imgui)
target_compile_definitions(engine PUBLIC PLATFORM_DESKTOP)
set(LOG_LEVELS TRACE DEBUG INFO WARN ERROR FATAL)
list(FIND LOG_LEVELS "${RAYLIBGAME_LOG_MIN}" LOG_MIN_INDEX)
if(LOG_MIN_INDEX LESS 0)
    message(FATAL_ERROR "RAYLIBGAME_LOG_MIN must be TRACE, DEBUG, INFO, WARN, ERROR or FATAL (got ${RAYLIBGAME_LOG_MIN})")
endif()
if(RAYLIBGAME_LOGGING)
    target_compile_definitions(engine PUBLIC RAYLIBGAME_LOGGING=1 RAYLIBGAME_LOG_MIN=${LOG_MIN_INDEX})
else()
    target_compile_definitions(engine PUBLIC RAYLIBGAME_LOGGING=0 RAYLIBGAME_LOG_MIN=${LOG_MIN_INDEX})
endif()
if(ENABLE_TRACY)
    target_link_libraries(engine PUBLIC Tracy::TracyClient) # also defines TRACY_ENABLE for everything
endif()
if(ENABLE_GNS)
    if(NOT EXISTS "${GNS_ROOT}/lib/libGameNetworkingSockets_s.a" OR NOT EXISTS "${GNS_ROOT}/lib/libprotobuf.a")
        message(FATAL_ERROR "ENABLE_GNS: GNS_ROOT (${GNS_ROOT}) needs lib/libGameNetworkingSockets_s.a and lib/libprotobuf.a: run tools/build_gns.sh first")
    endif()
    target_include_directories(engine PUBLIC "${GNS_ROOT}/include")
    target_compile_definitions(engine PUBLIC RAYLIBGAME_GNS STEAMNETWORKINGSOCKETS_STATIC_LINK)
    target_link_libraries(engine PUBLIC "${GNS_ROOT}/lib/libGameNetworkingSockets_s.a" "${GNS_ROOT}/lib/libprotobuf.a")
    if(WIN32)
        target_link_libraries(engine PUBLIC ws2_32 winmm bcrypt iphlpapi crypt32)
    else()
        # Linux: GNS uses OpenSSL's libcrypto (tools/build_gns.sh), linked statically so the server needs no libssl.
        set(OPENSSL_USE_STATIC_LIBS TRUE)
        find_package(OpenSSL REQUIRED COMPONENTS Crypto)
        target_link_libraries(engine PUBLIC OpenSSL::Crypto dl pthread)
    endif()
endif()
# Steam Audio (#1020, #1292): loaded at run time (src/engine/steamaudio.cpp), so building needs only its vendored headers.
# tools/fetch_steamaudio.sh puts phonon.dll in deps/steamaudio; when it is there (at configure time) the game build copies
# it next to the exe. Without it the game uses its own panning.
target_include_directories(engine PRIVATE third_party/steamaudio/include)
target_link_libraries(engine PUBLIC ${CMAKE_DL_LIBS})
set(STEAMAUDIO_DLL "${CMAKE_SOURCE_DIR}/deps/steamaudio/lib/windows-x64/phonon.dll")
if(RAYLIBGAME_HEADLESS)
    target_compile_definitions(engine PUBLIC RAYLIBGAME_HEADLESS_PLATFORM=1) # RaylibHeadless_GLMode() exists (gpucompute.cpp)
endif()
if(WIN32)
    target_link_libraries(engine PUBLIC ws2_32 iphlpapi) # (#341) the menu's Host checks its port is free (hostprocess.cpp), in every build
    target_link_libraries(engine PUBLIC dbghelp) # (#1526) the crash handler's minidump (logger.cpp)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(engine PRIVATE -Wall -Wno-missing-braces -Werror=format -Werror=format-security -ffp-contract=off) # no FMA contraction: the GPU terrain port (terrain.glsl) matches the CPU bit for bit
endif()

function(add_game target dir)
    file(GLOB GAME_SOURCES CONFIGURE_DEPENDS ${dir}/*.cpp)
    # Everything but main.cpp is an object library (#869), so the unit tests (game_tests) link the same objects
    # without compiling them twice. The exe still gets every object (object libraries are never dropped by the linker).
    list(FILTER GAME_SOURCES EXCLUDE REGEX "/main\\.cpp$")
    add_library(${target}_objs OBJECT ${GAME_SOURCES})
    raylibgame_pch(${target}_objs)
    target_include_directories(${target}_objs PUBLIC ${dir})
    target_link_libraries(${target}_objs PUBLIC engine)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target}_objs PRIVATE -Wall -Wno-missing-braces -Werror=format -Werror=format-security)
    endif()
    # Tier 2 scenarios (#870): tests/scenarios/*.cpp register themselves with SCENARIO(...), so they are compiled into the
    # exe itself (not the object library's archive, whose unreferenced members game_tests' linker would drop).
    file(GLOB SCENARIO_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/tests/scenarios/*.cpp)
    add_executable(${target} ${dir}/main.cpp ${SCENARIO_SOURCES})
    raylibgame_pch(${target})
    target_include_directories(${target} PRIVATE ${dir})
    target_link_libraries(${target} PRIVATE ${target}_objs engine)
    # Version for the window title (#501): branch + commit, regenerated at every build (cmake/game_version.cmake).
    set(version_src "${CMAKE_BINARY_DIR}/generated/${target}_version.cpp")
    add_custom_target(${target}_version
        COMMAND ${CMAKE_COMMAND} -DSRC=${CMAKE_SOURCE_DIR} -DOUT=${version_src} -P ${CMAKE_SOURCE_DIR}/cmake/game_version.cmake
        BYPRODUCTS ${version_src}
        COMMENT "" VERBATIM)
    target_sources(${target} PRIVATE ${version_src})
    add_dependencies(${target} ${target}_version)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -Wall -Wno-missing-braces -Werror=format -Werror=format-security)
    endif()
    if(MINGW)
        # Self-contained exe: no libstdc++/libgcc DLLs needed next to it.
        target_link_options(${target} PRIVATE -static)
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        # glibc stays shared (dlopen of libEGL needs it); the C++ runtime goes in, so any glibc 2.36+ distro runs it.
        target_link_options(${target} PRIVATE -static-libstdc++ -static-libgcc)
    endif()
    if(WIN32 AND EXISTS "${STEAMAUDIO_DLL}") # (#1292) Steam Audio next to the exe (it is found there first)
        add_custom_command(TARGET ${target} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${STEAMAUDIO_DLL}" "${GAME_OUTPUT_DIR}/phonon.dll" VERBATIM)
    endif()
    # Assets load relative to the working directory: put the exe at the repo root.
    set_target_properties(${target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${GAME_OUTPUT_DIR}"
        RUNTIME_OUTPUT_DIRECTORY_DEBUG "${GAME_OUTPUT_DIR}"
        RUNTIME_OUTPUT_DIRECTORY_RELEASE "${GAME_OUTPUT_DIR}"
        VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endfunction()

add_game(game src/space)        # space survival

# ---------------------------------------------------------------------------------------------
# Tier 1 unit tests (#869): GoogleTest, tests next to the code
# ---------------------------------------------------------------------------------------------
# A tested .cpp ends with its tests inside `#ifdef RAYLIBGAME_UNIT_TESTS`. The game never defines it (game.exe is
# unchanged and does not link GoogleTest). game_tests compiles only the files that contain such a block again, with the
# define, and links every other object from the game's own libraries: the linker takes a library member only for a
# symbol nobody defines yet, and the test build of a file defines everything its plain build does. Re-run cmake after
# adding the first test block to a file (the scan runs at configure time). Run: ctest --test-dir build/check, or
# build/check/bin/game_tests.exe (from the repo root: some tests read resources/).
option(RAYLIBGAME_UNIT_TESTS "Build the game_tests unit test exe (GoogleTest, fetched at a pinned version)" ON)
if(RAYLIBGAME_UNIT_TESTS)
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.15.2
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(googletest)

    file(GLOB_RECURSE UNIT_SCAN CONFIGURE_DEPENDS src/engine/*.cpp src/space/*.cpp)
    set(UNIT_TEST_SOURCES "")
    foreach(f ${UNIT_SCAN})
        file(STRINGS "${f}" hit REGEX "^#ifdef RAYLIBGAME_UNIT_TESTS" LIMIT_COUNT 1)
        if(hit)
            list(APPEND UNIT_TEST_SOURCES "${f}")
        endif()
    endforeach()
    list(LENGTH UNIT_TEST_SOURCES n_unit)
    message(STATUS "game_tests: ${n_unit} source files carry unit tests")

    add_library(game_space STATIC $<TARGET_OBJECTS:game_objs>) # the game's objects as a library (no recompile)
    add_executable(game_tests ${UNIT_TEST_SOURCES} tests/unit/unit_helpers.cpp "${CMAKE_BINARY_DIR}/generated/game_version.cpp")
    add_dependencies(game_tests game_version) # kGameVersion (menus.cpp)
    raylibgame_pch(game_tests)
    target_compile_definitions(game_tests PRIVATE RAYLIBGAME_UNIT_TESTS=1)
    target_include_directories(game_tests PRIVATE src/engine src/space tests/unit)
    # The game objects and the engine refer to each other: rescan both until every symbol resolves.
    target_link_libraries(game_tests PRIVATE "$<LINK_GROUP:RESCAN,game_space,engine>" GTest::gtest_main)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(game_tests PRIVATE -Wall -Wno-missing-braces -Werror=format -Werror=format-security -ffp-contract=off)
    endif()
    if(MINGW)
        target_link_options(game_tests PRIVATE -static)
    endif()
    set_target_properties(game_tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${GAME_OUTPUT_DIR}")
    enable_testing()
    include(GoogleTest)
    # One CTest entry per test, found at build time; run from the repo root (resources/ and the data files).
    gtest_discover_tests(game_tests WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}" DISCOVERY_TIMEOUT 60 PROPERTIES TIMEOUT 30)
endif()
