cmake_minimum_required(VERSION 3.24)
project(RaylibGame LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # for clangd / VS Code IntelliSense

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(USE_SYSTEM_RAYLIB "Link a prebuilt raylib at RAYLIB_PATH instead of fetching raylib 5.5" OFF)
set(RAYLIB_PATH "C:/raylib/raylib" CACHE PATH "raylib source tree containing src/libraylib.a (USE_SYSTEM_RAYLIB only)")
option(ENABLE_TRACY "Build with the Tracy profiler client (connect with tracy-profiler v0.11.1)" OFF)
set(GAME_OUTPUT_DIR "${CMAKE_SOURCE_DIR}" CACHE PATH "Where game.exe is written (default: repo root, which is also the working directory they expect)")

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
target_include_directories(engine PUBLIC src/engine)
target_link_libraries(engine PUBLIC raylib imgui)
target_compile_definitions(engine PUBLIC PLATFORM_DESKTOP)
if(ENABLE_TRACY)
    target_link_libraries(engine PUBLIC Tracy::TracyClient) # also defines TRACY_ENABLE for everything
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(engine PRIVATE -Wall -Wno-missing-braces)
endif()

function(add_game target dir)
    file(GLOB GAME_SOURCES CONFIGURE_DEPENDS ${dir}/*.cpp)
    add_executable(${target} ${GAME_SOURCES})
    target_include_directories(${target} PRIVATE ${dir})
    target_link_libraries(${target} PRIVATE engine)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -Wall -Wno-missing-braces)
    endif()
    if(MINGW)
        # Self-contained exe: no libstdc++/libgcc DLLs needed next to it.
        target_link_options(${target} PRIVATE -static)
    endif()
    # Assets load relative to the working directory: put the exe at the repo root.
    set_target_properties(${target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${GAME_OUTPUT_DIR}"
        RUNTIME_OUTPUT_DIRECTORY_DEBUG "${GAME_OUTPUT_DIR}"
        RUNTIME_OUTPUT_DIRECTORY_RELEASE "${GAME_OUTPUT_DIR}"
        VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endfunction()

add_game(game src/space)        # space survival
