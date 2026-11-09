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
# Game
# ---------------------------------------------------------------------------------------------
file(GLOB GAME_SOURCES CONFIGURE_DEPENDS src/*.cpp src/utils/*.cpp)
add_executable(game ${GAME_SOURCES})
target_include_directories(game PRIVATE src)
target_link_libraries(game PRIVATE raylib imgui)
target_compile_definitions(game PRIVATE PLATFORM_DESKTOP)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(game PRIVATE -Wall -Wno-missing-braces)
endif()
if(MINGW)
    # Self-contained exe: no libstdc++/libgcc DLLs needed next to it.
    target_link_options(game PRIVATE -static)
endif()

if(ENABLE_TRACY)
    set(TRACY_ENABLE ON CACHE BOOL "" FORCE)
    set(TRACY_ON_DEMAND ON CACHE BOOL "" FORCE)
    FetchContent_Declare(tracy
        GIT_REPOSITORY https://github.com/wolfpld/tracy.git
        GIT_TAG v0.11.1
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(tracy)
    target_link_libraries(game PRIVATE Tracy::TracyClient) # also defines TRACY_ENABLE
endif()

# The game loads resources/ and assets/ relative to the working directory, so put the exe
# at the repo root (where it always lived) and point IDE debuggers there too.
set_target_properties(game PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_SOURCE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}"
    VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
