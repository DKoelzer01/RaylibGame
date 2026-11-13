#ifndef MENUS_H
#define MENUS_H

// Main and pause menus: clickable buttons (raylib immediate mode), a save list, and settings.
//
// Main menu: New Game, Load Game (the saves in saves/, newest first), Settings, Exit.
// Pause menu (Esc in game): Resume, Save, Settings, Save & Quit to Menu; Esc resumes (or backs
// out of Settings). Settings (mouse sensitivity, field of view, volume, foliage, frame cap) are
// kept in settings.txt.

#include <string>

enum class MenuAction { None, NewGame, LoadGame, Resume, Save, SaveAndQuit, Exit };
struct MenuResult {
    MenuAction action = MenuAction::None;
    std::string path;            // LoadGame: the save file
};

struct GameSettings {
    float mouseSensitivity = 1.0f;  // multiplier on the default look speed
    float fov = 65.0f;              // degrees (vertical; narrow views are a motion-sickness trigger)
    float volume = 0.35f;
    bool foliage = true;
    int fpsCap = 0;                 // 0: the monitor's refresh rate
    bool vsync = true;              // no tearing (caps at the monitor's refresh rate)
    bool reduceMotion = false;      // eased walking starts/stops and a smoothed camera up (motion sickness)
    bool showHands = true;          // the hands and held tool (a fixed reference in view: motion sickness)
    int handBob = 0;                // 0 full, 1 reduced, 2 off
    int terrainDetail = 1;          // 0 Low, 1 Medium, 2 High, 3 Ultra: LOD factor 2 / 2.5 / 3 / 4 (#297)
    int difficulty = 1;             // new worlds: 0 peaceful (no threats, slower hunger, extra food), 1 normal, 2 hard (no Map, hungrier, fiercer nights)
};
extern GameSettings settings;
void Settings_Load();            // settings.txt (missing: defaults)
void Settings_Save();
void Settings_Apply();           // to the player, ship, camera, audio, foliage and frame cap

// Draw (between BeginDrawing/EndDrawing) and handle clicks. Esc is handled here too.
MenuResult Menu_Main();
MenuResult Menu_Pause();
void Menu_ResetPause();          // back to the pause menu's first page when it opens
// A new save's file: saves/world.sav, else the first free saves/world_N.sav.
std::string Menu_NewSavePath();

#endif
