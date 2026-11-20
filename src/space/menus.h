#ifndef MENUS_H
#define MENUS_H

// Main and pause menus: clickable buttons (raylib immediate mode), a save list, and settings.
//
// Main menu: New Game, Load Game (the saves in saves/, newest first), Multiplayer (#341: host a world, join by IP),
// Settings, Exit.
// Pause menu (Esc in game): Resume, Save, Settings, Save & Quit to Menu; Esc resumes (or backs
// out of Settings). Settings (mouse sensitivity, field of view, volume, foliage, frame cap) are
// kept in settings.txt.

#include "raylib.h"
#include <cstdint>
#include <string>
#include <vector>

enum class MenuAction { None, NewGame, LoadGame, Resume, Save, SaveAndQuit, Exit, ReportBug, HostGame, JoinGame }; // ReportBug: the bug report window (#169)
struct MenuResult {
    MenuAction action = MenuAction::None;
    std::string path;            // LoadGame: the save file; HostGame: the world to host ("": a new one)
    // Multiplayer (#341). JoinGame: the server ("ip[:port]"), our name, its password, our suit colour (0xRRGGBB, 0 the
    // default). HostGame: the port, the password, max players; then it joins its own server under name/suitColor.
    std::string address, name, password;
    int port = 27015, maxPlayers = 8;
    uint32_t suitColor = 0;
};

struct GameSettings {
    float mouseSensitivity = 1.0f;  // multiplier on the default look speed
    float fov = 65.0f;              // degrees (vertical; narrow views are a motion-sickness trigger)
    float volume = 0.35f;
    bool foliage = true;
    int fpsCap = 0;                 // 0: the monitor's refresh rate
    bool vsync = true;              // no tearing (caps at the monitor's refresh rate)
    bool showFps = false;           // the FPS counter, on every screen (#572)
    bool reduceMotion = false;      // eased walking starts/stops and a smoothed camera up (motion sickness)
    bool showHands = true;          // the hands and held tool (a fixed reference in view: motion sickness)
    int handBob = 0;                // 0 full, 1 reduced, 2 off
    int terrainDetail = 1;          // 0 Low, 1 Medium, 2 High, 3 Ultra: LOD factor 2 / 2.5 / 3 / 4 (#297)
    int cameraSmoothing = 2;        // 0 Off, 1 Low, 2 Default, 3 High: eye smoothing strength (#168)
    bool speedVignette = false;     // darken the screen edges when moving or turning fast (#168)
    bool fireFlicker = true;        // torches and campfires flicker their light (#558: off = a steady glow)
    float uiScale = 0.0f;            // (#821) the interface size: 0 auto (by the screen's height), else 0.75-1.5
    int keyHints = 0;                // (#827) the key hints line: 0 auto (fades after the first 2 hours played), 1 always, 2 off
    float xpRate = 1.0f;              // (#883) the last New World screen's skill XP rate (the next one's default)
    bool invertSlotWheel = false;     // (#862) the inventory's scroll wheel: false = down sends one to the other side, up pulls one in
    int tonemap = 1;                // #131: the HDR curve, as the Graphics tab lists them: 0 Clip, 1 Knee (the default: today's look), 2 ACES
    int difficulty = 1;             // new worlds: 0 peaceful (no threats, slower hunger, extra food), 1 normal, 2 hard (no Map, hungrier, fiercer nights)
    std::string logLevel;           // settings.txt "log_level trace|debug|info|warn|error|fatal|off" (#514; --log-level wins); empty: the build default
    // Multiplayer (#341): what the Host and Join pages last used (passwords are never kept).
    std::string mpName = "player", mpAddress = "127.0.0.1:27015";
    uint32_t mpSuit = 0;            // 0xRRGGBB; 0 the default suit
    std::vector<std::string> mpRecent; // servers joined, newest first (at most 5)
    int mpPort = 27015, mpMaxPlayers = 8;
    std::string playerUuid;          // (#346) this install's player id on every server (made once; a name is only what others see)
};
extern GameSettings settings;
void Settings_Load();            // settings.txt (missing: defaults)
void Settings_Save();
void Settings_Apply();           // to the player, ship, camera, audio, foliage and frame cap

// Draw (between BeginDrawing/EndDrawing) and handle clicks. Esc is handled here too.
MenuResult Menu_Main();
MenuResult Menu_Pause();
void Menu_ResetPause();          // back to the pause menu's first page when it opens
bool Menu_TestSaves();           // --test-saves (#279): copy, rename and delete a scratch save with regions and bake
void Menu_TestSettings(int tab, bool pause); // --shot-settings (#512): open the settings page on a tab (0 Graphics .. 3 Gameplay)
void Menu_SettingsLayout(float& contentBottom, Rectangle& back); // the last drawn page: its rows' bottom and the Back button
// A new save's file: saves/world.sav, else the first free saves/world_N.sav.
std::string Menu_NewSavePath();
// (#1110) A world the menus offer to load or host: a .sav that isn't the bench's and isn't a joined game's cache of a
// server's world (a .synced beside it; hosting one started a fresh world and player).
bool Menu_IsWorldSave(const std::string& path);
// Multiplayer (#341): back on the Join (or Host) page with why it failed (refused, timed out, the server didn't start).
void Menu_MultiplayerFailed(bool host, const std::string& why);
void Menu_RememberServer(const std::string& address); // a server joined: first in the recent list
// Tests (--test-mp-menu): open a multiplayer page (0 the Multiplayer page, 1 Join, 2 Host).
void Menu_TestMultiplayer(int page);
// The suit colours the Join page offers (0xRRGGBB; the first, 0, is the default suit).
int Menu_SuitColorCount();
uint32_t Menu_SuitColor(int i);

#endif
