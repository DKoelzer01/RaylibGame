#include "menus.h"
#include "framepacer.h"
#include "playerinput.h"
#include "magic.h" // (#1249) Magic_Enabled: the magic sight key only with magic on
#include "playerrecord.h"
#include "scenetarget.h"
#include "slideshow.h"
#include "player.h"
#include "ship.h"
#include "telescope.h"
#include "audio.h"
#include "combat.h" // (#1218) gDoubleTapDodge
#include "subtitles.h"
#include "uipalette.h"
#include "foliage.h"
#include "world.h"
#include "headless.h"
#include "bugreport.h"
#include "simclock.h"

#include <raylib.h>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

GameSettings settings;

namespace {
    constexpr float kDefaultMouse = 0.003f, kDefaultShipMouse = 0.0025f;
    enum class Page { Main, Load, Settings, Multiplayer, Join, Host };
    Page mainPage = Page::Main, pausePage = Page::Main;

    struct SaveEntry { std::string path, name, when; fs::file_time_type time; };
    std::vector<SaveEntry> saves;
    float saveScroll = 0.0f;

    const Color kPanel = { 18, 24, 34, 230 }, kButton = { 44, 58, 80, 255 }, kHover = { 70, 96, 132, 255 }, kText = { 230, 236, 245, 255 };

    float Scale() { return std::max(1.0f, GetScreenHeight() / 1080.0f); }

    bool Button(Rectangle r, const char* label, bool enabled = true) {
        const bool hover = enabled && CheckCollisionPointRec(GetMousePosition(), r);
        DrawRectangleRounded(r, 0.18f, 6, enabled ? (hover ? kHover : kButton) : Fade(kButton, 0.4f));
        if (hover) DrawRectangleRoundedLinesEx(r, 0.18f, 6, 2.0f, Fade(SKYBLUE, 0.8f));
        const int fs = (int)(r.height * 0.45f);
        DrawText(label, (int)(r.x + (r.width - MeasureText(label, fs)) * 0.5f), (int)(r.y + (r.height - fs) * 0.5f), fs, enabled ? kText : Fade(kText, 0.4f));
        return hover && IsMouseButtonReleased(MOUSE_BUTTON_LEFT);
    }
    // A slider row: label, bar, value text. Returns true while it changed.
    bool Slider(Rectangle r, const char* label, float& v, float lo, float hi, const char* fmt) {
        const int fs = (int)(r.height * 0.5f);
        DrawText(label, (int)r.x, (int)(r.y + (r.height - fs) * 0.5f), fs, kText);
        const Rectangle bar = { r.x + r.width * 0.38f, r.y + r.height * 0.4f, r.width * 0.44f, r.height * 0.2f };
        DrawRectangleRec(bar, kButton);
        const float f = (v - lo) / (hi - lo);
        DrawRectangle((int)bar.x, (int)bar.y, (int)(bar.width * f), (int)bar.height, SKYBLUE);
        DrawCircle((int)(bar.x + bar.width * f), (int)(bar.y + bar.height * 0.5f), r.height * 0.22f, kText);
        DrawText(TextFormat(fmt, v), (int)(bar.x + bar.width + 16), (int)(r.y + (r.height - fs) * 0.5f), fs, kText);
        const Rectangle hit = { bar.x - 10, r.y, bar.width + 20, r.height };
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) && CheckCollisionPointRec(GetMousePosition(), hit)) {
            const float nf = std::clamp((GetMousePosition().x - bar.x) / bar.width, 0.0f, 1.0f);
            v = lo + nf * (hi - lo);
            return true;
        }
        return false;
    }
    bool Toggle(Rectangle r, const char* label, bool& v) {
        const int fs = (int)(r.height * 0.5f);
        DrawText(label, (int)r.x, (int)(r.y + (r.height - fs) * 0.5f), fs, kText);
        const Rectangle b = { r.x + r.width * 0.38f, r.y + r.height * 0.1f, r.height * 2.2f, r.height * 0.8f };
        if (Button(b, v ? "On" : "Off")) { v = !v; return true; }
        return false;
    }

    void Backdrop(const char* title) {
        const float s = Scale();
        DrawRectangleGradientV(0, 0, GetScreenWidth(), GetScreenHeight(), { 8, 12, 24, 255 }, { 20, 34, 56, 255 });
        Slideshow_Draw(0.4f); // screenshots behind the menu (#386; none: the gradient shows)
        const int fs = (int)(64 * s);
        DrawText(title, (GetScreenWidth() - MeasureText(title, fs)) / 2, (int)(GetScreenHeight() * 0.18f), fs, kText);
    }
    Rectangle Row(int i, float width = 420.0f, float top = 0.36f) {
        const float s = Scale(), w = width * s, h = 58.0f * s;
        return { (GetScreenWidth() - w) * 0.5f, GetScreenHeight() * top + i * (h + 14.0f * s), w, h };
    }

    // Save management (#279), from the main menu only (no world loaded). A save is every entry of saves/ named
    // "<name>.sav" or "<name>.sav.<something>": the file, its .regions folder and its .erosion bake cache.
    std::vector<fs::path> SaveParts(const std::string& name) {
        std::vector<fs::path> out;
        std::error_code ec;
        const std::string base = name + ".sav";
        for (const auto& e : fs::directory_iterator("saves", ec)) {
            const std::string f = e.path().filename().string();
            if (f == base || f.rfind(base + ".", 0) == 0) out.push_back(e.path());
        }
        return out;
    }
    bool SaveNameOk(const std::string& n) {
        if (n.empty() || n.size() > 40 || n == "bench") return false;
        for (char c : n) if (!(isalnum((unsigned char)c) || c == ' ' || c == '_' || c == '-')) return false;
        return SaveParts(n).empty(); // not taken
    }
    bool SaveRename(const std::string& from, const std::string& to) {
        std::error_code ec;
        for (const fs::path& p : SaveParts(from)) {
            const std::string f = p.filename().string();
            fs::rename(p, fs::path("saves") / (to + f.substr(from.size())), ec);
            if (ec) return false;
        }
        return true;
    }
    bool SaveCopy(const std::string& from, const std::string& to) {
        std::error_code ec;
        for (const fs::path& p : SaveParts(from)) {
            const std::string f = p.filename().string();
            fs::copy(p, fs::path("saves") / (to + f.substr(from.size())), fs::copy_options::recursive, ec);
            if (ec) return false;
        }
        return true;
    }
    bool SaveDelete(const std::string& name) {
        std::error_code ec;
        for (const fs::path& p : SaveParts(name)) fs::remove_all(p, ec);
        return SaveParts(name).empty();
    }
    bool TestSaveOps() {
        std::error_code ec;
        fs::create_directories("saves/zz_test.sav.regions", ec);
        { FILE* f = fopen("saves/zz_test.sav", "w"); if (f) { fputs("test", f); fclose(f); } }
        { FILE* f = fopen("saves/zz_test.sav.regions/b1_0_0_0.rgn", "w"); if (f) { fputs("region", f); fclose(f); } }
        { FILE* f = fopen("saves/zz_test.sav.erosion", "w"); if (f) { fputs("bake", f); fclose(f); } }
        bool ok = true;
        auto check = [&](bool c, const char* what) { if (!c) { printf("bench: saves: %s: NO\n", what); ok = false; } };
        check(SaveParts("zz_test").size() == 3, "a save is its file, regions and bake");
        check(!SaveNameOk("zz_test") && !SaveNameOk("bad/name") && !SaveNameOk("bench") && SaveNameOk("zz_test2"), "names: taken, bad characters and bench refused");
        check(SaveCopy("zz_test", "zz_test2") && fs::exists("saves/zz_test2.sav.regions/b1_0_0_0.rgn", ec) && fs::exists("saves/zz_test2.sav.erosion", ec), "copy takes all three");
        check(SaveRename("zz_test2", "zz_test3") && SaveParts("zz_test2").empty() && SaveParts("zz_test3").size() == 3, "rename moves all three");
        check(SaveDelete("zz_test3") && SaveDelete("zz_test") && SaveParts("zz_test").empty(), "delete removes all of it");
        printf("bench: saves (#279): copy, rename, delete of a save with regions and bake -> %s\n", ok ? "PASS" : "FAIL");
        return ok;
    }

    // The Load page's row state: which save is being renamed (its text) or waits for a second click to delete.
    int renaming = -1, deleting = -1;
    std::string renameText;
    double deleteAsked = 0.0;
    std::string saveMessage;

    void RefreshSaves() {
        saves.clear();
        std::error_code ec;
        if (!fs::exists("saves", ec)) return;
        for (const auto& e : fs::directory_iterator("saves", ec)) {
            if (!e.is_regular_file() || !Menu_IsWorldSave(e.path().string())) continue;
            SaveEntry s;
            s.path = "saves/" + e.path().filename().string();
            s.name = e.path().stem().string();
            s.time = e.last_write_time(ec);
            // file_time -> time_t (C++17: via the clocks' current times).
            const auto sys = std::chrono::time_point_cast<std::chrono::system_clock::duration>(s.time - fs::file_time_type::clock::now() + std::chrono::system_clock::now());
            const std::time_t tt = std::chrono::system_clock::to_time_t(sys);
            char buf[64];
            std::strftime(buf, sizeof buf, "%Y-%m-%d %H:%M", std::localtime(&tt));
            s.when = buf;
            saves.push_back(s);
        }
        std::sort(saves.begin(), saves.end(), [](const SaveEntry& a, const SaveEntry& b) { return a.time > b.time; });
        saveScroll = 0.0f;
    }

    // A row of choice buttons: label, then one button per option (the chosen one in brackets). True when one is picked.
    bool Choice(Rectangle r, const char* label, const char* const* names, int n, int& v) {
        const int fs = (int)(r.height * 0.5f);
        DrawText(label, (int)r.x, (int)(r.y + (r.height - fs) * 0.5f), fs, kText);
        bool picked = false;
        for (int k = 0; k < n; ++k) {
            const Rectangle b = { r.x + r.width * 0.38f + k * (r.height * 2.4f), r.y + r.height * 0.1f, r.height * 2.2f, r.height * 0.8f };
            if (Button(b, v == k ? TextFormat("[%s]", names[k]) : names[k])) { v = k; picked = true; }
        }
        return picked;
    }

    // The settings page (shared by both menus), in tabs (#512: one long list ran its Back button off a 720p screen).
    // Back is pinned to the bottom edge; if a tab's rows still would not fit above it, they shrink. Returns true on Back.
    enum class SettingsTab { Graphics, Controls, Keys, Audio, Gameplay, Count };
    SettingsTab settingsTab = SettingsTab::Graphics;
    float settingsContentBottom = 0.0f; Rectangle settingsBackRect{}; // the last frame's layout (--shot-settings checks it)
    bool SettingsPage() {
        const float s = Scale();
        const float sh = (float)GetScreenHeight(), sw = (float)GetScreenWidth();
        const float backH = 58 * s, backY = sh - backH - 28 * s;
        // Tabs under the menu's title (the main menu's at 18%, the pause menu's at 20% of the height).
        const char* tabNames[] = { "Graphics", "Controls", "Keys", "Audio", "Gameplay" };
        const float tabW = 200 * s, tabH = 46 * s, tabGap = 10 * s, tabY = sh * 0.30f;
        const float tabsX = (sw - ((int)SettingsTab::Count * tabW + ((int)SettingsTab::Count - 1) * tabGap)) * 0.5f;
        for (int t = 0; t < (int)SettingsTab::Count; ++t) {
            const Rectangle r = { tabsX + t * (tabW + tabGap), tabY, tabW, tabH };
            const bool on = (int)settingsTab == t;
            if (on) DrawRectangleRounded({ r.x - 3, r.y - 3, r.width + 6, r.height + 6 }, 0.2f, 6, Fade(SKYBLUE, 0.7f));
            if (Button(r, tabNames[t])) settingsTab = (SettingsTab)t;
        }
        // The tallest tab holds 5 rows and a note: shrink the rows if that would reach the Back button.
        const float top = tabY + tabH + 24 * s;
        const float h = std::min(52.0f * s, (backY - 16 * s - top) / (5 * 1.12f + 0.5f)), gap = h * 0.12f;
        const float w = 760.0f * s, x = (sw - w) * 0.5f;
        float y = top;
        bool changed = false;
        auto next = [&] { y += h + gap; };
        switch (settingsTab) {
        case SettingsTab::Graphics: {
            changed |= Slider({ x, y, w, h }, "Field of view", settings.fov, 30.0f, 100.0f, "%.0f deg"); next();
            // Terrain detail (#297): how far out the terrain keeps its finer levels (the LOD factor). Flying drops to Low.
            // (#1610) "High (lighter)": LOD factor 2.75, ~5% faster than High (GUESS name; for the playtest).
            static const char* detail[] = { "Low", "Medium", "High (lighter)", "High", "Ultra" };
            changed |= Choice({ x, y, w, h }, "Terrain detail", detail, 5, settings.terrainDetail); next();
            changed |= Toggle({ x, y, w, h }, "Foliage", settings.foliage); next();
            changed |= Toggle({ x, y, w, h }, "VSync (no tearing)", settings.vsync); next();
            changed |= Toggle({ x, y, w, h }, "Show FPS", settings.showFps); next(); // on every screen (#572)
            static const int caps[] = { 0, 60, 144, 240 };
            static const char* capNames[] = { "Monitor", "60", "144", "240" };
            int cap = 0;
            for (int k = 0; k < 4; ++k) if (settings.fpsCap == caps[k]) cap = k;
            const int capBefore = settings.fpsCap;
            if (Choice({ x, y, w, h }, "Frame cap", capNames, 4, cap)) { settings.fpsCap = caps[cap]; changed |= capBefore != caps[cap]; }
            next();
            // #131: how bright light rolls off - Clip (hard white), Knee (today's look; highlights roll off), ACES (filmic).
            static const char* curves[] = { "Clip", "Knee", "ACES" };
            changed |= Choice({ x, y, w, h }, "Tonemap", curves, 3, settings.tonemap); next();
            changed |= Slider({ x, y, w, h }, "Brightness", settings.brightness, -1.0f, 1.0f, "%+.1f EV"); next(); // (#1166)
            // #302: Low also turns anti-aliasing (MSAA) off to save video memory; that takes a restart.
            static const int startDetail = settings.terrainDetail;
            if ((startDetail == 0) != (settings.terrainDetail == 0))
                DrawText("Anti-aliasing changes after a restart (Low: off)", (int)x, (int)y, (int)(h * 0.4f), Fade(kText, 0.7f));
            break;
        }
        case SettingsTab::Controls: {
            changed |= Slider({ x, y, w, h }, "Mouse sensitivity", settings.mouseSensitivity, 0.2f, 3.0f, "%.2fx"); next();
            // Camera smoothing (#168): how much the eye's small bumps (collision push-out on rough ground) are filtered.
            static const char* smoothing[] = { "Off", "Low", "Default", "High" };
            changed |= Choice({ x, y, w, h }, "Camera smoothing", smoothing, 4, settings.cameraSmoothing); next();
            // Interface size (#821): the game window and panels; Auto follows the screen's height (1080 p is 100%).
            static const float sizes[] = { 0.0f, 0.75f, 1.0f, 1.25f, 1.5f };
            static const char* sizeNames[] = { "Auto", "75%", "100%", "125%", "150%" };
            int sz = 0;
            for (int k = 0; k < 5; ++k) if (fabsf(settings.uiScale - sizes[k]) < 0.01f) sz = k;
            if (Choice({ x, y, w, h }, "Interface size", sizeNames, 5, sz)) { settings.uiScale = sizes[sz]; changed = true; }
            next();
            break;
        }
        case SettingsTab::Keys: {
            // Key bindings (#821): every action; click one, then press its new key or mouse button (Esc keeps the old one).
            static int waiting = -1;   // the action waiting for a key
            static int waitFrames = 0; // frames it has waited (the click that started it is no binding)
            static float scroll = 0.0f;
            static_assert((int)Act::MagicSight == (int)Act::Count - 1, "magic sight is the last action (hidden while magic is off)");
            const int n = (int)Act::Count - (Magic_Enabled() ? 0 : 1); // (#1249) no key for a feature that isn't there
            const float rowH = std::min(40.0f * s, h * 0.8f), colGap = 24 * s, colW = (std::min(sw - 64 * s, 1500.0f * s) - colGap) * 0.5f;
            const float x0 = (sw - (2 * colW + colGap)) * 0.5f, listTop = y, listBottom = backY - 70 * s;
            const int rowsPerCol = (n + 1) / 2;
            const float content = rowsPerCol * (rowH + 4 * s), view = listBottom - listTop;
            if (waiting < 0) scroll = std::clamp(scroll - GetMouseWheelMove() * rowH * 2.0f, 0.0f, std::max(0.0f, content - view));
            BeginScissorMode((int)x0, (int)listTop, (int)(2 * colW + colGap), (int)view);
            for (int a = 0; a < n; ++a) {
                const int col = a / rowsPerCol, row = a % rowsPerCol;
                const float rx = x0 + col * (colW + colGap), ry = listTop + row * (rowH + 4 * s) - scroll;
                if (ry + rowH < listTop || ry > listBottom) continue;
                const int fs = (int)(rowH * 0.45f);
                DrawText(PlayerInput_ActLabel((Act)a), (int)rx, (int)(ry + (rowH - fs) * 0.5f), fs, kText);
                const Rectangle kb = { rx + colW * 0.62f, ry, colW * 0.38f, rowH };
                const bool isWaiting = waiting == a;
                const int code = PlayerInput_Binding((Act)a);
                const bool changedKey = code != PlayerInput_DefaultBinding((Act)a);
                if (Button(kb, isWaiting ? "press a key..." : changedKey ? TextFormat("%s *", PlayerInput_KeyName(code)) : PlayerInput_KeyName(code)) && waiting < 0) {
                    waiting = a; waitFrames = 0;
                }
            }
            EndScissorMode();
            const Rectangle reset = { (sw - 360 * s) * 0.5f, listBottom + 12 * s, 360 * s, 44 * s };
            if (waiting < 0 && Button(reset, "Reset all to defaults")) { PlayerInput_ResetBindings(); changed = true; }
            y = listBottom;
            if (waiting >= 0) {
                // Waiting for the key: Esc keeps the old one (and stays on the page); the starting click is skipped.
                bool done = false;
                if (IsKeyPressed(KEY_ESCAPE)) done = true;
                else if (const int code = PlayerInput_PressedCode(); code >= 0 && waitFrames > 0) { PlayerInput_SetBinding((Act)waiting, code); changed = true; done = true; }
                ++waitFrames;
                if (done) waiting = -1;
                settingsContentBottom = y;
                settingsBackRect = { (sw - 300 * s) * 0.5f, backY, 300 * s, backH };
                if (changed) { Settings_Apply(); Settings_Save(); }
                DrawText("Press the new key or mouse button (Esc keeps the old one)", (int)x0, (int)(listBottom + 20 * s), (int)(22 * s), Fade(kText, 0.8f));
                return false;
            }
            break;
        }
        case SettingsTab::Audio:
            changed |= Slider({ x, y, w, h }, Audio_V2() ? "Master volume" : "Volume", settings.volume, 0.0f, 1.0f, "%.2f"); next();
            if (Audio_V2()) // (#1238) the buses; the old engine has only the master volume
                for (int b = 1; b < (int)AudioBus::Count; ++b) { changed |= Slider({ x, y, w, h }, AudioBus_Name((AudioBus)b), settings.busVolume[b], 0.0f, 1.0f, "%.2f"); next(); }
            if (Audio_V2()) { // (#1244) captions for sounds, with their direction
                static const char* const kSub[] = { "Off", "Important", "All" };
                changed |= Choice({ x, y, w, h }, "Subtitles", kSub, 3, settings.subtitles); next();
            }
            break;
        case SettingsTab::Gameplay: {
            // (Difficulty is chosen on the New World screen and kept by the world, #526: no setting here.)
            // Comfort (#168): motion-sickness options.
            changed |= Toggle({ x, y, w, h }, "Reduce motion", settings.reduceMotion); next();
            changed |= Toggle({ x, y, w, h }, "Show hands", settings.showHands); next();
            static const char* bob[] = { "Full", "Reduced", "Off" };
            changed |= Choice({ x, y, w, h }, "Hand bob", bob, 3, settings.handBob); next();
            changed |= Toggle({ x, y, w, h }, "Speed vignette", settings.speedVignette); next();
            changed |= Toggle({ x, y, w, h }, "Firelight flicker", settings.fireFlicker); next(); // torches and campfires (#558)
            changed |= Toggle({ x, y, w, h }, "Double-tap to dodge", settings.doubleTapDodge); next(); // (#1218)
            changed |= Toggle({ x, y, w, h }, "Reduce flashing", settings.reduceFlashing); next(); // lightning (#1162)
            static const char* vision[] = { "Normal", "Red-green", "Blue-yellow" }; // (#1432) colour-blind-safe meaning colours
            changed |= Choice({ x, y, w, h }, "Colour vision", vision, 3, settings.colorVision); next();
            static const char* hints[] = { "Auto", "Always", "Off" }; // (#827) Auto: fades after your first 2 hours of play
            changed |= Choice({ x, y, w, h }, "Key hints", hints, 3, settings.keyHints); next();
            changed |= Toggle({ x, y, w, h }, "Inventory scroll: up sends", settings.invertSlotWheel); next(); // (#862) off: down sends to the other side, up pulls in (Mouse Tweaks)
            break;
        }
        case SettingsTab::Count: break;
        }
        settingsContentBottom = y;
        settingsBackRect = { (sw - 300 * s) * 0.5f, backY, 300 * s, backH };
        if (changed) Settings_Apply();
        const bool back = Button(settingsBackRect, "Back") || IsKeyPressed(KEY_ESCAPE);
        if (back) Settings_Save();
        return back;
    }
}

void Settings_Load() {
    FILE* f = fopen("settings.txt", "r");
    if (f) {
        char key[64], val[64];
        while (fscanf(f, " %63s %63s", key, val) == 2) { // values are numbers, except log_level's word (#514)
            const float v = strtof(val, nullptr);
            if (!strcmp(key, "log_level")) { settings.logLevel = val; continue; }
            // Multiplayer (#341): words (a name, addresses), the suit colour in hex.
            if (!strcmp(key, "player_uuid")) { settings.playerUuid = val; continue; } // (#346)
            if (!strncmp(key, "bind_", 5)) { Act a; if (PlayerInput_ActFromId(key + 5, a)) PlayerInput_SetBinding(a, atoi(val)); continue; } // (#821) a key binding
            if (!strcmp(key, "uiscale")) { settings.uiScale = std::clamp(v, 0.0f, 1.5f); continue; } // (#821) 0 auto
            if (!strncmp(key, "tabwin_", 7) && key[7] >= '0' && key[7] <= '3' && !key[8]) { settings.tabWin[key[7] - '0'] = std::clamp(v, 0.0f, 1.0f); continue; } // (#1080)
            if (!strcmp(key, "subtitles")) { settings.subtitles = std::clamp((int)v, 0, 2); continue; } // (#1244)
            if (!strncmp(key, "volume_", 7)) { // (#1238) a bus volume
                for (int b = 1; b < (int)AudioBus::Count; ++b) if (!strcmp(key + 7, AudioBus_Key((AudioBus)b))) settings.busVolume[b] = std::clamp(v, 0.0f, 1.0f);
                continue;
            }
            if (!strcmp(key, "mp_name")) { settings.mpName = val; continue; }
            if (!strcmp(key, "mp_address")) { settings.mpAddress = val; continue; }
            if (!strcmp(key, "mp_recent")) { if (settings.mpRecent.size() < 5) settings.mpRecent.push_back(val); continue; }
            if (!strcmp(key, "mp_suit")) { settings.mpSuit = (uint32_t)strtoul(val, nullptr, 16) & 0xFFFFFFu; continue; }
            if (!strcmp(key, "mp_port")) { settings.mpPort = std::clamp((int)v, 1024, 65535); continue; }
            if (!strcmp(key, "mp_maxplayers")) { settings.mpMaxPlayers = std::clamp((int)v, 1, 64); continue; }
            if (!strcmp(key, "mouse")) settings.mouseSensitivity = std::clamp(v, 0.2f, 3.0f);
            else if (!strcmp(key, "fov")) settings.fov = std::clamp(v, 30.0f, 100.0f);
            else if (!strcmp(key, "volume")) settings.volume = std::clamp(v, 0.0f, 1.0f);
            else if (!strcmp(key, "foliage")) settings.foliage = v != 0.0f;
            else if (!strcmp(key, "fpscap")) settings.fpsCap = std::max(0, (int)v);
            else if (!strcmp(key, "vsync")) settings.vsync = v != 0.0f;
            else if (!strcmp(key, "showfps")) settings.showFps = v != 0.0f;
            else if (!strcmp(key, "reducemotion")) settings.reduceMotion = v != 0.0f;
            else if (!strcmp(key, "showhands")) settings.showHands = v != 0.0f;
            else if (!strcmp(key, "handbob")) settings.handBob = std::clamp((int)v, 0, 2);
            else if (!strcmp(key, "difficulty")) settings.difficulty = std::clamp((int)v, 0, 2);
            else if (!strcmp(key, "xprate")) settings.xpRate = std::clamp((float)v, 0.5f, 5.0f); // (#883)
            else if (!strcmp(key, "terraindetail")) { static const int kOld[] = { 0, 1, 3, 4 }; settings.terrainDetail = kOld[std::clamp((int)v, 0, 3)]; } // (the four-level numbering, before #1610)
            else if (!strcmp(key, "terraindetail5")) settings.terrainDetail = std::clamp((int)v, 0, 4); // (#1610) five levels (after the old key: it wins)
            else if (!strcmp(key, "camerasmoothing")) settings.cameraSmoothing = std::clamp((int)v, 0, 3);
            else if (!strcmp(key, "tonemap")) settings.tonemap = std::clamp((int)v, 0, 2); // (#131)
            else if (!strcmp(key, "brightness")) settings.brightness = std::clamp(v, -1.0f, 1.0f); // (#1166)
            else if (!strcmp(key, "speedvignette")) settings.speedVignette = v != 0.0f;
            else if (!strcmp(key, "fireflicker")) settings.fireFlicker = v != 0.0f;
            else if (!strcmp(key, "reduceflashing")) settings.reduceFlashing = v != 0.0f;
            else if (!strcmp(key, "colorvision")) settings.colorVision = std::clamp((int)v, 0, 2); // (#1432)
            else if (!strcmp(key, "keyhints")) settings.keyHints = std::clamp((int)v, 0, 2); // (#827)
            else if (!strcmp(key, "invertslotwheel")) settings.invertSlotWheel = v != 0; // (#862)
            else if (!strcmp(key, "doubletapdodge")) settings.doubleTapDodge = v != 0; // (#1218)
        }
        fclose(f);
    }
    if (!PlayerId_Valid(settings.playerUuid)) { settings.playerUuid = PlayerId_New(); Settings_Save(); } // (#346) made on the first start
    Settings_Apply();
}

void Settings_Save() {
    FILE* f = fopen("settings.txt", "w");
    if (!f) return;
    fprintf(f, "mouse %.3f\nfov %.1f\nvolume %.3f\nfoliage %d\nfpscap %d\nvsync %d\nreducemotion %d\n", settings.mouseSensitivity, settings.fov, settings.volume,
            settings.foliage ? 1 : 0, settings.fpsCap, settings.vsync ? 1 : 0, settings.reduceMotion ? 1 : 0);
    for (int b = 1; b < (int)AudioBus::Count; ++b) fprintf(f, "volume_%s %.3f\n", AudioBus_Key((AudioBus)b), settings.busVolume[b]); // (#1238)
    fprintf(f, "subtitles %d\n", settings.subtitles); // (#1244)
    fprintf(f, "showhands %d\nhandbob %d\n", settings.showHands ? 1 : 0, settings.handBob);
    fprintf(f, "difficulty %d\n", settings.difficulty);
    fprintf(f, "xprate %g\n", settings.xpRate); // (#883)
    fprintf(f, "showfps %d\n", settings.showFps ? 1 : 0);
    { static const int kOld[] = { 0, 1, 2, 2, 3 }; // (#1610) older builds read the four levels (High (lighter) -> High)
      fprintf(f, "terraindetail %d\ncamerasmoothing %d\nspeedvignette %d\n", kOld[std::clamp(settings.terrainDetail, 0, 4)], settings.cameraSmoothing, settings.speedVignette ? 1 : 0); } // (terrain detail was never saved, #297)
    fprintf(f, "terraindetail5 %d\n", settings.terrainDetail); // (#1610) the five levels; after terraindetail, so it wins on load
    fprintf(f, "tonemap %d\n", settings.tonemap); // (#131)
    fprintf(f, "brightness %.2f\n", settings.brightness); // (#1166)
    fprintf(f, "fireflicker %d\n", settings.fireFlicker ? 1 : 0);
    fprintf(f, "reduceflashing %d\n", settings.reduceFlashing ? 1 : 0);
    fprintf(f, "colorvision %d\n", settings.colorVision); // (#1432)
    fprintf(f, "uiscale %.2f\n", settings.uiScale); // (#821)
    for (int k = 0; k < 4; ++k) fprintf(f, "tabwin_%d %.4f\n", k, settings.tabWin[k]); // (#1080) the Tab window's place and size
    fprintf(f, "keyhints %d\n", settings.keyHints); // (#827)
    fprintf(f, "invertslotwheel %d\n", settings.invertSlotWheel ? 1 : 0); // (#862)
    fprintf(f, "doubletapdodge %d\n", settings.doubleTapDodge ? 1 : 0); // (#1218)
    for (int a = 0; a < (int)Act::Count; ++a) // (#821) only the key bindings changed from their defaults
        if (PlayerInput_Binding((Act)a) != PlayerInput_DefaultBinding((Act)a)) fprintf(f, "bind_%s %d\n", PlayerInput_ActId((Act)a), PlayerInput_Binding((Act)a));
    if (!settings.logLevel.empty()) fprintf(f, "log_level %s\n", settings.logLevel.c_str()); // set by hand in the file (#514)
    // Multiplayer (#341; never the password).
    if (!settings.mpName.empty()) fprintf(f, "mp_name %s\n", settings.mpName.c_str());
    if (!settings.mpAddress.empty()) fprintf(f, "mp_address %s\n", settings.mpAddress.c_str());
    for (const std::string& a : settings.mpRecent) fprintf(f, "mp_recent %s\n", a.c_str());
    fprintf(f, "mp_suit %06x\nmp_port %d\nmp_maxplayers %d\n", settings.mpSuit, settings.mpPort, settings.mpMaxPlayers);
    if (!settings.playerUuid.empty()) fprintf(f, "player_uuid %s\n", settings.playerUuid.c_str()); // (#346)
    fclose(f);
}

void Settings_Apply() {
    player.mouseSensitivity = kDefaultMouse * settings.mouseSensitivity;
    ship.mouseSensitivity = kDefaultShipMouse * settings.mouseSensitivity;
    telescope.baseFov = settings.fov; // the camera's fov is set by the game each frame it is not zoomed
    Audio_SetMasterVolume(settings.volume);
    static_assert(sizeof(settings.busVolume) / sizeof(settings.busVolume[0]) == (int)AudioBus::Count, "one volume per audio bus");
    for (int b = 1; b < (int)AudioBus::Count; ++b) Audio_SetBusVolume((AudioBus)b, settings.busVolume[b]); // (#1238)
    Subtitles_SetLevel((SubtitleLevel)std::clamp(settings.subtitles, 0, 2)); // (#1244)
    UiPalette_Set((ColorVision)std::clamp(settings.colorVision, 0, 2)); // (#1432)
    foliageEnabled = settings.foliage;
    { const float kDetail[] = { 2.0f, 2.5f, 2.75f, 3.0f, 4.0f }; lodSettings.detailFactor = kDetail[std::clamp(settings.terrainDetail, 0, 4)]; } // #297; (#1610) 2.75: High (lighter)
    lodSettings.preloadRadius = settings.terrainDetail == 0 ? 150.0f : 300.0f; // #302 (user 2026-11-15): Low keeps a smaller spawn ring (VRAM)
    if (Headless()) { FramePacer_SetTarget(kSimHz); SetFrameTarget(kSimHz); } // #326: a server (hidden window) ticks at the sim rate: vsync never throttles a hidden window, so "uncapped" spun a core
    else {
    if (settings.vsync) SetWindowState(FLAG_VSYNC_HINT); else ClearWindowState(FLAG_VSYNC_HINT);
    { const int fps = settings.fpsCap > 0 ? settings.fpsCap : (settings.vsync ? 0 : std::max(60, GetMonitorRefreshRate(GetCurrentMonitor())));
      FramePacer_SetTarget(fps); SetFrameTarget(fps, settings.vsync ? GetMonitorRefreshRate(GetCurrentMonitor()) : 0); } // (#545) the budget's target
    }
    { static const int kCurve[] = { 2, 0, 1 }; SceneTarget_SetTonemap(kCurve[std::clamp(settings.tonemap, 0, 2)]); } // (#131) Clip, Knee, ACES
    SceneTarget_SetBrightness(settings.brightness); // (#1166)
    player.smoothMotion = settings.reduceMotion;
    gDoubleTapDodge = settings.doubleTapDodge; // (#1218)
    { const float kScale[] = { 1.0f, 0.5f, 1.0f, 2.4f }; const int k = std::clamp(settings.cameraSmoothing, 0, 3); player.eyeSmoothing = k > 0; player.eyeSmoothScale = kScale[k]; } // #168
    player.speedVignette = settings.speedVignette;
}

bool Menu_IsWorldSave(const std::string& path) {
    const fs::path p(path);
    std::error_code ec;
    const std::string name = p.filename().string();
    const bool benchSave = name == "bench.sav" || name.rfind("bench-", 0) == 0; // (#1387) bench-<seed>.sav too
    return p.extension() == ".sav" && !benchSave && !fs::exists(path + ".synced", ec);
}

std::string Menu_NewSavePath() {
    std::error_code ec;
    fs::create_directories("saves", ec);
    if (!fs::exists("saves/world.sav", ec)) return "saves/world.sav";
    for (int n = 2;; ++n) {
        const std::string p = "saves/world_" + std::to_string(n) + ".sav";
        if (!fs::exists(p, ec)) return p;
    }
}

extern const char* kGameVersion; // "<branch> <commit>", generated at build time (#501)

// ---- Multiplayer (#341): host a game (others join while you play), join by IP ----
namespace {
    // Suit colours (0xRRGGBB): the default first, then ten that read apart from each other on a grey suit and at a distance.
    constexpr uint32_t kSuits[] = { 0x000000, 0xD8452F, 0xE8902A, 0xE3C93B, 0x58B548, 0x2FA39A, 0x3E86D8, 0x6B5BD6, 0xB04FC2, 0xE06FA4, 0xEDEDED };
    constexpr int kSuitCount = (int)(sizeof(kSuits) / sizeof(kSuits[0]));
    int focus = -1;                       // the text field being typed in
    std::string mpPassword, mpMessage;    // (the password is never saved)
    bool mpMessageBad = false;
    int hostWorld = -1;                   // the save to host: an index into `saves`, -1 a new world
    std::string portText, maxText;

    // One line of text: click to type; `allowed` the characters it takes.
    void TextField(Rectangle r, const char* label, std::string& value, int id, size_t maxLen, const char* allowed, bool hidden = false) {
        const float s = Scale();
        const int fs = (int)(r.height * 0.45f);
        DrawText(label, (int)r.x, (int)(r.y + (r.height - fs) * 0.5f), fs, kText);
        const Rectangle box = { r.x + r.width * 0.34f, r.y, r.width * 0.66f, r.height };
        const bool hover = CheckCollisionPointRec(GetMousePosition(), box);
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) { if (hover) focus = id; else if (focus == id) focus = -1; }
        DrawRectangleRounded(box, 0.18f, 6, focus == id ? kHover : kButton);
        if (focus == id) {
            DrawRectangleRoundedLinesEx(box, 0.18f, 6, 2.0f, Fade(SKYBLUE, 0.8f));
            for (int c = GetCharPressed(); c > 0; c = GetCharPressed())
                if (c < 128 && strchr(allowed, c) && value.size() < maxLen) value += (char)c;
            if ((IsKeyPressed(KEY_BACKSPACE) || IsKeyPressedRepeat(KEY_BACKSPACE)) && !value.empty()) value.pop_back();
        }
        const std::string shown = hidden ? std::string(value.size(), '*') : value;
        DrawText(TextFormat("%s%s", shown.c_str(), focus == id && fmodf((float)GetTime(), 1.0f) < 0.5f ? "_" : ""), (int)(box.x + 12 * s), (int)(box.y + (box.height - fs) * 0.5f), fs, kText);
    }
    // The suit colour swatches; returns true when one was picked.
    bool SuitPicker(Rectangle r, uint32_t& suit) {
        const int fs = (int)(r.height * 0.45f);
        DrawText("Suit colour", (int)r.x, (int)(r.y + (r.height - fs) * 0.5f), fs, kText);
        const float x0 = r.x + r.width * 0.34f, w = r.width * 0.66f, gap = 6.0f * Scale();
        const float side = std::min(r.height, (w - gap * (kSuitCount - 1)) / kSuitCount);
        bool picked = false;
        for (int i = 0; i < kSuitCount; ++i) {
            const Rectangle b = { x0 + i * (side + gap), r.y + (r.height - side) * 0.5f, side, side };
            const uint32_t c = kSuits[i];
            const Color col = c ? Color{ (unsigned char)(c >> 16), (unsigned char)(c >> 8), (unsigned char)c, 255 } : Color{ 150, 154, 160, 255 }; // (the default suit's grey)
            DrawRectangleRounded(b, 0.25f, 6, col);
            if (c == 0) DrawText("-", (int)(b.x + side * 0.38f), (int)(b.y + side * 0.18f), (int)(side * 0.6f), Fade(BLACK, 0.6f));
            const bool hover = CheckCollisionPointRec(GetMousePosition(), b);
            if (suit == c) DrawRectangleRoundedLinesEx(b, 0.25f, 6, 3.0f, WHITE);
            else if (hover) DrawRectangleRoundedLinesEx(b, 0.25f, 6, 2.0f, Fade(SKYBLUE, 0.8f));
            if (hover && IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) { suit = c; picked = true; }
        }
        return picked;
    }
    void Message(float y) {
        if (mpMessage.empty()) return;
        const int fs = (int)(24 * Scale());
        DrawText(mpMessage.c_str(), (GetScreenWidth() - MeasureText(mpMessage.c_str(), fs)) / 2, (int)y, fs, mpMessageBad ? Color{ 255, 140, 110, 255 } : Fade(kText, 0.85f));
    }
    const char* kNameChars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-";
    const char* kAddrChars = "0123456789.:abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ-";
    const char* kDigits = "0123456789";
    const char* kPassChars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-.!@#$%^&*+=?";

    void MultiplayerPage(MenuResult& r) {
        const float s = Scale();
        const float w = 760.0f;
        auto row = [&](int i) { return Row(i, w, 0.30f); };
        if (mainPage == Page::Multiplayer) {
            if (Button(Row(0), "Host game")) { RefreshSaves(); hostWorld = saves.empty() ? -1 : 0; portText = std::to_string(settings.mpPort); maxText = std::to_string(settings.mpMaxPlayers); mpMessage.clear(); mainPage = Page::Host; }
            if (Button(Row(1), "Join game")) { mpMessage.clear(); mainPage = Page::Join; }
            const int fs = (int)(20 * s);
            // Hosting here is not a dedicated server (#806): the world lives in this game and closes with it.
            const char* note = "Host game: others can join while you play; the world closes when you leave.";
            const char* note2 = "Players join by IP address; outside your network they need your port forwarded.";
            DrawText(note, (GetScreenWidth() - MeasureText(note, fs)) / 2, (int)(Row(2).y + 2 * s), fs, Fade(kText, 0.8f));
            DrawText(note2, (GetScreenWidth() - MeasureText(note2, fs)) / 2, (int)(Row(2).y + 2 * s + fs * 1.4f), fs, Fade(kText, 0.6f));
            if (Button(Row(3), "Back") || IsKeyPressed(KEY_ESCAPE)) mainPage = Page::Main;
            return;
        }
        const bool join = mainPage == Page::Join;
        int i = 0;
        TextField(row(i++), "Your name", settings.mpName, 1, 24, kNameChars);
        if (join) {
            TextField(row(i++), "Address", settings.mpAddress, 2, 64, kAddrChars);
            TextField(row(i++), "Password", mpPassword, 3, 32, kPassChars, true);
        } else {
            { // the world to host: the saves, or a new one
                const Rectangle rw = row(i++);
                const int fs = (int)(rw.height * 0.45f);
                DrawText("World", (int)rw.x, (int)(rw.y + (rw.height - fs) * 0.5f), fs, kText);
                const float bx = rw.x + rw.width * 0.34f, bw = rw.width * 0.66f, side = rw.height;
                if (Button({ bx, rw.y, side, side }, "<")) hostWorld = hostWorld < 0 ? (int)saves.size() - 1 : hostWorld - 1;
                if (Button({ bx + bw - side, rw.y, side, side }, ">")) hostWorld = hostWorld + 1 >= (int)saves.size() ? -1 : hostWorld + 1;
                const char* name = hostWorld >= 0 && hostWorld < (int)saves.size() ? saves[hostWorld].name.c_str() : "A new world";
                DrawText(name, (int)(bx + side + 12 * s), (int)(rw.y + (rw.height - fs) * 0.5f), fs, kText);
            }
            TextField(row(i++), "Port", portText, 4, 5, kDigits);
            TextField(row(i++), "Max players", maxText, 5, 2, kDigits);
            TextField(row(i++), "Password (optional)", mpPassword, 3, 32, kPassChars, true);
        }
        SuitPicker(row(i++), settings.mpSuit);
        if (join && !settings.mpRecent.empty()) { // recent addresses: click to fill it in
            const Rectangle rr = row(i++);
            const int fs = (int)(rr.height * 0.4f);
            DrawText("Recent", (int)rr.x, (int)(rr.y + (rr.height - fs) * 0.5f), fs, Fade(kText, 0.8f));
            float x = rr.x + rr.width * 0.34f;
            for (size_t k = 0; k < settings.mpRecent.size() && k < 3; ++k) {
                const std::string& a = settings.mpRecent[k];
                const float bw = std::min(rr.width * 0.66f / 3.0f - 8 * s, (float)MeasureText(a.c_str(), fs) + 24 * s);
                if (Button({ x, rr.y + rr.height * 0.1f, bw, rr.height * 0.8f }, a.c_str())) settings.mpAddress = a;
                x += bw + 8 * s;
            }
        }
        const Rectangle go = row(i + 1);
        Message(row(i).y + 8 * s);
        const bool nameOk = !settings.mpName.empty();
        const Rectangle bGo = { go.x, go.y, go.width * 0.49f, go.height }, bBack = { go.x + go.width * 0.51f, go.y, go.width * 0.49f, go.height };
        if (join) {
            const bool ok = nameOk && !settings.mpAddress.empty();
            if (Button(bGo, "Join", ok) || (ok && focus >= 0 && IsKeyPressed(KEY_ENTER))) {
                r.action = MenuAction::JoinGame; r.address = settings.mpAddress; r.name = settings.mpName; r.password = mpPassword; r.suitColor = settings.mpSuit;
                mpMessage = TextFormat("Joining %s ...", settings.mpAddress.c_str()); mpMessageBad = false; focus = -1;
                Settings_Save();
            }
        } else {
            const int port = atoi(portText.c_str()), maxp = atoi(maxText.c_str());
            const bool ok = nameOk && port >= 1024 && port <= 65535 && maxp >= 1 && maxp <= 64;
            if (!ok && !portText.empty() && !maxText.empty()) { mpMessage = "Port 1024-65535, 1-64 players, and a name"; mpMessageBad = true; }
            if (Button(bGo, "Host", ok)) {
                settings.mpPort = port; settings.mpMaxPlayers = maxp;
                r.action = MenuAction::HostGame; r.path = hostWorld >= 0 && hostWorld < (int)saves.size() ? saves[hostWorld].path : std::string();
                r.port = port; r.maxPlayers = maxp; r.password = mpPassword; r.name = settings.mpName; r.suitColor = settings.mpSuit;
                mpMessage = "Opening your world to others ..."; mpMessageBad = false; focus = -1;
                Settings_Save();
            }
        }
        if (Button(bBack, "Back") || (focus < 0 && IsKeyPressed(KEY_ESCAPE))) { mainPage = Page::Multiplayer; mpMessage.clear(); focus = -1; }
    }
}

void Menu_MultiplayerFailed(bool host, const std::string& why) { mainPage = host ? Page::Host : Page::Join; mpMessage = why; mpMessageBad = true; }
void Menu_RememberServer(const std::string& address) {
    auto& v = settings.mpRecent;
    v.erase(std::remove(v.begin(), v.end(), address), v.end());
    v.insert(v.begin(), address);
    if (v.size() > 5) v.resize(5);
    Settings_Save();
}
void Menu_TestMultiplayer(int page) { mainPage = page == 1 ? Page::Join : page == 2 ? Page::Host : Page::Multiplayer; if (page == 2) { RefreshSaves(); hostWorld = -1; portText = std::to_string(settings.mpPort); maxText = std::to_string(settings.mpMaxPlayers); } }
int Menu_SuitColorCount() { return kSuitCount; }
uint32_t Menu_SuitColor(int i) { return i >= 0 && i < kSuitCount ? kSuits[i] : 0; }

MenuResult Menu_Main() {
    MenuResult r;
    Backdrop("PLANETOID");
    const float s = Scale();
    { // #501: which build this is (branch + commit), bottom right, on every main-menu page
        const int vfs = (int)(18 * s), pad = (int)(12 * s);
        const int vx = GetScreenWidth() - MeasureText(kGameVersion, vfs) - pad, vy = GetScreenHeight() - vfs - pad;
        DrawText(kGameVersion, vx + 1, vy + 1, vfs, Fade(BLACK, 0.6f)); // shadow: readable over bright slideshow shots
        DrawText(kGameVersion, vx, vy, vfs, Fade(kText, 0.75f));
    }
    if (mainPage == Page::Main) {
        if (Button(Row(0), "New Game")) r.action = MenuAction::NewGame;
        if (Button(Row(1), "Load Game")) { RefreshSaves(); mainPage = Page::Load; }
        if (Button(Row(2), "Multiplayer")) mainPage = Page::Multiplayer; // (#341)
        if (Button(Row(3), "Settings")) mainPage = Page::Settings;
        if (Button(Row(4), "Exit")) r.action = MenuAction::Exit;
    } else if (mainPage == Page::Load) {
        const int fs = (int)(26 * s);
        if (saves.empty()) DrawText("No saves yet", (GetScreenWidth() - MeasureText("No saves yet", fs)) / 2, (int)(GetScreenHeight() * 0.38f), fs, Fade(kText, 0.7f));
        saveScroll = std::clamp(saveScroll - GetMouseWheelMove(), 0.0f, std::max(0.0f, (float)saves.size() - 7.0f));
        const int first = (int)saveScroll;
        bool changed = false;
        for (int i = first; i < (int)saves.size() && i < first + 7; ++i) {
            const SaveEntry& e = saves[i];
            // The row: the save (click to load), then Rename, Copy and Delete (#279).
            const Rectangle row = Row(i - first, 820.0f, 0.32f);
            const float sw = 92.0f * s, gap = 8.0f * s;
            const Rectangle main = { row.x, row.y, row.width - 3 * (sw + gap), row.height };
            const Rectangle bRen = { main.x + main.width + gap, row.y, sw, row.height }, bCopy = { bRen.x + sw + gap, row.y, sw, row.height }, bDel = { bCopy.x + sw + gap, row.y, sw, row.height };
            if (renaming == i) { // the name being typed
                DrawRectangleRounded(main, 0.18f, 6, kHover);
                const bool ok = SaveNameOk(renameText);
                const int tfs = (int)(main.height * 0.45f);
                DrawText(TextFormat("%s%s", renameText.c_str(), fmodf((float)GetTime(), 1.0f) < 0.5f ? "_" : ""), (int)(main.x + 14 * s), (int)(main.y + (main.height - tfs) * 0.5f), tfs, ok || renameText == e.name ? kText : ORANGE);
                for (int c = GetCharPressed(); c > 0; c = GetCharPressed()) if (c < 128 && renameText.size() < 40) renameText += (char)c;
                if (IsKeyPressed(KEY_BACKSPACE) && !renameText.empty()) renameText.pop_back();
                if ((IsKeyPressed(KEY_ENTER) || Button(bRen, "OK")) && ok) {
                    saveMessage = SaveRename(e.name, renameText) ? TextFormat("Renamed to %s", renameText.c_str()) : "Rename failed";
                    renaming = -1; changed = true;
                }
                if (Button(bCopy, "Cancel") || IsKeyPressed(KEY_ESCAPE)) renaming = -1;
                continue;
            }
            if (Button(main, TextFormat("%s   %s", e.name.c_str(), e.when.c_str()))) { r.action = MenuAction::LoadGame; r.path = e.path; mainPage = Page::Main; }
            if (Button(bRen, "Rename")) { renaming = i; renameText = e.name; deleting = -1; }
            if (Button(bCopy, "Copy")) {
                std::string to = e.name + " copy";
                for (int k = 2; !SaveNameOk(to) && k < 100; ++k) to = TextFormat("%s copy %d", e.name.c_str(), k);
                saveMessage = SaveCopy(e.name, to) ? TextFormat("Copied to %s", to.c_str()) : "Copy failed";
                changed = true;
            }
            const bool sure = deleting == i && GetTime() - deleteAsked < 4.0;
            if (Button(bDel, sure ? "Sure?" : "Delete")) {
                if (sure) { saveMessage = SaveDelete(e.name) ? TextFormat("Deleted %s", e.name.c_str()) : "Delete failed"; deleting = -1; changed = true; }
                else { deleting = i; deleteAsked = GetTime(); }
            }
        }
        if (!saveMessage.empty()) DrawText(saveMessage.c_str(), (GetScreenWidth() - MeasureText(saveMessage.c_str(), fs)) / 2, (int)(Row(7, 300.0f, 0.32f).y - 34 * s), fs, Fade(kText, 0.8f));
        if (changed) { const float keep = saveScroll; RefreshSaves(); saveScroll = std::min(keep, std::max(0.0f, (float)saves.size() - 7.0f)); }
        if (renaming < 0 && (Button(Row(7, 300.0f, 0.32f), "Back") || IsKeyPressed(KEY_ESCAPE))) { mainPage = Page::Main; deleting = -1; saveMessage.clear(); }
    } else if (mainPage == Page::Multiplayer || mainPage == Page::Join || mainPage == Page::Host) {
        MultiplayerPage(r);
    } else if (SettingsPage()) {
        mainPage = Page::Main;
    }
    return r;
}

void Menu_ResetPause() { pausePage = Page::Main; }
void Menu_TestSettings(int tab, bool pause) {
    settingsTab = (SettingsTab)std::clamp(tab, 0, (int)SettingsTab::Count - 1);
    if (pause) pausePage = Page::Settings; else mainPage = Page::Settings;
}
void Menu_SettingsLayout(float& contentBottom, Rectangle& back) { contentBottom = settingsContentBottom; back = settingsBackRect; }
bool Menu_TestSaves() { return TestSaveOps(); }

MenuResult Menu_Pause() {
    MenuResult r;
    DrawRectangle(0, 0, GetScreenWidth(), GetScreenHeight(), Fade(BLACK, 0.65f));
    const float s = Scale();
    const int fs = (int)(52 * s);
    DrawText("Paused", (GetScreenWidth() - MeasureText("Paused", fs)) / 2, (int)(GetScreenHeight() * 0.2f), fs, kText);
    if (pausePage == Page::Settings) {
        if (SettingsPage()) pausePage = Page::Main;
        return r;
    }
    const Rectangle panel = { Row(0).x - 30 * s, Row(0).y - 30 * s, Row(0).width + 60 * s, Row(4).y + Row(4).height - Row(0).y + 60 * s };
    DrawRectangleRounded(panel, 0.06f, 6, kPanel);
    if (Button(Row(0), "Resume") || IsKeyPressed(KEY_ESCAPE)) r.action = MenuAction::Resume;
    if (Button(Row(1), "Save")) r.action = MenuAction::Save;
    if (Button(Row(2), "Settings")) pausePage = Page::Settings;
    if (Button(Row(3), "Report a bug")) r.action = MenuAction::ReportBug; // #169
    if (const int n = BugReport_Unsent(); n > 0) { // (#617) beside the button: reports still waiting for the relay
        const char* t = TextFormat(n == 1 ? "1 unsent bug report" : "%d unsent bug reports", n);
        const int tfs = (int)(18 * s);
        DrawText(t, (int)(Row(3).x + Row(3).width + 16 * s), (int)(Row(3).y + (Row(3).height - tfs) * 0.5f), tfs, Color{ 255, 190, 90, 255 });
    }
    if (Button(Row(4), "Save & Quit to Menu")) r.action = MenuAction::SaveAndQuit;
    return r;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"
TEST(Menus, AJoinedGamesCacheIsNotAWorldToHost) { // (#1110: hosting the client's cache started a fresh world and player)
    const fs::path dir = fs::temp_directory_path() / "raylibgame_1110";
    std::error_code ec;
    fs::remove_all(dir, ec); fs::create_directories(dir, ec);
    auto touch = [&](const char* name) { if (FILE* f = fopen((dir / name).string().c_str(), "w")) fclose(f); };
    touch("world_25.sav"); touch("client_player.sav"); touch("client_player.sav.synced"); touch("bench.sav");
    EXPECT_TRUE(Menu_IsWorldSave((dir / "world_25.sav").string())) << "a world";
    EXPECT_FALSE(Menu_IsWorldSave((dir / "client_player.sav").string())) << "a joined game's cache (a .synced beside it)";
    EXPECT_FALSE(Menu_IsWorldSave((dir / "bench.sav").string())) << "the bench's";
    touch("bench-7.sav");
    EXPECT_FALSE(Menu_IsWorldSave((dir / "bench-7.sav").string())) << "(#1387) a bench save for another seed";
    EXPECT_FALSE(Menu_IsWorldSave((dir / "world_25.sav.synced").string())) << "not a .sav";
    fs::remove_all(dir, ec);
}
#endif
