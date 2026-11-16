#include "menus.h"
#include "slideshow.h"
#include "player.h"
#include "ship.h"
#include "telescope.h"
#include "audio.h"
#include "foliage.h"
#include "world.h"

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
    enum class Page { Main, Load, Settings };
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
            if (!e.is_regular_file() || e.path().extension() != ".sav" || e.path().filename() == "bench.sav") continue;
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

    // The settings page (shared by both menus). Returns true when Back is pressed.
    bool SettingsPage() {
        const float s = Scale();
        const float w = 760.0f * s, h = 52.0f * s;
        const float x = (GetScreenWidth() - w) * 0.5f;
        float y = GetScreenHeight() * 0.26f;
        bool changed = false;
        changed |= Slider({ x, y, w, h }, "Mouse sensitivity", settings.mouseSensitivity, 0.2f, 3.0f, "%.2fx"); y += h + 6 * s;
        changed |= Slider({ x, y, w, h }, "Field of view", settings.fov, 30.0f, 100.0f, "%.0f deg"); y += h + 6 * s;
        changed |= Slider({ x, y, w, h }, "Volume", settings.volume, 0.0f, 1.0f, "%.2f"); y += h + 6 * s;
        changed |= Toggle({ x, y, w, h }, "Foliage", settings.foliage); y += h + 6 * s;
        changed |= Toggle({ x, y, w, h }, "VSync (no tearing)", settings.vsync); y += h + 6 * s;
        changed |= Toggle({ x, y, w, h }, "Reduce motion", settings.reduceMotion); y += h + 6 * s;
        changed |= Toggle({ x, y, w, h }, "Show hands", settings.showHands); y += h + 6 * s;
        {
            const int fs = (int)(h * 0.5f);
            DrawText("Hand bob", (int)x, (int)(y + (h - fs) * 0.5f), fs, kText);
            const char* names[] = { "Full", "Reduced", "Off" };
            for (int k = 0; k < 3; ++k) {
                const Rectangle b = { x + w * 0.38f + k * (h * 2.4f), y + h * 0.1f, h * 2.2f, h * 0.8f };
                if (Button(b, settings.handBob == k ? TextFormat("[%s]", names[k]) : names[k])) { settings.handBob = k; changed = true; }
            }
            y += h + 6 * s;
        }
        {   // Camera smoothing (#168): how much the eye's small bumps (collision push-out on rough ground) are filtered.
            const int fs = (int)(h * 0.5f);
            DrawText("Camera smoothing", (int)x, (int)(y + (h - fs) * 0.5f), fs, kText);
            const char* names[] = { "Off", "Low", "Default", "High" };
            for (int k = 0; k < 4; ++k) {
                const Rectangle b = { x + w * 0.38f + k * (h * 2.4f), y + h * 0.1f, h * 2.2f, h * 0.8f };
                if (Button(b, settings.cameraSmoothing == k ? TextFormat("[%s]", names[k]) : names[k])) { settings.cameraSmoothing = k; changed = true; }
            }
            y += h + 6 * s;
        }
        changed |= Toggle({ x, y, w, h }, "Speed vignette", settings.speedVignette); y += h + 6 * s;
        {   // Terrain detail (#297): how far out the terrain keeps its finer levels (the LOD factor). Flying drops to Low.
            const int fs = (int)(h * 0.5f);
            DrawText("Terrain detail", (int)x, (int)(y + (h - fs) * 0.5f), fs, kText);
            const char* names[] = { "Low", "Medium", "High", "Ultra" };
            for (int k = 0; k < 4; ++k) {
                const Rectangle b = { x + w * 0.38f + k * (h * 2.4f), y + h * 0.1f, h * 2.2f, h * 0.8f };
                if (Button(b, settings.terrainDetail == k ? TextFormat("[%s]", names[k]) : names[k])) { settings.terrainDetail = k; changed = true; }
            }
            y += h + 6 * s;
            // #302: Low also turns anti-aliasing (MSAA) off to save video memory; that takes a restart.
            static const int startDetail = settings.terrainDetail;
            if ((startDetail == 0) != (settings.terrainDetail == 0)) {
                DrawText("Anti-aliasing changes after a restart (Low: off)", (int)x, (int)y, (int)(h * 0.4f), Fade(kText, 0.7f));
                y += h * 0.5f;
            }
        }
        {   // Difficulty for the next new world (each world keeps its own).
            const int fs = (int)(h * 0.5f);
            DrawText("New worlds", (int)x, (int)(y + (h - fs) * 0.5f), fs, kText);
            const char* names[] = { "Peaceful", "Normal", "Hard" };
            for (int k = 0; k < 3; ++k) {
                const Rectangle b = { x + w * 0.38f + k * (h * 2.4f), y + h * 0.1f, h * 2.2f, h * 0.8f };
                if (Button(b, settings.difficulty == k ? TextFormat("[%s]", names[k]) : names[k])) { settings.difficulty = k; changed = true; }
            }
            y += h + 6 * s;
        }
        {
            const int fs = (int)(h * 0.5f);
            DrawText("Frame cap", (int)x, (int)(y + (h - fs) * 0.5f), fs, kText);
            const int caps[] = { 0, 60, 144, 240 };
            const char* names[] = { "Monitor", "60", "144", "240" };
            for (int k = 0; k < 4; ++k) {
                const Rectangle b = { x + w * 0.38f + k * (h * 2.4f), y + h * 0.1f, h * 2.2f, h * 0.8f };
                if (Button(b, settings.fpsCap == caps[k] ? TextFormat("[%s]", names[k]) : names[k])) { settings.fpsCap = caps[k]; changed = true; }
            }
            y += h + 6 * s;
        }
        if (changed) Settings_Apply();
        const bool back = Button({ (GetScreenWidth() - 300 * s) * 0.5f, y + 30 * s, 300 * s, 58 * s }, "Back") || IsKeyPressed(KEY_ESCAPE);
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
            if (!strcmp(key, "mouse")) settings.mouseSensitivity = std::clamp(v, 0.2f, 3.0f);
            else if (!strcmp(key, "fov")) settings.fov = std::clamp(v, 30.0f, 100.0f);
            else if (!strcmp(key, "volume")) settings.volume = std::clamp(v, 0.0f, 1.0f);
            else if (!strcmp(key, "foliage")) settings.foliage = v != 0.0f;
            else if (!strcmp(key, "fpscap")) settings.fpsCap = std::max(0, (int)v);
            else if (!strcmp(key, "vsync")) settings.vsync = v != 0.0f;
            else if (!strcmp(key, "reducemotion")) settings.reduceMotion = v != 0.0f;
            else if (!strcmp(key, "showhands")) settings.showHands = v != 0.0f;
            else if (!strcmp(key, "handbob")) settings.handBob = std::clamp((int)v, 0, 2);
            else if (!strcmp(key, "difficulty")) settings.difficulty = std::clamp((int)v, 0, 2);
            else if (!strcmp(key, "terraindetail")) settings.terrainDetail = std::clamp((int)v, 0, 3);
            else if (!strcmp(key, "camerasmoothing")) settings.cameraSmoothing = std::clamp((int)v, 0, 3);
            else if (!strcmp(key, "speedvignette")) settings.speedVignette = v != 0.0f;
        }
        fclose(f);
    }
    Settings_Apply();
}

void Settings_Save() {
    FILE* f = fopen("settings.txt", "w");
    if (!f) return;
    fprintf(f, "mouse %.3f\nfov %.1f\nvolume %.3f\nfoliage %d\nfpscap %d\nvsync %d\nreducemotion %d\n", settings.mouseSensitivity, settings.fov, settings.volume,
            settings.foliage ? 1 : 0, settings.fpsCap, settings.vsync ? 1 : 0, settings.reduceMotion ? 1 : 0);
    fprintf(f, "showhands %d\nhandbob %d\n", settings.showHands ? 1 : 0, settings.handBob);
    fprintf(f, "difficulty %d\n", settings.difficulty);
    fprintf(f, "terraindetail %d\ncamerasmoothing %d\nspeedvignette %d\n", settings.terrainDetail, settings.cameraSmoothing, settings.speedVignette ? 1 : 0); // (terrain detail was never saved, #297)
    if (!settings.logLevel.empty()) fprintf(f, "log_level %s\n", settings.logLevel.c_str()); // set by hand in the file (#514)
    fclose(f);
}

void Settings_Apply() {
    player.mouseSensitivity = kDefaultMouse * settings.mouseSensitivity;
    ship.mouseSensitivity = kDefaultShipMouse * settings.mouseSensitivity;
    telescope.baseFov = settings.fov; // the camera's fov is set by the game each frame it is not zoomed
    Audio_SetMasterVolume(settings.volume);
    foliageEnabled = settings.foliage;
    { const float kDetail[] = { 2.0f, 2.5f, 3.0f, 4.0f }; lodSettings.detailFactor = kDetail[std::clamp(settings.terrainDetail, 0, 3)]; } // #297
    lodSettings.preloadRadius = settings.terrainDetail == 0 ? 150.0f : 300.0f; // #302 (user 2026-11-15): Low keeps a smaller spawn ring (VRAM)
    if (settings.vsync) SetWindowState(FLAG_VSYNC_HINT); else ClearWindowState(FLAG_VSYNC_HINT);
    SetTargetFPS(settings.fpsCap > 0 ? settings.fpsCap : (settings.vsync ? 0 : std::max(60, GetMonitorRefreshRate(GetCurrentMonitor()))));
    player.smoothMotion = settings.reduceMotion;
    { const float kScale[] = { 1.0f, 0.5f, 1.0f, 2.4f }; const int k = std::clamp(settings.cameraSmoothing, 0, 3); player.eyeSmoothing = k > 0; player.eyeSmoothScale = kScale[k]; } // #168
    player.speedVignette = settings.speedVignette;
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
        if (Button(Row(2), "Settings")) mainPage = Page::Settings;
        if (Button(Row(3), "Exit")) r.action = MenuAction::Exit;
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
    } else if (SettingsPage()) {
        mainPage = Page::Main;
    }
    return r;
}

void Menu_ResetPause() { pausePage = Page::Main; }
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
    if (Button(Row(4), "Save & Quit to Menu")) r.action = MenuAction::SaveAndQuit;
    return r;
}
