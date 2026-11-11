#include "menus.h"
#include "player.h"
#include "ship.h"
#include "telescope.h"
#include "audio.h"
#include "foliage.h"

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
        const int fs = (int)(64 * s);
        DrawText(title, (GetScreenWidth() - MeasureText(title, fs)) / 2, (int)(GetScreenHeight() * 0.18f), fs, kText);
    }
    Rectangle Row(int i, float width = 420.0f, float top = 0.36f) {
        const float s = Scale(), w = width * s, h = 58.0f * s;
        return { (GetScreenWidth() - w) * 0.5f, GetScreenHeight() * top + i * (h + 14.0f * s), w, h };
    }

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
        float y = GetScreenHeight() * 0.32f;
        bool changed = false;
        changed |= Slider({ x, y, w, h }, "Mouse sensitivity", settings.mouseSensitivity, 0.2f, 3.0f, "%.2fx"); y += h + 10 * s;
        changed |= Slider({ x, y, w, h }, "Field of view", settings.fov, 30.0f, 100.0f, "%.0f deg"); y += h + 10 * s;
        changed |= Slider({ x, y, w, h }, "Volume", settings.volume, 0.0f, 1.0f, "%.2f"); y += h + 10 * s;
        changed |= Toggle({ x, y, w, h }, "Foliage", settings.foliage); y += h + 10 * s;
        changed |= Toggle({ x, y, w, h }, "VSync (no tearing)", settings.vsync); y += h + 10 * s;
        changed |= Toggle({ x, y, w, h }, "Reduce motion", settings.reduceMotion); y += h + 10 * s;
        {
            const int fs = (int)(h * 0.5f);
            DrawText("Frame cap", (int)x, (int)(y + (h - fs) * 0.5f), fs, kText);
            const int caps[] = { 0, 60, 144, 240 };
            const char* names[] = { "Monitor", "60", "144", "240" };
            for (int k = 0; k < 4; ++k) {
                const Rectangle b = { x + w * 0.38f + k * (h * 2.4f), y + h * 0.1f, h * 2.2f, h * 0.8f };
                if (Button(b, settings.fpsCap == caps[k] ? TextFormat("[%s]", names[k]) : names[k])) { settings.fpsCap = caps[k]; changed = true; }
            }
            y += h + 10 * s;
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
        char key[64];
        float v;
        while (fscanf(f, " %63s %f", key, &v) == 2) {
            if (!strcmp(key, "mouse")) settings.mouseSensitivity = std::clamp(v, 0.2f, 3.0f);
            else if (!strcmp(key, "fov")) settings.fov = std::clamp(v, 30.0f, 100.0f);
            else if (!strcmp(key, "volume")) settings.volume = std::clamp(v, 0.0f, 1.0f);
            else if (!strcmp(key, "foliage")) settings.foliage = v != 0.0f;
            else if (!strcmp(key, "fpscap")) settings.fpsCap = std::max(0, (int)v);
            else if (!strcmp(key, "vsync")) settings.vsync = v != 0.0f;
            else if (!strcmp(key, "reducemotion")) settings.reduceMotion = v != 0.0f;
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
    fclose(f);
}

void Settings_Apply() {
    player.mouseSensitivity = kDefaultMouse * settings.mouseSensitivity;
    ship.mouseSensitivity = kDefaultShipMouse * settings.mouseSensitivity;
    telescope.baseFov = settings.fov; // the camera's fov is set by the game each frame it is not zoomed
    Audio_SetMasterVolume(settings.volume);
    foliageEnabled = settings.foliage;
    if (settings.vsync) SetWindowState(FLAG_VSYNC_HINT); else ClearWindowState(FLAG_VSYNC_HINT);
    SetTargetFPS(settings.fpsCap > 0 ? settings.fpsCap : (settings.vsync ? 0 : std::max(60, GetMonitorRefreshRate(GetCurrentMonitor()))));
    player.smoothMotion = settings.reduceMotion;
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

MenuResult Menu_Main() {
    MenuResult r;
    Backdrop("PLANETOID");
    const float s = Scale();
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
        for (int i = first; i < (int)saves.size() && i < first + 7; ++i) {
            const SaveEntry& e = saves[i];
            if (Button(Row(i - first, 620.0f, 0.32f), TextFormat("%s   %s", e.name.c_str(), e.when.c_str()))) { r.action = MenuAction::LoadGame; r.path = e.path; mainPage = Page::Main; }
        }
        if (Button(Row(7, 300.0f, 0.32f), "Back") || IsKeyPressed(KEY_ESCAPE)) mainPage = Page::Main;
    } else if (SettingsPage()) {
        mainPage = Page::Main;
    }
    return r;
}

void Menu_ResetPause() { pausePage = Page::Main; }

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
    const Rectangle panel = { Row(0).x - 30 * s, Row(0).y - 30 * s, Row(0).width + 60 * s, Row(3).y + Row(3).height - Row(0).y + 60 * s };
    DrawRectangleRounded(panel, 0.06f, 6, kPanel);
    if (Button(Row(0), "Resume") || IsKeyPressed(KEY_ESCAPE)) r.action = MenuAction::Resume;
    if (Button(Row(1), "Save")) r.action = MenuAction::Save;
    if (Button(Row(2), "Settings")) pausePage = Page::Settings;
    if (Button(Row(3), "Save & Quit to Menu")) r.action = MenuAction::SaveAndQuit;
    return r;
}
