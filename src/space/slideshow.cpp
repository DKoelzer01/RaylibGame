#include "slideshow.h"
#include "menus.h"

#include <raylib.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <future>
#include <map>
#include <random>
#include <string>
#include <vector>

namespace {
    constexpr double kShow = 9.0, kFade = 1.5; // seconds per slide (including the fade in), crossfade
    const char* kDir = "resources/menu/slides";

    struct Slide { std::string path, caption; };
    struct Shown { Texture2D tex{}; int slide = -1; double start = 0.0; float panX = 0.0f, panY = 0.0f; bool zoomIn = true; };

    bool scanned = false;
    std::vector<Slide> slides;
    Shown cur, prev;                       // the slide showing and the one fading out under it
    std::future<Image> loading;            // the next slide, decoded on a worker
    int loadingSlide = -1;
    std::mt19937 rng{ std::random_device{}() };

    void Scan() {
        scanned = true;
        std::map<std::string, std::string> captions;
        if (FILE* f = fopen((std::string(kDir) + "/slides.txt").c_str(), "r")) {
            char line[512];
            while (fgets(line, sizeof line, f)) {
                std::string s = line;
                while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
                if (s.empty() || s[0] == '#') continue;
                const size_t bar = s.find('|');
                auto trim = [](std::string t) { const size_t a = t.find_first_not_of(" \t"), b = t.find_last_not_of(" \t"); return a == std::string::npos ? std::string() : t.substr(a, b - a + 1); };
                if (bar != std::string::npos) captions[trim(s.substr(0, bar))] = trim(s.substr(bar + 1));
            }
            fclose(f);
        }
        std::error_code ec;
        for (const auto& e : std::filesystem::directory_iterator(kDir, ec)) {
            if (!e.is_regular_file(ec)) continue;
            std::string ext = e.path().extension().string();
            for (char& c : ext) c = (char)tolower(c);
            if (ext != ".png" && ext != ".jpg" && ext != ".jpeg") continue;
            const std::string name = e.path().filename().string();
            slides.push_back({ e.path().generic_string(), captions.count(name) ? captions[name] : std::string() });
        }
        std::sort(slides.begin(), slides.end(), [](const Slide& a, const Slide& b) { return a.path < b.path; });
    }

    int PickNext(int notThis) {
        if (slides.size() <= 1) return slides.empty() ? -1 : 0;
        std::uniform_int_distribution<int> d(0, (int)slides.size() - 1);
        int k = d(rng);
        while (k == notThis) k = d(rng); // no immediate repeat
        return k;
    }

    void StartLoad(int k) {
        if (k < 0) return;
        loadingSlide = k;
        const std::string path = slides[k].path;
        loading = std::async(std::launch::async, [path] { return LoadImage(path.c_str()); }); // decode off the main thread
    }

    // A finished decode becomes the showing slide (the old one fades out under it).
    void TakeLoaded(double now) {
        if (!loading.valid() || loading.wait_for(std::chrono::seconds(0)) != std::future_status::ready) return;
        Image img = loading.get();
        if (!img.data) { loadingSlide = -1; return; } // unreadable: skip it (the next pick tries another)
        if (prev.tex.id) UnloadTexture(prev.tex);
        prev = cur;
        cur = Shown{};
        cur.tex = LoadTextureFromImage(img); // the upload: one texture, on the main thread
        UnloadImage(img);
        SetTextureFilter(cur.tex, TEXTURE_FILTER_BILINEAR);
        cur.slide = loadingSlide;
        cur.start = now;
        std::uniform_real_distribution<float> u(-1.0f, 1.0f);
        cur.panX = u(rng); cur.panY = u(rng) * 0.5f; cur.zoomIn = u(rng) > 0.0f;
        loadingSlide = -1;
    }

    void DrawCover(const Shown& s, double now, float alpha) {
        if (!s.tex.id || alpha <= 0.0f) return;
        const float sw = (float)GetScreenWidth(), sh = (float)GetScreenHeight();
        // Cover the window (cropping the image's longer side), then the slow pan and zoom (Ken Burns), off with Reduce motion.
        const float t = (float)std::clamp((now - s.start) / (kShow + kFade), 0.0, 1.0);
        const bool still = settings.reduceMotion;
        const float zoom = still ? 1.0f : (s.zoomIn ? 1.0f + 0.08f * t : 1.08f - 0.08f * t);
        const float k = std::max(sw / s.tex.width, sh / s.tex.height) * zoom;
        const float dw = s.tex.width * k, dh = s.tex.height * k;
        const float spareX = dw - sw, spareY = dh - sh;
        const float px = still ? 0.5f : 0.5f + 0.5f * s.panX * (t - 0.5f), py = still ? 0.5f : 0.5f + 0.5f * s.panY * (t - 0.5f);
        DrawTexturePro(s.tex, { 0, 0, (float)s.tex.width, (float)s.tex.height }, { -spareX * px, -spareY * py, dw, dh }, { 0, 0 }, 0.0f, Fade(WHITE, alpha));
    }
}

bool Slideshow_Draw(float dim) {
    if (!scanned) Scan();
    if (slides.empty()) return false;
    const double now = GetTime();
    TakeLoaded(now);
    // The first slide, or the next one due a fade before this one ends.
    if (!loading.valid() && loadingSlide < 0 && (!cur.tex.id || now - cur.start > kShow - 0.5)) StartLoad(PickNext(cur.slide));
    if (!cur.tex.id) return false; // the first decode is still running: the caller's backdrop shows
    const float fade = (float)std::clamp((now - cur.start) / kFade, 0.0, 1.0);
    if (fade < 1.0f) DrawCover(prev, now, 1.0f);
    DrawCover(cur, now, fade);
    if (fade >= 1.0f && prev.tex.id) { UnloadTexture(prev.tex); prev = Shown{}; }
    // Readability: darker in the middle band where the menu and the loading bar sit, and at the bottom.
    const int w = GetScreenWidth(), h = GetScreenHeight();
    const Color clear = Fade(BLACK, 0.15f), dark = Fade(BLACK, std::clamp(dim, 0.0f, 1.0f));
    DrawRectangleGradientV(0, 0, w, h / 2, clear, dark);
    DrawRectangleGradientV(0, h / 2, w, h - h / 2, dark, Fade(BLACK, std::clamp(dim + 0.15f, 0.0f, 1.0f)));
    // The caption, bottom left.
    const std::string& cap = slides[cur.slide].caption;
    if (!cap.empty()) DrawText(cap.c_str(), 24, h - 40, 20, Fade(RAYWHITE, 0.75f * fade));
    return true;
}

void Slideshow_Unload() {
    if (loading.valid()) { Image img = loading.get(); if (img.data) UnloadImage(img); }
    if (cur.tex.id) UnloadTexture(cur.tex);
    if (prev.tex.id) UnloadTexture(prev.tex);
    cur = prev = Shown{};
    loadingSlide = -1;
}
