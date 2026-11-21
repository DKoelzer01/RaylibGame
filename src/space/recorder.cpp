#include "recorder.h"
#include "world.h"
#include "gameui.h" // (#1119) the start and stop toasts

#include <raylib.h>
#include "external/msf_gif.h" // raylib bundles msf_gif (its own GIF recording) and exports it

#include <algorithm>
#include <climits>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <deque>
#include <filesystem>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace {
    // Settings (command line, or the F10 defaults).
    std::string dir;
    int framesWanted = 120, every = 1;
    float startDelay = 0.0f, fixedDt = 0.0f;
    bool hudFrames = true, keepRunning = false, wantGif = false, wantSheet = false, cli = false;
    bool onCue = false, cued = false; // --record-on-cue: hold the first frame until the scenario calls Recorder_Cue()
    std::string compareDir;      // --record-compare: an earlier capture to diff each frame against
    double diffSum = 0.0, diffMax = 0.0; int diffFrames = 0, diffWorst = -1;

    // Run state.
    bool active = false, waiting = false, finished = false;
    int captured = 0, rendered = 0;
    double waitStart = 0.0;
    uint64_t firstHash = 0, lastHash = 0;
    struct Thumb { Image img; int frame; float t; };
    std::vector<Thumb> thumbs;     // for the contact sheet: every thumbStride-th frame, 240 px wide
    // (#1119) F10 records until F10 again: the sheet's thumbnails thin out to every 2nd, 4th... frame at kMaxThumbs (still
    // evenly spaced, first included), and the GIF covers the first gifMaxFrames frames (0: all, --record).
    constexpr int kMaxThumbs = 256;
    int thumbStride = 1, gifMaxFrames = 0;
    constexpr double kF10MaxSeconds = 30.0 * 60.0;  // a safety cap on one F10 recording
    constexpr double kF10MinFreeBytes = 2.0e9;     // stop before the disk fills
    double FreeBytes(const std::string& path) {    // free space on the drive or share of path (its nearest existing folder)
        std::error_code ec;
        std::filesystem::path p(path);
        while (!p.empty() && !std::filesystem::exists(p, ec)) p = p.parent_path();
        const auto s = std::filesystem::space(p.empty() ? std::filesystem::path(".") : p, ec);
        return ec ? 1e18 : (double)s.available; // unknown: don't stop a recording over it
    }

    // The writer thread takes captured frames in order (hash, compare, contact-sheet thumb, GIF) and hands their PNGs to
    // a small pool (ExportImage is the slow part). The main thread only reads the screen back (#537: hashing, two resizes
    // and a single PNG writer behind a full queue on the main thread slowed F10 play to a crawl). F10 drops a frame
    // when the writer is behind instead of waiting; --record waits (its frames are numbered).
    struct Job { Image img; std::string path; bool gif = false; bool end = false; int centis = 2; bool frame = false; int index = 0; float t = 0.0f; };
    std::atomic<int> dropped{0};
    std::deque<Job> pngQueue;
    std::mutex pngM;
    std::condition_variable pngCv, pngRoom;
    std::vector<std::thread> pngPool;
    bool pngEnd = false;
    void PngLoop() {
        for (;;) {
            Job j;
            {
                std::unique_lock<std::mutex> lk(pngM);
                pngCv.wait(lk, [] { return !pngQueue.empty() || pngEnd; });
                if (pngQueue.empty()) return;
                j = pngQueue.front();
                pngQueue.pop_front();
            }
            pngRoom.notify_all();
            ExportImage(j.img, j.path.c_str());
            UnloadImage(j.img);
        }
    }
    void PushPng(Job j) {
        std::unique_lock<std::mutex> lk(pngM);
        if (!cli && pngQueue.size() >= 24) { lk.unlock(); UnloadImage(j.img); ++dropped; return; } // F10: skip this PNG
        pngRoom.wait(lk, [] { return pngQueue.size() < 24; });
        pngQueue.push_back(j);
        lk.unlock();
        pngCv.notify_one();
    }
    std::deque<Job> queue;
    std::mutex m;
    std::condition_variable cv, cvRoom;
    std::thread writer;
    bool writerRunning = false;
    MsfGifState gifState{};
    bool gifStarted = false;
    std::string gifPath;
    int gifCentis = 2;
    // GIF pacing (#375: the GIF played ~1.2x fast): delays are whole centiseconds, at least 2 (viewers clamp shorter ones),
    // so one GIF frame is held until the next capture at least 2 cs later, for the time up to it, the rounding carried on:
    // the GIF lasts exactly as long as the capture.
    Image gifPending{};
    double gifPendingAt = 0.0, gifEmittedCs = 0.0;
    // Frame times: the fixed steps with --fixed-dt (--record); in F10 play the real clock since the first frame (the sheet
    // labelled 0.5 s apart frames 1/60 s apart, #375).
    double clockStart = 0.0;
    // F10 finishes in the background (draining PNGs to the share, the GIF, the sheet, pruning: a freeze, #375).
    std::thread finisher;
    std::mutex lastDirM;
    std::string lastDir; // the last finished recording (bug reports, #169)

    uint64_t Hash(const Image& img) {
        uint64_t h = 1469598103934665603ull;
        const unsigned char* p = (const unsigned char*)img.data;
        const size_t n = (size_t)img.width * img.height * 4;
        for (size_t i = 0; i < n; ++i) { h ^= p[i]; h *= 1099511628211ull; }
        return h;
    }

    void EmitGif(const Image& img, int centis) {
        if (!gifStarted) { msf_gif_begin(&gifState, img.width, img.height); gifStarted = true; }
        msf_gif_frame(&gifState, (uint8_t*)img.data, centis, 16, img.width * 4);
    }

    // One captured frame, in order, on the writer thread.
    void ProcessFrame(Job& j) {
        Image& img = j.img;
        const uint64_t h = Hash(img);
        if (!compareDir.empty()) { // the share of pixels that differ (any channel by more than 2) from the same frame of the earlier capture
            char oldPath[1024]; snprintf(oldPath, sizeof oldPath, "%s/frame_%04d.png", compareDir.c_str(), j.index); // (TextFormat is not thread-safe)
            Image old = LoadImage(oldPath);
            if (old.data && old.width == img.width && old.height == img.height) {
                ImageFormat(&old, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
                const unsigned char* a = (const unsigned char*)img.data; const unsigned char* b = (const unsigned char*)old.data;
                size_t bad = 0; const size_t px = (size_t)img.width * img.height;
                for (size_t q = 0; q < px; ++q) if (abs(a[q * 4] - b[q * 4]) > 2 || abs(a[q * 4 + 1] - b[q * 4 + 1]) > 2 || abs(a[q * 4 + 2] - b[q * 4 + 2]) > 2) ++bad;
                const double share = 100.0 * (double)bad / (double)px;
                diffSum += share; ++diffFrames;
                if (share > diffMax) { diffMax = share; diffWorst = j.index; }
            }
            if (old.data) UnloadImage(old);
        }
        if (j.index == 0) firstHash = h;
        lastHash = h;
        if (wantSheet && j.index % thumbStride == 0) {
            Image s = ImageCopy(img); ImageResize(&s, 240, std::max(1, img.height * 240 / std::max(1, img.width))); thumbs.push_back({ s, j.index, j.t });
            if ((int)thumbs.size() >= kMaxThumbs) { // keep every other one: half as many, twice as far apart (#1119)
                size_t w = 0;
                for (size_t k = 0; k < thumbs.size(); ++k) { if (k % 2 == 0) thumbs[w++] = thumbs[k]; else UnloadImage(thumbs[k].img); }
                thumbs.resize(w);
                thumbStride *= 2;
            }
        }
        if (wantGif && (gifMaxFrames <= 0 || j.index < gifMaxFrames)) {
            const double nowCs = 100.0 * (double)j.t;          // this capture's time in the recording
            if (!gifPending.data) { gifEmittedCs = nowCs; }
            else if (nowCs - gifEmittedCs >= 2.0) {            // the held frame lasts until now
                const int d = (int)std::lround(nowCs - gifEmittedCs);
                EmitGif(gifPending, d); UnloadImage(gifPending);
                gifEmittedCs += d; gifPending = Image{};
            }
            if (!gifPending.data) { Image g = ImageCopy(img); ImageResize(&g, 640, std::max(1, img.height * 640 / std::max(1, img.width))); gifPending = g; gifPendingAt = nowCs; }
        }
        PushPng(Job{ img, j.path });
    }

    void WriterLoop() {
        pngEnd = false;
        for (int k = 0; k < 3; ++k) pngPool.emplace_back(PngLoop);
        for (;;) {
            Job j;
            {
                std::unique_lock<std::mutex> lk(m);
                cv.wait(lk, [] { return !queue.empty(); });
                j = queue.front();
                queue.pop_front();
            }
            cvRoom.notify_all();
            if (j.frame) { ProcessFrame(j); continue; }
            if (j.end) {
                if (gifPending.data) { EmitGif(gifPending, std::max(2, gifCentis)); UnloadImage(gifPending); gifPending = Image{}; } // the last frame: one interval
                { std::lock_guard<std::mutex> lk(pngM); pngEnd = true; }
                pngCv.notify_all();
                for (std::thread& t : pngPool) t.join();
                pngPool.clear();
                if (gifStarted) {
                    MsfGifResult r = msf_gif_end(&gifState);
                    if (r.data) { if (FILE* f = fopen(gifPath.c_str(), "wb")) { fwrite(r.data, r.dataSize, 1, f); fclose(f); } }
                    msf_gif_free(r);
                    gifStarted = false;
                }
                return;
            }
            if (j.gif) {
                if (!gifStarted) { msf_gif_begin(&gifState, j.img.width, j.img.height); gifStarted = true; }
                msf_gif_frame(&gifState, (uint8_t*)j.img.data, j.centis, 16, j.img.width * 4);
            } else {
                ExportImage(j.img, j.path.c_str());
            }
            UnloadImage(j.img);
        }
    }

    // False when F10 dropped the frame (the writer is behind).
    bool Push(Job j) {
        std::unique_lock<std::mutex> lk(m);
        if (!cli && j.frame && queue.size() >= 8) return false;
        cvRoom.wait(lk, [] { return queue.size() < 24; }); // bounded: a slow disk slows the capture, not memory
        queue.push_back(j);
        lk.unlock();
        cv.notify_one();
        return true;
    }

    // Captures live outside the repo, one folder per calendar day (user 2026-11-14): <root>/<YYYY-MM-DD>/<name>. The root is
    // RAYLIBGAME_CAPTURE_ROOT (e.g. \\DEWEYSERVER\files\raylib-captures) or %USERPROFILE%\RaylibCaptures; an absolute
    // --record path is used as it is.
    std::string CaptureRoot() {
        if (const char* r = getenv("RAYLIBGAME_CAPTURE_ROOT"); r && r[0]) return r;
        // The Claude media store on deweyserver's 8 TB drive (user 2026-11-14), when reachable (checked once a run: an
        // offline server can take a few seconds to time out); else a local folder.
        static int server = -1;
        static const char* kServer = "//DEWEYSERVER/files/claude-media/raylib-captures";
        if (server < 0) { std::error_code ec; server = std::filesystem::is_directory(kServer, ec) ? 1 : 0; }
        if (server == 1) return kServer;
        const char* home = getenv("USERPROFILE");
        return std::string(home && home[0] ? home : ".") + "/RaylibCaptures";
    }
    std::string Today(const char* fmt) {
        char s[32];
        const std::time_t now = std::time(nullptr);
        std::strftime(s, sizeof s, fmt, std::localtime(&now));
        return s;
    }
    std::string Resolve(const std::string& d) {
        const std::filesystem::path p(d);
        // A UNC path (\\server\share\... or //server/...) is absolute too: this libstdc++ does not say so, and one was joined
        // under the root as <root>/<day>/DEWEYSERVER/files/... (#628).
        if (p.is_absolute() || (d.size() > 2 && (d[0] == '\\' || d[0] == '/') && (d[1] == '\\' || d[1] == '/'))) return d;
        return CaptureRoot() + "/" + Today("%Y-%m-%d") + "/" + d;
    }

    // Keep the captures under RAYLIBGAME_CAPTURE_MAX_GB (default 10 GB): delete whole captures, oldest day first and within
    // a day the oldest capture first, never the one just written.
    void Prune(const std::string& keep) {
        namespace fs = std::filesystem;
        std::error_code ec;
        const fs::path root(CaptureRoot());
        if (!fs::is_directory(root, ec)) return;
        const char* g = getenv("RAYLIBGAME_CAPTURE_MAX_GB");
        const double maxBytes = (g && atof(g) > 0.0 ? atof(g) : 10.0) * 1024.0 * 1024.0 * 1024.0;
        struct Cap { fs::path path; std::string day; fs::file_time_type time; uintmax_t bytes; };
        std::vector<Cap> caps;
        uintmax_t total = 0;
        for (const auto& day : fs::directory_iterator(root, ec)) {
            if (!day.is_directory(ec)) continue;
            for (const auto& c : fs::directory_iterator(day.path(), ec)) {
                if (!c.is_directory(ec)) continue;
                uintmax_t b = 0;
                for (const auto& f : fs::recursive_directory_iterator(c.path(), ec)) if (f.is_regular_file(ec)) b += f.file_size(ec);
                caps.push_back({ c.path(), day.path().filename().string(), c.last_write_time(ec), b });
                total += b;
            }
        }
        std::sort(caps.begin(), caps.end(), [](const Cap& a, const Cap& b) { return a.day != b.day ? a.day < b.day : a.time < b.time; });
        const fs::path keepPath = fs::weakly_canonical(fs::path(keep), ec);
        int removed = 0;
        for (const Cap& c : caps) {
            if ((double)total <= maxBytes) break;
            if (fs::weakly_canonical(c.path, ec) == keepPath) continue;
            fs::remove_all(c.path, ec);
            if (!ec) { total -= c.bytes; ++removed; }
        }
        for (const auto& day : fs::directory_iterator(root, ec)) if (day.is_directory(ec) && fs::is_empty(day.path(), ec)) fs::remove(day.path(), ec);
        if (removed > 0) printf("record: removed the %d oldest captures to stay under %.1f GB\n", removed, maxBytes / (1024.0 * 1024.0 * 1024.0));
    }

    // (#1119) a capture folder as the player knows it: "<date>/<HHMMSS>" (the last two parts of the path)
    std::string FolderName(const std::string& d) {
        const std::filesystem::path p(d);
        return (p.parent_path().filename() / p.filename()).generic_string();
    }
    void Begin(const std::string& d) {
        if (finisher.joinable()) finisher.join(); // the last F10 recording is still being written (shares the writer queue)
        dir = d;
        clockStart = GetTime();
        std::error_code ec;
        std::filesystem::create_directories(dir, ec);
        active = true; finished = false; captured = 0; rendered = 0; dropped = 0;
        waiting = startDelay > 0.0f || cli; waitStart = GetTime();
        firstHash = lastHash = 0;
        thumbs.clear(); thumbStride = 1;
        const float dt = fixedDt > 0.0f ? fixedDt : 1.0f / 60.0f;
        gifCentis = std::max(2, (int)lroundf(100.0f * dt * (float)every));
        gifPath = dir + "/capture.gif";
        writer = std::thread(WriterLoop);
        writerRunning = true;
    }

    // Finish: drain the writer (and the GIF), then the contact sheet (CPU images; ImageDrawText reads the default font's
    // glyph images) and pruning. On this thread for --record (the run exits after it); for F10 on a finisher thread.
    void Finish();
    void End(bool background = false) {
        if (!active) return;
        active = false; // (the writer emits the last GIF frame when it reaches the end job)
        if (!cli) GameUI_Toast(TextFormat("Recording saved: %s", FolderName(dir).c_str())); // (#1119) F10's, by the key or when full
        if (background) { finisher = std::thread(Finish); return; }
        Finish();
    }
    void Finish() {
        if (writerRunning) { Push(Job{ Image{}, "", false, true }); writer.join(); writerRunning = false; }
        std::string sheetPath = "-";
        if (wantSheet && !thumbs.empty()) {
            const int n = std::min<int>(25, (int)thumbs.size());
            const int cols = n <= 4 ? n : n <= 9 ? 3 : n <= 16 ? 4 : 5, rows = (n + cols - 1) / cols;
            const int tw = thumbs[0].img.width, th = thumbs[0].img.height, pad = 6, label = 22;
            Image sheet = GenImageColor(cols * (tw + pad) + pad, rows * (th + label + pad) + pad, { 18, 18, 22, 255 });
            for (int k = 0; k < n; ++k) {
                const int idx = n == 1 ? 0 : (int)lroundf((float)k * (float)(thumbs.size() - 1) / (float)(n - 1)); // evenly spaced, first and last included
                const Thumb& t = thumbs[idx];
                const int x = pad + (k % cols) * (tw + pad), y = pad + (k / cols) * (th + label + pad);
                ImageDraw(&sheet, t.img, { 0, 0, (float)tw, (float)th }, { (float)x, (float)(y + label), (float)tw, (float)th }, WHITE);
                ImageDrawText(&sheet, TextFormat("#%d  %.2f s", t.frame, t.t), x + 2, y + 3, 16, RAYWHITE);
            }
            sheetPath = dir + "/sheet.png";
            ExportImage(sheet, sheetPath.c_str());
            UnloadImage(sheet);
        }
        for (Thumb& t : thumbs) UnloadImage(t.img);
        thumbs.clear();
        printf("record: %s frames=%d gif=%s sheet=%s hash=%016llx..%016llx\n", dir.c_str(), captured, wantGif && captured > 0 ? gifPath.c_str() : "-",
               sheetPath.c_str(), (unsigned long long)firstHash, (unsigned long long)lastHash);
        if (dropped > 0) printf("record: skipped %d frames / PNGs while the writer was behind (F10 keeps the game running, #537)\n", dropped.load());
        fflush(stdout);
        if (diffFrames > 0) printf("record: compared with %s: %d frames, %.3f%% of pixels differ on average, worst %.3f%% (frame %d)\n", compareDir.c_str(), diffFrames, diffSum / diffFrames, diffMax, diffWorst);
        Prune(dir);
        { std::lock_guard<std::mutex> lk(lastDirM); lastDir = dir; }
        if (cli && !keepRunning) finished = true;
    }
}

bool Recorder_ParseArg(int argc, char** argv, int& i) {
    const std::string a = argv[i];
    auto next = [&](const char*& v) { if (i + 1 >= argc) return false; v = argv[++i]; return true; };
    const char* v = nullptr;
    if (a == "--record" && next(v)) { dir = v; cli = true; if (fixedDt <= 0.0f) fixedDt = 1.0f / 60.0f; return true; }
    if (a == "--frames" && next(v)) { framesWanted = std::max(1, atoi(v)); return true; }
    if (a == "--record-every" && next(v)) { every = std::max(1, atoi(v)); return true; }
    if (a == "--record-start" && next(v)) { startDelay = (float)atof(v); return true; }
    if (a == "--record-hud" && next(v)) { hudFrames = atoi(v) != 0; return true; }
    if (a == "--record-continue") { keepRunning = true; return true; }
    if (a == "--record-on-cue") { onCue = true; return true; }
    if (a == "--gif") { wantGif = true; return true; }
    if (a == "--sheet") { wantSheet = true; return true; }
    if (a == "--record-compare" && next(v)) { compareDir = v; return true; }
    if (a == "--fixed-dt" && next(v)) { fixedDt = std::max(0.0f, (float)atof(v)); return true; }
    return false;
}

void Recorder_Start() {
    gFixedFrameTime = fixedDt; // GameFrameTime() (world.h): the simulation's frame time while recording
    if (worldGen.seed != 0) SetRandomSeed(worldGen.seed); // raylib's GetRandomValue, as the seed says
    if (!compareDir.empty()) compareDir = Resolve(compareDir); // a name: today's capture of that name
    if (cli) Begin(Resolve(dir));
}

namespace {
    // F2 (#579): one screenshot, HUD included, into <root>/<YYYY-MM-DD>/shot-<HHMMSS>/ beside the recordings
    // (screenshot.png, and sheet.png for the media tools, which show a capture folder by its sheet), written on its own
    // thread so the frame doesn't wait for the share.
    std::atomic<bool> shotWanted{ false };
    bool RecordingFrame(bool hud) { return active && !waiting && hud == hudFrames && rendered % every == 0; }
    std::function<std::string(const Image&)> shotExtra; // (#910) the LOD dump written beside the screenshot
    void TakeShot() {
        Image img = LoadImageFromScreen();
        ImageFormat(&img, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
        const std::string d = Resolve("shot-" + Today("%H%M%S"));
        std::string extra = shotExtra ? shotExtra(img) : std::string();
        std::thread([img, d, extra]() {
            std::error_code ec;
            std::filesystem::create_directories(d, ec);
            if (!extra.empty()) if (FILE* f = fopen((d + "/lod.txt").c_str(), "wb")) { fwrite(extra.data(), 1, extra.size(), f); fclose(f); }
            const bool ok = ExportImage(img, (d + "/screenshot.png").c_str()) && ExportImage(img, (d + "/sheet.png").c_str());
            UnloadImage(img);
            printf("screenshot: %s %s\n", d.c_str(), ok ? "saved" : "FAILED");
            Prune(d);
        }).detach();
    }
}

void Recorder_Screenshot() { shotWanted = true; }

// (#1119) While an F10 recording runs: a blinking red dot, REC, the time so far and what is left, top right. Drawn after the
// frame was captured (main.cpp calls it after Recorder_Capture), so it never shows in the recording itself.
void Recorder_DrawIndicator() {
    if (!active || cli) return;
    const int t = std::max(0, (int)(GetTime() - clockStart)); // (#1119) the time so far only: it records until F10 again
    const char* text = TextFormat("REC %d:%02d", t / 60, t % 60);
    const int fs = 20, w = MeasureText(text, fs), x = GetScreenWidth() - w - 46, y = 14;
    DrawRectangle(x - 34, y - 6, w + 46, fs + 12, Fade(BLACK, 0.55f));
    if (fmod(GetTime(), 1.0) < 0.6) DrawCircle(x - 17, y + fs / 2, 7.0f, RED);
    DrawText(text, x, y, fs, RAYWHITE);
}

void Recorder_SetShotExtra(std::function<std::string(const Image&)> f) { shotExtra = std::move(f); }
void Recorder_Cue() { cued = true; }

void Recorder_Update(bool worldReady) {
    if (IsKeyPressed(KEY_F2)) shotWanted = true; // a screenshot into the captures store (#579)
    if (IsKeyPressed(KEY_F10) && !cli) { // normal play: start / stop a recording for a bug report
        if (active) End(true);
        else {
            char stamp[32];
            const std::time_t now = std::time(nullptr);
            std::strftime(stamp, sizeof stamp, "%Y%m%d-%H%M%S", std::localtime(&now));
            // (#1119) until F10 again (or the safety cap, or a nearly full disk): every 2nd frame; the GIF is its first 30 s
            framesWanted = INT_MAX; every = 2; startDelay = 0.0f; hudFrames = true; wantGif = true; wantSheet = true; gifMaxFrames = 900;
            const std::string where = Resolve(Today("%H%M%S"));
            if (FreeBytes(where) < kF10MinFreeBytes) { GameUI_Toast("Not recording: the disk is nearly full"); return; }
            Begin(where);
            waiting = false;
            GameUI_Toast(TextFormat("Recording to %s (F10 stops)", FolderName(dir).c_str())); // (#1119)
        }
    }
    if (active && waiting && GetTime() - waitStart >= startDelay && worldReady && (!onCue || cued)) { waiting = false; gGameClock = 0.0; clockStart = GetTime(); } // animation clocks start with the capture
}

bool Recorder_WantsFrame(bool hud) {
    if (hud && shotWanted) return true; // F2
    return RecordingFrame(hud);
}

void Recorder_Capture() {
    if (shotWanted.exchange(false)) { TakeShot(); if (!RecordingFrame(true)) return; } // F2 (#579); a recording may want this frame too
    // Called where WantsFrame said yes; counts rendered frames on the HUD pass (every frame reaches it).
    // Only the readback here; the writer thread does the rest (#537). F10 skips the frame when the writer is behind.
    if (!cli) { std::lock_guard<std::mutex> lk(m); if (queue.size() >= 8) { ++dropped; return; } }
    Image img = LoadImageFromScreen();
    ImageFormat(&img, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    const float t = fixedDt > 0.0f ? (float)(captured * every) * fixedDt : (float)(GetTime() - clockStart);
    Job j{ img, TextFormat("%s/frame_%04d.png", dir.c_str(), captured) };
    j.frame = true; j.index = captured; j.t = t;
    if (!Push(j)) { UnloadImage(img); ++dropped; return; }
    ++captured;
    if (!cli) { // (#1119) F10: the safety cap, and a disk check every 300 frames
        const bool tooLong = GetTime() - clockStart >= kF10MaxSeconds;
        const bool diskFull = captured % 300 == 0 && FreeBytes(dir) < kF10MinFreeBytes;
        if (tooLong || diskFull) { GameUI_Toast(tooLong ? "Recording stopped at 30 minutes" : "Recording stopped: the disk is nearly full"); End(true); return; }
    }
    if (captured >= framesWanted) End(!cli);
}

void Recorder_FrameDone() { if (active && !waiting) ++rendered; }

bool Recorder_Finished() { return finished; }
bool Recorder_Active() { return active; }
int Recorder_Captured() { return active && !waiting ? captured : -1; } // (#545) the next frame_NNNN index
void Recorder_Shutdown() { End(); if (finisher.joinable()) finisher.join(); }
std::string Recorder_LastDir() { std::lock_guard<std::mutex> lk(lastDirM); return lastDir; }
