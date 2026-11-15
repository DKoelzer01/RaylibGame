#include "recorder.h"
#include "world.h"

#include <raylib.h>
#include "external/msf_gif.h" // raylib bundles msf_gif (its own GIF recording) and exports it

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstdio>
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
    std::string compareDir;      // --record-compare: an earlier capture to diff each frame against
    double diffSum = 0.0, diffMax = 0.0; int diffFrames = 0, diffWorst = -1;

    // Run state.
    bool active = false, waiting = false, finished = false;
    int captured = 0, rendered = 0;
    double waitStart = 0.0;
    uint64_t firstHash = 0, lastHash = 0;
    struct Thumb { Image img; int frame; float t; };
    std::vector<Thumb> thumbs;     // for the contact sheet: every frame, 240 px wide

    // The writer thread: PNGs (ExportImage is slow) and the GIF, in order.
    struct Job { Image img; std::string path; bool gif = false; bool end = false; int centis = 2; };
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

    void WriterLoop() {
        for (;;) {
            Job j;
            {
                std::unique_lock<std::mutex> lk(m);
                cv.wait(lk, [] { return !queue.empty(); });
                j = queue.front();
                queue.pop_front();
            }
            cvRoom.notify_all();
            if (j.end) {
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

    void Push(Job j) {
        std::unique_lock<std::mutex> lk(m);
        cvRoom.wait(lk, [] { return queue.size() < 24; }); // bounded: a slow disk slows the capture, not memory
        queue.push_back(j);
        lk.unlock();
        cv.notify_one();
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
        if (p.is_absolute()) return d;
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

    void Begin(const std::string& d) {
        if (finisher.joinable()) finisher.join(); // the last F10 recording is still being written (shares the writer queue)
        dir = d;
        clockStart = GetTime();
        std::error_code ec;
        std::filesystem::create_directories(dir, ec);
        active = true; finished = false; captured = 0; rendered = 0;
        waiting = startDelay > 0.0f || cli; waitStart = GetTime();
        firstHash = lastHash = 0;
        thumbs.clear();
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
        active = false;
        if (gifPending.data) { // the last GIF frame: one capture interval
            Job j{ gifPending, "", true, false }; j.centis = std::max(2, gifCentis); gifPending = Image{};
            if (!background) Push(j);
            else { std::lock_guard<std::mutex> lk(m); queue.push_back(j); cv.notify_one(); } // never wait for room here
        }
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

void Recorder_Update(bool worldReady) {
    if (IsKeyPressed(KEY_F10) && !cli) { // normal play: start / stop a recording for a bug report
        if (active) End(true);
        else {
            char stamp[32];
            const std::time_t now = std::time(nullptr);
            std::strftime(stamp, sizeof stamp, "%Y%m%d-%H%M%S", std::localtime(&now));
            framesWanted = 600; every = 2; startDelay = 0.0f; hudFrames = true; wantGif = true; wantSheet = true;
            Begin(Resolve(Today("%H%M%S")));
            waiting = false;
        }
    }
    if (active && waiting && GetTime() - waitStart >= startDelay && worldReady) { waiting = false; gGameClock = 0.0; clockStart = GetTime(); } // animation clocks start with the capture
}

bool Recorder_WantsFrame(bool hud) {
    if (!active || waiting || hud != hudFrames) return false;
    return rendered % every == 0;
}

void Recorder_Capture() {
    // Called where WantsFrame said yes; counts rendered frames on the HUD pass (every frame reaches it).
    Image img = LoadImageFromScreen();
    ImageFormat(&img, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    const uint64_t h = Hash(img);
    if (!compareDir.empty()) { // the share of pixels that differ (any channel by more than 2) from the same frame of the earlier capture
        Image old = LoadImage(TextFormat("%s/frame_%04d.png", compareDir.c_str(), captured));
        if (old.data && old.width == img.width && old.height == img.height) {
            ImageFormat(&old, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
            const unsigned char* a = (const unsigned char*)img.data; const unsigned char* b = (const unsigned char*)old.data;
            size_t bad = 0; const size_t px = (size_t)img.width * img.height;
            for (size_t q = 0; q < px; ++q) if (abs(a[q * 4] - b[q * 4]) > 2 || abs(a[q * 4 + 1] - b[q * 4 + 1]) > 2 || abs(a[q * 4 + 2] - b[q * 4 + 2]) > 2) ++bad;
            const double share = 100.0 * (double)bad / (double)px;
            diffSum += share; ++diffFrames;
            if (share > diffMax) { diffMax = share; diffWorst = captured; }
        }
        if (old.data) UnloadImage(old);
    }
    if (captured == 0) firstHash = h;
    lastHash = h;
    const float t = fixedDt > 0.0f ? (float)(captured * every) * fixedDt : (float)(GetTime() - clockStart);
    if (wantSheet) { Image s = ImageCopy(img); ImageResize(&s, 240, std::max(1, img.height * 240 / std::max(1, img.width))); thumbs.push_back({ s, captured, t }); }
    if (wantGif) {
        const double nowCs = 100.0 * (double)t;            // this capture's time in the recording
        if (!gifPending.data) { gifEmittedCs = nowCs; }
        else if (nowCs - gifEmittedCs >= 2.0) {            // the held frame lasts until now
            const int d = (int)std::lround(nowCs - gifEmittedCs);
            Job j{ gifPending, "", true, false }; j.centis = d; Push(j);
            gifEmittedCs += d; gifPending = Image{};
        }
        if (!gifPending.data) { Image g = ImageCopy(img); ImageResize(&g, 640, std::max(1, img.height * 640 / std::max(1, img.width))); gifPending = g; gifPendingAt = nowCs; }
    }
    Push(Job{ img, TextFormat("%s/frame_%04d.png", dir.c_str(), captured), false, false });
    ++captured;
    if (captured >= framesWanted) End(!cli);
}

void Recorder_FrameDone() { if (active && !waiting) ++rendered; }

bool Recorder_Finished() { return finished; }
bool Recorder_Active() { return active; }
void Recorder_Shutdown() { End(); if (finisher.joinable()) finisher.join(); }
std::string Recorder_LastDir() { std::lock_guard<std::mutex> lk(lastDirM); return lastDir; }
