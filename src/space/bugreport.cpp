#include "bugreport.h"
#include "gameui.h"
#include "recorder.h"
#include "player.h"
#include "world.h"
#include "utils/logger.h"

#include "imgui.h"
#include <raylib.h>
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace fs = std::filesystem;

namespace {
    // The tags offered (each becomes a label on the issue; the relay adds "bug report").
    const char* kTags[] = { "rendering", "terrain", "water", "ui", "performance", "crash", "gameplay", "building", "machines", "animals", "audio" };
    constexpr int kTagCount = (int)(sizeof(kTags) / sizeof(kTags[0]));
    const char* kSeverity[] = { "minor", "major", "blocker" };

    bool open = false;
    char title[160] = {};
    char body[4096] = {};
    bool tags[kTagCount] = {};
    int severity = 0;
    bool attachShot = true, attachLog = true, attachCapture = true, attachSave = false;
    Image shot{};                       // the world as the game paused
    Scene* world = nullptr;


    std::atomic<bool> sending{ false }, done{ false };
    std::mutex resultM;
    std::string result;                 // the relay's answer or the error, for the toast ("" after a quiet retry)
    std::string logLine;                // the same for log.txt, in full (the folder path too: #617)
    std::string notice;                 // (#617) a failed send: shown until dismissed
    int unsent = -1;                    // reports in the outbox (-1: not counted yet)

    // (#617) Unsent reports wait in an outbox under saves/ (not %TEMP%, which gets cleaned): one folder each, with a
    // fields.txt listing what curl sends; they are sent again at start-up and every few minutes, and removed once sent.
    fs::path OutboxDir() { return fs::path("saves") / "bug-outbox"; }
    int CountOutbox() {
        std::error_code ec;
        int n = 0;
        if (fs::exists(OutboxDir(), ec))
            for (const auto& e : fs::directory_iterator(OutboxDir(), ec)) n += e.is_directory(ec) && fs::exists(e.path() / "fields.txt", ec) ? 1 : 0;
        return n;
    }

    std::string RelayUrl() {
        const char* e = getenv("RAYLIBGAME_BUG_RELAY");
        return e && e[0] ? e : "http://deweyserver:8765/report";
    }

    bool WriteText(const fs::path& p, const std::string& s) {
        FILE* f = fopen(p.string().c_str(), "wb");
        if (!f) return false;
        fwrite(s.data(), 1, s.size(), f);
        fclose(f);
        return true;
    }

    // Everything the relay gets, as files in a fresh folder (text goes in files too: no quoting problems on the
    // command line), then curl -F name=@file / name=<file.
    void Send(fs::path dir, std::vector<std::pair<std::string, std::string>> fields, bool retry) {
#ifdef _WIN32
        std::string cmd = "curl.exe -s -S -m 60";
#else
        std::string cmd = "curl -s -S -m 60";
#define _popen popen
#define _pclose pclose
#endif
        for (const auto& [name, spec] : fields) cmd += " -F \"" + name + "=" + spec + "\"";
        cmd += " \"" + RelayUrl() + "\" 2>&1";
        std::string out;
        if (FILE* p = _popen(cmd.c_str(), "r")) {
            char buf[512];
            while (fgets(buf, sizeof buf, p)) out += buf;
            const int rc = _pclose(p);
            if (rc != 0 && out.empty()) out = "curl failed (" + std::to_string(rc) + ")";
        } else out = "could not run curl.exe";
        // The relay answers {"url": "...", "number": N}; anything else is an error to show.
        std::string msg, log, note;
        const size_t u = out.find("\"url\"");
        if (u != std::string::npos) {
            const size_t a = out.find('"', out.find(':', u) + 1), b = a == std::string::npos ? a : out.find('"', a + 1);
            const std::string link = b != std::string::npos ? out.substr(a + 1, b - a - 1) : std::string("(no link)");
            msg = (retry ? "An unsent bug report was sent: " : "Bug reported: ") + link;
            log = msg + " (from " + dir.string() + ")";
            std::error_code ec;
            fs::remove_all(dir, ec); // sent: out of the outbox
        } else {
            while (!out.empty() && (out.back() == '\n' || out.back() == '\r')) out.pop_back();
            log = "Bug report not sent: " + out.substr(0, 400) + " (kept in " + fs::absolute(dir).string() + ")"; // (the path in full: #617)
            if (!retry) { // (a retry that fails again stays quiet: the first failure said it)
                msg = "Bug report not sent: the relay is offline";
                note = "Not sent: the bug relay is offline. Kept in the outbox; it will be sent when the relay is back.";
            }
        }
        { std::lock_guard<std::mutex> lk(resultM); result = msg; logLine = log; if (!note.empty()) notice = note; }
        done = true;
        sending = false;
    }

    // (#617) The oldest report in the outbox, sent again.
    bool ResendOne() {
        std::error_code ec;
        if (!fs::exists(OutboxDir(), ec)) return false;
        std::vector<fs::path> dirs;
        for (const auto& e : fs::directory_iterator(OutboxDir(), ec)) if (e.is_directory(ec) && fs::exists(e.path() / "fields.txt", ec)) dirs.push_back(e.path());
        if (dirs.empty()) return false;
        std::sort(dirs.begin(), dirs.end()); // (folders are named by time)
        std::vector<std::pair<std::string, std::string>> fields;
        if (FILE* f = fopen((dirs.front() / "fields.txt").string().c_str(), "rb")) {
            char line[2048];
            while (fgets(line, sizeof line, f)) {
                std::string s = line;
                while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
                const size_t tab = s.find('\t');
                if (tab != std::string::npos) fields.push_back({ s.substr(0, tab), s.substr(tab + 1) });
            }
            fclose(f);
        }
        if (fields.empty()) return false;
        sending = true;
        done = false;
        std::thread(Send, dirs.front(), fields, true).detach();
        return true;
    }

    void Submit() {
        char stamp[32];
        const std::time_t now = std::time(nullptr);
        std::strftime(stamp, sizeof stamp, "%Y%m%d-%H%M%S", std::localtime(&now));
        std::error_code ec;
        const fs::path dir = fs::absolute(OutboxDir() / stamp, ec); // (#617: kept until sent; not %TEMP%)
        fs::create_directories(dir, ec);
        std::vector<std::pair<std::string, std::string>> fields;
        WriteText(dir / "title.txt", title);
        WriteText(dir / "body.txt", body);
        std::string labels;
        for (int i = 0; i < kTagCount; ++i) if (tags[i]) labels += std::string(labels.empty() ? "" : ",") + kTags[i];
        WriteText(dir / "labels.txt", labels);
        WriteText(dir / "severity.txt", kSeverity[severity]);
        // What the game knows that helps a triage: build, seed, save, where the player is, the frame rate.
        std::string meta = TextFormat("build: %s %s\nseed: %u\nsave: %s\nfps: %d\nwindow: %dx%d\n", __DATE__, __TIME__, worldGen.seed, WorldSavePath(), GetFPS(), GetScreenWidth(), GetScreenHeight());
        if (world && player.nearest) {
            const Vector3 local = player.nearest->toLocal((player.position - world->renderOrigin).toFloat());
            meta += TextFormat("body: %s\nposition (body-local, m): %.2f %.2f %.2f\n", player.nearest->name.c_str(), local.x, local.y, local.z);
            const float d = Vector3Length(local);
            if (d > 1.0f) meta += TextFormat("spawn dir (RAYLIBGAME_SPAWN_DIR): %.6f,%.6f,%.6f\n", local.x / d, local.y / d, local.z / d);
        }
        WriteText(dir / "meta.txt", meta);
        for (const char* t : { "title", "body", "labels", "severity", "meta" }) fields.push_back({ t, "<" + (dir / (std::string(t) + ".txt")).string() });
        if (attachShot && shot.data) { const fs::path p = dir / "screenshot.png"; ExportImage(shot, p.string().c_str()); fields.push_back({ "screenshot", "@" + p.string() }); }
        if (attachLog) { // flushed first (it buffers), then read: the logger holds it open
            Log_Flush();
            std::string text;
            if (FILE* f = fopen("log.txt", "rb")) { char buf[8192]; size_t k; while ((k = fread(buf, 1, sizeof buf, f)) > 0) text.append(buf, k); fclose(f); }
            if (!text.empty() && WriteText(dir / "log.txt", text)) fields.push_back({ "log", "@" + (dir / "log.txt").string() });
        }
        if (attachCapture) { // (#617) copied in: the outbox folder holds everything a later resend needs
            const std::string cap = Recorder_LastDir();
            for (const char* f : { "sheet.png", "capture.gif" })
                if (!cap.empty() && fs::exists(fs::path(cap) / f, ec) && fs::copy_file(fs::path(cap) / f, dir / f, fs::copy_options::overwrite_existing, ec))
                    fields.push_back({ std::string("capture_") + (f[0] == 's' ? "sheet" : "gif"), "@" + (dir / f).string() });
        }
        if (attachSave && fs::exists(WorldSavePath(), ec)) {
            fs::copy_file(WorldSavePath(), dir / "world.sav", fs::copy_options::overwrite_existing, ec);
            fields.push_back({ "save", "@" + (dir / "world.sav").string() });
        }
        { // (#617) what curl sends, for a resend: name TAB spec, one a line
            std::string list;
            for (const auto& [name, spec] : fields) list += name + "\t" + spec + "\n";
            WriteText(dir / "fields.txt", list);
        }
        unsent = -1;

        sending = true;
        done = false;
        std::thread(Send, dir, fields, false).detach(); // (detached: a send still running at exit must not abort the game)
        GameUI_Toast("Sending the bug report...");
    }
}

void BugReport_CaptureScene() {
    if (shot.data) UnloadImage(shot);
    rlDrawRenderBatchActive();
    shot = LoadImageFromScreen();
}

void BugReport_Open(Scene& w) {
    world = &w;
    open = true;
    title[0] = body[0] = 0;
    for (bool& t : tags) t = false;
    severity = 0;
    attachShot = shot.data != nullptr;
    attachCapture = !Recorder_LastDir().empty();
}

bool BugReport_IsOpen() { return open; }

void BugReport_Draw() {
    const ImGuiIO& io = ImGui::GetIO();
    { // (#617) a failed send's notice: it stays until dismissed (a 2.5 s toast went unnoticed)
        std::string text;
        { std::lock_guard<std::mutex> lk(resultM); text = notice; }
        if (!text.empty()) {
            ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.12f), ImGuiCond_Always, ImVec2(0.5f, 0.0f));
            ImGui::SetNextWindowBgAlpha(0.92f);
            if (ImGui::Begin("##bugnotice", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing)) {
                ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.35f, 1.0f), "%s", text.c_str());
                if (ImGui::Button("OK")) { std::lock_guard<std::mutex> lk(resultM); notice.clear(); }
            }
            ImGui::End();
        }
    }
    if (!open) return;
    ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(620, 0), ImGuiCond_Appearing);
    if (ImGui::Begin("Report a bug", &open, ImGuiWindowFlags_NoCollapse)) {
        ImGui::TextDisabled("Goes to the project's GitHub issues (through the relay on deweyserver).");
        ImGui::InputText("Title", title, sizeof title);
        ImGui::InputTextMultiline("##what", body, sizeof body, ImVec2(-1, 160));
        ImGui::TextDisabled("What happened, what you expected, how to make it happen again.");
        ImGui::Separator();
        ImGui::TextUnformatted("Tags");
        for (int i = 0; i < kTagCount; ++i) { if (i % 6) ImGui::SameLine(); ImGui::Checkbox(kTags[i], &tags[i]); }
        ImGui::Combo("Severity", &severity, kSeverity, 3);
        ImGui::Separator();
        ImGui::TextUnformatted("Attach");
        ImGui::BeginDisabled(!shot.data); ImGui::Checkbox("Screenshot (as the game paused)", &attachShot); ImGui::EndDisabled();
        ImGui::Checkbox("log.txt", &attachLog);
        ImGui::BeginDisabled(Recorder_LastDir().empty()); ImGui::Checkbox("Last F10 capture (contact sheet and GIF)", &attachCapture); ImGui::EndDisabled();
        ImGui::Checkbox("The save (larger)", &attachSave);
        ImGui::Separator();
        ImGui::BeginDisabled(sending || title[0] == 0);
        if (ImGui::Button("Send", ImVec2(120, 0))) { Submit(); open = false; }
        ImGui::EndDisabled();
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0)) || ImGui::IsKeyPressed(ImGuiKey_Escape)) open = false;
        if (sending) { ImGui::SameLine(); ImGui::TextDisabled("sending the last one..."); }
    }
    ImGui::End();
}

void BugReport_TestSend(Scene& w) {
    BugReport_Open(w);
    snprintf(title, sizeof title, "Test report from --test-bug-report");
    snprintf(body, sizeof body, "A scripted report (#169): the relay should store the files and answer with an issue link.");
    tags[0] = tags[4] = true; // rendering, performance
    severity = 0;
    attachShot = shot.data != nullptr; attachLog = true; attachCapture = false; attachSave = false;
    Submit();
    open = false;
}
bool BugReport_Sending() { return sending; }

void BugReport_Update() {
    // (#617) The outbox: one report sent again at start-up (the first call) and every 5 minutes while any wait.
    static double nextRetry = 0.0;
    if (!sending && GetTime() >= nextRetry) {
        nextRetry = GetTime() + 300.0;
        ResendOne();
    }
    if (!done.exchange(false)) return;
    unsent = -1; // (recounted)
    std::lock_guard<std::mutex> lk(resultM);
    if (!result.empty()) GameUI_Toast(result.c_str());
    printf("bug report: %s\n", logLine.c_str()); // stdout too: --test-bugreport reads it
    if (logLine.rfind("Bug report not sent", 0) == 0) LOG_WARN("BugReport", "%s", logLine.c_str()); else LOG_INFO("BugReport", "%s", logLine.c_str());
}

int BugReport_Unsent() { // (#617) for the pause menu
    if (unsent < 0) unsent = CountOutbox();
    return unsent;
}
