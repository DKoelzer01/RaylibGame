#include "bugreport.h"
#include "gameui.h"
#include "recorder.h"
#include "player.h"
#include "world.h"
#include "utils/logger.h"

#include "imgui.h"
#include <raylib.h>
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
    std::string result;                 // the relay's answer or the error, for the toast

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
    void Send(fs::path dir, std::vector<std::pair<std::string, std::string>> fields) {
        std::string cmd = "curl.exe -s -S -m 60";
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
        std::string msg;
        const size_t u = out.find("\"url\"");
        if (u != std::string::npos) {
            const size_t a = out.find('"', out.find(':', u) + 1), b = a == std::string::npos ? a : out.find('"', a + 1);
            msg = "Bug reported: " + (b != std::string::npos ? out.substr(a + 1, b - a - 1) : std::string("(no link)"));
        } else {
            msg = "Bug report not sent: " + out.substr(0, 160) + " (kept in " + dir.string() + ")";
        }
        { std::lock_guard<std::mutex> lk(resultM); result = msg; }
        done = true;
        sending = false;
    }

    void Submit() {
        char stamp[32];
        const std::time_t now = std::time(nullptr);
        std::strftime(stamp, sizeof stamp, "%Y%m%d-%H%M%S", std::localtime(&now));
        const char* tmp = getenv("TEMP");
        const fs::path dir = fs::path(tmp && tmp[0] ? tmp : ".") / "raylibgame-bug" / stamp;
        std::error_code ec;
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
            logger.flush();
            std::string text;
            if (FILE* f = fopen("log.txt", "rb")) { char buf[8192]; size_t k; while ((k = fread(buf, 1, sizeof buf, f)) > 0) text.append(buf, k); fclose(f); }
            if (!text.empty() && WriteText(dir / "log.txt", text)) fields.push_back({ "log", "@" + (dir / "log.txt").string() });
        }
        if (attachCapture) {
            const std::string cap = Recorder_LastDir();
            for (const char* f : { "sheet.png", "capture.gif" })
                if (!cap.empty() && fs::exists(fs::path(cap) / f, ec)) fields.push_back({ std::string("capture_") + (f[0] == 's' ? "sheet" : "gif"), "@" + (fs::path(cap) / f).string() });
        }
        if (attachSave && fs::exists(WorldSavePath(), ec)) {
            fs::copy_file(WorldSavePath(), dir / "world.sav", fs::copy_options::overwrite_existing, ec);
            fields.push_back({ "save", "@" + (dir / "world.sav").string() });
        }

        sending = true;
        done = false;
        std::thread(Send, dir, fields).detach(); // (detached: a send still running at exit must not abort the game)
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
    if (!open) return;
    const ImGuiIO& io = ImGui::GetIO();
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
    if (!done.exchange(false)) return;
    std::lock_guard<std::mutex> lk(resultM);
    GameUI_Toast(result.c_str());
    printf("bug report: %s\n", result.c_str());
}
