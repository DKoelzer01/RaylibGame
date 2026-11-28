// Tier 2 scenario runner (#870); see scenario.h.
#include "scenario.h"
#include "cliptelemetry.h"
#include "survival.h"
#include "fluid.h"
#include "player.h"
#include "playersession.h"
#include "recorder.h"
#include "scene.h"
#include "simclock.h"
#include "voxel.h"
#include "world.h"
#include "profiler.h"
#include "gputimer.h"

#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <map>
#include <atomic>
#include <mutex>
#include <set>
#include <thread>

extern const char* kGameVersion;
int ScenarioProc_Run(const std::string& cmdline, double timeoutSec);
// (#1054) A scenario's timeout is sim time (fixed dt): a loaded machine ran fluid/flow-meter at a third of real time and it
// timed out on the wall clock while simulating fine. Wall time is only a backstop for a run that stops advancing: this
// many times the timeout (RAYLIBGAME_SCENARIO_WALL_X).
static double WallFactor() { static const double f = getenv("RAYLIBGAME_SCENARIO_WALL_X") ? std::max(1.0, atof(getenv("RAYLIBGAME_SCENARIO_WALL_X"))) : 4.0; return f; }
std::string ScenarioProc_Exe();
int ScenarioProc_Pid();
int ScenarioReport_Main(int argc, char** argv); // scenarioreport.cpp
int ScenarioAB_Main(int argc, char** argv);     // scenarioab.cpp (#871)

namespace fs = std::filesystem;
void Scenario_LegacySetup();  // (below: the legacyDrives mode)
void Scenario_LegacyFinish();

namespace {

struct Entry { std::string name, file; ScenarioBuild build; };
std::vector<Entry>& Registry() { static std::vector<Entry> r; return r; }

ScenarioDef Built(const Entry& e) {
    ScenarioDef d;
    d.name = e.name; d.file = e.file;
    e.build(d);
    return d;
}
std::vector<ScenarioDef> AllDefs() {
    std::vector<ScenarioDef> out;
    for (const Entry& e : Registry()) out.push_back(Built(e));
    std::sort(out.begin(), out.end(), [](const ScenarioDef& a, const ScenarioDef& b) { return a.name < b.name; });
    return out;
}
// sel: a name, tag:TAG, all, or several of those separated by commas.
std::vector<ScenarioDef> Matching(const std::string& sel) {
    std::vector<std::string> parts;
    for (size_t a = 0; a <= sel.size();) { const size_t b = std::min(sel.find(',', a), sel.size()); if (b > a) parts.push_back(sel.substr(a, b - a)); a = b + 1; }
    std::vector<ScenarioDef> out;
    for (ScenarioDef& d : AllDefs())
        for (const std::string& p : parts)
            if (p == "all" || p == d.name || (p.rfind("tag:", 0) == 0 && d.hasTag(p.substr(4)))) { out.push_back(std::move(d)); break; }
    return out;
}

std::string JsonEsc(const std::string& s) {
    std::string o;
    for (unsigned char c : s) {
        if (c == '"' || c == '\\') { o += '\\'; o += (char)c; }
        else if (c == '\n') o += "\\n";
        else if (c == '\r') o += "\\r";
        else if (c == '\t') o += "\\t";
        else if (c < 0x20) { char b[8]; snprintf(b, sizeof b, "\\u%04x", c); o += b; }
        else o += (char)c;
    }
    return o;
}
std::string Num(double v) {
    if (!std::isfinite(v)) return "null";
    char b[64]; snprintf(b, sizeof b, "%.9g", v); return b;
}
std::string NowIso() {
    const std::time_t t = std::time(nullptr);
    char b[32]; std::strftime(b, sizeof b, "%Y-%m-%dT%H:%M:%S", std::localtime(&t)); return b;
}
std::string Stamp() {
    const std::time_t t = std::time(nullptr);
    char b[32]; std::strftime(b, sizeof b, "%Y%m%d-%H%M%S", std::localtime(&t)); return b;
}
std::string Quote(const std::string& s) { return s.find_first_of(" \t\"") == std::string::npos ? s : "\"" + s + "\""; }

// The active run (one scenario per process).
struct Run {
    bool active = false, finished = false, keep = false, record = true;
    ScenarioDef def;
    std::string outDir, runDir, scratch, capture, verdict = "RUNNING", reason, keepReason, waiting, args;
    std::vector<ScenarioCheck> checks;
    std::vector<ScenarioMetric> metrics;
    std::vector<std::string> notes;
    struct StepTime { std::string name; double sim = 0.0, wall = 0.0; bool done = false; };
    std::vector<StepTime> stepTimes;
    int step = 0;
    double simStart = -1.0, stepSim = 0.0, wallStart = 0.0, stepWall = 0.0, lastBeat = 0.0, startedAt = 0.0;
    long frames = 0;
    std::string started;
    std::string profile; // (#871) the profiler's summary, taken when the scenario finished
} R;

double Wall() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

void Finish(const std::string& verdict, const std::string& reason);

std::string ResultsJson() {
    std::string j = "{\n";
    j += "  \"schema\": 1,\n";
    j += "  \"scenario\": \"" + JsonEsc(R.def.name) + "\",\n";
    j += "  \"file\": \"" + JsonEsc(R.def.file) + "\",\n";
    j += "  \"tags\": [";
    for (size_t i = 0; i < R.def.tags.size(); ++i) j += (i ? ", \"" : "\"") + JsonEsc(R.def.tags[i]) + "\"";
    j += "],\n";
    j += "  \"legacy_flag\": \"" + JsonEsc(R.def.legacyFlag) + "\",\n";
    j += "  \"verdict\": \"" + R.verdict + "\",\n";
    j += "  \"reason\": \"" + JsonEsc(R.reason) + "\",\n";
    j += "  \"build\": \"" + JsonEsc(kGameVersion) + "\",\n";
    j += "  \"started\": \"" + R.started + "\",\n";
    j += "  \"seed\": " + std::to_string(R.def.seed) + ",\n";
    j += "  \"world\": \"" + std::string(R.def.world == ScenarioWorld::Sandbox ? "sandbox" : "planet") + "\",\n";
    j += "  \"args\": \"" + JsonEsc(R.args) + "\",\n";
    j += "  \"wall_seconds\": " + Num(R.wallStart > 0.0 ? Wall() - R.wallStart : 0.0) + ",\n";
    j += "  \"sim_seconds\": " + Num(R.simStart >= 0.0 ? SimTime() - R.simStart : 0.0) + ",\n";
    j += "  \"load_seconds\": " + Num(R.wallStart > 0.0 ? R.wallStart - R.startedAt : 0.0) + ",\n";
    j += "  \"steps\": [";
    for (size_t i = 0; i < R.stepTimes.size(); ++i)
        j += std::string(i ? ",\n    " : "\n    ") + "{\"name\": \"" + JsonEsc(R.stepTimes[i].name) + "\", \"sim_seconds\": " + Num(R.stepTimes[i].sim) +
             ", \"wall_seconds\": " + Num(R.stepTimes[i].wall) + ", \"done\": " + (R.stepTimes[i].done ? "true" : "false") + "}";
    j += R.stepTimes.empty() ? "],\n" : "\n  ],\n";
    j += "  \"checks\": [";
    for (size_t i = 0; i < R.checks.size(); ++i) {
        const ScenarioCheck& c = R.checks[i];
        j += std::string(i ? ",\n    " : "\n    ") + "{\"name\": \"" + JsonEsc(c.name) + "\", \"value\": " + Num(c.value) + ", \"op\": \"" + JsonEsc(c.op) +
             "\", \"expected\": \"" + JsonEsc(c.expected) + "\", \"pass\": " + (c.pass ? "true" : "false") + ", \"detail\": \"" + JsonEsc(c.detail) + "\"}";
    }
    j += R.checks.empty() ? "],\n" : "\n  ],\n";
    j += "  \"metrics\": [";
    for (size_t i = 0; i < R.metrics.size(); ++i) {
        const ScenarioMetric& m = R.metrics[i];
        j += std::string(i ? ",\n    " : "\n    ") + "{\"name\": \"" + JsonEsc(m.name) + "\", \"value\": " + Num(m.value) + ", \"unit\": \"" + JsonEsc(m.unit) +
             "\", \"lower_is_better\": " + (m.lowerIsBetter ? "true" : "false") + "}";
    }
    j += R.metrics.empty() ? "],\n" : "\n  ],\n";
    j += "  \"notes\": [";
    for (size_t i = 0; i < R.notes.size(); ++i) j += std::string(i ? ",\n    " : "\n    ") + "\"" + JsonEsc(R.notes[i]) + "\"";
    j += R.notes.empty() ? "],\n" : "\n  ],\n";
    const bool kept = R.record && (R.verdict != "PASS" || R.keep);
    j += "  \"capture\": {\"recorded\": " + std::string(R.record ? "true" : "false") + ", \"kept\": " + (kept ? "true" : "false") +
         ", \"dir\": \"" + JsonEsc(kept ? R.capture : "") + "\", \"sheet\": \"" + JsonEsc(kept && fs::exists(R.capture + "/sheet.png") ? R.capture + "/sheet.png" : "") +
         "\", \"reason\": \"" + JsonEsc(R.keepReason) + "\"}";
    const std::string prof = R.finished ? R.profile : Profiler_SummaryJson(); // (#871) the run summary; "" when off
    if (!prof.empty()) j += ",\n  \"profile\": " + prof;
    j += "\n}\n";
    return j;
}
void WriteResults() {
    if (R.runDir.empty()) return;
    std::error_code ec;
    fs::create_directories(R.runDir, ec);
    const std::string path = R.runDir + "/results.json", tmp = path + ".tmp";
    if (FILE* f = fopen(tmp.c_str(), "wb")) { const std::string j = ResultsJson(); fwrite(j.data(), 1, j.size(), f); fclose(f); fs::rename(tmp, path, ec); }
}

void Finish(const std::string& verdict, const std::string& reason) {
    if (R.finished) return;
    if (R.step < (int)R.stepTimes.size()) { R.stepTimes[R.step].sim = SimTime() - R.stepSim; R.stepTimes[R.step].wall = Wall() - R.stepWall; }
    {   // (#1442) frames whose view was not normal gameplay: the camera or the player inside solid
        const ClipStats& k = ClipTelemetry_Stats();
        R.metrics.push_back({ "camera_inside_solid_frames", "frames", (double)k.cameraFrames, true });
        R.metrics.push_back({ "player_inside_solid_frames", "frames", (double)k.playerFrames, true });
        if (k.cameraFrames) R.notes.push_back(TextFormat("the camera was inside solid for %ld frames (%.2f s, t=%.2f-%.2f s, %.1f m deep at most)",
                                                         k.cameraFrames, k.cameraSeconds, k.cameraFirst, k.cameraLast, k.cameraDepthMax));
        if (k.playerFrames) R.notes.push_back(TextFormat("the player was inside solid for %ld frames (%.2f s, t=%.2f-%.2f s, %.1f m deep at most)",
                                                         k.playerFrames, k.playerSeconds, k.playerFirst, k.playerLast, k.playerDepthMax));
        if (R.def.expectNoClip && k.cameraFrames + k.playerFrames > 0 && verdict == "PASS") {
            R.checks.push_back({ "no frame with the camera or the player inside solid", "==", "0", std::to_string(k.cameraFrames) + " camera, " + std::to_string(k.playerFrames) + " player frames",
                                 (double)(k.cameraFrames + k.playerFrames), false });
            R.finished = false; Finish("FAIL", "the camera or the player was inside solid (expectNoClip)"); return;
        }
    }
    R.finished = true;
    R.verdict = verdict; R.reason = reason;
    int failed = 0; for (const auto& c : R.checks) failed += c.pass ? 0 : 1;
    printf("scenario: %s -> %s%s%s (%zu checks, %d failed; %.1f s sim, %.1f s wall)\n", R.def.name.c_str(), verdict.c_str(), reason.empty() ? "" : ": ", reason.c_str(),
           R.checks.size(), failed, SimTime() - R.simStart, Wall() - R.wallStart);
    // (#871) The profile covers the scenario, not the bench frames after it.
    R.profile = Profiler_SummaryJson();
    Profiler_PrintSummary();
    Profiler_Enable(false);
    WriteResults();
}

// (#1422) A hang the frame-driven backstops can't see: the main loop blocked inside a step, or the run held in a loading
// screen (Scenario_Frame runs only out of it). The main thread touches this clock each frame (Scenario_Frame,
// Scenario_Loading); a watchdog thread ends the run with a clear TIMEOUT when it stops for the scenario's wall backstop (at
// least 120 s), naming the step and the last note, instead of leaving the parent to kill a silent process.
std::atomic<double> gLastTouch{ 0.0 };
std::mutex gWatchMx; // the watchdog's message vs the main thread's notes
void StartWatchdog() {
    static bool started = false;
    if (started) return;
    started = true;
    gLastTouch = Wall();
    std::thread([] {
        const double limit = std::max(120.0, R.def.timeout * WallFactor());
        for (;;) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            if (R.finished) return;
            const double idle = Wall() - gLastTouch.load();
            if (idle < limit) continue;
            std::lock_guard<std::mutex> lock(gWatchMx);
            const std::string step = R.step < (int)R.def.steps.size() ? R.def.steps[R.step].name : std::string("?");
            const std::string last = R.notes.empty() ? std::string("none") : R.notes.back();
            Finish("TIMEOUT", "the main loop stopped in step '" + step + "' (no frame for " + Num(idle) + " s; last note: " + last + ")");
            fflush(stdout);
            std::_Exit(3);
        }
    }).detach();
}

void SetEnv(const std::string& k, const std::string& v) {
#ifdef _WIN32
    _putenv_s(k.c_str(), v.c_str());
#else
    setenv(k.c_str(), v.c_str(), 1);
#endif
}

// argv storage for the rewritten command line (lives for the process).
std::vector<std::string> gArgStore;
std::vector<char*> gArgPtrs;

int RunMany(const std::vector<ScenarioDef>& defs, const std::string& outDir, const std::vector<std::string>& pass) {
    const std::string exe = ScenarioProc_Exe();
    std::error_code ec;
    fs::create_directories(outDir, ec);
    int passed = 0, failed = 0;
    std::string combined = "{\n  \"schema\": 1,\n  \"started\": \"" + NowIso() + "\",\n  \"build\": \"" + JsonEsc(kGameVersion) + "\",\n  \"runs\": [";
    bool first = true;
    for (const ScenarioDef& d : defs) {
        std::string cmd = Quote(exe) + " --scenario " + Quote(d.name) + " --scenario-out " + Quote(outDir);
        for (const std::string& a : pass) cmd += " " + Quote(a);
        printf("scenario: === %s\n", d.name.c_str());
        fflush(stdout);
        const auto t0 = std::chrono::steady_clock::now();
        // The scenario's own wall backstop covers its run (WallFactor x its sim timeout, #1054); loading (a planet's first bake
        // ~80 s) and the exit get 240 s more.
        const int rc = ScenarioProc_Run(cmd, d.timeout * WallFactor() + 240.0);
        const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::string safe = d.name;
        for (char& c : safe) if (c == ':') c = '_'; // (the name only: the out dir has a drive letter)
        const std::string dir = outDir + "/" + safe;
        std::string j;
        if (FILE* f = fopen((dir + "/results.json").c_str(), "rb")) { char b[4096]; size_t n; while ((n = fread(b, 1, sizeof b, f)) > 0) j.append(b, n); fclose(f); }
        const bool ok = j.find("\"verdict\": \"PASS\"") != std::string::npos;
        if (j.empty()) {
            // No results: the process died or was ended before the scenario finished (a crash, a hang in loading).
            j = "{\"schema\": 1, \"scenario\": \"" + JsonEsc(d.name) + "\", \"verdict\": \"" + (rc == -2 ? "TIMEOUT" : "ERROR") + "\", \"reason\": \"" +
                (rc == -2 ? "the process was ended after " + Num(secs) + " s" : "no results.json (exit code " + std::to_string(rc) + ")") +
                "\", \"tags\": [], \"checks\": [], \"metrics\": [], \"notes\": [], \"capture\": {\"kept\": false}}";
        }
        while (!j.empty() && (j.back() == '\n' || j.back() == ' ')) j.pop_back();
        combined += (first ? "\n" : ",\n") + j;
        first = false;
        (ok ? passed : failed)++;
        printf("scenario: === %s %s in %.0f s\n", d.name.c_str(), ok ? "PASS" : "FAIL", secs);
    }
    combined += "\n  ],\n  \"passed\": " + std::to_string(passed) + ",\n  \"failed\": " + std::to_string(failed) + "\n}\n";
    const std::string path = outDir + "/results.json";
    if (FILE* f = fopen(path.c_str(), "wb")) { fwrite(combined.data(), 1, combined.size(), f); fclose(f); }
    printf("scenario: %d passed, %d failed; results %s\n", passed, failed, path.c_str());
    return failed ? 1 : 0;
}

// --scenario-parity NAME [TICKS] (#968): the fluid's CPU and GPU paths, the same scenario run twice in lockstep (GPU, then
// FLUID_CPU=1), each writing a per-tick hash of every brick's cells (RAYLIBGAME_FLUID_TRACE_FILE) and ending after TICKS
// traced ticks; PASS when the two traces are identical tick for tick, else the first tick they part. A lockstep scenario only
// (one whose s.env leaves RAYLIBGAME_FLUID_LOCKSTEP on).
int RunParity(const ScenarioDef& d, const std::string& outDir, const std::vector<std::string>& pass, long long ticks) {
    const std::string exe = ScenarioProc_Exe();
    std::error_code ec;
    fs::create_directories(outDir, ec);
    const char* modes[2] = { "gpu", "cpu" };
    std::vector<std::string> lines[2];
    for (int m = 0; m < 2; ++m) {
        const std::string trace = outDir + "/parity-" + modes[m] + ".txt";
        fs::remove(trace, ec);
        SetEnv("FLUID_CPU", m ? "1" : "0");
        SetEnv("RAYLIBGAME_FLUID_TRACE", "1");
        SetEnv("RAYLIBGAME_FLUID_TRACE_FILE", trace);
        SetEnv("RAYLIBGAME_FLUID_TRACE_STOP", std::to_string(ticks));
        std::string cmd = Quote(exe) + " --scenario " + Quote(d.name) + " --scenario-out " + Quote(outDir + "/parity-" + modes[m]) + " --scenario-no-record";
        for (const std::string& a : pass) cmd += " " + Quote(a);
        printf("scenario: === parity %s (%s)\n", d.name.c_str(), modes[m]);
        fflush(stdout);
        const int rc = ScenarioProc_Run(cmd, d.timeout * 6.0 + 240.0); // lockstep on the CPU path runs well under real time
        if (FILE* f = fopen(trace.c_str(), "rb")) {
            char b[256];
            while (fgets(b, sizeof b, f)) { std::string s = b; while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back(); lines[m].push_back(s); }
            fclose(f);
        }
        printf("scenario: parity %s (%s): %zu ticks traced (exit %d)\n", d.name.c_str(), modes[m], lines[m].size(), rc);
    }
    SetEnv("FLUID_CPU", "0"); SetEnv("RAYLIBGAME_FLUID_TRACE", "0"); SetEnv("RAYLIBGAME_FLUID_TRACE_FILE", ""); SetEnv("RAYLIBGAME_FLUID_TRACE_STOP", "0");
    const size_t n = std::min(lines[0].size(), lines[1].size());
    size_t diff = n;
    for (size_t i = 0; i < n; ++i) if (lines[0][i] != lines[1][i]) { diff = i; break; }
    std::string verdict, reason;
    if (n == 0) { verdict = "ERROR"; reason = "no ticks traced (is the scenario lockstep, and did both runs reach the fluid?)"; }
    else if (diff < n) { verdict = "FAIL"; reason = "the CPU and GPU paths part at traced tick " + std::to_string(diff + 1) + " of " + std::to_string(n) + ": gpu '" + lines[0][diff] + "', cpu '" + lines[1][diff] + "'"; }
    else if (n < (size_t)ticks) { verdict = "FAIL"; reason = "identical for all " + std::to_string(n) + " common ticks, but a run ended early (gpu " + std::to_string(lines[0].size()) + ", cpu " + std::to_string(lines[1].size()) + " of " + std::to_string(ticks) + ")"; }
    else { verdict = "PASS"; reason = "identical for all " + std::to_string(n) + " traced ticks"; }
    printf("scenario: parity %s -> %s: %s\n", d.name.c_str(), verdict.c_str(), reason.c_str());
    const std::string j = "{\n  \"schema\": 1,\n  \"scenario\": \"parity:" + JsonEsc(d.name) + "\",\n  \"verdict\": \"" + verdict + "\",\n  \"reason\": \"" + JsonEsc(reason) +
                          "\",\n  \"ticks\": " + std::to_string(n) + ",\n  \"first_difference\": " + (diff < n ? std::to_string(diff + 1) : std::string("null")) + "\n}\n";
    if (FILE* f = fopen((outDir + "/parity-results.json").c_str(), "wb")) { fwrite(j.data(), 1, j.size(), f); fclose(f); }
    return verdict == "PASS" ? 0 : 1;
}

} // namespace

ScenarioRegistrar::ScenarioRegistrar(const char* name, const char* file, ScenarioBuild build) {
    std::string f = file ? file : "";
    const size_t at = f.find("tests/scenarios/") != std::string::npos ? f.find("tests/scenarios/") : f.find("tests\\scenarios\\");
    if (at != std::string::npos) f = f.substr(at);
    Registry().push_back({ name, f, build });
}

bool Scenario_Active() { return R.active; }

void Scenario_PrepareArgs(int& argc, char**& argv) {
    std::string sel, outDir;
    bool list = false, parity = false;
    long long parityTicks = 400;
    std::vector<std::string> rest; // everything else, passed on
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--scenario-report") { std::exit(ScenarioReport_Main(argc - i - 1, argv + i + 1)); }
        else if (a == "--tier-report" && i + 1 < argc) { std::exit(QuestBook_TierReport(argv[i + 1])); } // (#1199)
        else if (a == "--scenario-ab") { std::exit(ScenarioAB_Main(argc - i - 1, argv + i + 1)); } // (#871) scenarioab.cpp
        else if (a == "--scenario" && i + 1 < argc) sel = argv[++i];
        else if (a == "--scenario-parity" && i + 1 < argc) { // (#968) NAME [TICKS]: the fluid's CPU and GPU paths compared
            parity = true; sel = argv[++i];
            if (i + 1 < argc && argv[i + 1][0] >= '0' && argv[i + 1][0] <= '9') parityTicks = std::max(1ll, atoll(argv[++i]));
        }
        else if (a == "--scenario-out" && i + 1 < argc) outDir = argv[++i];
        else if (a == "--scenario-list") { list = true; if (i + 1 < argc && argv[i + 1][0] != '-') sel = argv[++i]; }
        else if (a == "--scenario-keep") { R.keep = true; R.keepReason = "--scenario-keep"; rest.push_back(a); }
        else if (a == "--scenario-no-record") { R.record = false; rest.push_back(a); }
        else rest.push_back(a);
    }
    if (list) {
        for (const ScenarioDef& d : sel.empty() ? AllDefs() : Matching(sel)) {
            std::string tags;
            for (const auto& t : d.tags) tags += (tags.empty() ? "" : ",") + t;
            printf("%-28s %-8s seed %-3u timeout %4.0f s  tags %-22s %s%s\n", d.name.c_str(), d.world == ScenarioWorld::Sandbox ? "sandbox" : "planet", d.seed, d.timeout,
                   tags.c_str(), d.file.c_str(), d.legacyFlag.empty() ? "" : ("  (was " + d.legacyFlag + ")").c_str());
        }
        std::exit(0);
    }
    if (sel.empty()) return;
    if (outDir.empty()) {
        std::error_code ec;
        outDir = (fs::temp_directory_path(ec) / "raylibgame_scenarios" / Stamp()).generic_string();
    }
    {
        std::error_code ec;
        outDir = fs::absolute(outDir, ec).generic_string(); // the recorder takes an absolute path as is (a relative one goes under the share)
    }
    setvbuf(stdout, nullptr, _IONBF, 0); // a run ended by a timeout still shows every line
    const std::vector<ScenarioDef> defs = Matching(sel);
    if (defs.empty()) { printf("scenario: no scenario matches '%s' (--scenario-list lists them)\n", sel.c_str()); std::exit(2); }
    if (parity) std::exit(RunParity(defs[0], outDir, rest, parityTicks)); // (#968) one scenario, its CPU and GPU runs compared
    if (sel.rfind("tag:", 0) == 0 || sel == "all" || defs.size() > 1) std::exit(RunMany(defs, outDir, rest));

    // One scenario, in this process.
    R.def = defs[0];
    if (!R.def.record) R.record = false; // (#871) perf scenarios: no recording
    R.active = true;
    R.outDir = outDir;
    std::string safe = R.def.name;
    for (char& c : safe) if (c == ':') c = '_';
    R.runDir = outDir + "/" + safe;
    R.scratch = R.runDir + "/scratch";
    R.capture = R.runDir + "/capture";
    R.started = NowIso();
    R.startedAt = Wall();
    for (const ScenarioStep& st : R.def.steps) R.stepTimes.push_back({ st.name, 0.0, 0.0, false });
    std::error_code ec;
    fs::remove_all(R.runDir, ec);
    fs::create_directories(R.scratch, ec);
    // Deterministic by default: lockstep fluids (#680/#782: the GPU's dispatches are read back as they end), no bench warmup.
    SetEnv("RAYLIBGAME_FLUID_LOCKSTEP", "1");
    SetEnv("RAYLIBGAME_BENCH_WARMUP", "0");
    for (const auto& kv : R.def.env) SetEnv(kv.first, kv.second);

    // RAYLIBGAME_SCENARIO_SEED: another world (a place a playtest saw, its save's seed); results then aren't comparable.
    const char* seedEnv = getenv("RAYLIBGAME_SCENARIO_SEED");
    const std::string seed = seedEnv && *seedEnv ? std::string(seedEnv) : std::to_string(R.def.seed);
    std::vector<std::string> a = { argv[0], "--bench", "0.1", "--seed", seed, "--preload", "0", "--cache-verts", "1", "--window", "960", "540" };
    if (R.def.world == ScenarioWorld::Sandbox) a.push_back("--sandbox");
    if (R.record) {
        // Low-res and every 4th frame, to the run's own folder (local temp): deleted on PASS (scenario.h). Implies a fixed
        // 1/60 s frame time; the contact sheet is what the report shows.
        const int frames = std::max(30, (int)(R.def.timeout * 60.0 / 4.0));
        for (const char* x : { "--record", R.capture.c_str(), "--record-every", "4", "--record-continue", "--sheet", "--record-hud", "0" }) a.push_back(x);
        a.push_back("--frames"); a.push_back(std::to_string(frames));
    } else { a.push_back("--fixed-dt"); a.push_back("0.0166666667"); }
    for (const std::string& x : R.def.args) a.push_back(x);
    for (const std::string& x : rest) a.push_back(x);
    for (size_t i = 1; i < a.size(); ++i) R.args += (i > 1 ? " " : "") + a[i];
    gArgStore = a;
    gArgPtrs.clear();
    for (std::string& s : gArgStore) gArgPtrs.push_back(s.data());
    gArgPtrs.push_back(nullptr);
    argc = (int)gArgStore.size();
    argv = gArgPtrs.data();
    printf("scenario: %s (%s) -> %s\n", R.def.name.c_str(), R.def.file.c_str(), R.runDir.c_str());
    printf("scenario: args %s\n", R.args.c_str());
    // (#871) The profiler records every scenario (RAYLIBGAME_PROFILE=0: off, the A/B overhead check's other side).
    if (const char* e = getenv("RAYLIBGAME_PROFILE"); !(e && e[0] == '0')) Profiler_Enable(true);
    WriteResults(); // RUNNING: a crash leaves this behind
}

void Scenario_Setup() {
    if (!R.active) return;
    // The save in the run's scratch dir (a sandbox or bench save shared between processes made parallel runs flaky). A
    // planet's erosion bake is cached beside its save: the bench save's cache is borrowed when there is one (a stale one
    // is rebaked: the cache is keyed on its parameters).
    const std::string save = R.scratch + "/world.sav";
    std::error_code ec;
    const char* shared = R.def.world == ScenarioWorld::Sandbox ? "saves/sandbox.sav.erosion" : "saves/bench.sav.erosion";
    if (fs::exists(shared, ec)) fs::copy_file(shared, save + ".erosion", fs::copy_options::overwrite_existing, ec);
    SetWorldSavePath(save.c_str());
    Scenario_LegacySetup(); // (legacyDrives: which bench*.png exist before the run)
}

bool Scenario_Frame(Scene& world) {
    if (!R.active || R.finished) return false;
    const double now = Wall();
    gLastTouch = now; // (#1422) the watchdog's clock
    if (R.simStart < 0.0) {
        R.simStart = SimTime(); R.stepSim = SimTime(); R.wallStart = now; R.stepWall = now; R.lastBeat = now;
        printf("scenario: %s: world ready after %.1f s; %zu steps\n", R.def.name.c_str(), now - R.startedAt, R.def.steps.size());
        ClipTelemetry_Reset(); // (#1442) the scenario's own frames, not the loading screen's
        StartWatchdog(); // (#1422)
        if (Profiler_On()) { GpuTimer_KeepOn(true); Profiler_Reset(); } // (#871) the passes as cpu:/gpu: zones; the summary from here on (not the loading)
    }
    ScenarioCtx c{ world };
    c.body = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get());
    c.field = c.body ? dynamic_cast<const PlanetField*>(c.body->field.get()) : nullptr;
    c.T = SimTime() - R.simStart;
    c.t = SimTime() - R.stepSim;
    c.wall = now - R.stepWall;
    c.clock = R.def.realTime ? std::min(c.t, c.wall) : c.t;
    c.frame = R.frames++;
    c.stepIndex = R.step;
    if (!c.body || (!c.field && R.def.needsPlanetField)) { Finish("ERROR", c.body ? "no planet field (body 0)" : "no body 0"); return false; }
    if (R.def.legacyDrives) { // a legacy bench drives: never hold it; its shots are compared when the game ends
        if (now - R.wallStart > R.def.timeout) Finish("TIMEOUT", "the legacy bench (" + R.args + ") ran past the scenario's " + Num(R.def.timeout) + " s");
        return false;
    }
    const double runClock = R.def.realTime ? std::min(c.T, now - R.wallStart) : c.T; // (#1972) as c.clock
    if (runClock > R.def.timeout || now - R.wallStart > R.def.timeout * WallFactor()) { // sim time; wall time a backstop (#1054)
        const bool wall = runClock <= R.def.timeout;
        Finish("TIMEOUT", "step '" + R.def.steps[R.step].name + "' still running after the scenario's " + Num(R.def.timeout) + " s" +
                              (wall ? " x " + Num(WallFactor()) + " of wall time (sim only " + Num(c.T) + " s: the run stopped advancing)" : " of sim time") +
                              (R.waiting.empty() ? "" : " (waiting for " + R.waiting + ")"));
        return false;
    }
    if (now - R.lastBeat >= 5.0) {
        R.lastBeat = now;
        printf("scenario: %s heartbeat: step %d/%zu '%s' t=%.1f s (sim %.1f s, wall %.0f s)%s\n", R.def.name.c_str(), R.step + 1, R.def.steps.size(),
               R.def.steps[R.step].name.c_str(), c.t, c.T, now - R.wallStart, R.waiting.empty() ? "" : (", waiting for " + R.waiting).c_str());
    }
    if (Profiler_On()) { // (#871) named counters, once a frame
        const Fluid::Stats fs = Fluid::GetStats();
        Profiler_Counter("fluid bricks awake", fs.activeBricks);
        Profiler_Counter("fluid volume m^3", fs.volume);
        Profiler_Counter("fluid step ms", fs.stepMs);
        const WorldStats& ws = GetWorldStats();
        Profiler_Counter("chunks alive", ws.nodes);
        Profiler_Counter("chunk jobs queued", ws.queuedJobs);
    }
    const ScenarioStep& st = R.def.steps[R.step];
    const bool done = st.run(c);
    if (R.finished) return false; // c.fail()
    if (done) {
        R.stepTimes[R.step] = { st.name, SimTime() - R.stepSim, now - R.stepWall, true };
        printf("scenario: %s: step '%s' done after %.1f s sim\n", R.def.name.c_str(), st.name.c_str(), SimTime() - R.stepSim);
        R.waiting.clear();
        R.stepSim = SimTime(); R.stepWall = now;
        if (++R.step >= (int)R.def.steps.size()) {
            int failed = 0; for (const auto& k : R.checks) failed += k.pass ? 0 : 1;
            Finish(failed || R.checks.empty() ? "FAIL" : "PASS", failed ? std::to_string(failed) + " of " + std::to_string(R.checks.size()) + " checks failed" : R.checks.empty() ? "no checks" : "");
            return false;
        }
    } else if (c.clock > st.deadline) {
        Finish("FAIL", "step '" + st.name + "' passed its deadline of " + Num(st.deadline) + (R.def.realTime ? " s of both sim and wall time" : " s sim") +
                           (R.waiting.empty() ? "" : " (waiting for " + R.waiting + ")"));
        return false;
    }
    return true;
}

// (#1422) Each frame of a loading screen once the world was ready: a heartbeat, and the wall backstop (Scenario_Frame doesn't
// run then, so a run held in a loading screen ended only when the parent killed it, with nothing said).
void Scenario_Loading() {
    if (!R.active || R.finished || R.simStart < 0.0) return;
    const double now = Wall();
    gLastTouch = now;
    // A legacy scenario (a bench drives it) has no steps: don't index past the list.
    const std::string stepName = R.step >= 0 && R.step < (int)R.def.steps.size() ? R.def.steps[R.step].name : std::string("(legacy bench)");
    if (now - R.lastBeat >= 5.0) {
        R.lastBeat = now;
        printf("scenario: %s heartbeat: a loading screen during step %d/%zu '%s' (wall %.0f s)\n", R.def.name.c_str(), R.step + 1, R.def.steps.size(),
               stepName.c_str(), now - R.wallStart);
    }
    if (now - R.wallStart > R.def.timeout * WallFactor()) {
        Finish("TIMEOUT", "a loading screen held step '" + stepName + "' past the scenario's wall backstop (" + Num(R.def.timeout * WallFactor()) + " s)");
        fflush(stdout);
        std::_Exit(3);
    }
}

void Scenario_Shutdown() {
    if (!R.active) return;
    Scenario_LegacyFinish(); // (legacyDrives: the shots against their baselines)
    if (!R.finished) Finish("ERROR", "the game ended before the scenario finished");
    std::error_code ec;
    // Give the planet's erosion bake back to the bench save, so the next run doesn't bake again.
    const std::string bake = R.scratch + "/world.sav.erosion";
    const char* shared = R.def.world == ScenarioWorld::Sandbox ? "saves/sandbox.sav.erosion" : "saves/bench.sav.erosion";
    if (fs::exists(bake, ec) && !fs::exists(shared, ec)) { fs::create_directories("saves", ec); fs::copy_file(bake, shared, ec); }
    fs::remove_all(R.scratch, ec);
    const bool keep = R.record && (R.verdict != "PASS" || R.keep);
    if (R.record && !keep) fs::remove_all(R.capture, ec);
    WriteResults();
    printf("scenario: results %s/results.json%s\n", R.runDir.c_str(), keep ? (", capture kept: " + R.capture).c_str() : "");
}

// ---- ScenarioCtx ----
namespace {
bool Record(const std::string& name, double v, const std::string& op, const std::string& expected, bool pass, const std::string& detail = "") {
    R.checks.push_back({ name, op, expected, detail, v, pass });
    printf("scenario: check %s = %s (want %s %s) %s%s%s\n", name.c_str(), Num(v).c_str(), op.c_str(), expected.c_str(), pass ? "ok" : "FAIL",
           detail.empty() ? "" : ": ", detail.c_str());
    return pass;
}
}
bool ScenarioCtx::check(const std::string& name, bool ok, const std::string& detail) { return Record(name, ok ? 1.0 : 0.0, "is", "true", ok, detail); }
bool ScenarioCtx::eq(const std::string& name, double v, double want, double tol) {
    return Record(name, v, "==", tol > 0.0 ? Num(want) + " +- " + Num(tol) : Num(want), std::fabs(v - want) <= tol);
}
bool ScenarioCtx::near(const std::string& name, double v, double want, double tol) { return Record(name, v, "~", Num(want) + " +- " + Num(tol), std::fabs(v - want) <= tol); }
bool ScenarioCtx::lt(const std::string& name, double v, double bound) { return Record(name, v, "<", Num(bound), v < bound); }
bool ScenarioCtx::timingLt(const std::string& name, double v, double bound) {
    // (#1034) Timing is judged on the quiet PC only: where RAYLIBGAME_TIMING_GATES=0 (the shared GPU runner) the value is
    // reported, with the bound in the note, and never fails the scenario.
    static const bool gates = !(getenv("RAYLIBGAME_TIMING_GATES") && atoi(getenv("RAYLIBGAME_TIMING_GATES")) == 0);
    if (gates) return lt(name, v, bound);
    metric(name, v, "", true);
    note(name + ": " + Num(v) + " (gate < " + Num(bound) + " not judged here: RAYLIBGAME_TIMING_GATES=0)");
    return true;
}
bool ScenarioCtx::gt(const std::string& name, double v, double bound) { return Record(name, v, ">", Num(bound), v > bound); }
bool ScenarioCtx::between(const std::string& name, double v, double lo, double hi) { return Record(name, v, "in", "[" + Num(lo) + ", " + Num(hi) + "]", v >= lo && v <= hi); }
void ScenarioCtx::metric(const std::string& name, double v, const std::string& unit, bool lowerIsBetter) {
    R.metrics.push_back({ name, unit, v, lowerIsBetter });
    printf("scenario: metric %s = %s %s\n", name.c_str(), Num(v).c_str(), unit.c_str());
}
void ScenarioCtx::note(const std::string& text) { { std::lock_guard<std::mutex> lock(gWatchMx); R.notes.push_back(text); } printf("scenario: %s\n", text.c_str()); }
void ScenarioCtx::fail(const std::string& reason) { Finish("FAIL", reason); }
void ScenarioCtx::waitingFor(const std::string& what) { R.waiting = what; }
void ScenarioCtx::keepCapture(const std::string& reason) { R.keep = true; R.keepReason = reason; }
void ScenarioCtx::cue() { Recorder_Cue(); }

Vector3 ScenarioCtx::spawnDirLocal() const {
    Vector3 d{};
    if (const char* e = getenv("RAYLIBGAME_SPAWN_DIR"); e && sscanf(e, "%f,%f,%f", &d.x, &d.y, &d.z) == 3 && Vector3Length(d) > 1e-6f) return Vector3Normalize(d);
    return field && body ? field->findSpawnDir(body->seed) : Vector3{ 0, 1, 0 };
}
Vector3 ScenarioCtx::spawnDirWorld() const { return body ? body->dirToWorld(spawnDirLocal()) : Vector3{ 0, 1, 0 }; }
void ScenarioCtx::camera(Vector3 eyeLocal, Vector3 targetLocal, Vector3 upLocal) {
    if (!body) return;
    world.camera.position = body->toWorld(eyeLocal);
    world.camera.target = body->toWorld(targetLocal);
    world.camera.up = body->dirToWorld(upLocal);
    player.flyMode = true;
    player.position = world.renderOrigin + world.camera.position;
}

// ---- Visual checks against baselines (#870 visual tier) ----
std::string Scenario_RunDir() { return R.runDir; }
std::string Scenario_BaselineDir() {
    const char* e = getenv("RAYLIBGAME_BASELINES");
    return e && *e ? e : "//DEWEYSERVER/files/claude-media/raylib-baselines";
}
namespace {
// Compares the image file `shotPath` with the scenario's baseline NAME.png (both shrunk to 320x180): records a check (the %
// of pixels that differ by more than 24 in a channel, under `tolPct`) and a metric (the mean difference). No baseline (or
// RAYLIBGAME_BASELINE_UPDATE=1): the shot becomes it. On a failure a difference image goes next to the shot.
bool CompareShot(const std::string& name, const std::string& shotPath, double tolPct) {
    std::error_code ec;
    std::string scen = R.def.name;
    for (char& ch : scen) if (ch == '/' || ch == ':' || ch == '\\') ch = '_';
    const std::string dir = Scenario_BaselineDir() + "/" + scen, base = dir + "/" + name;
    const bool update = getenv("RAYLIBGAME_BASELINE_UPDATE") && atoi(getenv("RAYLIBGAME_BASELINE_UPDATE")) != 0;
    if (!fs::exists(shotPath, ec)) return Record("shot " + name, 0.0, "is", "taken", false, "no image at " + shotPath);
    if (update || !fs::exists(base, ec)) {
        fs::create_directories(dir, ec);
        fs::copy_file(shotPath, base, fs::copy_options::overwrite_existing, ec);
        const std::string why = update ? "updated (RAYLIBGAME_BASELINE_UPDATE)" : "recorded (there was none)";
        R.notes.push_back("baseline " + name + " " + why + ": " + base);
        printf("scenario: baseline %s %s: %s\n", name.c_str(), why.c_str(), base.c_str());
        return Record("shot " + name, 0.0, "<=", Num(tolPct) + " % differing", !ec, ec ? "could not write the baseline: " + ec.message() : "baseline " + why);
    }
    Image a = LoadImage(shotPath.c_str()), b = LoadImage(base.c_str());
    if (!a.data || !b.data) { if (a.data) UnloadImage(a); if (b.data) UnloadImage(b); return Record("shot " + name, 100.0, "<=", Num(tolPct) + " % differing", false, "could not read the shot or its baseline"); }
    constexpr int W = 320, H = 180;
    ImageFormat(&a, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8); ImageFormat(&b, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    ImageResize(&a, W, H); ImageResize(&b, W, H);
    const Color* pa = (const Color*)a.data; const Color* pb = (const Color*)b.data;
    int differ = 0; double sum = 0.0;
    Image d = GenImageColor(W, H, BLACK);
    Color* pd = (Color*)d.data;
    for (int i = 0; i < W * H; ++i) {
        const int dr = abs(pa[i].r - pb[i].r), dg = abs(pa[i].g - pb[i].g), db = abs(pa[i].b - pb[i].b);
        const int m = std::max(dr, std::max(dg, db));
        sum += (dr + dg + db) / 3.0;
        if (m > 24) ++differ;
        pd[i] = Color{ (unsigned char)std::min(255, dr * 4), (unsigned char)std::min(255, dg * 4), (unsigned char)std::min(255, db * 4), 255 };
    }
    const double pct = 100.0 * differ / (W * H);
    R.metrics.push_back({ "shot " + name + " mean difference", "/255", sum / (W * H), true });
    const bool ok = Record("shot " + name, pct, "<=", Num(tolPct) + " % differing", pct <= tolPct, "against " + base);
    if (!ok) {
        const std::string diffPath = shotPath.substr(0, shotPath.size() - 4) + "_diff.png";
        ExportImage(d, diffPath.c_str());
        R.keep = true; R.keepReason = "shot " + name + " differs from its baseline (" + diffPath + ")";
    }
    UnloadImage(a); UnloadImage(b); UnloadImage(d);
    return ok;
}
// (legacyDrives) the bench_*.png files in the working directory and their write times when the run began
std::map<std::string, fs::file_time_type> pngsBefore;
void SnapshotPngs() {
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(fs::current_path(ec), ec)) {
        const std::string n = e.path().filename().string();
        if (n.rfind("bench", 0) == 0 && e.path().extension() == ".png") pngsBefore[n] = fs::last_write_time(e.path(), ec);
    }
}
void CompareLegacyShots() {
    std::error_code ec;
    std::vector<std::string> names = R.def.shots;
    if (names.empty())
        for (const auto& e : fs::directory_iterator(fs::current_path(ec), ec)) {
            const std::string n = e.path().filename().string();
            if (n.rfind("bench", 0) != 0 || e.path().extension() != ".png") continue;
            auto it = pngsBefore.find(n);
            if (it == pngsBefore.end() || fs::last_write_time(e.path(), ec) > it->second) names.push_back(n);
        }
    std::sort(names.begin(), names.end());
    fs::create_directories(R.runDir + "/shots", ec);
    for (const std::string& n : names) {
        const std::string kept = R.runDir + "/shots/" + n;
        fs::copy_file(n, kept, fs::copy_options::overwrite_existing, ec);
        CompareShot(n, ec ? n : kept, R.def.shotTolerance);
    }
    if (names.empty()) Record("screenshots written", 0.0, ">", "0", false, "the legacy bench wrote no bench*.png");
}
}
bool ScenarioCtx::shot(const std::string& name, double tolerancePct) {
    std::error_code ec;
    fs::create_directories(R.runDir + "/shots", ec);
    const std::string file = name.size() > 4 && name.substr(name.size() - 4) == ".png" ? name : name + ".png";
    const std::string path = R.runDir + "/shots/" + file;
    Image im = LoadImageFromScreen();
    ExportImage(im, path.c_str());
    UnloadImage(im);
    return CompareShot(file, path, tolerancePct);
}
void ScenarioCtx::grab(const std::string& name) {
    const char* dir = getenv("RAYLIBGAME_SLIDE_SHOTS");
    if (!dir || !*dir) return;
    static std::set<std::string> taken; // once per name: steps call it every frame
    if (!taken.insert(name).second) return;
    std::string base = R.def.name;
    for (char& ch : base) if (ch == '/' || ch == ':') ch = '_';
    std::error_code ec;
    fs::create_directories(dir, ec);
    const std::string path = std::string(dir) + "/" + base + "-" + name + ".png";
    Image im = LoadImageFromScreen();
    const bool ok = ExportImage(im, path.c_str());
    UnloadImage(im);
    note("slide shot " + path + (ok ? "" : " FAILED"));
}
void Scenario_LegacySetup() { if (R.active && R.def.legacyDrives) SnapshotPngs(); }
void Scenario_LegacyFinish() {
    if (!R.active || !R.def.legacyDrives || R.finished) return;
    CompareLegacyShots();
    int failed = 0; for (const auto& k : R.checks) failed += k.pass ? 0 : 1;
    Finish(failed || R.checks.empty() ? "FAIL" : "PASS", failed ? std::to_string(failed) + " of " + std::to_string(R.checks.size()) + " shots differ from their baselines" : "");
}
