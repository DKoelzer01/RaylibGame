// The A/B runner (#871): one scenario run 2N times, A and B interleaved (A B A B ...: drift in the machine's load hits both
// alike), each in its own process with its side's environment. Compares every metric the runs report (the perf scenarios'
// frame times, the fluid scenarios' numbers) as mean +- standard deviation and the change from A to B.
//   game.exe --scenario-ab NAME --ab-a "K=V[;K=V]" --ab-b "K=V[;K=V]" [--ab-runs N] [--ab-max-delta PCT [--ab-metric NAME]]
//            [--scenario-out DIR] [other flags passed to every run]
// --ab-max-delta: exit 1 when |B - A| / A of the metric (default "frame mean") exceeds PCT: the profiler's overhead check is
//   --scenario-ab perf/planet-view --ab-a RAYLIBGAME_PROFILE=0 --ab-b RAYLIBGAME_PROFILE=1 --ab-runs 3 --ab-max-delta 1
// Timing comparisons need a quiet machine (CLAUDE.md: ask the coordinator for a quiet window).
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

int ScenarioProc_Run(const std::string& cmdline, double timeoutSec);
std::string ScenarioProc_Exe();

namespace {
std::string Q(const std::string& s) { return s.find_first_of(" \t\"") == std::string::npos ? s : "\"" + s + "\""; }
std::string Slurp(const std::string& p) {
    std::string s;
    if (FILE* f = fopen(p.c_str(), "rb")) { char b[8192]; size_t n; while ((n = fread(b, 1, sizeof b, f)) > 0) s.append(b, n); fclose(f); }
    return s;
}
// Every {"name": "...", "value": X} under "metrics" (results.json as scenario.cpp writes it).
std::map<std::string, double> Metrics(const std::string& j) {
    std::map<std::string, double> out;
    const size_t m = j.find("\"metrics\": [");
    if (m == std::string::npos) return out;
    const size_t end = j.find("\n  ]", m);
    size_t at = m;
    while ((at = j.find("{\"name\": \"", at)) != std::string::npos && (end == std::string::npos || at < end)) {
        at += 10;
        const size_t q = j.find('"', at);
        const std::string name = j.substr(at, q - at);
        const size_t v = j.find("\"value\": ", q);
        if (v == std::string::npos) break;
        const char* p = j.c_str() + v + 9;
        char* e = nullptr;
        const double val = strtod(p, &e);
        if (e != p) out[name] = val;
        at = q;
    }
    return out;
}
void SetEnv(const std::string& k, const std::string& v) {
#ifdef _WIN32
    _putenv_s(k.c_str(), v.c_str());
#else
    setenv(k.c_str(), v.c_str(), 1);
#endif
}
std::vector<std::pair<std::string, std::string>> ParseEnv(const std::string& s) {
    std::vector<std::pair<std::string, std::string>> out;
    size_t a = 0;
    while (a < s.size()) {
        size_t b = s.find(';', a); if (b == std::string::npos) b = s.size();
        const std::string kv = s.substr(a, b - a);
        const size_t eq = kv.find('=');
        if (eq != std::string::npos) out.push_back({ kv.substr(0, eq), kv.substr(eq + 1) });
        a = b + 1;
    }
    return out;
}
}

int ScenarioAB_Main(int argc, char** argv) {
    std::string name, envA, envB, outDir, metricName = "frame mean";
    int runs = 3;
    double maxDelta = -1.0;
    std::vector<std::string> pass;
    for (int i = 0; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--ab-a" && i + 1 < argc) envA = argv[++i];
        else if (a == "--ab-b" && i + 1 < argc) envB = argv[++i];
        else if (a == "--ab-runs" && i + 1 < argc) runs = std::max(1, atoi(argv[++i]));
        else if (a == "--ab-max-delta" && i + 1 < argc) maxDelta = atof(argv[++i]);
        else if (a == "--ab-metric" && i + 1 < argc) metricName = argv[++i];
        else if (a == "--scenario-out" && i + 1 < argc) outDir = argv[++i];
        else if (name.empty() && a[0] != '-') name = a;
        else pass.push_back(a);
    }
    if (name.empty()) { printf("usage: --scenario-ab NAME --ab-a K=V --ab-b K=V [--ab-runs N] [--ab-max-delta PCT] [--ab-metric NAME]\n"); return 2; }
    std::error_code ec;
    if (outDir.empty()) {
        const std::time_t t = std::time(nullptr);
        char b[32]; std::strftime(b, sizeof b, "%Y%m%d-%H%M%S", std::localtime(&t));
        outDir = (std::filesystem::temp_directory_path(ec) / "raylibgame_scenarios" / (std::string("ab-") + b)).generic_string();
    }
    outDir = std::filesystem::absolute(outDir, ec).generic_string();
    std::filesystem::create_directories(outDir, ec);
    setvbuf(stdout, nullptr, _IONBF, 0);
    const auto a = ParseEnv(envA), b = ParseEnv(envB);
    std::map<std::string, std::vector<double>> va, vb;
    const std::string exe = ScenarioProc_Exe();
    int bad = 0;
    for (int r = 0; r < runs * 2; ++r) {
        const bool isA = r % 2 == 0;
        // The other side's keys are cleared to empty first (a key only one side names would leak into the next run).
        for (const auto& kv : a) SetEnv(kv.first, "");
        for (const auto& kv : b) SetEnv(kv.first, "");
        for (const auto& kv : isA ? a : b) SetEnv(kv.first, kv.second);
        const std::string dir = outDir + "/" + (isA ? "A" : "B") + std::to_string(r / 2 + 1);
        std::string cmd = Q(exe) + " --scenario " + Q(name) + " --scenario-out " + Q(dir) + " --scenario-no-record";
        for (const std::string& p : pass) cmd += " " + Q(p);
        printf("ab: run %d of %d (%s: %s)\n", r + 1, runs * 2, isA ? "A" : "B", isA ? envA.c_str() : envB.c_str());
        const int rc = ScenarioProc_Run(cmd, 900.0);
        std::string safe = name; for (char& c : safe) if (c == ':') c = '_';
        const std::string j = Slurp(dir + "/" + safe + "/results.json");
        if (j.find("\"verdict\": \"PASS\"") == std::string::npos) { ++bad; printf("ab: run %d did not pass (exit %d): its numbers are still used if present\n", r + 1, rc); }
        for (const auto& kv : Metrics(j)) (isA ? va : vb)[kv.first].push_back(kv.second);
    }
    for (const auto& kv : a) SetEnv(kv.first, "");
    for (const auto& kv : b) SetEnv(kv.first, "");
    auto stats = [](const std::vector<double>& v, double& mean, double& sd) {
        mean = sd = 0.0;
        if (v.empty()) return;
        for (double x : v) mean += x;
        mean /= v.size();
        for (double x : v) sd += (x - mean) * (x - mean);
        sd = v.size() > 1 ? std::sqrt(sd / (v.size() - 1)) : 0.0;
    };
    std::string json = "{\n  \"scenario\": \"" + name + "\", \"a\": \"" + envA + "\", \"b\": \"" + envB + "\", \"runs\": " + std::to_string(runs) + ",\n  \"metrics\": [";
    printf("ab: %s, %d runs a side (A: %s | B: %s)\n", name.c_str(), runs, envA.c_str(), envB.c_str());
    printf("ab: %-34s %22s %22s %9s\n", "metric", "A mean +- sd", "B mean +- sd", "B vs A");
    double deltaChecked = 0.0; bool haveChecked = false;
    bool first = true;
    for (const auto& kv : va) {
        auto it = vb.find(kv.first);
        if (it == vb.end()) continue;
        double ma, sa, mb, sb;
        stats(kv.second, ma, sa); stats(it->second, mb, sb);
        const double pct = ma != 0.0 ? 100.0 * (mb - ma) / std::fabs(ma) : 0.0;
        printf("ab: %-34s %12.4g +- %-7.3g %12.4g +- %-7.3g %+8.2f%%\n", kv.first.c_str(), ma, sa, mb, sb, pct);
        char line[512];
        snprintf(line, sizeof line, "%s\n    {\"name\": \"%s\", \"a_mean\": %.6g, \"a_sd\": %.6g, \"b_mean\": %.6g, \"b_sd\": %.6g, \"change_pct\": %.4g}", first ? "" : ",", kv.first.c_str(), ma, sa, mb, sb, pct);
        json += line; first = false;
        if (kv.first == metricName) { deltaChecked = pct; haveChecked = true; }
    }
    json += "\n  ]\n}\n";
    if (FILE* f = fopen((outDir + "/ab.json").c_str(), "wb")) { fwrite(json.data(), 1, json.size(), f); fclose(f); }
    printf("ab: %s/ab.json%s\n", outDir.c_str(), bad ? (" (" + std::to_string(bad) + " runs did not pass)").c_str() : "");
    if (maxDelta >= 0.0) {
        if (!haveChecked) { printf("ab: no '%s' in both sides -> FAIL\n", metricName.c_str()); return 1; }
        const bool ok = std::fabs(deltaChecked) <= maxDelta;
        printf("ab: %s changed %+.2f%% (limit %.2f%%) -> %s\n", metricName.c_str(), deltaChecked, maxDelta, ok ? "PASS" : "FAIL");
        return ok ? 0 : 1;
    }
    return 0;
}
