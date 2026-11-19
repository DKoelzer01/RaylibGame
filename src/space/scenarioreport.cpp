// Report v1 (#870): results.json (one scenario's, or a tag run's "runs") -> one HTML page, and a line per run appended to the
// history (JSONL), whose last 20 runs per scenario give the pass rates (flaky = under 100%) and the perf baselines.
//   game.exe --scenario-report <results.json | dir> [out.html] [--history <file.jsonl> | --no-history]
// A dir is searched for results.json files (a tag run's combined one wins). No window is opened.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

// ---- a small JSON reader (objects, arrays, strings, numbers, true/false/null) ----
struct J {
    enum Kind { Null, Bool, Num, Str, Arr, Obj } k = Null;
    bool b = false; double n = 0.0; std::string s;
    std::vector<J> a;
    std::vector<std::pair<std::string, J>> o;
    const J& operator[](const char* key) const {
        static const J none;
        for (const auto& kv : o) if (kv.first == key) return kv.second;
        return none;
    }
    std::string str(const char* key, const char* def = "") const { const J& v = (*this)[key]; return v.k == Str ? v.s : def; }
    double num(const char* key, double def = 0.0) const { const J& v = (*this)[key]; return v.k == Num ? v.n : def; }
    bool flag(const char* key) const { const J& v = (*this)[key]; return v.k == Bool && v.b; }
};
struct Parser {
    const char* p; const char* e; bool ok = true;
    void ws() { while (p < e && (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t')) ++p; }
    J value() {
        ws();
        J v;
        if (p >= e) { ok = false; return v; }
        if (*p == '{') {
            v.k = J::Obj; ++p; ws();
            if (p < e && *p == '}') { ++p; return v; }
            while (ok && p < e) {
                ws(); J key = value(); ws();
                if (key.k != J::Str || p >= e || *p != ':') { ok = false; break; }
                ++p;
                J val = value();
                v.o.push_back({ key.s, std::move(val) });
                ws();
                if (p < e && *p == ',') { ++p; continue; }
                if (p < e && *p == '}') { ++p; break; }
                ok = false;
            }
        } else if (*p == '[') {
            v.k = J::Arr; ++p; ws();
            if (p < e && *p == ']') { ++p; return v; }
            while (ok && p < e) {
                v.a.push_back(value()); ws();
                if (p < e && *p == ',') { ++p; continue; }
                if (p < e && *p == ']') { ++p; break; }
                ok = false;
            }
        } else if (*p == '"') {
            v.k = J::Str; ++p;
            while (p < e && *p != '"') {
                if (*p == '\\' && p + 1 < e) {
                    ++p;
                    switch (*p) {
                        case 'n': v.s += '\n'; break; case 't': v.s += '\t'; break; case 'r': v.s += '\r'; break;
                        case 'u': if (p + 4 < e) { v.s += (char)strtol(std::string(p + 1, p + 5).c_str(), nullptr, 16); p += 4; } break;
                        default: v.s += *p;
                    }
                } else v.s += *p;
                ++p;
            }
            if (p < e) ++p; else ok = false;
        } else if (!strncmp(p, "true", 4)) { v.k = J::Bool; v.b = true; p += 4; }
        else if (!strncmp(p, "false", 5)) { v.k = J::Bool; p += 5; }
        else if (!strncmp(p, "null", 4)) { p += 4; }
        else { char* end = nullptr; v.n = strtod(p, &end); if (end == p) { ok = false; return v; } v.k = J::Num; p = end; }
        return v;
    }
};
bool ParseJson(const std::string& text, J& out) { Parser ps{ text.data(), text.data() + text.size() }; out = ps.value(); return ps.ok; }
std::string Slurp(const std::string& path) {
    std::string s;
    if (FILE* f = fopen(path.c_str(), "rb")) { char b[8192]; size_t n; while ((n = fread(b, 1, sizeof b, f)) > 0) s.append(b, n); fclose(f); }
    return s;
}
std::string Html(const std::string& s) {
    std::string o;
    for (char c : s) { if (c == '<') o += "&lt;"; else if (c == '>') o += "&gt;"; else if (c == '&') o += "&amp;"; else if (c == '"') o += "&quot;"; else o += c; }
    return o;
}
std::string Fmt(double v) { char b[64]; snprintf(b, sizeof b, std::fabs(v) >= 1000 || v == (long long)v ? "%.0f" : "%.3g", v); return b; }
std::string JsonEsc(const std::string& s) {
    std::string o;
    for (unsigned char c : s) { if (c == '"' || c == '\\') { o += '\\'; o += (char)c; } else if (c < 0x20) o += ' '; else o += (char)c; }
    return o;
}
std::string FileUrl(const std::string& path) {
    std::string p = path;
    for (char& c : p) if (c == '\\') c = '/';
    if (p.rfind("//", 0) == 0) return "file:" + p;        // a share: file://server/...
    return "file:///" + (p.size() > 1 && p[1] == ':' ? p : fs::absolute(p).generic_string());
}

struct HistRun { std::string when, verdict, build; std::map<std::string, double> metrics; };

} // namespace

int ScenarioReport_Main(int argc, char** argv) {
    std::string in, out, history;
    bool useHistory = true;
    for (int i = 0; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--history" && i + 1 < argc) history = argv[++i];
        else if (a == "--no-history") useHistory = false;
        else if (in.empty()) in = a;
        else if (out.empty()) out = a;
    }
    if (in.empty()) { printf("usage: --scenario-report <results.json | dir> [out.html] [--history file.jsonl | --no-history]\n"); return 2; }
    // The runs: a tag run's combined file, else every results.json under the dir.
    std::vector<J> runs;
    std::vector<std::string> files;
    std::error_code ec;
    if (fs::is_directory(in, ec)) {
        if (fs::exists(fs::path(in) / "results.json", ec)) files.push_back((fs::path(in) / "results.json").string());
        else for (const auto& e : fs::recursive_directory_iterator(in, ec)) if (e.path().filename() == "results.json") files.push_back(e.path().string());
    } else files.push_back(in);
    // A tag run's combined file holds its scenarios' own results: those under its folder are not read again.
    std::vector<std::string> combinedDirs;
    std::vector<std::pair<std::string, J>> parsed;
    for (const std::string& f : files) {
        J j;
        if (!ParseJson(Slurp(f), j)) { printf("report: %s is not valid JSON\n", f.c_str()); continue; }
        if (j["runs"].k == J::Arr) combinedDirs.push_back(fs::path(f).parent_path().generic_string() + "/");
        parsed.push_back({ fs::path(f).generic_string(), std::move(j) });
    }
    for (auto& pj : parsed) {
        const bool combined = pj.second["runs"].k == J::Arr;
        bool inside = false;
        for (const std::string& d : combinedDirs) inside = inside || (!combined && pj.first.rfind(d, 0) == 0);
        if (inside) continue;
        if (combined) for (const J& r : pj.second["runs"].a) runs.push_back(r); else runs.push_back(pj.second);
    }
    if (runs.empty()) { printf("report: no runs found in %s\n", in.c_str()); return 2; }
    if (out.empty()) out = (fs::is_directory(in, ec) ? fs::path(in) : fs::path(in).parent_path()).append("report.html").string();

    // History: read (for pass rates and baselines), then append this report's runs.
    if (history.empty()) {
        const char* e = getenv("RAYLIBGAME_TEST_HISTORY");
        std::string dir = e && *e ? e : "//DEWEYSERVER/files/claude-media/raylib-test-history";
        if (!fs::is_directory(dir, ec)) { const char* up = getenv("USERPROFILE"); dir = std::string(up ? up : ".") + "/RaylibTestHistory"; fs::create_directories(dir, ec); }
        history = dir + "/history.jsonl";
    }
    std::map<std::string, std::vector<HistRun>> past;
    if (useHistory) {
        if (FILE* f = fopen(history.c_str(), "rb")) {
            char line[1 << 16];
            while (fgets(line, sizeof line, f)) {
                J h;
                if (!ParseJson(line, h) || h.k != J::Obj) continue;
                HistRun r{ h.str("when"), h.str("verdict"), h.str("build"), {} };
                for (const auto& kv : h["metrics"].o) if (kv.second.k == J::Num) r.metrics[kv.first] = kv.second.n;
                past[h.str("scenario")].push_back(r);
            }
            fclose(f);
        }
    }
    const std::time_t now = std::time(nullptr);
    char when[32]; std::strftime(when, sizeof when, "%Y-%m-%dT%H:%M:%S", std::localtime(&now));

    // Pass rates over the last 20 runs, this one included.
    auto rate = [&](const std::string& name, const std::string& verdict, int& n) {
        std::vector<HistRun> v = past[name];
        v.push_back({ when, verdict, "", {} });
        if (v.size() > 20) v.erase(v.begin(), v.end() - 20);
        int pass = 0; for (const HistRun& r : v) pass += r.verdict == "PASS";
        n = (int)v.size();
        return 100.0 * pass / std::max(1, n);
    };
    // The baseline of a metric: the median of the last 10 passing runs' values.
    auto baseline = [&](const std::string& name, const std::string& metric, double& base) {
        std::vector<double> vals;
        const auto& v = past[name];
        for (auto it = v.rbegin(); it != v.rend() && vals.size() < 10; ++it)
            if (it->verdict == "PASS") { auto m = it->metrics.find(metric); if (m != it->metrics.end()) vals.push_back(m->second); }
        if (vals.empty()) return false;
        std::sort(vals.begin(), vals.end());
        base = vals[vals.size() / 2];
        return true;
    };

    std::vector<const J*> failed, passed;
    for (const J& r : runs) (r.str("verdict") == "PASS" ? passed : failed).push_back(&r);

    std::string h;
    h += "<!doctype html><html><head><meta charset=utf-8><meta name=viewport content=\"width=device-width,initial-scale=1\"><title>Test report</title><style>"
         ":root{--bg:#f6f7f5;--panel:#fff;--fg:#1d2321;--muted:#5d6763;--line:#d8ddd9;--ok:#2f6f5e;--bad:#a3342c;--warn:#9a5b12}"
         "@media (prefers-color-scheme:dark){:root{--bg:#141917;--panel:#1b2220;--fg:#e4e9e6;--muted:#9aa5a0;--line:#2d3633;--ok:#7cc4ad;--bad:#ec8f86;--warn:#e0b070}}"
         "body{background:var(--bg);color:var(--fg);font:14px/1.5 system-ui,sans-serif;margin:0}main{max-width:1000px;margin:0 auto;padding:24px 16px}"
         "h1{font-size:24px;margin:0 0 4px}h2{font-size:18px;margin:28px 0 8px;border-top:1px solid var(--line);padding-top:12px}"
         ".card{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:12px;margin:10px 0}"
         "table{border-collapse:collapse;width:100%;font-size:13px}th,td{text-align:left;padding:4px 8px;border-bottom:1px solid var(--line);vertical-align:top}"
         "th{color:var(--muted);font-weight:600}.ok{color:var(--ok)}.bad{color:var(--bad);font-weight:600}.warn{color:var(--warn)}.muted{color:var(--muted)}"
         "img{max-width:100%;border-radius:6px;margin-top:8px}code{font-size:12px}.tbl{overflow-x:auto}</style></head><body><main>";
    h += "<h1>Test report</h1><p class=muted>" + std::string(when) + " &middot; " + std::to_string(runs.size()) + " runs: <span class=ok>" + std::to_string(passed.size()) +
         " passed</span>, <span class=" + (failed.empty() ? "ok" : "bad") + ">" + std::to_string(failed.size()) + " failed</span> &middot; build <code>" +
         Html(runs[0].str("build")) + "</code></p>";

    auto checksTable = [&](const J& r, bool onlyFailed) {
        std::string t = "<div class=tbl><table><tr><th>Check</th><th>Value</th><th>Wanted</th><th></th></tr>";
        for (const J& c : r["checks"].a) {
            if (onlyFailed && c.flag("pass")) continue;
            t += "<tr><td>" + Html(c.str("name")) + "</td><td>" + Fmt(c.num("value")) + "</td><td>" + Html(c.str("op")) + " " + Html(c.str("expected")) + "</td><td class=" +
                 (c.flag("pass") ? "ok>ok" : "bad>FAIL") + (c.str("detail").empty() ? "" : " <span class=muted>" + Html(c.str("detail")) + "</span>") + "</td></tr>";
        }
        return t + "</table></div>";
    };

    h += "<h2>Failures</h2>";
    if (failed.empty()) h += "<p class=ok>None.</p>";
    for (const J* r : failed) {
        h += "<div class=card><b class=bad>" + Html(r->str("verdict")) + "</b> <b>" + Html(r->str("scenario")) + "</b> <span class=muted>" + Html(r->str("file")) + "</span>";
        if (!r->str("reason").empty()) h += "<p>" + Html(r->str("reason")) + "</p>";
        h += checksTable(*r, false);
        for (const J& n : (*r)["notes"].a) h += "<div class=muted>" + Html(n.s) + "</div>";
        const J& cap = (*r)["capture"];
        if (!cap.str("sheet").empty()) h += "<a href=\"" + FileUrl(cap.str("dir")) + "\">capture</a><br><img src=\"" + FileUrl(cap.str("sheet")) + "\" alt=\"contact sheet\">";
        else if (!cap.str("dir").empty()) h += "<a href=\"" + FileUrl(cap.str("dir")) + "\">capture</a>";
        h += "<div class=muted>args <code>" + Html(r->str("args")) + "</code></div></div>";
    }

    // Perf: every metric against its baseline (the median of the last 10 passing runs). Reported always; the policy's
    // failures (+15% frame time, a new hitch over 30 ms) are checks in the perf scenarios themselves.
    h += "<h2>Perf against the last passing runs</h2><div class=tbl><table><tr><th>Scenario</th><th>Metric</th><th>Now</th><th>Baseline</th><th>Change</th></tr>";
    int perfRows = 0;
    for (const J& r : runs) {
        std::vector<std::pair<std::string, std::pair<double, std::string>>> ms;
        for (const J& m : r["metrics"].a) ms.push_back({ m.str("name"), { m.num("value"), m.str("unit") } });
        const J& prof = r["profile"];
        if (prof.k == J::Obj) {
            if (prof["frame"].k == J::Obj) { ms.push_back({ "frame p50", { prof["frame"].num("p50_ms"), "ms" } }); ms.push_back({ "frame p99", { prof["frame"].num("p99_ms"), "ms" } }); }
        }
        for (const auto& m : ms) {
            double base = 0.0;
            const bool has = baseline(r.str("scenario"), m.first, base);
            const double pct = has && base != 0.0 ? 100.0 * (m.second.first - base) / std::fabs(base) : 0.0;
            h += "<tr><td>" + Html(r.str("scenario")) + "</td><td>" + Html(m.first) + "</td><td>" + Fmt(m.second.first) + " " + Html(m.second.second) + "</td><td>" +
                 (has ? Fmt(base) : "<span class=muted>none yet</span>") + "</td><td class=" + (!has ? "muted>-" : std::fabs(pct) >= 15.0 ? "warn>" + Fmt(pct) + "%" : "muted>" + Fmt(pct) + "%") + "</td></tr>";
            ++perfRows;
        }
    }
    h += "</table></div>";
    if (!perfRows) h += "<p class=muted>No metrics in these runs.</p>";

    h += "<h2>Flaky (pass rate over the last 20 runs)</h2>";
    std::string flaky;
    for (const J& r : runs) {
        int n = 0;
        const double pr = rate(r.str("scenario"), r.str("verdict"), n);
        if (pr < 100.0 && n >= 2) flaky += "<tr><td>" + Html(r.str("scenario")) + "</td><td class=warn>" + Fmt(pr) + "%</td><td>" + std::to_string(n) + " runs</td></tr>";
    }
    h += flaky.empty() ? "<p class=ok>None (a scenario needs two runs in the history to count).</p>" : "<div class=tbl><table><tr><th>Scenario</th><th>Pass rate</th><th></th></tr>" + flaky + "</table></div>";

    h += "<h2>Passed</h2>";
    for (const J* r : passed) {
        h += "<details class=card><summary><b class=ok>PASS</b> " + Html(r->str("scenario")) + " <span class=muted>" + Fmt(r->num("sim_seconds")) + " s sim, " +
             Fmt(r->num("wall_seconds")) + " s wall" + (r->str("legacy_flag").empty() ? "" : ", was " + Html(r->str("legacy_flag"))) + "</span></summary>" + checksTable(*r, false) + "</details>";
    }
    h += "<p class=muted>History: " + Html(useHistory ? history : "off") + "</p></main></body></html>\n";

    if (FILE* f = fopen(out.c_str(), "wb")) { fwrite(h.data(), 1, h.size(), f); fclose(f); }
    else { printf("report: cannot write %s\n", out.c_str()); return 1; }

    if (useHistory) {
        fs::create_directories(fs::path(history).parent_path(), ec);
        if (FILE* f = fopen(history.c_str(), "ab")) {
            for (const J& r : runs) {
                std::string line = "{\"when\": \"" + std::string(when) + "\", \"scenario\": \"" + JsonEsc(r.str("scenario")) + "\", \"verdict\": \"" + JsonEsc(r.str("verdict")) +
                                   "\", \"build\": \"" + JsonEsc(r.str("build")) + "\", \"wall_seconds\": " + Fmt(r.num("wall_seconds")) + ", \"metrics\": {";
                bool first = true;
                for (const J& m : r["metrics"].a) { line += std::string(first ? "" : ", ") + "\"" + JsonEsc(m.str("name")) + "\": " + Fmt(m.num("value")); first = false; }
                const J& prof = r["profile"];
                if (prof.k == J::Obj && prof["frame"].k == J::Obj) {
                    line += std::string(first ? "" : ", ") + "\"frame p50\": " + Fmt(prof["frame"].num("p50_ms")) + ", \"frame p99\": " + Fmt(prof["frame"].num("p99_ms"));
                }
                line += "}}\n";
                fputs(line.c_str(), f);
            }
            fclose(f);
        } else printf("report: cannot append to the history %s\n", history.c_str());
    }
    printf("report: %zu runs (%zu failed) -> %s%s\n", runs.size(), failed.size(), out.c_str(), useHistory ? (" (history " + history + ")").c_str() : "");
    return failed.empty() ? 0 : 1;
}
