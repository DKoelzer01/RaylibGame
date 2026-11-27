#include "savemigrate.h"
#include "saveversion.h"
#include "plots.h"
#include <cstdio>
#include <iterator>

namespace {
// The world file's steps, in order. A step turns format `from` into `from + 1`; each new version adds a golden save
// (tests/saves/v<N>, made by saves/golden-write). A step may rename tags or drop what no longer means anything
// (ctx.drop("what"): logged and toasted).

// 13 -> 14 (#1881): part types were numbered by the enum's order; 14 names them in a palette line after the camera line.
// A 13 file gets the order it was written in (frozen in plots.cpp), so a type added mid-enum since still maps right.
bool Step13to14(std::vector<std::string>& lines, MigrateContext&) {
    for (size_t i = 0; i < lines.size(); ++i)
        if (lines[i].rfind("camera ", 0) == 0) { lines.insert(lines.begin() + (long)i + 1, PartPalette_V13Line()); return true; }
    return false; // (no header camera line: not a world file this reader knows)
}
const SaveMigration kWorldStepTable[] = { { 13, "part types by name (a part_palette line)", Step13to14 } };
const SaveMigration* const kWorldSteps = kWorldStepTable;
constexpr int kWorldStepCount = (int)std::size(kWorldStepTable);

const SaveMigration* const kPlayerSteps = nullptr;
constexpr int kPlayerStepCount = 0;

bool Run(const SaveMigration* steps, int count, int from, int to, const char* header, std::vector<std::string>& lines, MigrateContext& ctx, std::string& error) {
    ctx.from = from;
    ctx.to = to;
    for (int v = from; v < to; ++v) {
        const SaveMigration* m = nullptr;
        for (int i = 0; i < count; ++i) if (steps[i].from == v) m = &steps[i];
        if (!m) { error = "no migration from version " + std::to_string(v) + " to " + std::to_string(v + 1); return false; }
        if (!m->step(lines, ctx)) { error = std::string("migration from version ") + std::to_string(v) + " failed (" + m->what + ")"; return false; }
        if (!lines.empty()) lines[0] = std::string(header) + " " + std::to_string(v + 1);
    }
    return true;
}

int gLoadingVersion = kSaveVersion;
}

bool Save_VersionSupported(int version, std::string& why) {
    if (version > kSaveVersion) {
        why = "This save was made by a newer build (version " + std::to_string(version) + "; this build reads up to " + std::to_string(kSaveVersion) + ").";
        return false;
    }
    if (version < kOldestSaveVersion) {
        why = "This save is from an older build (version " + std::to_string(version) + "); saves from version " + std::to_string(kOldestSaveVersion) + " on load.";
        return false;
    }
    return true;
}

bool Save_Migrate(int fromVersion, std::vector<std::string>& lines, MigrateContext& ctx, std::string& error) {
    return Run(kWorldSteps, kWorldStepCount, fromVersion, kSaveVersion, "RAYLIBGAME_SAVE", lines, ctx, error);
}

bool PlayerRecord_Migrate(int fromVersion, std::vector<std::string>& lines, MigrateContext& ctx, std::string& error) {
    return Run(kPlayerSteps, kPlayerStepCount, fromVersion, kPlayerRecordVersion, "RAYLIBGAME_PLAYER", lines, ctx, error);
}

std::string Save_UpgradeNote(const MigrateContext& ctx) {
    std::string s = "Save upgraded from version " + std::to_string(ctx.from);
    for (size_t i = 0; i < ctx.dropped.size(); ++i) s += (i ? ", " : "; dropped: ") + ctx.dropped[i];
    return s;
}

bool Save_ReadLines(const std::string& path, std::vector<std::string>& lines) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    std::string cur;
    char buf[4096];
    size_t n;
    while ((n = fread(buf, 1, sizeof buf, f)) > 0)
        for (size_t i = 0; i < n; ++i) {
            if (buf[i] == '\n') { if (!cur.empty() && cur.back() == '\r') cur.pop_back(); lines.push_back(cur); cur.clear(); }
            else cur += buf[i];
        }
    if (!cur.empty()) lines.push_back(cur);
    fclose(f);
    return true;
}

bool Save_WriteLines(const std::string& path, const std::vector<std::string>& lines) {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) return false;
    for (const std::string& l : lines) { fwrite(l.data(), 1, l.size(), f); fputc('\n', f); }
    return fclose(f) == 0;
}

int Save_LoadingVersion() { return gLoadingVersion; }
void Save_SetLoadingVersion(int version) { gLoadingVersion = version; }
