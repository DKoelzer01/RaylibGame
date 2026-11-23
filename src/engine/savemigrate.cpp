#include "savemigrate.h"
#include "saveversion.h"
#include <cstdio>

namespace {
// The world file's steps, in order. A step turns format `from` into `from + 1`. Empty today (13 is the oldest and the newest).
//
// Example of what the 13 -> 14 step would look like when a change renames a line tag and drops a section:
//
//   bool Step13to14(std::vector<std::string>& lines, MigrateContext& ctx) {
//       std::vector<std::string> out;
//       for (std::string& l : lines) {
//           if (l.rfind("old_tag ", 0) == 0) l.replace(0, 7, "new_tag");        // renamed
//           if (l.rfind("retired_thing ", 0) == 0) { ctx.drop("retired things"); continue; } // meaningless now: dropped, told
//           out.push_back(l);
//       }
//       lines.swap(out);
//       return true;
//   }
//   const SaveMigration kWorldSteps[] = { { 13, "old_tag renamed, retired things dropped", Step13to14 } };
//
// and kSaveVersion becomes 14. Add a golden save (tests/saves/v14) made by the new build.
const SaveMigration* const kWorldSteps = nullptr;
constexpr int kWorldStepCount = 0;

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
