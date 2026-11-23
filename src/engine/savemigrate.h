#ifndef SAVEMIGRATE_H
#define SAVEMIGRATE_H

// Save migrations (#1693): an older save is brought to the current format as text lines in memory (the file on disk is
// untouched; the upgraded save is written by the next normal save), then parsed by the one reader that knows only the
// current format. One SaveMigration per version bump, run in order.
#include <string>
#include <vector>

struct MigrateContext {
    int from = 0, to = 0;                    // the version the file was written with, and the one it is brought to
    std::vector<std::string> dropped;        // data a step made meaningless and removed ("what", logged and toasted)
    void drop(const std::string& what) { dropped.push_back(what); }
};

// One step: rewrites `lines` (the world file's text lines, without line ends, the first being the version line) from format
// `from` to `from + 1`. Returns false (and logs why via ctx.drop or its own LOG) only when the save cannot be brought on.
struct SaveMigration {
    int from;
    const char* what;
    bool (*step)(std::vector<std::string>& lines, MigrateContext& ctx);
};

// Whether this build opens a world file of that version; else `why` is the message for the log and the screen.
bool Save_VersionSupported(int version, std::string& why);
// Runs the steps from `fromVersion` to kSaveVersion over `lines`, updating the version line. False with `error` on failure.
bool Save_Migrate(int fromVersion, std::vector<std::string>& lines, MigrateContext& ctx, std::string& error);
// The player record's chain (its own number, kPlayerRecordVersion): same shape, no steps yet.
bool PlayerRecord_Migrate(int fromVersion, std::vector<std::string>& lines, MigrateContext& ctx, std::string& error);
// "Save upgraded from version N; dropped: a, b" (the dropped part only when something was).
std::string Save_UpgradeNote(const MigrateContext& ctx);

// Text file <-> lines (line ends stripped; written back with '\n').
bool Save_ReadLines(const std::string& path, std::vector<std::string>& lines);
bool Save_WriteLines(const std::string& path, const std::vector<std::string>& lines);

// The format version the readers are parsing right now (kSaveVersion once any migration ran; set by the world loader).
// The file's own version is MigrateContext::from.
int Save_LoadingVersion();
void Save_SetLoadingVersion(int version);

#endif
