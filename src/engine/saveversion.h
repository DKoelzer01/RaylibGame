#ifndef SAVEVERSION_H
#define SAVEVERSION_H

// The save's one version (#1693). The world file's first line is "RAYLIBGAME_SAVE <kSaveVersion>"; every section the save
// holds (extras, plots, machines, fluid, rings) is in that format. A save-breaking change bumps this and adds one step to
// the migration table (savemigrate.cpp); the readers only ever know the current format. See docs/notes/architecture.md,
// "Saves and migrations (#1693)".
constexpr int kSaveVersion = 13;
// The oldest world file this build still opens (older ones are refused with a message, not half-loaded).
constexpr int kOldestSaveVersion = 13;
// The player record's own number ("RAYLIBGAME_PLAYER <n>", <save>.players/<id>.player); it has its own migration chain.
constexpr int kPlayerRecordVersion = 1;

#endif
