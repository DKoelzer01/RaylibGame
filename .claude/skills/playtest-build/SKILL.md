---
name: playtest-build
description: Coordinator only, and only when the user says to. Install a playtest build into the user's root game.exe - hold the machine, back up the current exe, build the requested branch tip in release (with multiplayer transport), smoke-test it, install it, and hold the slot as "playtest" so no session runs the game while they play. Use when the user says they want to playtest.
argument-hint: "[branch, default m84-integration]"
---

# Playtest build (coordinator, on the user's word only)

The root `C:\Users\dewey\RaylibGame\game.exe` is the user's. Never touch it without their explicit "go" for this build.

1. **Hold the machine** for the user: `SLOT_NAME=playtest bash tools/slot.sh take 600 quiet "user playtest"`
   (waits for nobody: if a session holds it, ask it to finish first). Quiet = no session builds or runs meanwhile.
2. **Build** the branch tip in a clean worktree build dir configured with the GNS release preset (see
   docs/notes/architecture.md, Multiplayer): configure `local-release`, then set ENABLE_GNS=ON, GNS_ROOT, and
   `GAME_OUTPUT_DIR=<that dir>/bin` in the cache; check CMakeCache before building. Read warnings.
3. **Smoke test** the new exe from its own folder: one short scenario that loads the planet (e.g. `--scenario world/spawn`
   or `--test-caves`), 0 crashes, no frame over 100 ms after loading.
4. **Install:** close nothing the user has open; copy the root `game.exe` to `game.exe.<old short sha>-backup`, copy the
   new exe in, and move the root checkout's resources to the same commit (`git -C C:/Users/dewey/RaylibGame checkout
   <sha>` only if it's clean; else tell the user). Write a checklist of what to test from the playtest board.
5. **Tell the user** the commit, the backup name and the checklist path. Keep the `playtest` hold until they say
   they're done; then `SLOT_NAME=playtest bash tools/slot.sh free`.
