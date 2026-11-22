---
name: new-worktree
description: Set up a fresh RaylibGame worktree's agent build - copy CMakeUserPresets.json, configure build/check with its own GAME_OUTPUT_DIR (never the user's root game.exe), retry the antivirus a.exe failure, and build game + game_tests at low priority. Use right after creating or entering a new worktree, before the first build.
argument-hint: "[--no-build]"
---

# New worktree setup

```bash
bash .claude/skills/new-worktree/scripts/setup.sh            # configure + build (low priority, -j 4, waits out quiet windows)
bash .claude/skills/new-worktree/scripts/setup.sh --no-build # configure only
```

It refuses to run in the main checkout (whose game.exe is the user's playtest build) and checks that
GAME_OUTPUT_DIR points into this worktree before building. If it fails, read `build/configure.log` or `build/build.log`.
For a GNS (multiplayer transport) build, see docs/notes/architecture.md (Multiplayer) instead.
