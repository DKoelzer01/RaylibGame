---
name: capture
description: Record the running game to check something visual. Builds the agent build if needed, runs game.exe with a scenario and in-game recording, and returns a GIF for the user plus a contact sheet and frames for you to inspect. Use before claiming any visual fix (rendering, shaders, UI, camera, animation, terrain, water), and when the user reports something that looks wrong.
argument-hint: "<name> [game flags...]  e.g. walk-forest --day --bench 10 --test-walk"
---

# Capture: record the game and look at it

You can read images, not video. This skill records frames **inside the game** (issue #375): numbered PNGs, a GIF for the user and a contact sheet (a grid of frames with numbers and times) for you.

## Run it
```
powershell -ExecutionPolicy Bypass -File "${CLAUDE_SKILL_DIR}/scripts/capture.ps1" -Name <name> -GameArgs "<scenario flags>" [-Frames 60] [-Every 5] [-NoHud]
```
- `<name>`: a short label; outputs go to `captures/<name>-<timestamp>/` in the repo (gitignored).
- `-GameArgs`: the scenario. Pick the one that shows the problem:
  - a scripted bench or test (`--seed 42 --bench 10 --test-walk`, `--test-ledge`, `--shot-forest`, `--day`, ...; see CLAUDE.md "Testing");
  - a saved world (`--load --at-save`, read-only copy as in CLAUDE.md);
  - a fixed camera (`--cam dx dy dz`).
  Always pass `--seed N` so runs are reproducible (the script adds `--fixed-dt` through `--record`).
- `-Frames` (default 60) and `-Every` (default 5): 60 frames, one every 5 rendered frames, gives ~5 s at 60 FPS. Use `-Every 1` for flicker or jitter (frame-to-frame changes).
- `-NoHud` to capture the world without HUD/UI.

The script builds `build/check` first (stops if the build fails), runs the game, waits for it, and prints the output paths. If the build dir doesn't exist yet, configure it as in CLAUDE.md (`-DGAME_OUTPUT_DIR=build/check/bin`).

## Look at it
1. Read `sheet.png` with the Read tool: it shows the whole run at a glance.
2. Read the individual `frame_NNNN.png` files where something changes or looks wrong (zoom in by reading a single frame).
3. For before/after: capture the same scenario and seed on the old build (or with the fix disabled) and compare sheets.
4. Describe what you see in concrete terms (what, where in the frame, which frames) and tie it to the code.

## Give it to the user
Report the path of `capture.gif` (and the sheet) so the user can open it. Don't paste images into the chat.

## Rules
- Never while the user is playtesting (ask first). Never overwrite `./game.exe`.
- Captures can take a minute: run in the background and check the output 15-45 s in (CLAUDE.md).
- Until #375's `--record` flags exist in the build, the script falls back to the bench screenshot (`bench.png`, `--shot-*` files) and says so.
