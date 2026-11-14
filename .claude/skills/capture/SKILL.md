---
name: capture
description: Record the running game to check something visual. Builds the agent build if needed, runs game.exe with a scenario and the in-game recorder, and returns a GIF for the user plus a contact sheet and frames for you to inspect; can diff against an earlier capture. Use before claiming any visual fix (rendering, shaders, UI, camera, animation, terrain, water), and when the user reports something that looks wrong.
argument-hint: "<name> [game flags...]  e.g. walk-forest --seed 42 --day --bench 10 --test-walk"
---

# Capture: record the game and look at it

You can read images, not video. The game records itself (#375, `src/space/recorder.*`): numbered PNG frames, a GIF for the user and a contact sheet (up to 25 evenly spaced frames labelled `#<frame>  <t> s`) for you.

## Run it
```
powershell -ExecutionPolicy Bypass -File "${CLAUDE_SKILL_DIR}/scripts/capture.ps1" -Name <name> -GameArgs "<scenario flags>" [-Frames 120] [-Every 1] [-Start 0] [-NoHud] [-Compare <earlier capture>] [-Continue]
```
- `-Name`: a short label. The script appends the time; the game files it under `<root>/<YYYY-MM-DD>/<name>-<HHMMSS>/`.
- `-GameArgs`: the scenario that shows the problem. Always include `--seed N`. Options:
  - a scripted bench or test (`--bench 10 --test-walk`, `--test-ledge`, `--shot-forest`, `--day`, ...; see CLAUDE.md "Testing");
  - a saved world (`--load --at-save`, a read-only copy as in CLAUDE.md);
  - a fixed camera (`--cam dx dy dz`).
  The recorder waits for the loading screen and `WorldIdle()` before frame 0, and recording implies `--fixed-dt 1/60` with gameplay randomness seeded from `--seed`, so runs are reproducible.
- `-Frames` (120) and `-Every` (1): 120 frames, one per rendered frame (2 s at 60 FPS). Use `-Every 1` for flicker or jitter; `-Every 5` to cover a longer action.
- `-Start S`: wait S more seconds after idle before frame 0 (e.g. to let an animation start).
- `-NoHud`: frames are grabbed before the UI/HUD (the 3D scene + viewmodel only).
- `-Compare <name or path>`: diff each frame against the same frame of an earlier capture (before/after a fix). It prints the % of pixels that differ by more than 2/255 (average, worst frame). Compare captures with this, **not** with the printed hashes (tiny ±1 noise changes hashes between identical runs).
- The game exits after the last frame (`-Continue` keeps it running). If the `--bench` is shorter than the capture, the capture is finished at exit with what it has.

The script builds `build/check` first (it stops if the build fails), runs the game, waits, and prints the recorder's `record:` lines and the resolved paths (`dir`, `sheet`, `gif`, `log`).

## Where captures go
The game resolves the root: `RAYLIBGAME_CAPTURE_ROOT` if set; else the Claude media store on deweyserver, `\\DEWEYSERVER\files\claude-media\raylib-captures`, when reachable; else `%USERPROFILE%\RaylibCaptures`. The root is kept under 10 GB (`RAYLIBGAME_CAPTURE_MAX_GB`): the oldest captures are deleted first. Don't rely on old captures existing; re-capture a baseline when you need one. Captures never go in the repo. The user can also press **F10** in normal play to record (`<root>/<date>/<HHMMSS>/`, every 2 frames, up to 600, GIF + sheet).

## Look at it
1. Read `sheet.png` with the Read tool (UNC paths work: `\\DEWEYSERVER\files\...`): the whole run at a glance.
2. Read individual `frame_NNNN.png` files where something changes or looks wrong.
3. Before/after: capture the old build (or the fix disabled via its A/B switch) first, then the fix with `-Compare <first capture>`, and read both sheets.
4. Describe what you see concretely (what, where in the frame, which frames) and tie it to the code.

## Give it to the user
Report the `capture.gif` path (and the sheet). Don't paste images into chat.

## Rules
- Never while the user is playtesting (ask first). Never overwrite `./game.exe` (the script uses `build/check`).
- Captures take a while (loading + frames + writing; ~1 s per 30 frames extra over SMB): run in the background and check the output 15-45 s in (CLAUDE.md).
- If the build has no recorder yet (#375 not merged into the branch), the script falls back to the bench screenshots (`bench*.png`) and says so.
