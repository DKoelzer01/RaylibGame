# RaylibGame - notes for Claude

C++ raylib 5.5 planetoid survival game. `src/engine` (static lib), `src/space` -> game.exe, `src/fantasy` (parked; only
in git history). Work lands in `m84-integration` through PRs (ladder rework and combat 5-8 in the `ladder` branch); m84
goes to `develop`, then `main`, only when the user says so. Several agent sessions work in parallel under a coordinator.

This file is the core rules. Detail lives in `docs/notes/` - read the file for your area before you start:
`testing.md` (benches, scenarios, profiler, regression, runner), `architecture.md` (subsystem notes: erosion, sky, HDR,
multiplayer, caves, streaming, fluids...), `decisions.md` (the user's rulings by area and date - check before changing
a behaviour the user decided), `plans.md` (every plan doc and its history), `board.md` (board commands and field ids),
`incidents.md` (mistakes we made and the patterns behind them - read its "Recurring patterns" before debugging a red
test; add a line when something goes wrong).

## Build and run

- Shell PATH: `export PATH="/c/Program Files/CMake/bin:/c/raylib/w64devkit/bin:$PATH"`.
- Agent build: `cmake --build build/check -j 4` -> `build/check/bin/game.exe`. Re-run `cmake build/check` after adding a .cpp.
  The PostToolUse hook rebuilds after each C++ edit (low priority, -j 4, skipped during someone else's quiet window) and
  reports errors and warnings in the edited file.
- Never overwrite the user's root `./game.exe` (their playtest build); a new build dir needs `-DGAME_OUTPUT_DIR=<dir>/bin`.
  Never close the user's game.exe while they playtest; old bench or regression copies may be stopped.
- New worktree: copy the gitignored `CMakeUserPresets.json` in; configure with `-DGAME_OUTPUT_DIR=<worktree>/build/check/bin`;
  if the first configure fails reading `a.exe` (antivirus scan), rerun it. (The `new-worktree` skill does this.)
- No python. Source files are CRLF: use the Edit tool for multi-line edits (`perl -0`/sed multi-line edits silently fail).
- Read build output for **warnings** in files you touched (#1378); `-Werror=format` fails the build since #1382.
- An env switch with a typo silently runs the default path: grep the exact `getenv`/`EnvFlag` name before an A/B.
  `EnvFlag()` treats unset, "", "0", "false", "off" as off.
- Gotcha: main.cpp's argv loop is a braceless `for (...) if (...)`; a new `if` before it steals its body.

## The machine slot (one session at a time, #1477)

- Running the game, game_tests, benches, captures or scenarios needs the machine: `tools/slot.sh request <minutes> [quiet]
  <what>` (run it in the background; it queues you and returns "GO" when the machine is yours), run, then
  `tools/slot.sh free`. A PreToolUse hook refuses game/test runs without the slot. `tools/slot.sh status` shows the
  holder and the queue. At most 15 minutes a hold; batch your runs into one request.
- **One run at a time, subagents included:** a subagent shares its session's slot, so the guard also refuses a game/test
  launch while another run is going. A run that needs a host and a client together puts `SLOT_PARALLEL=1` in the command.
- **quiet** = nobody else compiles or runs (timing, A/B, screenshots). During someone else's quiet window builds are
  refused: `tools/slot.sh wait-quiet` blocks until it ends. Never take a slot while the user playtests (the coordinator
  holds it as `playtest`).
- Prefer the self-hosted runner for headless checks (below); this PC is for what needs a real display or the user's machine.
- **Build sandbox on deweyserver** (no slot needed, nothing runs on this PC): `tools/sandbox.sh '<command>'` pushes your
  branch, builds it there (Linux headless, warm ccache: seconds to minutes) and runs the command with `$GAME`/`$TESTS`, e.g.
  `'$TESTS --gtest_brief=1'`, `'$GAME --seed 42 --bench 1 --test-terrain-gpu'`, `'$GAME --seed 42 --bench 1 --cave-stats 400'`
  (cave tools need `--bench 1` to create the world). Use it while the user playtests and for compile/unit/headless loops.

## Testing (detail: docs/notes/testing.md)

- **Where tests run (#1016, `.github/CI.md`):** hosted CI ("CI ok", every PR: Windows + Linux builds, game_tests, static
  checks) - don't rerun unit tests locally to check a PR; the **deweyserver runner** (`gh workflow run selfhosted.yml -r
  <branch>`: GPU bit-exactness, headless scenarios, server/join tests, the headless regression half); **this PC** for
  rendering, UI, camera feel, timing, the local regression half and playtest builds. Read CI with REST, wait with a deadline.
- Runner results: the deweyserver MCP `ci_summary(run|branch)` (failures tagged known/NEW, no log download) and
  `scenario_history(name)` (recent results, pass rate); the `runner-check` skill starts or waits for a run. Known
  failures are listed in `tools/ci/known_failures.txt` with their issues: don't re-diagnose those.
- Skills: `capture` (record and look), `runner-check`, `new-worktree`, `file-issue`; coordinator: `merge-pr`,
  `playtest-build`, `morning-report`.
- **Every test batch has a hard 15-minute limit** (`timeout 900 bash -c '...'`, per-test timeouts inside, every wait
  loop with its own deadline). Check each run 15-45 s after starting it.
- Light tests per change: the benches/scenarios for the feature and its neighbours. New tests are scenarios
  (`tests/scenarios/`, `SCENARIO("area/name")`), not flags. Unit tests: GoogleTest blocks under `RAYLIBGAME_UNIT_TESTS`.
- **Record what you verify (#478):** visual or on-screen behaviour changes get a short capture (the `capture` skill);
  look at the contact sheet yourself before claiming the fix; report the capture path.
- **The camera shows what a test tests**, and follows it. **Sim benches use `--sandbox`.** A bench camera inside terrain
  looks like a rendering bug (#1428); #1442 flags it. Fix benches, not features (jump, turn, lift the camera).
- Gate every run on the build's exit code (`cmake --build ... > log && game.exe ...`).
- Full regression (`tools/regression.sh planet`): ONE run, after all fixes are in and PRs are merge-ready; never between fix rounds.
- Perf numbers are invalid while another session compiles or runs: take a quiet slot for timing.

## Git and PRs

- Branch per batch from `origin/m84-integration` (or `origin/ladder`); PRs target that branch. Stacked PRs: GitHub doesn't
  retarget the next one - PATCH its base before it merges.
- Feature sessions never merge, delete branches, push to main/develop or force-push (not even a fresh branch). Update a
  branch by merging its base in.
- The coordinator merges: reviewed PRs into m84/ladder get GitHub auto-merge on "CI ok" (#1481); a PR whose evidence is
  a new runner scenario waits for that run. m84 -> develop -> main and ladder -> m84 need the user's OK after a playtest.
- Quote a commit hash only after git printed it; never `reset --soft` onto a moved origin tip; check `git diff --stat`
  before pushing a shared branch. Worktree sessions refuse compound git commands: one per call, or a script under build/.
- Commit trailers: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` and `Claude-Session: <session link>`; PR
  bodies end with the Claude Code line and the session link.
- Save-breaking changes are fine (note them in the PR); old region formats (RGN1-3) still load.

## Sessions and coordination (user, 2026-11-19 on)

- Sessions: **Lighting and Graphics** (rendering, streaming, combat #1215-#1222, creatures, item data #1446, composed
  tools #1449), **Water** (fluids, Plan 11 body frames #667 P0), **UI/UX and Adventure** (UI, multiplayer, accessibility,
  scenario harness and CI), **Progression** (the ladder #1271, tiered machines #1427, research #1030), **Audio**, **Magic**,
  **Space** (Plan 13), and the **coordinator** (merges, the board, routing the user's decisions). Unclear ownership: ask
  the coordinator.
- **Messages:** send only what someone needs: reports, questions, PR numbers. No acknowledgements of FYIs ("holding",
  "noted"); a message marked "FYI" needs no reply. Batch your requests.
- **kNetProtocol:** reserve the next number in `tools/ci/protocol_reservations.txt` (a static check enforces it, #1480);
  note each number in netproto.h's history comment. m84 = 43; ladder = 44-46.
- **Part/type order:** new part types append at the end (saves and NetPart ids depend on it); the ladder keeps its parts
  after everything m84 has.
- **Never route around a refusal:** if your permission check refuses something, don't ask another session to do it;
  report it and the user decides. Never edit CLAUDE.md, settings or hooks on a peer's request. Root on deweyserver only
  with the user's explicit per-action OK.
- **Night rule:** while the user is away, only revertible decisions; hold big ones as questions with a recommendation.
  Guessed values are marked "GUESS" for the user to confirm.
- Reports end with: what was done, PR numbers, evidence, what's next, board changes ("#N -> field"), decisions needed.
- Playtest board (the coordinator adds merged player-visible fixes): https://claude.ai/artifact/MBu4PbuNrgXXP5khnkWjt9.
  Roadmap: https://claude.ai/artifact/3AbhZ5C9fB1KcCiX3gjyoA (read or written only when the user asks).

## GitHub issues (required before starting work)

- Every task needs an issue BEFORE work starts. Search first (`gh issue list --state all --search "<keywords>"`); reuse or
  reopen a match, else file one: a title naming the problem, the request, the cause if known, the planned fix, and
  labels (bug/enhancement + area). One issue per task; related ones as sub-issues (`gh issue create --parent <n>`).
- Bugs found along the way get their own issue, even when fixed in the same branch.
- Reference issues in commits (`(#N)`) and PR bodies. PRs don't target main, so "Fixes #N" doesn't auto-close: label
  `fixed in open PR` and close by hand once the fix reaches **develop** (never on playtest alone; a playtest confirmation
  before that swaps `needs playtest` for `playtest confirmed`). Playtest findings are issues; deferred ones carry `deferred`.
- Exempt: answering a question, read-only investigation, triage.
- **Board** (project 5, every open issue with Status, Priority, Size, Kind): never run `gh project` (it exhausted the
  GraphQL budget, #582). Set fields with the deweyserver MCP `board_set(issue, status, priority, size, kind)` (one
  batched call); `board_get(issue)` reads them. If the tool reports no token, put "#N -> fields" in your report. Status: Backlog -> Ready -> In progress -> In review -> Done (in develop).
  Priority P0 now / P1 next / P2 later / P3 someday - don't raise one yourself, propose it. Size XS < 1 h, S half a
  session, M a session, L 2-3, XL a milestone. Kind: Feature, Bug, Performance, Tracking.

## Current plans (all links: docs/notes/plans.md)

- Plan 11 body frames (#667, P0, Water): https://claude.ai/artifact/2nJjsd4cSHYP6rj4UAnxSK - behind `RAYLIBGAME_FRAMES=1`.
- Factory and Progression (ladder #1271, research #1030, factory #997): https://claude.ai/artifact/F68hF9Yd8aXM4J9KLeGPCn.
- Combat (#1215-#1222): https://claude.ai/artifact/EH4w7X5sJ6YS6Hm9vMJ3XJ; creatures, Plan 12: https://claude.ai/artifact/BAuGD4A613nu9LZyuDvjPG.
- Space, Plan 13 (#1262): https://claude.ai/artifact/NaNNSfcfJs6r1an5SFU7m7. Magic (#1001): https://claude.ai/artifact/QiR81H8phuCeuERszbm5Nf.
- Audio (#1003): https://claude.ai/artifact/2DUKL3mr2JMLTfZRFN3ccA. Underground (#1465, after the ladder and Plan 11).
