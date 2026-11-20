# RaylibGame - notes for Claude

C++ raylib 5.5 planetoid survival game. `src/engine` (static lib), `src/space` -> game.exe,
`src/fantasy` (parked). Work lands on `develop` through PRs; `main` lags.

## Build and run

- Shell PATH: `export PATH="/c/Program Files/CMake/bin:/c/raylib/w64devkit/bin:$PATH"`.
- Agent build: `cmake --build build/check -j 8` -> `build/check/bin/game.exe`. Re-run
  `cmake build/check` after adding a .cpp (GLOB with CONFIGURE_DEPENDS usually catches it).
- Never overwrite the user's root `./game.exe` (their playtest build). Any new build dir must be
  configured with `-DGAME_OUTPUT_DIR=<dir>/bin`.
- While the user is playtesting, never close their game.exe. Other game processes (old bench or
  regression copies) may be stopped.
- New worktree: copy the gitignored `CMakeUserPresets.json` in; configure with
  `-DGAME_OUTPUT_DIR=<worktree>/build/check/bin`; if the first configure fails reading `a.exe` (antivirus scan), rerun it.
- No python on this machine. Source files are CRLF: use the Edit tool for multi-line edits; `perl -0`/sed
  multi-line substitutions and perl replacements containing `\n` silently fail or insert literal newlines.
- Gotcha: main.cpp's argv loop is a braceless `for (...) if (...)`; inserting a new `if` before it steals its body.

## Testing

- **Where tests run (user, 2026-11-19, #1016; details in `.github/CI.md`).** Run each check where it belongs; this PC is
  only for what can't run elsewhere or has to be validated here:
  - **GitHub-hosted CI, every PR (required check "CI ok"):** the Windows build (w64devkit), the Linux server build,
    game_tests, static checks (NetMsg ids / kNetProtocol, bake-version reminder). Don't rebuild and rerun unit tests
    locally just to check a PR: push and read CI. Never ask for a merge while "CI ok" is red or pending.
  - **Self-hosted on deweyserver (GPU, isolated from the LAN):** `--test-terrain-gpu` bit-exactness, headless scenarios
    (fluids, craft, skills, world), the multiplayer server and join tests, the Linux server package, and the headless
    half of the regression (`tools/regression.sh planet --only headless`, #1028; nightly, or add `-f regression=headless`).
    Runs on pushes to m84-integration/develop and nightly; for your own branch, dispatch it: `gh workflow run
    selfhosted.yml -r <branch>` (one GPU run at a time, so batch your changes). Use it instead of running those tests here.
  - **This PC only:** what needs the Windows client with a real display/GPU or the user's machine: rendering, lighting,
    UI and any visual change (captures/screenshots), camera and movement feel, performance and timing (quiet windows),
    the local half of the regression (`tools/regression.sh planet --only local`), playtest builds, and final validation
    of anything the user will see.
  - Read CI with REST (`gh api repos/DKoelzer01/RaylibGame/commits/<sha>/check-runs`, `gh run view <id> --log-failed`);
    wait with a deadline, never an open-ended poll loop.
- Light tests per change: the scripted benches for the feature touched and its neighbours
  (`game.exe --seed 42 --bench 3 --test-<name>`), plus a screenshot (bench.png) when visuals change.
- **Unit tests (Tier 1, #869):** GoogleTest, at the end of the .cpp they test inside `#ifdef RAYLIBGAME_UNIT_TESTS` (the game
  never defines it). `cmake --build build/check --target game_tests` then `ctest --test-dir build/check -j 4` (or
  `build/check/bin/game_tests.exe` from the repo root); the suite stays under 30 s. Re-run `cmake build/check` after adding
  the first test block to a file. No window or GL in them (`UnitTest::FieldOnlyScene()` for code that reads only renderOrigin).
- **Scenarios (Tier 2, #870):** one file per scenario in `tests/scenarios/` (`SCENARIO("area/name") { s.step(...); }`, see
  src/space/scenario.h). `game.exe --scenario NAME` (or `tag:fluid`, `all`, `a,b`; `--scenario-list`) runs with fixed dt, lockstep
  fluids, seed 42, a scratch save dir, a heartbeat every 5 s and a low-res recording (deleted on PASS, kept on FAIL); results go
  to `<--scenario-out DIR>/<name>/results.json`. `game.exe --scenario-report DIR [out.html]` makes the HTML report and appends to
  the history on the share. Several names or a tag run each scenario in its own process. New tests are scenarios, not flags.
- **Profiler (#871):** without Tracy, `PROFILE_SCOPE` records into a per-thread ring (~20 ns a zone on, ~1 ns off); scenarios
  turn it on (`RAYLIBGAME_PROFILE=1` for any run) and results.json gets its summary (zones p50/p95/p99, frame histogram, hitches
  over 20 ms with their top zones, gputimer's cpu:/gpu: passes, counters, peak commit/VRAM). A/B timing:
  `game.exe --scenario-ab perf/planet-view --ab-a K=V --ab-b K=V --ab-runs 3` (interleaved; quiet window only). Perf policy: the
  report always shows the numbers and fails only on +15% mean frame time or a new hitch over 30 ms against recent passes.
  In play (not a scenario), `RAYLIBGAME_PROFILE=1` writes the session's summary on quit to saves/logs/profile-<date-time>.txt and
  .json (#886): zones, GPU passes, the frame histogram, and every hitch with its game time, top zones and what was going on.
- **Record what you verify (user, 2026-11-15, #478):** when verifying a graphical change (rendering, shaders,
  lighting, UI, animation, terrain or water look) or a behaviour change that shows on screen (movement, camera,
  physics, animals, machines, tools), record the test run with the `capture` skill (the in-game recorder, `--record`)
  and look at it yourself before claiming the fix, so visual problems are caught without the user. Keep it cheap:
  - short captures: 30-120 frames, `--record-every` 2-5 unless you're hunting flicker or jitter (then every 1);
  - look at the contact sheet first, one image via `mcp__deweyserver__media_get(path=<capture>)` (default 1280 px
    wide; use a smaller `max_width` when that's enough); open single frames (`frame=N`) only where something
    changes or looks wrong; don't read every frame;
  - for before/after, use `-Compare` (`--record-compare`) and read its % difference instead of comparing frames by eye;
  - skip it for changes with nothing to see (data, saves, pure logic covered by a scripted test), and when a
    recording would cost more than the check is worth (say so in the report);
  - never while the user playtests, and report the capture path (and the GIF for the user) with the result.
- **The camera shows what a test tests (user, 2026-11-16):** every test keeps its subject in view - the water, the
  structure being built, the pre-placed resource being mined, the machine running - and follows it when it moves (a stream
  running off, a dam break spreading, a player walking). Set the camera directly each frame (`world.camera.*`, as
  `--test-dam` does); placing the player every frame keeps the loading screen up. Fix a test that doesn't when you touch it.
- **Sim benches use the sandbox (user, 2026-11-16, #662):** `--sandbox` is a flat plane (no bake, rivers, sea, caves, rocks,
  trees; world ready < 1 s) that a test shapes with edits and parts. Use it for simulation work (fluids, machines,
  building); real world generation only for tests that need a real environment. Feature suites keep their real-world
  tests; full regressions wait for major milestones. Fluid perf: `--sandbox --test-fluid-load S` (an S x S m basin).
- Full regression (`tools/regression.sh planet`, `EXE=... SNAPSHOT=1`): ONE run only, after ALL fixes are in and
  the PRs are merge-ready (user, 2026-11-14), never between fix rounds; re-check failures with targeted runs
  first. It is split (#1028): each test is tagged headless or local (`tools/regression.sh planet --list`);
  `--only headless` runs on the deweyserver runner (nightly, or `gh workflow run selfhosted.yml -r <branch> -f
  regression=headless`), `--only local` (rendering, camera feel, perf and ms gates, menus) is what this PC runs. Erosion is the default terrain (#362, user 2026-11-15; `RAYLIBGAME_EROSION=0` = the old noise terrain for A/B): the first bench on a seed bakes ~80 s, cached with its save after that.
  Minor targeted tests are fine any time.
- Known flaky: `--test-progression` (power chain, #372), `--test-animals` flee distance (#373). Re-run before
  treating a failure as a regression.
- Perf numbers are invalid while another agent compiles or benches; compare A/B under the same conditions.
  `RAYLIBGAME_GPU_VERIFY_ALL=1` checks every GPU density block against the CPU (slow: ~20 FPS).
- Bench gotchas: `player.placeAt` sets the velocity itself (passing body-local coordinates to
  `surfaceVelocity` killed the player); start a bench clock after the loading screen (set start = -1, take
  `GetTime()` on the next frame); `BuildTool_Update` needs `player.nearest` set.
- Gate any test run on the build's exit code (`cmake --build ... > log && game.exe ...`): a failed build
  plus a stale exe runs a plain game window. Never start benches while the user playtests.
- Perf benches: `--bench 20 --bench-fly 60` (streaming; prints "bench: streaming:" counters, frame
  percentiles and camera drift), `--bench-walk`, `--day` (mid-morning sun; seed 42 starts at night). A
  plain planet `--bench` starts 40 m above the spawn looking along the ground (`--cam` = old offset).
  Other tests added 2026-11-13: --test-ledge, --test-tor (seed 1), --test-machines, --test-fall,
  --test-telescope-close (all in tools/regression.sh).
- Streaming benches (#262): `--preload R` (m, 0 = off), `--cache-verts M` (millions), `--preload-spawn` (time the
  spawn preload first, as a new game does: prints "bench: spawn preload:"), `--bench-fly-back` (fly out, then back:
  the cache should spare the return). "bench: preload:" prints ring/cache counts. The planet is the default
  bodies.txt (bodies_planet.txt is gone: a missing --bodies file falls back to the old asteroid system).
- Leave space/asteroid tests (--test-rocket, --test-ring, asteroid scenes) out unless the change is
  about them.
- Stop stale regression runs by command line (`*regression.sh*`, game.exe under `\Temp\tmp.*`);
  TaskStop does not kill the child scripts.
- Check on every test run 15-45 s after starting it (user, 2026-11-14): run it in the background and look at its
  output early, so an error or a stall is caught then, not after the whole timeout.
- **Every batch of tests has a hard limit of 15 minutes** (user, 2026-11-17: a wait loop hung for 8 hours), except a
  large-scale regression (`tools/regression.sh`). Wrap the whole batch, not only each test: `timeout 900 bash -c '...'`
  (per-test `timeout` inside it as well), and give every wait or poll loop its own deadline (`until <check> || [ $SECONDS
  -gt 900 ]`), never an open-ended `until grep ...`. Split longer work into several batches.
- `--test-planet-walk M` (#363, user 2026-11-14): M metres as 100 m legs at places drawn from the seed (teleports
  between legs, hovering while each place loads; every other leg prefers ground above 300 m). PASS = most legs walked
  and never under the surface or inside solid; a leg lost to a death or a 45 s timeout on rough ground is not a failure.
  Don't tailor the picks to walkable ground: the test is meant to sample the terrain.
- Slow frames: any bench frame over 20 ms prints "bench: slow frame" with its CPU passes over 1 ms and the time.
- Bench helpers: per-pass GPU + CPU timings print at the end of every --bench run (gputimer.*);
  `--native-window` (the game's borderless monitor-size window), `--no-msaa`,
  `--load` / `--at-save` (bench saves/world.sav from its own camera; copy a save into a temp dir
  with `resources`/`assets` junctions to bench it read-only), `--test-paint`, `--test-dig`
  (casts down from the player: the bench camera is not placed yet when early tests run).

## Git and PRs

- Branch per batch; stacked PRs merge only via the merge-async API, and GitHub does NOT retarget
  the next stacked PR - PATCH its base to develop before merging it.
- Merging needs the user's explicit OK, and only after their playtest of the branch, when they say so
  (user, 2026-11-14). Commit/push only when the work calls for it.
- Commit trailers: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` and
  `Claude-Session: <session link>`; PR bodies end with the Claude Code line and the session link.
- Save-breaking changes are fine (early development); note them in the PR. No compatibility shims
  needed, but old region formats (RGN1-3) still load.

## GitHub issues (required before starting work)

- Every task needs a GitHub issue BEFORE any work starts (code, shaders, data, tests, or a plan doc
  section). First search for a matching issue: `gh issue list --state all --search "<keywords>"`.
  If one matches, use it, and reopen it if it was closed but the problem is back. If none matches, file one
  with `gh issue create`: a title naming the problem or feature, a body with the problem/request, the cause
  if known, and the planned fix, plus labels (bug/enhancement plus area labels: rendering, shadows,
  terrain, water, physics, camera, ui, performance, tests, ...).
- A batch of several tasks gets one issue per task. Group related ones as sub-issues
  (`gh issue create --parent <n>`) under a tracking issue (e.g. #128 Graphics Pipeline Overhaul,
  #143 Plan 6).
- Reference the issue in commit messages (`(#N)`) and in PR bodies (a "Fixes" list with the numbers).
  PRs target develop, not the default branch (main), so "Fixes #N" does NOT auto-close: close the
  issues by hand (`gh issue close N --comment "Fixed in #PR"`) once the PR merges. Until then, label them
  `fixed in open PR`.
- Close an issue only when its fix is merged into develop (then board Done), never on playtest alone
  (user, 2026-11-14). A playtest confirmation before the merge removes `needs playtest`, adds
  `playtest confirmed`, and keeps the issue In review.
- Playtest findings are filed as GitHub issues (CurrentIssues.md was retired 2026-11-15, #481).
  Fixed-but-unconfirmed items carry `needs playtest`; deferred ones carry `deferred`.
- Bugs found along the way (a failing test, a regression) get their own issue, even when fixed in
  the same branch.
- Exempt: answering a question, read-only investigation, and triage itself (which comments on the
  issue it triages).

### Project board (keep it current)

Board: "RaylibGame Roadmap", https://github.com/users/DKoelzer01/projects/5 (project number 5,
owner DKoelzer01, id `PVT_kwHOAe7Hf84BrZzd`). Every open issue belongs on it, with all four fields set. The
board's auto-add workflow adds new issues itself, but sets no fields, so:
- **New issue:** right after `gh issue create`, set Status, Priority, Size and Kind. `item-add` either adds it or fails
  with "Content already exists" (auto-add was first); then look the item id up instead:
  `gh project item-add 5 --owner DKoelzer01 --url <issue url> --format json --jq .id` gives the item id, then
  `gh project item-edit --project-id PVT_kwHOAe7Hf84BrZzd --id <item id> --field-id <field> --single-select-option-id <option>`
  once per field. Find an existing item's id with
  `gh project item-list 5 --owner DKoelzer01 --limit 500 --format json --jq '.items[] | select(.content.number==N) | .id'`.
- **Status as work moves:** Backlog (not started) -> Ready (requirements clear, unblocked) -> In progress
  (starting work) -> In review (fixed/built in an open PR; also label `fixed in open PR`) -> Done (merged
  into develop, issue closed).
- **Priority:** P0 now (current plan / blocking), P1 next, P2 later, P3 someday. Follow the user's order in
  #158 (features) and #128 (graphics). Don't raise a priority on your own; propose it in the summary instead.
- **Size:** XS < 1 h, S ~ half a session, M ~ 1 session, L 2-3 sessions, XL a milestone of several
  (usually a tracking issue). Re-estimate when the scope turns out different.
- **Kind:** Feature, Bug, Performance, Tracking (an issue with sub-issues).
- Field and option ids (single-select):

  | Field | Field id | Options |
  | --- | --- | --- |
  | Status | `PVTSSF_lAHOAe7Hf84BrZzdzhpssRE` | Backlog `bf767b79`, Ready `2c8b23de`, In progress `4d99a4b6`, In review `ecd2f2aa`, Done `e280656f` |
  | Priority | `PVTSSF_lAHOAe7Hf84BrZzdzhpssSs` | P0 `e93b5bb8`, P1 `026eb06a`, P2 `2906a46d`, P3 `c36f5ea9` |
  | Size | `PVTSSF_lAHOAe7Hf84BrZzdzhpssSw` | XS `39a448a0`, S `032ae8bd`, M `6bd53109`, L `d7a23c47`, XL `ee6d351f` |
  | Kind | `PVTSSF_lAHOAe7Hf84BrZzdzhpssZg` | Feature `dcf904b0`, Bug `fae89327`, Performance `8172895f`, Tracking `d260109d` |

  If an id stops working, re-read them with `gh project field-list 5 --owner DKoelzer01 --format json`.
  ("Type" is reserved by GitHub, hence "Kind".)
- The gh token needs the `project` scope (added 2026-11-12). If a project command fails on scopes, ask the
  user to run `gh auth refresh -h github.com -s project`.

## Plans and trackers

- Plan 2 (early game, done): https://claude.ai/artifact/RL22kZ7VBMugAFgGRPS3cE
- Plan 3 (Terra tiers, space): https://claude.ai/artifact/XVw45AC6vwYitbdhaeP2xL - M21-M24 done (M24 on PR #28).
- Plan 4 (flowing water): https://claude.ai/artifact/27xL5j1zkirSe8pDZgPz59
- **Current plans:** Plan 7 (erosion-shaped terrain + all-voxel water, #351):
  https://claude.ai/artifact/EXix8zmVSS2WRUSLAsovJg - merged into develop (#371) on by default since the #362 gate (RAYLIBGAME_EROSION=0: off).
  Plan 8: https://claude.ai/artifact/LFnRaLaCPdgEPsuX7L6dJb (PR #374, branch m51-plan8, waits on the user's playtest).
  Roadmap trackers: #158 features, #128 graphics, #303 multiplayer; playtest findings #187, #250.
- Plan 7 details: with erosion (the default) the 3D band is off (RAYLIBGAME_3D=0/1 overrides). Gullies (ErosionGullies,
  RAYLIBGAME_GULLIES=0 A/B); glacial (RAYLIBGAME_GLACIAL); rivers on the bake grid, carve-only water (--water-stats raises
  nothing). Tests: --test-flow-sea, --test-pond, --test-pump-hand, --test-hover. Open questions for the user are on #351.
- Plan 6 (GPU terrain generation + GPU-driven drawing, approved 2026-11-12): https://claude.ai/artifact/VJdKtUZEpBg5uefnqewRyZ -
  issues #143 (steps #144-#150), under #128. Done: steps 1-4 + 7 (PRs #153, #157, #179); step 5 deferred (#148), step 6
  not needed (#149: GPU density is bit-identical, collision keeps the CPU field). One path on GL 4.3 (base-instance records).
- Plan 5 (art, life, world variety): https://claude.ai/artifact/5igt2VrsTyNhZ3BZ94teHX - A animal models/AI, B biomes,
  C 3D terrain, D stylized PBR terrain: all built and merged.

## Architecture notes

- Erosion (Plan 7): the bake runs on a thread behind a progress screen (the `gWorldGenWait` hook). Cache
  `<save>.erosion` (ERO6) keyed on `ErosionParams` + `kErosionBakeVersion`: bump the version on ANY bake code change,
  or stale bakes load. The 3D band's `Column3D` (mask, f) must match `terrain.glsl` / `terrain_density.comp`
  bit-exactly (the grid is 2 vec4 per node; the second is unused since the cliff terraces and lips went, #543);
  `--test-terrain-gpu` checks it. Plan 9 (#586): steep terrain comes from formations in the bake input, not post-bake steps.
- Sky light: `fragSky` from depth below the surface, times a short-ray openness probe (#369: terrain via chunk
  samples, parts via plot cells); `RAYLIBGAME_SKY_PROBES=0` for A/B. Moonlight scale defaults to 0.12 (#897; was 0.25)
  (`RAYLIBGAME_MOONLIGHT`). Night light = that moonlit-sky ambient alone: the fixed 1/255 base ambient fades out under the
  air (#897: it was a ~10/255 floor that hid the moonlight setting; `RAYLIBGAME_AMBIENT_BASE=1` = before). Within 40 m
  under the surface, cave mouths and shafts get at least 0.6 x (1 - depth/40 m) of their open sky (#898; `RAYLIBGAME_CAVE_MOUTH_SKY=0` = off;
  `RAYLIBGAME_SHOW_SKY=1` draws fragSky as grey).
  Captures at a spot: `--load --shot-cave-mouth` with `RAYLIBGAME_MOUTH_AT="lat lon h"` and `RAYLIBGAME_MOUTH_VIEW=in|out`.
- HDR (#131, scenetarget.cpp): the scene target is RGBA16F; world shaders still write display-encoded colour (pow 1/2.2),
  now unclamped. The composite decodes, applies auto-exposure (log-average of a 64x36 sample grid, never brightens:
  nights stay as tuned), tonemaps (default "knee": today's look up to 0.8, highlights roll off; RAYLIBGAME_TONEMAP=aces|clip)
  and re-encodes. RAYLIBGAME_HDR=0 = the old RGBA8 target; RAYLIBGAME_EXPOSURE_LOG=1 prints exposure/luminance.
- Physical sky (#142, skyphys.*): Hillaire LUTs over an Earth-sized air shell (terra's 25 km of air is too short a path for red
  sunsets); the sky and the haze colour (not its amount: ~35 km visibility) come from it, 35% toward the old palette. On by
  default (user, 2026-11-18); RAYLIBGAME_SKY_PHYS=0 = the old analytic sky and haze. Tint/exposure: RAYLIBGAME_SKY_PHYS_TINT / _EXPOSURE.
- Multiplayer phase 1 (#304, PRs #435, #438, #440, #442): gameplay runs in fixed 60 Hz steps (`simclock.h`:
  `SimTime()`/`SimDt()`, never `GetTime()` in sim code; the camera interpolates the player by `SimAlpha()`; mouse look is
  per frame). Gameplay reads `PlayerInput` (`playerinput.h`; `PlayerInput_Step()` in steps, `PlayerInput_Frame()` for the
  per-frame tools), never `IsKey*`/`IsMouse*`. Per-player state lives in `PlayerSession` (`player`, `inventory`, ... are
  references to `LocalSession()`'s members). Replication refers to `NetId`s (`netid.h`): bodies by creation order, ring
  rocks by name hash + index, animals/vehicles/stations allocated and saved. Benches can find a frame with no sim step:
  wait a little after a teleport before reading `grounded`. Transport (#310): `src/engine/net.*` over GameNetworkingSockets,
  only with `-DENABLE_GNS=ON` after `tools/build_gns.sh` (~2-7 min) builds protobuf 3.21 + GNS 1.4.1 (pinned, w64devkit, BCrypt) into the gitignored `deps/gns/install` (the default GNS_ROOT;
  +4 MB exe), else stubs; `--test-net` (loopback echo) SKIPs without it. `gNetRole` (net.h): Offline/Host roll gameplay
  and world events (meteor showers are one world-wide event, #316), a Client applies what the host sends. Dedicated server
  (#311): `--host <save> [--port N] [--seed N] [--autosave-minutes N]` loads the save (keeping its seed) or creates it (the
  --seed, else a printed random seed), autosaves every 30 min (0 = off) and on shutdown, runs `--headless` and answers "join" with
  the JoinInfo + shower (src/space/hostserver.*); `--test-host` (loopback client; SKIP without GNS), `--test-shower`.
- Inventory: 45 grid cells + 10 hotbar cells (indices 45-54). `Hotbar_Bind` moves a stack (items must be carried);
  `Inventory::preferred` = the selected hotbar cell is used first.
- Caves (Plan 10, #615): all cave carving is `PlanetField::caveDistance` (voxel.cpp) and its bit-exact twin in terrain.glsl; change both
  together. `--cave-stats R [D]` (counts, the trench metric, a density hash that refactors must keep), `--test-caves` (fly into the nearest
  mouth). Graph tunnels (step 2 #995, on by default since step 3 #1022; `RAYLIBGAME_CAVES=0` = the old noise caves): per-256 m-cell passages,
  per-block segment lists (6 vec4 each, to the GPU at binding 4; <= 32 a block, else the CPU path); `RAYLIBGAME_CAVE_ABUNDANCE` 0-2 (the New
  World slider to come). Shapes (#1022): the segment's bore is the skeleton; heavy 3D noise only grows the wall inward (none inside half the
  bore), looked up squeezed radially (x0.1) so there are no voids behind the faces; ~9 m rock forms with ~3.2 m hard-edged creases
  (`RAYLIBGAME_CAVE_WALLS=0` = the first merged ~3 m / ~1 m walls). Mouths only in real faces (5 m+ rise within 10 m, no tors in
  front). Karst (#1052): dolines in karst country are sinks (shaft, river passage down the drainage, keyholes, mazes, resurgences);
  `RAYLIBGAME_CAVES_KARST=1 --test-caves` finds one, `RAYLIBGAME_CAVES_INSIDE=1` starts in a passage. Noise factors include each noise's own base frequency (1/70-1/700: forgetting it makes every term constant). Captures:
  `RAYLIBGAME_CAMERA_LIGHT=1` (a light at the eye), `RAYLIBGAME_CAVES_VIEW=front` (eye-level frontal approach), `RAYLIBGAME_CAVES_MOUTH="x y z"`
  (that mouth, at any abundance). In GLSL, store every a*b+c handed to a function (xdiv, xsqrt, fractal) in a `precise` variable first.
- Terrain internals still worth knowing: `PlanetField::rockFeature` (tors, boulders) and `terrain3D` (the 3D band);
  SSE2 `SimplexNoise::noise4` is bit-identical to `noise()` (`--profile-density` checks it); static water feeds only
  dug cells, and ponds are taken over by the fluid sim (`PlanetField::claimPond`) and become finite.

## Decisions (playtest 2026-11-11 and after)

Player and building (2026-11-14)
- Player collision capsule 1.7 m tall (5 spheres at most 0.23 m apart; `RAYLIBGAME_PLAYER_HEIGHT` for A/B); 2 m
  doorways and tunnels must pass.
- Parts face automatically (from the hit face and view) until the player presses R/V.
- The Help tab is a Factoriopedia-style database; its "Discovered only" choice is saved per world.

Controls and tools
- Hotbar starts empty; an empty slot is bare hands (no "Hands" item). Loose fillers (gravel, sand,
  regolith) go into the hand from the hotbar and RMB places them; dirt/stone/etc. are build blocks.
- LMB breaks terrain and built parts; parts need a 0.45 s hold (also outside build mode). A part in
  front of the ground blocks the terrain tool (no digging/tilling/placing through structures).
- One swing per cooldown (`1/strokesPerSecond`, chops 0.3 s); the cooldown survives button releases,
  the hand animation lasts exactly the cooldown, and a bar under the crosshair shows it (amber when
  breaking a part).
- Digs: the carve sphere sinks half its radius along the view; the highlight is a depth-tested circle
  on the surface (the sphere's mouth), never a see-through sphere.
- Multi-cell parts placed against a face slide out until no cell is behind that face.
- Tools are held upright; hand bob = one side swing per 1.8 m walked.

Materials and progression
- Placed ground keeps its material: delta chunks carry a material channel (`DeltaChunk::m`), regions
  save as RGN4. Placing uses the filler in hand first.
- Stone tools, their repair and the tier-0 quest use stone; loose stones (tree-system variant 9) are
  picked up with LMB for 1-2 stone.
- Gems are rare: 1 in 3000 rock strokes (quartz 1 in ~1900); veins come later. Tests use
  `terrainTool.gemBoost`.

World and rendering
- Lakes = whole filled basins of the river flood, level held 3 m under the land sampled around them;
  water region = basin cells + 60 m margin; a capped (4 m) bank in the outer margin. Rivers in a lake
  build no banks. Chunks find lakes by BFS distance; dry water-lattice points sit under their ground.
- Ring rocks are cut at 50 km true distance (the mid pass's 60 km far plane cut them by depth).
- Clouds: a shell 1.5 km up in atmosphere.fs, in the planet frame, drifting with the wind; cover
  never below ~0.22 so a few are always out.
- Perf: chunks draw nearest-first (early-Z), CPU copies of chunk meshes are freed after upload.
  Scene::KickGpu (glFlush) at pass boundaries: without it CPU and GPU ran serialized (+25% FPS). GPU timer = elapsed queries.
  Reverse-Z (#135, scenetarget.*): the 3D world draws into an offscreen float-depth target (window sample count, glClipControl
  0..1, infinite far, GL_GREATER) in ONE pass (no far/mid/near depth ranges or clears); RAYLIBGAME_REVZ=0 = old passes.
  Terrain arena (#144, terrainarena.*): one multi-draw-indirect per body per pass; RAYLIBGAME_ARENA=0 = per-chunk draws.
  (Old note "remaining time is present/composition" was wrong: it was CPU/GPU serialization.)
- OpenGL 4.3 is required (user, 2026-11-13, #185): raylib is built for a 4.3 core context, so without it there is no window
  (main.cpp shows an error box). There is no 3.3 renderer; RAYLIBGAME_ARENA=0 / RAYLIBGAME_REVZ=0 / RAYLIBGAME_GPU_TERRAIN=0
  are A/B and debugging switches on a 4.3 context, not fallbacks for older hardware.
- Terrain lattice 0.25 m (kVoxelSize; RAYLIBGAME_VOXEL_1M builds the old 1 m one): digs carve exactly the
  previewed sphere (--test-dig-shape IoU ~0.85, pass bar 0.83: #425). Delta chunks hold sparse 4^3 bricks; regions are RGN5.
- Shadows (2026-11-12): two cascades (near +-32 m, far +-200 m, 2048 each), bilinear PCF, normal-offset
  lookup + slope bias; trees and bushes draw into both (tree vertex shader + depth.fs). Grass casts none
  by choice (blades thinner than a near texel).
- Foliage (user, 2026-11-12): grass = 3 instanced meshes (tuft, 14-blade clump, 24-blade dense clump;
  variant rides in scale + 10*variant), 7 per m^2, darker lower half, drawn to 100 m; trees to 2 km,
  bushes 300 m. The user trades a few FPS for view distance (walk bench ~220 FPS). Draw distances are
  runtime values (--grass-dist, --tree-dist); the terrain canopy tint follows treeDistance.

GL features and streaming (2026-11-13, stack #287-#293)
- The driver hands out a 4.3 core context with newer features as extensions; src/engine/glcaps.* (`GL_Caps()`)
  detects them at runtime and each path falls back when missing. `RAYLIBGAME_GL_LIMIT=43` forces the 4.3 paths.
- A/B switches: RAYLIBGAME_STREAM_OLD=1 (old streaming), RAYLIBGAME_GPU_READBACK_COPY=1 (old readback),
  RAYLIBGAME_UPLOAD_STAGING=0 (direct uploads), RAYLIBGAME_TREE_CULL=0 (no GPU tree culling),
  RAYLIBGAME_INDIRECT_COUNT=1 (opt-in multi-draw-indirect-count with chunk_compact.comp: measured no gain).
- Streaming: <= 96 requests per frame, frees capped at 256/frame, nearest-first uploads, telescope keeps the
  zoom-1 selection resident. GPU readback = persistent mapped ring; uploads go through a 16 MB persistent write
  ring with per-upload fences.
- Preload ring + chunk cache (#262: #294-#296, 2026-11-13): every node the LOD would pick from a camera within
  300 m of the player along the ground (0-8 m over it; nodes touching +-12 m of the height map) is built and kept
  resident, not drawn. A new or loaded game shows "Preparing spawn" until it is built (~8 s, ~38k nodes,
  ~27 M vertices = ~2 GB VRAM on seed 42); afterwards it re-centres every 75 m in spare time only (the view's
  unused request budget, a short queue, 1 ms of uploads, one 8-block GPU batch a frame) and pauses above
  15 m/s. A ball-shaped ring reached 300 m down into caves (415k nodes, 80 M vertices) - keep it on the ground.
  Unwanted built nodes go to an LRU cache (8 M vertices / 40k nodes, ring and view not counted); the arena is
  grown once behind the spawn screen so it never regrows (GB copy) mid-walk. WorldIdle() ignores ring builds.
  The cache also evicts while the arena is over 85% full, and a ring install reserves nodes x 800 + cache + 4 M
  vertices up front (#300: a late regrow was a 30 ms frame closing the telescope).
- Terrain detail (#297): a settings choice Low/Medium/High/Ultra = LOD factor 2 / 2.5 / 3 / 4, default Medium. It
  only sharpens mid-range terrain (the planet is selected to the horizon already). Above 15 m/s against the nearest
  body (back under 10) the factor drops to 2: at 3, flying 60 m/s cost 4.5 ms GPU streaming a frame (115 FPS).
  Seed 42 at Medium: still 267 / walk 250 / fly 262 FPS. Benches: `--lod F` sets the detail, `--lod-fast F` the
  flying factor (99 = no drop); benches do not read settings.txt.
- Tree shadows: caster sets per cascade chosen by distance to the sun axis (1.26 -> 0.09 ms); main view
  culled on the GPU (tree_cull.comp + glDrawArraysIndirect).
- Fly camera follows the body whose sphere of influence it is in, rotating with it within 2 radii (#127).
- Gameplay in the same stack: boilers need fuel + water (pumps fill them), generators are demand driven;
  fall damage safe below 6 m/s (8.5 with servo legs); animals walk on the density surface with ledge
  reach/detours. Ravines are generated (voxel.cpp), not streaming holes. There is no start ship.

Water (Plan 4, user answers)
- Sea, rivers and large lakes (> ~1 km^2) are infinite; rivers can be diverted; small lakes/ponds
  are finite once disturbed. The sea can flood a channel dug inland below sea level.
- Multi-level (cave) water must work. Buckets AND pipes/pumps. All of F1-F6 before the next playtest.
- Decided (2026-11-11): F2 is a sparse 3D cell simulation, Terraria-style, 1 m voxels (fill amount
  per cell; down, then sideways; a slow levelling pass for U-bends / cave pressure), surface meshed
  with marching cubes on the fill amounts. Replaces the per-column span design.
- Fluid cells are 0.25 m by default (user, 2026-11-17, #714); `RAYLIBGAME_FLUID_CELL=1` is the old 1 m grid (opt-in).
  A save records its cell size; sites saved at another size are dropped on load (their simulated water starts over).

Deferred by the user: moon mining launch (#116), motion sickness (#117; comfort options are #168). The history of
the earlier playtest batches is in git history and the closed issues.

## Art pass (branch art-pass-5, 2026-11-12, Plan 5 light items)

- Shaders take `#include "file"` (HotShader expands it, relative to the shader, and watches the file).
  `resources/shadow.glsl` = the two-cascade sun shadow lookup + headlamp for lighting/trees/grass/parts/water.
  Shadow maps use hardware compare (sampler2DShadow, linear): never sample them as sampler2D (debug
  previews call SetShadowCompare(id, false) first). Near cascade blends into the far one over 8% of its box.
- Parts receive shadows + headlamp (ctx.bindSunLighting binds slots 6-7 before plot and water draws).
- Part surface patterns: id in the part mesh's vertex-colour alpha (plots.cpp PartPatternOf), drawn in
  section-local plot metres by parts.fs (stone, brick, planks, metal, glass, tiles, marble, gems, logs,
  lamp glow, machine casing, conveyor, cable, farmland, hull).
- Item icons: src/space/itemicons.* paints every item at startup (one atlas); --shot-icons.
- Particles: src/engine/particles.* (body-local, pooled, rlgl billboards); emitters in terraintool
  (dig/place/chop/bucket), rocket.cpp Exhaust, src/space/ambience.* (rain splashes, fireflies, machine
  smoke, campfire flame); --shot-particles.
- Lighting review leftovers for the user: tonemapping (ACES), hemisphere/SH ambient + SSAO, GGX specular.
  Staggering the far cascade was rejected (planet spin makes a stale map shimmer).
- Don't run screenshot benches while another game instance runs (e.g. a background bisect): chunk
  loading slows and shots show a half-loaded world.

## Plan 5 answers (user, 2026-11-12) - all approved, order B -> C -> A -> machine models -> D
- A: procedural code-built animal models (glTF slot later); predators only chase (no kills yet);
  machine models after the animals.
- B: all 10 biomes (jungle, savanna, desert dunes, badlands/mesas, temperate forest, meadow, swamp,
  taiga, tundra, alpine); regional size (2-5 km, several near spawn); save-breaking is fine.
- C: dramatic terrain tied to biomes, with a good variety of flat building spots as well as
  interesting ones; floating islands rare; ravines and cave mouths in the same pass.
- D: some PBR; stylized look - higher resolution than Minecraft but fairly unrealistic.

## Generic fluids (user, 2026-11-12) - on m25-biomes / PR #30
- One fluid model, not per-fluid copies: kinds come from resources/data/fluids.txt (water = 0, crude_oil;
  room for 8, e.g. lava later: colour, opacity, glow, viscosity, min spread, bucket item). The kind rides in
  cell flag bits 5-7; kinds never mix (a fluid enters only empty cells or its own kind). Meshes per kind,
  drawn by water.fs (fluidColor.a > 0 = a non-water kind). Saves: fill + 256 x kind; sources on a "springs" line.
- Sources (Fluid::AddSource): a rate (m^3/s) and a head; fill the column's surface up to the head, then rest
  (equilibrium, bricks sleep) or overflow into a stream. Slow sources add in drops of >= 8/255 (thinner
  specks soak away: kEvaporate). Place them on the real ground (probe density from ~8 m up: the mesh can
  stand metres above the height map).
- Natural sources (ambience.cpp kSourceTypes): oil seeps in warm biomes (carved basin, head under the rim:
  pools) and water springs in wet hill biomes (head 1.2 m over the ground: streams). --test-oil, --test-spring.
- Follow-up ideas from the user: express rivers/lakes (and later magma) as sources too; springs show as
  thin films on steep scree (consider a carved channel or a lower film limit for springs).
