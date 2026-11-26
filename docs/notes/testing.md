# Testing: full notes

Detail behind CLAUDE.md's Testing section (#1478). Read the parts for what you're testing.

- **Where tests run (user, 2026-11-19, #1016; details in `.github/CI.md`).** Run each check where it belongs; this PC is
  only for what can't run elsewhere or has to be validated here:
  - **GitHub-hosted CI, every PR (required check "CI ok"):** the Windows build (w64devkit), the Linux server build,
    game_tests, static checks (NetMsg ids / kNetProtocol, bake-version reminder). Don't rebuild and rerun unit tests
    locally just to check a PR: push and read CI. Never ask for a merge while "CI ok" is red or pending.
  - **Self-hosted on deweyserver (GPU, isolated from the LAN):** `--test-terrain-gpu` bit-exactness, headless scenarios
    (fluids, craft, skills, world), the multiplayer server and join tests, the Linux server package, and the headless
    half of the regression (`tools/regression.sh planet --only headless`, #1028; nightly, or add `-f regression=headless`).
    Runs on pushes to staging/develop and nightly; for your own branch, dispatch it: `gh workflow run
    selfhosted.yml -r <branch>` (one GPU run at a time, so batch your changes). Use it instead of running those tests here.
    Three ephemeral runners (raylibgame-runner@1-3, 4 CPU/10 GB each) share one `ci-persist` volume; each job takes a
    flock'd slot (`ci-persist/slotN`) with its own out/ and run/ (#1411: shared dirs made join tests stall each other).
    The GPU gate (#1394) runs the GPU checks only when shader, terrain, fluid, render or build files change. Static
    checks also enforce scenario placement (check 5) and the renderOrigin ratchet (check 4) (#1420, #1423).
    The runner is slower than this PC: fluid ticks lag sim time there, so fluid scenarios measure in fluid time, not
    sim seconds (#1431); a cold shader cache costs ~37 s on it (#1047). Superseded runs may be cancelled (user's standing
    OK, 2026-11-21); a PR whose evidence is a new runner scenario waits for that run's result before merging (#1165).
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
- **Perf budgets (#1690):** four reference scenes tagged `budget` (`perf/planet-view`, `perf/creatures-60`, `perf/water-shore`,
  `perf/factory-base`: 5 s warm-up, ~20 s on a fixed camera path; metrics "frame mean", "frame p95" and the scene's own cost).
  The runner's "Perf budgets" step runs them 3 times under `gpu_lock.sh exclusive` and `tools/ci/perf_budgets.sh` judges the
  medians against `tools/ci/perf_budgets.txt` (runner-only limits: today's median x1.25 for means, x1.35 for p95; 0 = not judged
  yet). Up to 20% over a limit only flags (a warning); more than 20% over blocks once `PERF_BUDGETS_ENFORCE=1` (user ruling
  2026-11-23: the first week is report-only, and the step is `continue-on-error`). A lock timeout reports "not judged". A PR may
  raise a limit in the same diff with a reason. Not done yet: rerunning only the failing scene before failing.
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
- A bench camera inside terrain looks like a rendering bug (black tubes, missing faces: #1428 was a view 26 m inside a
  hill). Before blaming the renderer, check the eye is in open air (#1442 adds an inside-solid flag to benches/scenarios).
- Bench and test fixes go into the bench, not the feature (adapt the test: jump, turn, lift the camera).
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
