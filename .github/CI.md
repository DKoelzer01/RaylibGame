# CI (#1016)

Three places run checks, each for what it can do safely.

| Where | Workflow | Triggers | What |
| --- | --- | --- | --- |
| GitHub-hosted runners | `workflows/pr-checks.yml` | every PR into `m84-integration` / `develop` / `main` (forks too), pushes to `m84-integration` / `develop`, manual | static checks, Windows build (w64devkit) + `game_tests`, Linux headless build + `game_tests` + a dedicated-server smoke test |
| Self-hosted runner on deweyserver (labels `raylibgame-linux`, `gpu`) | `workflows/selfhosted.yml` | pushes to `develop` / `m84-integration`, nightly 05:00 UTC, manual. **Never** `pull_request` / `pull_request_target` | Linux server build, `game_tests` and the server package (their own job), the GPU bit-exact check (`--test-terrain-gpu` on EGL), headless scenarios (fluid, craft, skills, world), `tools/test_server.sh` + `tools/test_join.sh`, the server tarball as an artifact; nightly (and on dispatch with `regression=headless`) the headless half of the regression, `tools/regression.sh planet --only headless` |
| The user's PC | none (agents and the user) | by hand, the midnight regression | rendering, visual baselines and captures, perf timing (quiet windows), the local half of the regression (`tools/regression.sh planet --only local`; the full run when asked) |

Why the split: the repo is public, so hosted runners are free and safe for code from anyone (read-only token, no
secrets, throwaway VMs). The self-hosted runner sits on the home server, so it only runs code that is already on a
protected branch (or that the owner starts by hand). Rendering and timing need the real GPU and display and a quiet
machine, which only the PC has.

## PR checks (hosted)

- **Static checks** (seconds, `tools/ci/static_checks.sh`): NetMsg ids unique and 1..255; a PR that adds, removes or
  renumbers a NetMsg or bumps `kNetProtocol` says so in the job summary (a NetMsg change without a bump is a warning);
  erosion bake code changed without a `kErosionBakeVersion` bump is a warning; a new `gh project` command in `tools/` is
  an error; a `kNetProtocol` not reserved in `tools/ci/protocol_reservations.txt` (or one marked "(unused)") is an error
  (#1480). It also decides whether the builds run: a change touching only Markdown, `docs/`, root images or `.claude/`
  skips them.
- **Windows build + unit tests (w64devkit)**: w64devkit 2.0.0 (GCC 14.2, the shipping toolchain, SHA-256 pinned),
  MinGW Makefiles, `-DENABLE_GNS=ON`, targets `game` and `game_tests`, then `game_tests.exe` from the repo root.
- **Linux build + unit tests (headless server, GNS)**: Ubuntu 24.04, `RAYLIBGAME_HEADLESS=ON` (the server's raylib
  backend), `-DENABLE_GNS=ON`, `game` + `game_tests`, then `game_tests` and `tools/test_server.sh` on the null GL.
- **CI ok**: passes when every job above passed or was skipped. **Make this the required check** (branch protection on
  `m84-integration` and `develop`: "Require status checks to pass", check name `CI ok`). It is reported on every PR,
  also when the builds are skipped, so docs-only PRs are not left waiting.

Caches (actions/cache, 10 GB per repo, least recently used goes first): ccache per branch (a PR falls back to its base
branch's, which the push runs keep warm), the toolchain downloads, the raylib / googletest sources
(`tools/ci/fetch_deps.sh`, keyed on the tags in CMakeLists.txt) and the GameNetworkingSockets install (keyed on
`tools/build_gns.sh`). A cold run builds GNS (protobuf) once per OS; later runs only compile what changed.
`concurrency` cancels a PR's older run when it is pushed again.

Timings (PR #1024, 2026-11-20): cold (nothing cached) Windows ~11.5 min (GNS 3.3 min + build 7.4 min), Linux ~5.5 min
(GNS 1.5 min + build 2.8 min); warm (all cached, nothing to recompile) Windows 35 s, Linux 51 s. A PR touching a few
files recompiles only those plus the links.

## Self-hosted (deweyserver)

Off until the repository variable is set: Settings > Secrets and variables > Actions > Variables > New repository
variable `SELF_HOSTED_ENABLED` = `true`. While it is unset, the job is skipped at once (never queued for a missing
runner). `NIGHTLY_BRANCH` (optional) picks the branch the nightly run builds, default `m84-integration`; a schedule only
fires from the default branch's (`main`'s) copy of the workflow, so the nightly starts once this file reaches `main`.

The runners (three ephemeral GPU runners in Docker on deweyserver, all with the labels below; a job takes whichever is
free) need: the packages of `tools/linux/Dockerfile.build` (build-essential, cmake, git,
curl, libssl-dev, libegl1) plus ccache, the NVIDIA container runtime with `NVIDIA_VISIBLE_DEVICES=all` and
`NVIDIA_DRIVER_CAPABILITIES=graphics,compute,utility` (the GPU step fails loudly on the null GL), and the labels
`raylibgame-linux` and `gpu`. Builds, the erosion bake cache and ccache persist in `<runner work>/RaylibGame/ci-persist`.
**Slots (#1410):** jobs run at once on several runners, which may share that volume, so each job claims a slot,
`ci-persist/slotN` with its own `out/` (the build and binary) and `run/` (saves, the bake, test scratch), held by a
`flock` for the whole job (`tools/ci/selfhosted_env.sh`); ccache is shared. A job's log starts with its runner, the
volume's device and inode, and its slot. (Before, concurrent jobs shared one `out/` and `run/` and overwrote each
other's binary and saves mid-run: the failures of 2026-11-21 from ~15:30 UTC, #1407, #1406.)
After the tests, the GPU job's summary lists every `-> FAIL` line of its logs marked **known (#N)** from
`tools/ci/known_failures.txt` or **NEW** (`tools/ci/mark_failures.sh`, #1479; the runner-check skill's matching).
The runner has no LAN access: the server package reaches the Pterodactyl server through the host-side script
`tools/ci/deweyserver-update-server.sh` (see `docs/ci-deweyserver.md`).

Three jobs, each after the one before (#1028), each on any free runner: `package` (build, `game_tests`, the server
tarball as an artifact), `linux-gpu` (the GPU bit-exact check, the headless scenarios in two steps of 45 min each (#1367: fluid, frames and factory; then the rest), `test_server.sh` / `test_join.sh`,
whose host and client logs go into the `selfhosted-test-results` artifact) and `regression` (below; nightly or on
request). **GPU steps only when needed (#1392):** on a push, `tools/ci/gpu_changes.sh` diffs the push against its base;
the GPU bit-exact check (`--test-terrain-gpu`) and the fluid CPU = GPU parity run only when a shader, the terrain
generation (voxel, noise, erosion, rivers, strata, the GPU terrain, meshing, world-gen data), fluid, render-target or
build file changed. The headless scenarios and the server and join tests always run; a schedule or a manual run runs
everything, as does a push whose base isn't in the checkout's last 50 commits. The step summary says which and why. The
host-side fetch script takes the newest push run whose `package` job succeeded, so a failing test elsewhere in the run
no longer holds a server build back. Each later job rebuilds incrementally first (a no-op unless another branch's run
used the build folder in between).

## The regression, split (#1028)

`tools/regression.sh` tags every test: **H** (headless) checks game state the Linux headless/EGL build computes the same
way (simulation, fluids, building, saves, inventory and crafting, multiplayer, terrain logic); **L** (local) needs the
Windows client on the PC: rendering and judged screenshots (`--crack-check`, torch and lamp luminance), camera and
movement feel (`--test-stand`, `--test-flyfollow`, `--test-ledge`, `--test-tor`, `--test-viewmodel`), the perf benches,
verdicts that gate on milliseconds (`--test-persist`, `--test-factory`, `--test-animals`, `--test-flow` and
`--test-flow-sea`'s step-time gate, `--test-telescope-close`) and the menus. `tools/regression.sh planet --list` prints
the table; a comment on each L line says why. Tests still timed on the wall clock (phases and timeouts on `GetTime()`)
are H with a "wall clock" note: a much slower machine could change them. Four tests gave other verdicts on the runner
than on the PC because the runner is slower, and are L until they have sim-time versions: `--test-suffocate` and
`--test-smallbody` (wall-clock deadlines), `--test-structures` (a screenshot tour waiting for streaming) and
`--test-leak` (the flag-mode fluid job runs at the machine's speed; scenarios use `RAYLIBGAME_FLUID_LOCKSTEP=1`).

- `tools/regression.sh planet --only headless`: the runner's half (the `regression` job of `selfhosted.yml`, after the
  checks; 71 tests, ~85 min on the PC, ~2.5 h on the runner, which runs at about half the PC's frame rate). It runs
  as 4 jobs (`--shard 1/4` .. `4/4`, every 4th test; ~35-40 min each, 90 min timeout), so queued runs get a GPU
  runner in between; never on pushes.
- `tools/regression.sh planet --only local`: the PC's half (30 tests, ~30 min instead of ~2 h).
- No `--only`: every test, as before the split.

The job runs in the persistent run folder (`RUNDIR`, so the seed-42 erosion bake is cached), writes a table of every test
(verdict, seconds) to the job summary (`tools/ci/regression_summary.sh`) and uploads `regression.log` as the artifact
`regression-headless-log-<shard>`. **Known failures** are listed in `tools/regression_known.txt` (the test's arguments and its
issue): the summary shows new failures first, then the known ones with their issue, and names a listed test that passed
(fixed? take its line out). The job fails only on new failures. The PC's half can use the same summary:
`tools/regression.sh planet --only local | tee r.log; tools/ci/regression_summary.sh r.log tools/regression_known.txt`.
Start it by hand on a branch:
`gh workflow run selfhosted.yml -r <branch> -f regression=headless` (`-f ref=<branch>` builds another branch with this
branch's workflow; `-f regression_sets="planet multiplayer"` adds the multiplayer set).

**The nightly trigger.** A `schedule` fires only from the default branch's (`main`'s) copy of the workflow, so the
`cron` in this file does nothing until `selfhosted.yml` is on `main`. Until then the nightly is started from outside:
a dispatch (`gh workflow run selfhosted.yml -r m84-integration -f regression=headless`, or the REST call
`POST /repos/DKoelzer01/RaylibGame/actions/workflows/selfhosted.yml/dispatches` with
`{"ref":"m84-integration","inputs":{"regression":"headless"}}`) needs a token with "Actions: write" on this repository.

## Rules for workflow changes

- Pin every action by full commit SHA with the version in a comment; `permissions` stay minimal (`contents: read`).
- Never add `pull_request` or `pull_request_target` to a workflow that runs on the self-hosted runner.
- Lint with [actionlint](https://github.com/rhysd/actionlint) (`.github/actionlint.yaml` lists the runner labels and
  variables). `.github/` is in CODEOWNERS: workflow changes need the owner's review.
