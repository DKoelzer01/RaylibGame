# CI (#1016)

Three places run checks, each for what it can do safely.

| Where | Workflow | Triggers | What |
| --- | --- | --- | --- |
| GitHub-hosted runners | `workflows/pr-checks.yml` | every PR into `m84-integration` / `develop` / `main` (forks too), pushes to `m84-integration` / `develop`, manual | static checks, Windows build (w64devkit) + `game_tests`, Linux headless build + `game_tests` + a dedicated-server smoke test |
| Self-hosted runner on deweyserver (labels `raylibgame-linux`, `gpu`) | `workflows/selfhosted.yml` | pushes to `develop` / `m84-integration`, nightly 05:00 UTC, manual. **Never** `pull_request` / `pull_request_target` | Linux server build, `game_tests`, the GPU bit-exact check (`--test-terrain-gpu` on EGL), headless scenarios (fluid, craft, skills, world), `tools/test_server.sh` + `tools/test_join.sh`, the server tarball as an artifact |
| The user's PC | none (agents and the user) | by hand | rendering, visual baselines and captures, perf timing (quiet windows), the full `tools/regression.sh` |

Why the split: the repo is public, so hosted runners are free and safe for code from anyone (read-only token, no
secrets, throwaway VMs). The self-hosted runner sits on the home server, so it only runs code that is already on a
protected branch (or that the owner starts by hand). Rendering and timing need the real GPU and display and a quiet
machine, which only the PC has.

## PR checks (hosted)

- **Static checks** (seconds, `tools/ci/static_checks.sh`): NetMsg ids unique and 1..255; a PR that adds, removes or
  renumbers a NetMsg or bumps `kNetProtocol` says so in the job summary (a NetMsg change without a bump is a warning);
  erosion bake code changed without a `kErosionBakeVersion` bump is a warning; a new `gh project` command in `tools/` is
  an error. It also decides whether the builds run: a change touching only Markdown, `docs/`, root images or `.claude/`
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

The runner (Docker, on deweyserver) needs: the packages of `tools/linux/Dockerfile.build` (build-essential, cmake, git,
curl, libssl-dev, libegl1) plus ccache, the NVIDIA container runtime with `NVIDIA_VISIBLE_DEVICES=all` and
`NVIDIA_DRIVER_CAPABILITIES=graphics,compute,utility` (the GPU step fails loudly on the null GL), and the labels
`raylibgame-linux` and `gpu`. Builds, the erosion bake cache and ccache persist in `<runner work>/RaylibGame/ci-persist`.
The runner has no LAN access: the server package reaches the Pterodactyl server through the host-side script
`tools/ci/deweyserver-update-server.sh` (see `docs/ci-deweyserver.md`).

## Rules for workflow changes

- Pin every action by full commit SHA with the version in a comment; `permissions` stay minimal (`contents: read`).
- Never add `pull_request` or `pull_request_target` to a workflow that runs on the self-hosted runner.
- Lint with [actionlint](https://github.com/rhysd/actionlint) (`.github/actionlint.yaml` lists the runner labels and
  variables). `.github/` is in CODEOWNERS: workflow changes need the owner's review.
