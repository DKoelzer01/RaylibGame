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
- No python on this machine. Use the Edit tool (not sed) for strings containing `\n`; `perl -0`
  multi-line substitutions silently fail on CRLF files - prefer Edit.

## Testing

- Light tests per change: the scripted benches for the feature touched and its neighbours
  (`game.exe --seed 42 --bench 3 --test-<name>`), plus a screenshot (bench.png) when visuals change.
- Full regression (`tools/regression.sh planet`, `EXE=... SNAPSHOT=1`) only at milestones/PRs - and,
  per the user (2026-11-11), NOT until the current playtest feature set (incl. flowing water) is done.
  Minor targeted tests for PRs are fine.
- Leave space/asteroid tests (--test-rocket, --test-ring, asteroid scenes) out unless the change is
  about them.
- Stop stale regression runs by command line (`*regression.sh*`, game.exe under `\Temp\tmp.*`);
  TaskStop does not kill the child scripts.
- Bench helpers: per-pass GPU + CPU timings print at the end of every --bench run (gputimer.*);
  `--native-window` (the game's borderless monitor-size window), `--no-msaa`,
  `--load` / `--at-save` (bench saves/world.sav from its own camera; copy a save into a temp dir
  with `resources`/`assets` junctions to bench it read-only), `--test-paint`, `--test-dig`
  (casts down from the player: the bench camera is not placed yet when early tests run).

## Git and PRs

- Branch per batch; stacked PRs merge only via the merge-async API, and GitHub does NOT retarget
  the next stacked PR - PATCH its base to develop before merging it.
- Merging needs the user's explicit OK. Commit/push only when the work calls for it.
- Commit trailers: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` and
  `Claude-Session: <session link>`; PR bodies end with the Claude Code line and the session link.
- Save-breaking changes are fine (early development); note them in the PR. No compatibility shims
  needed, but old region formats (RGN1-3) still load.

## Plans and trackers

- Plan 2 (early game, done): https://claude.ai/artifact/RL22kZ7VBMugAFgGRPS3cE
- Plan 3 (Terra tiers, space): https://claude.ai/artifact/XVw45AC6vwYitbdhaeP2xL - next: M24 moon/ring tiers.
- Plan 4 (flowing water): https://claude.ai/artifact/27xL5j1zkirSe8pDZgPz59
- `CurrentIssues.md` (untracked, the user's playtest list): erase lines once fully resolved; mark
  deferred or changed-needs-confirmation items inline.

## Decisions (playtest 2026-11-11 and after)

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
  Remaining frame time is mostly present/composition, not our passes.

Water (Plan 4, user answers)
- Sea, rivers and large lakes (> ~1 km^2) are infinite; rivers can be diverted; small lakes/ponds
  are finite once disturbed. The sea can flood a channel dug inland below sea level.
- Multi-level (cave) water must work. Buckets AND pipes/pumps. All of F1-F6 before the next playtest.
- Decided (2026-11-11): F2 is a sparse 3D cell simulation, Terraria-style, 1 m voxels (fill amount
  per cell; down, then sideways; a slow levelling pass for U-bends / cave pressure), surface meshed
  with marching cubes on the fill amounts. Replaces the per-column span design.

## Remaining playtest batches (branch playtest-fixes-4)

1. Controls: camera jitter on steep slopes - done (eye smoothing, d0ecdf3), confirm in play.
2. UI: right-click container/machine inventories like the player's; crafting as an icon grid with
   hover tooltips (name, cost, description) - done (30111b8; --shot-tabs checks a chest slot take).
3. World: 60/20/20 tree mix per area + more wood types, recipes take any planks (log -> planks of its
   wood), placeable log/stone blocks, rocks/boulders/cliffs/overhangs, less far-LOD texture tiling -
   done (c75af72): five woods (WoodSpeciesInfo), any_planks/any_log wildcards (ItemMatches), tors and
   boulders in PlanetField::rockFeature, two-scale terrain texturing; --shot-forest.
4. Flowing water F1-F6 (3D voxel sim) - done: src/engine/fluid.*, src/space/waterworks.*; tests --test-flow,
   --test-pond, --test-bucket, --test-pump. Static water feeds only dug cells; ponds are taken over
   (PlanetField::claimPond) and finite; static chunk water shows only over the original bed near by.
5. Parallelism and GPU compute (user, 2026-11-11): audit main-thread work for the job system;
   raylib built with GRAPHICS_API_OPENGL_43 exposes compute shaders + SSBOs (rlgl) - candidates
   are the water sim, chunk density/meshing, foliage placement; keep a GL 3.3 fallback.
Deferred by the user: moon mining launch, motion sickness.
