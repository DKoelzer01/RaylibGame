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
- Playtest items: each CurrentIssues.md line has a matching issue. Fixed-but-unconfirmed items carry
  `needs playtest`; deferred ones carry `deferred`. Keep the issue and the CurrentIssues.md line in sync.
- Bugs found along the way (a failing test, a regression) get their own issue, even when fixed in
  the same branch.
- Exempt: answering a question, read-only investigation, and triage itself (which comments on the
  issue it triages).

## Plans and trackers

- Plan 2 (early game, done): https://claude.ai/artifact/RL22kZ7VBMugAFgGRPS3cE
- Plan 3 (Terra tiers, space): https://claude.ai/artifact/XVw45AC6vwYitbdhaeP2xL - M21-M24 done (M24 on PR #28).
- Plan 4 (flowing water): https://claude.ai/artifact/27xL5j1zkirSe8pDZgPz59
- Plan 6 (GPU terrain generation + GPU-driven drawing, approved 2026-11-12: whole plan, GL 4.6 gl_DrawID with
  the 4.3 fallback, CPU terrain as a tested fallback): https://claude.ai/artifact/VJdKtUZEpBg5uefnqewRyZ - issues #143
  (steps #144-#150), under #128 Graphics Pipeline Overhaul.
- Plan 5 (art, life, world variety): https://claude.ai/artifact/5igt2VrsTyNhZ3BZ94teHX - A animal models/AI, B biomes,
  C Minecraft-style 3D terrain, D authored textures await the user's answers; the light items are done on art-pass-5.
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
  Scene::KickGpu (glFlush) at pass boundaries: without it CPU and GPU ran serialized (+25% FPS). GPU timer = elapsed queries.
  Reverse-Z (#135, scenetarget.*): the 3D world draws into an offscreen float-depth target (window sample count, glClipControl
  0..1, infinite far, GL_GREATER) in ONE pass (no far/mid/near depth ranges or clears); RAYLIBGAME_REVZ=0 = old passes.
  Terrain arena (#144, terrainarena.*): one multi-draw-indirect per body per pass; RAYLIBGAME_ARENA=0 = per-chunk draws.
  (Old note "remaining time is present/composition" was wrong: it was CPU/GPU serialization.)
- Terrain lattice 0.25 m (kVoxelSize; RAYLIBGAME_VOXEL_1M builds the old 1 m one): digs carve exactly the
  previewed sphere (--test-dig-shape IoU 0.85). Delta chunks hold sparse 4^3 bricks; regions are RGN5.
- Shadows (2026-11-12): two cascades (near +-32 m, far +-200 m, 2048 each), bilinear PCF, normal-offset
  lookup + slope bias; trees and bushes draw into both (tree vertex shader + depth.fs). Grass casts none
  by choice (blades thinner than a near texel).
- Foliage (user, 2026-11-12): grass = 3 instanced meshes (tuft, 14-blade clump, 24-blade dense clump;
  variant rides in scale + 10*variant), 7 per m^2, darker lower half, drawn to 100 m; trees to 2 km,
  bushes 300 m. The user trades a few FPS for view distance (walk bench ~220 FPS). Draw distances are
  runtime values (--grass-dist, --tree-dist); the terrain canopy tint follows treeDistance.

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
   Done (da1cde5): shader hot-reload checks throttled to 2/s (frame 2.93 -> 2.55 ms), water brick
   classification on the job pool (worst step 26.8 -> 11.9 ms). GPU compute not adopted: it needs a GL 4.3
   raylib build and a shader copy of the density field that the CPU (collisions, queries) would still
   need too; revisit only if a CPU system becomes the bottleneck.
5b. M24 moon and ring tiers - done (15fae76 + the ladder gate: --test-progression ev, PASS).
6. 3D terrain (user, 2026-11-12): full 3D noise for hills and cliffs with natural overhangs and
   alcoves, on top of the height map and the tors/boulders (PlanetField::rockFeature). Plan: a 3D
   displacement band around the surface (~-15..+25 m), strength from a regional "rugged" mask and
   the slope (gentle plains stay smooth), domain-warped ridged noise for cliffs and undercuts;
   material = rock where the 3D term carves or builds; removes floating bits (no solid fully above
   open air beyond a few m). Systems that assume a height map must use density near the surface:
   tree/stone placement (already checks density), static water bed, spawn, map, materialSky depth.
   Coarse LOD keeps the low octaves so far silhouettes match. Gate: --shot-terrain (cliff and
   overhang screenshots), --crack-check, --test-planet-walk, chunk build time within +30%.
   Built: PlanetField::terrain3D (smooth warped noise stretched upward + a little ridged, raised more
   than carved, dips capped at 3 m; mask from rugged noise and elevation); materialSky depth follows it;
   --shot-terrain. Tuning history: all-ridged made paper-thin shelves, deep carving made pits.
7. Chunk generation speed (2026-11-12): builds already run in parallel on the job pool (one job per
   chunk, ~7 ms each, 3-5 ms of it density sampling). Next: SIMD noise (4-8 samples at once, the
   biggest win), a per-region height-map cache shared across chunks and LOD levels, cheap
   all-rock/all-air rejection before full sampling; split a single chunk across threads only for
   urgent rebuilds. Do this with or right after 6 (3D noise adds sampling cost).
   Built: SimplexNoise::noise4 (SSE2, bit-identical to noise(), checked by --profile-density) used by
   fractal/fractalLimited 4 octaves at a time (density sample ~1.7x faster); river pieces carry a
   bounding sphere for nearest(). --profile-density prints ns per call by part.
Deferred by the user: moon mining launch, motion sickness.
Follow-ups done (2026-11-12, after the issue list): minimap mask at draw time (f4afeb0), tor strata kept to
the walls (5e19f49; --shot-forest --outcrop-high), planet spin as a wrapped double (d187759), trees cast
shadows (d591eee), sparse 4^3 bricks in delta chunks (66f2dbf; --test-persist prints delta memory), fuller
grass + longer view distances (8201a55; --shot-shadow adds bench_grass.png, --shot-far-forest). Full planet
regression 47/47 on 040f1eb. All CurrentIssues items are "changed - confirm in play".

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
