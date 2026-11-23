# Architecture notes

Subsystem notes moved out of CLAUDE.md (#1478). Read the entries for the code you touch.

- Erosion (Plan 7): the bake runs on a thread behind a progress screen (the `gWorldGenWait` hook). Cache
  `<save>.erosion` (ERO6) keyed on `ErosionParams` + `kErosionBakeVersion`: bump the version on ANY bake code change,
  or stale bakes load. The 3D band's `Column3D` (mask, f) must match `terrain.glsl` / `terrain_density.comp`
  bit-exactly (the grid is 2 vec4 per node; the second is unused since the cliff terraces and lips went, #543);
  `--test-terrain-gpu` checks it. Plan 9 (#586): steep terrain comes from formations in the bake input, not post-bake steps.
- Sky light: `fragSky` from depth below the surface, times a short-ray openness probe (#369: terrain via chunk
  samples, parts via plot cells); `RAYLIBGAME_SKY_PROBES=0` for A/B. Moonlight scale defaults to 0.06 (#1166; 0.09 #1117, 0.12 #897, was 0.25; Settings > Brightness offsets the exposure, -1..+1 EV)
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
- Inventory: `Inventory::kGrid` (50) grid cells drawn 10 wide, then 10 hotbar cells (50-59), then the off hand (60)
  (#495; `kOldGrid` = 45 for old saves). Chests have as many cells as the grid. `Hotbar_Bind` moves a stack (items must
  be carried); `Inventory::preferred` = the selected hotbar cell is used first.
- Streaming budget (#1202): chunk work per frame is capped by a time budget (`RAYLIBGAME_STREAM_BUDGET=0` = off,
  `_BUSY=0` = the old frame-time measure). Anything the player must see at once bypasses it: closing the telescope
  hit 300 ms frames until zoom changes went urgent (#1424/#1441).
- Shader compiles (#1386): compute shaders compile asynchronously where the driver supports it
  (`RAYLIBGAME_ASYNC_SHADERS=0` = off; off by default on the headless Linux platform, where NVIDIA's driver does the
  work on the first status read anyway). `RAYLIBGAME_COLD_SHADERS=1` forces a cold compile to measure it: a cold cache
  is ~2.9 s for fluid_step.comp + 2.5 s for terrain_density.comp here, ~37 s on the runner (#1047).
- Fluids: the U-bend pressure pass (#1168/#1395) is CPU-only and opt-in (`RAYLIBGAME_FLUID_PRESSURE=1`) until
  fluid_step.comp has it - CPU and GPU kernels must agree (`RAYLIBGAME_FLUID_TRACE=1` = the parity trace, #968). Rivers
  are infinite carved water with no dry reaches: a dry river bed is a bug (#1438). Shore stripes while a channel fills
  (#1091) are a moving flow pattern held by two deliberate rules (giver-not-resting, minSpread), parked for a spreading
  design.
- Fly camera (#1436): on the frame its follow mode flips (spinning frame <-> centred at 2 radii) or its reference body
  changes (SOI), it is carried with the body and mode it was recorded in; skipping that carry left it a frame of orbital
  motion behind (3 m hops at the moon).
- Factory (#997): machines sleep when idle (`RAYLIBGAME_FACTORY_SLEEP=0` = off), distant factories run as a summary
  sim (`RAYLIBGAME_FACTORY_SUMMARY=0` = off), worker threads <= 4 and <= half the cores (`RAYLIBGAME_FACTORY_WORKERS`).
- Weather captures: `RAYLIBGAME_WEATHER_FORCE=0..3` (clear .. storm).
- Caves (Plan 10, #615): all cave carving is `PlanetField::caveDistance` (voxel.cpp) and its bit-exact twin in terrain.glsl; change both
  together. `--cave-stats R [D]` (counts, the trench metric, a density hash that refactors must keep), `--test-caves` (fly into the nearest
  mouth). Graph tunnels (step 2 #995, on by default since step 3 #1022; `RAYLIBGAME_CAVES=0` = the old noise caves): per-256 m-cell passages,
  per-block segment lists (6 vec4 each, to the GPU at binding 4; <= 32 a block, else the CPU path); `RAYLIBGAME_CAVE_ABUNDANCE` 0-2 (the New
  World slider to come). Shapes (#1022): the segment's bore is the skeleton; heavy 3D noise only grows the wall inward (none inside half the
  bore), looked up squeezed radially (x0.1) so there are no voids behind the faces; ~9 m rock forms with ~3.2 m hard-edged creases
  (`RAYLIBGAME_CAVE_WALLS=0` = the first merged ~3 m / ~1 m walls). Mouths only in real faces (5 m+ rise within 10 m, no tors in
  front). Karst (#1052): dolines in karst country are sinks (shaft, river passage down the drainage, keyholes, mazes, resurgences);
  `RAYLIBGAME_CAVES_KARST=1 --test-caves` finds one, `RAYLIBGAME_CAVES_INSIDE=1` starts in a passage. Chambers (#1058): at deep passage
  vertices (segment kind 3: a level floor), with solid boulders and columns; graph mode has no noise caverns. Lava tubes (#1063): on
  volcano flanks, kind-3 tubes 5-25 m down with ledges, skylights and an exit; `RAYLIBGAME_CAVES_TUBE=1 --test-caves` finds one. Noise factors include each noise's own base frequency (1/70-1/700: forgetting it makes every term constant). Captures:
  `RAYLIBGAME_CAMERA_LIGHT=1` (a light at the eye), `RAYLIBGAME_CAVES_VIEW=front` (eye-level frontal approach), `RAYLIBGAME_CAVES_MOUTH="x y z"`
  (that mouth, at any abundance). In GLSL, store every a*b+c handed to a function (xdiv, xsqrt, fractal) in a `precise` variable first.
- Terrain internals still worth knowing: `PlanetField::rockFeature` (tors, boulders) and `terrain3D` (the 3D band);
  SSE2 `SimplexNoise::noise4` is bit-identical to `noise()` (`--profile-density` checks it); static water feeds only
  dug cells, and ponds are taken over by the fluid sim (`PlanetField::claimPond`) and become finite.


## Saves and migrations (#1693)

- **One version**: `kSaveVersion` in `src/engine/saveversion.h`, written as the world file's first line (`RAYLIBGAME_SAVE 13`); the
  player record has its own `kPlayerRecordVersion` (`RAYLIBGAME_PLAYER 1`). `kOldestSaveVersion` (13) is the oldest world file
  that opens; older ones are refused ("This save is from an older build (version N); saves from version 13 on load."), as is a
  newer one ("made by a newer build ... reads up to 13"). Both refuse before anything is touched; `LastLoadError()` has the
  message (the menu's load toast shows it).
- **The readers know only the current format.** Section tags (`plots5`, `planetoid`, `extras`...) stay as checks, not versions.
  `Save_LoadingVersion()` is the format being parsed (after migration, the current one).
- **Loading an older save** (`LoadWorldSave`): the existing backup (`<save>.backups/`) is taken first, the file's lines go through
  `Save_Migrate` (`src/engine/savemigrate.cpp`, the table `kWorldSteps`) into `<save>.migrating`, which the one reader parses and
  then deletes. The file on disk is untouched; the upgraded save is written by the next normal save. Data a step makes meaningless
  is removed with `ctx.drop("what")`: logged, and one toast after the load ("Save upgraded from version N; dropped: ...").
  A multiplayer client (`gNetRole == Client`) never migrates and never autosaves (the autosave in main.cpp is gated on it): the host's
  save is the only one. Player records: `PlayerRecord_Migrate`, same shape (`ReadRecordInto` runs it for a record below the
  current number; a record from a newer build is not read and the next save keeps it as a `.bak`).
- **To bump the version**: change the format; add a step `{ N, "what it does", Step_N_to_N1 }` to `kWorldSteps` (the commented example
  in savemigrate.cpp), a function over the text lines that turns format N into N+1 (rewrite lines, `ctx.drop` what no longer means
  anything); set `kSaveVersion` to N+1; reserve nothing else. Binary region files keep their own magic (RGN5; older ones skipped
  with a warning): a region change adds its step beside the world file's. Then make a golden save with the new build:
  `$GAME --seed 42 --scenario saves/golden-write --scenario-out <dir>` and commit `<dir>/saves/golden-write/golden-vN/` as
  `tests/saves/vN/`. `saves/golden-load` (scenario, runner or local) loads every `tests/saves/v*/` through the chain and checks the
  facts in each `facts.txt`; keep the old directories, they are the test that old saves still open.
