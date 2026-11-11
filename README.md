# RaylibGame
A space open-world survival prototype on a from-scratch raylib voxel planet engine (marching
cubes, Simplex fBm, octree LOD, shadowed lighting). You start on a 60 km home planet with a moon and an
asteroid ring overhead.

## Layout
- `src/engine/` - the reusable engine library (`engine` target): density fields, meshing/LOD, job
  system, player controller, inventory, audio, scenes, debug UI.
- `src/space/` - the space survival game (`game` target): ship, survival, save/load, HUD, `main.cpp`.
- `resources/` - shaders and data (`resources/data/bodies.txt`, `items.txt`, `recipes.txt`).

## Build (CMake)
Needs CMake 3.24+ and a GCC/MinGW toolchain on PATH (e.g. `C:/raylib/w64devkit/bin`).
raylib 5.5 is fetched and built automatically on first configure.

```sh
cmake --preset release          # or: debug, profile (Tracy)
cmake --build --preset release
./game.exe                      # the exe is written to the repo root (assets are loaded relative to it)
```

Options: `-DUSE_SYSTEM_RAYLIB=ON -DRAYLIB_PATH=C:/raylib/raylib` links a prebuilt raylib instead of
fetching; `-DENABLE_TRACY=ON` (the `profile` preset) builds the Tracy v0.11.1 client; connect with
the matching `tracy-profiler` GUI.

VS Code: the CMake Tools extension picks up `CMakePresets.json`; `compile_commands.json` is
generated in the build directory for IntelliSense/clangd.

## Controls
| Key | Action |
| --- | --- |
| Enter | start (main menu) |
| Mouse / arrows | look |
| WASD, Space, Ctrl | move (Shift/Alt = faster) |
| LMB / RMB (hold) | dig / build terrain at the crosshair (build costs stone or regolith) |
| Tab | inventory + crafting |
| L | headlamp (on foot: suit power; piloting: ship nose light) |
| F | toggle free-fly camera (WASD along the view, Space/Ctrl along camera up, mouse turns, Q/E roll, Shift/Alt faster) |
| Space | jump; hold for jetpack |
| E | board / exit the ship |
| 1 / 2 | use oxygen canister / power cell |
| Mouse wheel | brush size (zoom while the telescope is up) |
| T | raise / lower the telescope |
| F5 / F9 | save / load (`saves/world.sav`) |
| Esc | pause; R resumes |
| F1 | debug panel (frame graph, streaming stats, sun/shadow/texture tweaks, brush, world regen by seed) |
| F3 | debug overlay (light, shadow box, shadow map, chunk bounds) |

Shaders in `resources/` hot-reload when saved.

## How the terrain works
- Each body has an immutable density field (`src/engine/voxel.*`): asteroids a `PlanetoidField` (seeded fBm density with a radial
  falloff), planets and moons a `PlanetField` (below); each also gives the material at a point.
- Terrain edits (`src/engine/terraindeltas.*`) are spherical density brushes folded straight into a lattice of
  density deltas (1-voxel spacing, 32^3 chunks of 16-bit values, copy-on-write snapshots shared with the workers).
  Density = base field + the trilinear delta, so builds cost the same after a million strokes. Saves keep the
  deltas in region files next to the save (`<save>.regions/`, 16^3 chunks per file, run-length coded), written
  incrementally: only regions changed since the last save. The game autosaves every 5 minutes. Older saves (brush
  lists) are folded into deltas on load.
- Chunks (32^3 cells) are built on a worker pool (`src/engine/jobsystem.*`): sample density, marching cubes,
  normals from the density gradient, material weights in vertex colours. The main thread only
  uploads finished meshes, time-boxed per frame. Edits re-mesh just the touched chunks.
- Each planetoid is an octree of LOD nodes (`UpdateLod` in `src/engine/world.cpp`): a node at level L covers
  32*2^L voxels with 32 cells, so near terrain is fine and far terrain costs a handful of draw calls.
  Nodes overlap their +x/+y/+z neighbours by one cell to hide cracks between levels, and old nodes
  stay on screen until their replacements are built (no holes while moving).
- World positions are doubles (`DVec3`); the scene keeps a floating `renderOrigin` near the camera,
  so rendering stays precise millions of units out.


## Gameplay (space survival prototype)
- **Player** (`src/engine/player.*`): walks on planetoids with gravity toward the nearest body; "up" follows
  the surface so you can walk all the way around. Collision is resolved directly against the density
  field, so tunnels you dig are walkable. Jump (buffered, with coyote time), WASD air control while airborne
  in gravity, jetpack on Space/Ctrl, 6-DOF suit thrusters in open
  space, and you ride along on spinning bodies. **F** toggles the free-fly dev camera.
- **Items & crafting** (`src/engine/inventory.*`, `resources/data/items.txt`, `recipes.txt`): mining yields
  items by voxel material (rock, regolith, ice, ore); building consumes stone/regolith/lunacrete.
  **Tab** opens inventory + crafting. Add items/recipes by editing the text files.
- **Survival** (`src/space/survival.*`): suit oxygen drains outside the ship's life-support field and a planet's breathable air (~5 min a
  tank), suit power pays for recharging the jetpack and the headlamp, health drops when suffocating or from hard landings and
  crashes. Death respawns you at the ship. **1** uses an oxygen canister, **2** a power cell (charges
  the ship instead while piloting). An objectives list guides the loop: dig stone, find iron ore, craft a power cell,
  board the ship, fly to the moon, visit the asteroid ring.
- **Ship** (`src/space/ship.*`): **E** to board/exit. 6-DOF flight (WASD, Space/Ctrl, mouse, Z/C roll,
  Shift boost, X brake) with flight assist, chase camera, terrain collision, crash damage, and an
  energy store. Once at rest on the ground it anchors to that body and rides its spin (position and
  heading) until you thrust off, keeping the surface's velocity. It feels gravity: piloted with energy,
  flight assist hovers against it (costing energy like thrust); unpiloted, it falls, lands and parks.
- **Gravity** (`GravityAt` in `src/engine/world.*`): summed over every body, inverse-square outside a
  body's mean surface and linear to zero at its centre inside. "Up" is against it, and is held where
  gravity is too weak to define it (deep space, the balance point between bodies).
- **Star system** (`src/engine/bodies.*`, `resources/data/bodies.txt`): one body (asteroid, planet, moon)
  or procedural ring per line, with radius, position or orbit, spin or day length, gravity and seed. The
  default is the 60 km home planet, its moon and a thin asteroid ring 250 km up (members on their own
  circular orbits); the old asteroid start is `bodies_asteroids.txt`. Saves embed the system they were
  made in, so editing the file never breaks a save. Errors are reported with line numbers in `log.txt`.
- **Audio** (`src/engine/audio.*`): procedural sound effects, no assets.
- **Planets and moons** (`PlanetField` in `src/engine/voxel.*`): a height map over the sphere (continents and
  ocean floor, ridged mountains, hills, river valleys), worm caves, ore veins that get richer with depth, snow
  at high latitudes and on peaks, and an unmineable mantle 1.5% of the radius below sea level (900 m on the
  60 km planet). Coarse LOD chunks leave out detail finer than their sample spacing; nodes over the
  horizon, above the highest terrain or inside the mantle are never built. A sea fills the planet to sea level
  (meshed per chunk where it meets open air; polar seas are frozen into ice shelves), and the planet has an
  atmosphere: sky colour by sun angle and altitude, dusk glow, stars at night and in space, haze with
  distance, and dense blue fog underwater. It is the default start system: you start on temperate
  grassland (picked from the seed), where the air is breathable below ~15 km (the suit refills its oxygen).
- **Distance sphere** (`src/engine/farsphere.*`): beyond 3 radii (telescope zoom included) planets and moons draw as a
  sphere with colour and normal maps baked on the workers from the height map and biome rules (sea by depth, forest
  canopy), lit by the sun with the night side and a sea glint, under the same haze as the terrain.
- **Foliage** (`src/engine/foliage.*`, `trees.*`): instanced grass tufts on grass terrain (from each finest chunk's mesh,
  to 70 m), and trees and bushes built from prisms, placed on a 64 m cube-sphere tile grid from the seed (a
  forest field, grass or dirt, gentle slopes, solid ground, below the tree line; broadleaf in the temperate band,
  conifers further north and up high), drawn instanced to 1 km (simpler meshes past 200 m). Trunks collide.
  LMB on a tree chops it (wood), on a bush harvests it (fibre, sometimes berries); plants regrow from a sapling
  after 10 minutes and are full size 20 minutes later; removals are saved. Bench runs save to `saves/bench.sav`.
  Beyond the trees' 1 km range, the terrain paints in forest canopy from the same forest field (forests to the
  horizon and from orbit); forests are densest in the tropics and thin out toward the tundra. The start has a
  14 m clearing for the player and the parked ship.
- **Telescope** (`src/space/telescope.*`): you start with one (craft more from glass and iron plate). T raises it, the
  wheel zooms 4x-64x; terrain in view refines as if that much closer (`SetLodZoom`). Bodies more than 20 km away
  draw in a separate far pass with their own depth range, and haze only counts the air the view ray passes
  through (an 8 km layer), so the moon and the ring are visible from the ground.
- Asteroid density is a displaced sphere (mean surface 0.72 R) with a uniform 0.25/voxel gradient: crisp surfaces and
  predictable dig depth.
## Benchmark / test mode
```sh
./game.exe --bench 10                       # 1920x1080 uncapped; waits for the world to load, prints frame stats
./game.exe --bench 5 --seed 42 --size 200 --cam 400 200 -400
./game.exe --bench 2 --seed 1553044913 --cam 60 45 -50 --test-dig 10 --test-saveload
```
Writes `bench.txt` and `bench.png`. `--cam` = camera offset from the first planetoid, `--count` =
planetoid count, `--system R` = belt radius, `--lod F` = LOD factor, `--first x y z` = first body position
(`--size`, `--count`, `--system` and `--first` override the first body and first ring of `bodies.txt`),
`--bodies file` uses another system file,
`--dump-system` prints every body, `--load` loads `saves/world.sav` at start, `--fps N` caps the frame rate,
`--bench-walk` walks the player through the measured window (frame time while terrain streams), `--crack-check`
clears to magenta without a skybox and reports the share of pixels showing through the terrain (look down at it),
`--ui` opens the debug panel, `--overlay` the F3 overlay, `--nocull` disables
frustum culling, `--test-dig N` carves a tunnel at the crosshair, `--test-saveload` saves and reloads,
`--test-walk S` / `--test-stand S` spawn the player and walk/stand (reports grounding + penetration),
`--test-planet-walk M` walks M metres on body 0 and prints PASS/FAIL (milestone 1 gate: 1000),
`--test-hole` (walk into a fresh pit, build a pillar underfoot), `--test-trees` (walk into a trunk), `--test-chop` (fell a tree, save/reload, regrow), `--test-persist N`
(milestone 2 gate: N strokes, flat builds, save and load under a second), `--bench-fly V`
/ `--bench-fly-alt A` (camera flight over body 0, LOD streaming; negative A: underwater), `--no-foliage`,
`--test-mine N` digs N strokes at the feet, crafts everything possible and builds 3 strokes,
`--test-fly S` autopilots the ship to the nearest other asteroid, `--test-suffocate` checks death/respawn,
`--test-moon` (moon on its analytic orbit; standing and jumping on it), `--time-scale k` (bodies spin and orbit
k times faster), `--test-mantle` (planet terrain: height band, mantle refuses digging, ore by depth, caves; on the
default system), `--test-ride` (standing/jumping on a spinning body; `--ride-spin D` sets its rate), `--test-orbit S` (jetpack orbit stability),
`--test-park` (parked ship rides a spinning body), `--test-gravity` (field values, continuity between
bodies, an unpiloted ship dropping and parking, a piloted ship hovering), `--test-precision` (float vs double
placement error standing on body 0; run it with a large body via `--bodies`), `--spin-test` (writes `spin_*.png` frames),
`--sun x y z` sets the sun direction, `--headlamp` / `--inventory` start with the lamp on / inventory
open, `--look-down` tilts the view to the ground (lighting checks).

## Third-party
See `third_party/README.md` (Dear ImGui v1.92.7, rlImGui tag Raylib_5_5).
