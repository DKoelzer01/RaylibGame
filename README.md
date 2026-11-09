# RaylibGame
A space open-world survival prototype on a from-scratch raylib voxel planetoid engine (marching
cubes, Simplex fBm, octree LOD, shadowed lighting).

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
| F | toggle free-fly camera |
| Space | jump; hold for jetpack |
| E | board / exit the ship |
| 1 / 2 | use oxygen canister / power cell |
| Mouse wheel | brush size |
| F5 / F9 | save / load (`saves/world.sav`) |
| Esc | pause; R resumes |
| F1 | debug panel (frame graph, streaming stats, sun/shadow/texture tweaks, brush, world regen by seed) |
| F3 | debug overlay (light, shadow box, shadow map, chunk bounds) |

Shaders in `resources/` hot-reload when saved.

## How the terrain works
- Each planetoid has an immutable `PlanetoidField` (`src/engine/voxel.*`): seeded fBm density with a radial
  falloff, plus a material function (rock, regolith, ice, ore).
- Terrain edits are spherical density brushes stored per planetoid. The field + edit list fully
  define the terrain, so a save file is just the generation settings and the edits.
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
  field, so tunnels you dig are walkable. Jump, jetpack in low gravity, 6-DOF suit thrusters in open
  space, and you ride along on spinning bodies. **F** toggles the free-fly dev camera.
- **Items & crafting** (`src/engine/inventory.*`, `resources/data/items.txt`, `recipes.txt`): mining yields
  items by voxel material (rock, regolith, ice, ore); building consumes stone/regolith/lunacrete.
  **Tab** opens inventory + crafting. Add items/recipes by editing the text files.
- **Survival** (`src/space/survival.*`): suit oxygen drains outside the ship's life-support field (~5 min a
  tank), suit power drains with jetpack use, health drops when suffocating or from hard landings and
  crashes. Death respawns you at the ship. **1** uses an oxygen canister, **2** a power cell (charges
  the ship instead while piloting). An objectives list guides the loop: mine ice, craft oxygen, board
  the ship, fly to another asteroid, refine ore into power cells.
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
  default is the start asteroid plus a ring of 7 more within ~1.5 km. Saves embed the system they were
  made in, so editing the file never breaks a save. Errors are reported with line numbers in `log.txt`.
- **Audio** (`src/engine/audio.*`): procedural sound effects, no assets.
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
`--ui` opens the debug panel, `--overlay` the F3 overlay, `--nocull` disables
frustum culling, `--test-dig N` carves a tunnel at the crosshair, `--test-saveload` saves and reloads,
`--test-walk S` / `--test-stand S` spawn the player and walk/stand (reports grounding + penetration),
`--test-mine N` digs N strokes at the feet, crafts everything possible and builds 3 strokes,
`--test-fly S` autopilots the ship to the nearest other asteroid, `--test-suffocate` checks death/respawn,
`--test-ride` (standing/jumping on a spinning body; `--ride-spin D` sets its rate), `--test-orbit S` (jetpack orbit stability),
`--test-park` (parked ship rides a spinning body), `--test-gravity` (field values, continuity between
bodies, an unpiloted ship dropping and parking, a piloted ship hovering), `--spin-test` (writes `spin_*.png` frames),
`--sun x y z` sets the sun direction, `--headlamp` / `--inventory` start with the lamp on / inventory
open, `--look-down` tilts the view to the ground (lighting checks).

## Third-party
See `third_party/README.md` (Dear ImGui v1.92.7, rlImGui tag Raylib_5_5).
