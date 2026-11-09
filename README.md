# RaylibGame
First attempt at Raylib: a voxel planetoid engine (marching cubes, Simplex fBm, shadowed lighting).

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
- Each planetoid has an immutable `PlanetoidField` (`src/voxel.*`): seeded fBm density with a radial
  falloff, plus a material function (rock, regolith, ice, ore).
- Terrain edits are spherical density brushes stored per planetoid. The field + edit list fully
  define the terrain, so a save file is just the generation settings and the edits.
- Chunks (32^3 cells) are built on a worker pool (`src/jobsystem.*`): sample density, marching cubes,
  normals from the density gradient, material weights in vertex colours. The main thread only
  uploads finished meshes, time-boxed per frame. Edits re-mesh just the touched chunks.
- Each planetoid is an octree of LOD nodes (`UpdateLod` in `src/world.cpp`): a node at level L covers
  32*2^L voxels with 32 cells, so near terrain is fine and far terrain costs a handful of draw calls.
  Nodes overlap their +x/+y/+z neighbours by one cell to hide cracks between levels, and old nodes
  stay on screen until their replacements are built (no holes while moving).
- World positions are doubles (`DVec3`); the scene keeps a floating `renderOrigin` near the camera,
  so rendering stays precise millions of units out.


## Gameplay (space survival prototype)
- **Player** (`src/player.*`): walks on planetoids with gravity toward the nearest body; "up" follows
  the surface so you can walk all the way around. Collision is resolved directly against the density
  field, so tunnels you dig are walkable. Jump, jetpack in low gravity, 6-DOF suit thrusters in open
  space, and you ride along on spinning bodies. **F** toggles the free-fly dev camera.
- **Items & crafting** (`src/inventory.*`, `resources/data/items.txt`, `recipes.txt`): mining yields
  items by voxel material (rock, regolith, ice, ore); building consumes stone/regolith/lunacrete.
  **Tab** opens inventory + crafting. Add items/recipes by editing the text files.
- **Survival** (`src/survival.*`): suit oxygen drains outside the ship's life-support field (~5 min a
  tank), suit power drains with jetpack use, health drops when suffocating or from hard landings and
  crashes. Death respawns you at the ship. **1** uses an oxygen canister, **2** a power cell (charges
  the ship instead while piloting). An objectives list guides the loop: mine ice, craft oxygen, board
  the ship, fly to another asteroid, refine ore into power cells.
- **Ship** (`src/ship.*`): **E** to board/exit. 6-DOF flight (WASD, Space/Ctrl, mouse, Z/C roll,
  Shift boost, X brake) with flight assist, chase camera, terrain collision, crash damage, and an
  energy store. The start system has 8 asteroids within ~1.5 km.
- **Audio** (`src/audio.*`): procedural sound effects, no assets.
- Asteroid density is a displaced sphere with a uniform 0.25/voxel gradient: crisp surfaces and
  predictable dig depth.
## Benchmark / test mode
```sh
./game.exe --bench 10                       # 1920x1080 uncapped; waits for the world to load, prints frame stats
./game.exe --bench 5 --seed 42 --size 200 --cam 400 200 -400
./game.exe --bench 2 --seed 1553044913 --cam 60 45 -50 --test-dig 10 --test-saveload
```
Writes `bench.txt` and `bench.png`. `--cam` = camera offset from the first planetoid, `--count` =
planetoid count, `--system R` = belt radius, `--lod F` = LOD factor, `--first x y z` = first body position,
`--ui` opens the debug panel, `--overlay` the F3 overlay, `--nocull` disables
frustum culling, `--test-dig N` carves a tunnel at the crosshair, `--test-saveload` saves and reloads,
`--test-walk S` / `--test-stand S` spawn the player and walk/stand (reports grounding + penetration),
`--test-mine N` digs N strokes at the feet, crafts everything possible and builds 3 strokes,
`--test-fly S` autopilots the ship to the nearest other asteroid, `--test-suffocate` checks death/respawn.

## Third-party
See `third_party/README.md` (Dear ImGui v1.92.7, rlImGui tag Raylib_5_5).
