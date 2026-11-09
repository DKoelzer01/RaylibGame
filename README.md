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
| LMB / RMB (hold) | dig / add terrain at the crosshair |
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

## Benchmark / test mode
```sh
./game.exe --bench 10                       # 1920x1080 uncapped; waits for the world to load, prints frame stats
./game.exe --bench 5 --seed 42 --size 200 --cam 400 200 -400
./game.exe --bench 2 --seed 1553044913 --cam 60 45 -50 --test-dig 10 --test-saveload
```
Writes `bench.txt` and `bench.png`. `--cam` = camera offset from the first planetoid, `--count` =
planetoid count, `--ui` opens the debug panel, `--overlay` the F3 overlay, `--nocull` disables
frustum culling, `--test-dig N` carves a tunnel at the crosshair, `--test-saveload` saves and reloads.

## Third-party
See `third_party/README.md` (Dear ImGui v1.92.7, rlImGui tag Raylib_5_5).
