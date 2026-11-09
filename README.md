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
| Esc | pause; R resumes |
| F1 | debug panel (frame graph, sun/shadow/texture tweaks, world regen by seed) |
| F3 | debug overlay (light, shadow box, shadow map, chunk bounds) |

Shaders in `resources/` hot-reload when saved.

## Benchmark mode
```sh
./game.exe --bench 10                      # 1920x1080 uncapped, prints avg/worst frame time
./game.exe --bench 5 --seed 1553044913 --cam 60 45 -50 --overlay --ui
```
Writes `bench.txt` and `bench.png`. `--cam` is the camera offset from the first planetoid,
`--nocull` disables frustum culling.

## Third-party
See `third_party/README.md` (Dear ImGui v1.92.7, rlImGui tag Raylib_5_5).
