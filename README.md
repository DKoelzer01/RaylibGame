# RaylibGame
First attempt at Raylib: a voxel planetoid engine (marching cubes, Simplex fBm, shadowed lighting).

## Build (Windows, w64devkit)
Requires raylib 5.5 built at `C:/raylib/raylib` (override with `RAYLIB_PATH=...`).

```sh
export PATH=/c/raylib/w64devkit/bin:$PATH   # or run from the w64devkit shell
mingw32-make -j8
./game.exe
```

Optional Tracy profiler (tracy v0.11.x checkout; delete `*.o` when switching it on/off):
```sh
mingw32-make -j8 TRACY_PATH=C:/path/to/tracy
```

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
