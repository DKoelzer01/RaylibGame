Vendored third-party sources (unmodified).

- imgui: https://github.com/ocornut/imgui tag v1.92.7 (MIT; matches the rlImGui tag below)
- rlImGui: https://github.com/raylib-extras/rlImGui tag Raylib_5_5 (zlib)
- steamaudio: https://github.com/ValveSoftware/steam-audio release v4.8.1, headers only (include/phonon.h,
  phonon_version.h; Apache-2.0: LICENSE.md, THIRDPARTY.md). The library itself is loaded at run time;
  tools/fetch_steamaudio.sh downloads it into deps/steamaudio (#1292).

Tracy (optional profiler) is not vendored: CMake fetches v0.11.1 when ENABLE_TRACY=ON.
