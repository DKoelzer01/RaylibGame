#!/bin/bash
# Builds the Linux dedicated server (#792, #326) and packs a self-contained runtime folder:
#
#   tools/build_linux.sh [--out DIR] [--no-gns] [--jobs N]            on a Linux host with Docker: builds in a
#                                                                     Debian 12 container (tools/linux/Dockerfile.build)
#   tools/build_linux.sh --inside [--out DIR] [--no-gns] [--jobs N]   in that container, or on any Debian 12 / Ubuntu
#                                                                     22.04+ host with the packages the Dockerfile lists
#
# Result: DIR/raylibgame-server-linux-x86_64.tar.gz (DIR defaults to <repo>/build/linux), which unpacks to:
#   game                    the server (raylib's headless backend: no X11/Wayland; OpenGL 4.3 from EGL when a GPU is
#                           in the container, else the null GL and the CPU paths), libstdc++/libgcc linked in
#   resources/ assets/      the game data
#   lib/                    libglvnd's libEGL.so.1 + libGLdispatch.so.0 and egl_vendor.d/ (NVIDIA, Mesa), for images
#                           that have no libEGL: the driver's own libEGL_nvidia.so.0 comes from the NVIDIA runtime
#   server.cfg.example      run ./game --server server.cfg (a missing server.cfg is written with the defaults)
#
# Built in Debian 12 (glibc 2.36), so it runs on Debian 12+ and Ubuntu 22.04+ (and the Pterodactyl yolks:debian image).
# GameNetworkingSockets (tools/build_gns.sh, with OpenSSL's static libcrypto) is cached in DIR/deps; without it
# (--no-gns) the server has no network transport and is only good for --test-* runs.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$REPO/build/linux"
INSIDE=0; GNS=1; JOBS="${JOBS:-4}"
while [ $# -gt 0 ]; do
    case "$1" in
        --inside) INSIDE=1 ;;
        --out) OUT="$2"; shift ;;
        --no-gns) GNS=0 ;;
        --jobs) JOBS="$2"; shift ;;
        -h|--help) sed -n '2,24p' "$0"; exit 0 ;;
        *) echo "build_linux: unknown argument $1"; exit 2 ;;
    esac
    shift
done
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"

if [ "$INSIDE" = 0 ]; then
    command -v docker >/dev/null || { echo "build_linux: no docker: run with --inside on a Debian 12 host with the packages in tools/linux/Dockerfile.build"; exit 1; }
    IMAGE=raylibgame-build:bookworm
    docker build -q -t "$IMAGE" -f "$REPO/tools/linux/Dockerfile.build" "$REPO/tools/linux" >/dev/null
    VERSION="$(git -C "$REPO" rev-parse --abbrev-ref HEAD 2>/dev/null || echo linux) $(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo unknown)"
    [ -f "$REPO/.raylibgame-version" ] && VERSION="$(cat "$REPO/.raylibgame-version")"
    ARGS=(--inside --out /out --jobs "$JOBS"); [ "$GNS" = 0 ] && ARGS+=(--no-gns)
    # As the calling user (files in OUT stay theirs), the source read-only.
    exec docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -e RAYLIBGAME_VERSION="$VERSION" \
        -v "$REPO":/src:ro -v "$OUT":/out raylibgame-build:bookworm bash /src/tools/build_linux.sh "${ARGS[@]}"
fi

# ---- inside: configure, build, stage, pack ----
for t in cmake gcc g++ make git curl tar; do command -v "$t" >/dev/null || { echo "build_linux: '$t' missing"; exit 1; }; done
export RAYLIBGAME_VERSION="${RAYLIBGAME_VERSION:-linux}"
BUILD="$OUT/build"; STAGE="$OUT/stage/raylibgame-server"
CMAKE_ARGS=(-DCMAKE_BUILD_TYPE=Release -DRAYLIBGAME_HEADLESS=ON -DGAME_OUTPUT_DIR="$OUT/bin")
if [ "$GNS" = 1 ]; then
    JOBS="$JOBS" nice -n 10 bash "$REPO/tools/build_gns.sh" "$OUT/deps/gns"
    CMAKE_ARGS+=(-DENABLE_GNS=ON -DGNS_ROOT="$OUT/deps/gns/install")
fi
# FetchContent keeps its raylib checkout in the build folder: a re-run does not download it again.
cmake -S "$REPO" -B "$BUILD" -G "Unix Makefiles" "${CMAKE_ARGS[@]}"
nice -n 10 cmake --build "$BUILD" -j "$JOBS" --target game

rm -rf "$OUT/stage"; mkdir -p "$STAGE/lib/egl_vendor.d"
cp "$OUT/bin/game" "$STAGE/game"
strip "$STAGE/game"
cp -r "$REPO/resources" "$STAGE/resources"
cp -r "$REPO/assets" "$STAGE/assets"
# libglvnd's dispatcher (tiny, vendor-neutral): only used when the image has no libEGL.so.1 of its own.
for l in libEGL.so.1 libGLdispatch.so.0; do
    src="$(ls /usr/lib/x86_64-linux-gnu/$l 2>/dev/null || true)"
    [ -n "$src" ] && cp -L "$src" "$STAGE/lib/$l" || echo "build_linux: $l not found (install libegl1): the runtime will need the image's own"
done
printf '{ "file_format_version": "1.0.0", "ICD": { "library_path": "libEGL_nvidia.so.0" } }\n' > "$STAGE/lib/egl_vendor.d/10_nvidia.json"
printf '{ "file_format_version": "1.0.0", "ICD": { "library_path": "libEGL_mesa.so.0" } }\n' > "$STAGE/lib/egl_vendor.d/50_mesa.json"
cat > "$STAGE/server.cfg.example" <<'EOF'
# RaylibGame dedicated server (#326). key = value; '#' starts a comment. Restart the server to apply changes.
save = server/world.sav
port = 27015
max_players = 8
password =
autosave_minutes = 30
motd =
seed = 0
log = server/server.log
EOF
echo "$RAYLIBGAME_VERSION" > "$STAGE/VERSION"
tar -C "$OUT/stage" -czf "$OUT/raylibgame-server-linux-x86_64.tar.gz" raylibgame-server
echo "build_linux: $(du -h "$OUT/raylibgame-server-linux-x86_64.tar.gz" | cut -f1) -> $OUT/raylibgame-server-linux-x86_64.tar.gz ($RAYLIBGAME_VERSION)"
