#!/bin/bash
# Builds Valve's GameNetworkingSockets (the multiplayer transport, #310) and its one dependency, protobuf, from
# source with the game's own compiler (w64devkit MinGW GCC), statically, into a gitignored folder:
#
#   tools/build_gns.sh [OUT]        OUT defaults to <repo>/deps/gns
#
# Result: OUT/install/{include/steam/*.h, lib/libGameNetworkingSockets_s.a, lib/libprotobuf.a}, which is the
# default GNS_ROOT of `cmake -DENABLE_GNS=ON`. Takes ~7 min the first time (protobuf ~5, GNS ~2); a re-run with
# the install present and the same versions does nothing (delete OUT, or pass --force, to rebuild).
#
# Versions are pinned here (and the downloads checked by SHA-256). No binaries are checked in.
# Needs: cmake, w64devkit (gcc/g++/mingw32-make) on PATH, curl, tar.
set -euo pipefail

PROTOBUF_VERSION=3.21.12   # the last protobuf without the abseil dependency
PROTOBUF_URL="https://github.com/protocolbuffers/protobuf/releases/download/v21.12/protobuf-cpp-${PROTOBUF_VERSION}.tar.gz"
PROTOBUF_SHA256=4eab9b524aa5913c6fffb20b2a8abf5ef7f95a80bc0701f3a6dbb4c607f73460
GNS_VERSION=1.4.1
GNS_URL="https://github.com/ValveSoftware/GameNetworkingSockets/archive/refs/tags/v${GNS_VERSION}.tar.gz"
GNS_SHA256=1cfb2bf79c51a08ae4e8b7ff5e9c1266b43cfff6f53ecd3e7bc5e3fcb2a22503
STAMP="protobuf-${PROTOBUF_VERSION} gns-${GNS_VERSION}"

FORCE=0
OUT=""
for a in "$@"; do
    case "$a" in
        --force) FORCE=1 ;;
        -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
        *) OUT="$a" ;;
    esac
done
REPO="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${OUT:-$REPO/deps/gns}"
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
INSTALL="$OUT/install"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"

if [ "$FORCE" = 0 ] && [ -f "$INSTALL/stamp.txt" ] && [ "$(cat "$INSTALL/stamp.txt")" = "$STAMP" ] \
   && [ -f "$INSTALL/lib/libGameNetworkingSockets_s.a" ] && [ -f "$INSTALL/lib/libprotobuf.a" ]; then
    echo "build_gns: $STAMP already built in $INSTALL"
    exit 0
fi

for t in cmake gcc g++ mingw32-make curl tar; do
    command -v "$t" >/dev/null || { echo "build_gns: '$t' not on PATH (export PATH=\"/c/Program Files/CMake/bin:/c/raylib/w64devkit/bin:\$PATH\")"; exit 1; }
done

sha_of() { sha256sum "$1" | cut -d' ' -f1; }
fetch() { # url sha256 file
    local url=$1 sum=$2 file=$3
    if [ ! -f "$file" ] || [ "$(sha_of "$file")" != "$sum" ]; then
        echo "build_gns: downloading $url"
        curl -fL --retry 3 -o "$file.part" "$url"
        mv "$file.part" "$file"
    fi
    [ "$(sha_of "$file")" = "$sum" ] || { echo "build_gns: SHA-256 mismatch for $file"; exit 1; }
}

# The first MinGW configure sometimes fails reading the compiler probe's a.exe while antivirus scans it:
# configure once more before giving up.
configure() {
    cmake "$@" || { echo "build_gns: configure failed, retrying once"; cmake "$@"; }
}

DL="$OUT/dl"; SRC="$OUT/src"; BLD="$OUT/build"
mkdir -p "$DL" "$SRC" "$BLD"
fetch "$PROTOBUF_URL" "$PROTOBUF_SHA256" "$DL/protobuf-cpp-${PROTOBUF_VERSION}.tar.gz"
fetch "$GNS_URL" "$GNS_SHA256" "$DL/gns-${GNS_VERSION}.tar.gz"

rm -rf "$SRC"/* "$BLD"/* "$INSTALL"
tar -xzf "$DL/protobuf-cpp-${PROTOBUF_VERSION}.tar.gz" -C "$SRC"
tar --exclude='*/.github' -xzf "$DL/gns-${GNS_VERSION}.tar.gz" -C "$SRC"  # .github holds a symlink busybox tar can't make
PB_SRC="$SRC/protobuf-${PROTOBUF_VERSION}"
GNS_SRC="$SRC/GameNetworkingSockets-${GNS_VERSION}"
PB_INSTALL="$BLD/pb-install"

echo "build_gns: protobuf $PROTOBUF_VERSION ($JOBS jobs)"
configure -S "$PB_SRC" -B "$BLD/pb" -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX="$PB_INSTALL" \
    -Dprotobuf_BUILD_TESTS=OFF -Dprotobuf_BUILD_SHARED_LIBS=OFF -Dprotobuf_WITH_ZLIB=OFF
cmake --build "$BLD/pb" -j "$JOBS"
cmake --install "$BLD/pb"

echo "build_gns: GameNetworkingSockets $GNS_VERSION"
# BCrypt = Windows' own crypto: no OpenSSL or libsodium.
configure -S "$GNS_SRC" -B "$BLD/gns" -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_PREFIX_PATH="$PB_INSTALL" \
    -DUSE_CRYPTO=BCrypt -DBUILD_STATIC_LIB=ON -DBUILD_SHARED_LIB=OFF -DProtobuf_USE_STATIC_LIBS=ON \
    -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DBUILD_TOOLS=OFF
cmake --build "$BLD/gns" -j "$JOBS" --target GameNetworkingSockets_s

mkdir -p "$INSTALL/include/steam" "$INSTALL/lib"
cp "$GNS_SRC"/include/steam/*.h "$INSTALL/include/steam/"
cp "$(find "$BLD/gns" -name libGameNetworkingSockets_s.a | head -1)" "$INSTALL/lib/"
cp "$PB_INSTALL/lib/libprotobuf.a" "$INSTALL/lib/"
echo "$STAMP" > "$INSTALL/stamp.txt"
echo "build_gns: done -> $INSTALL"
echo "build_gns: configure the game with -DENABLE_GNS=ON (GNS_ROOT defaults to <repo>/deps/gns/install)"
