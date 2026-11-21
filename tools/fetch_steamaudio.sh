#!/bin/bash
# Fetches Valve's Steam Audio SDK (#1020, #1292: HRTF, later occlusion and reverb) into a gitignored folder:
#
#   tools/fetch_steamaudio.sh [OUT]   OUT defaults to <repo>/deps/steamaudio
#
# Result: OUT/lib/windows-x64/phonon.dll and OUT/lib/linux-x64/libphonon.so (the game loads them at run time:
# src/engine/steamaudio.cpp; the headers are vendored in third_party/steamaudio). Re-run cmake afterwards: the build
# then copies phonon.dll next to game.exe. Without it the game falls back to its own panning (no HRTF).
# The version is pinned and the download checked by SHA-256 (~181 MB zip; phonon.dll is 53 MB). Apache-2.0
# (third_party/steamaudio/LICENSE.md, THIRDPARTY.md). Needs: curl, unzip, sha256sum.
set -euo pipefail

SA_VERSION=4.8.1
SA_URL="https://github.com/ValveSoftware/steam-audio/releases/download/v${SA_VERSION}/steamaudio_${SA_VERSION}.zip"
SA_SHA256=4a0aa5ec1176f38f0b0993a37c2259d9e86f27e22d5e24f83ec4c3cb9a1d5449

case "${1:-}" in -h|--help) sed -n '2,10p' "$0"; exit 0 ;; esac
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${1:-$ROOT/deps/steamaudio}"
if [ -f "$OUT/VERSION" ] && [ "$(cat "$OUT/VERSION")" = "$SA_VERSION" ] && [ -f "$OUT/lib/windows-x64/phonon.dll" ]; then
    echo "steam audio $SA_VERSION already in $OUT"
    exit 0
fi
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
echo "downloading steam audio $SA_VERSION ..."
curl -sSfL -o "$TMP/sa.zip" "$SA_URL"
echo "$SA_SHA256  $TMP/sa.zip" | sha256sum -c -
unzip -q "$TMP/sa.zip" 'steamaudio/lib/windows-x64/phonon.dll' 'steamaudio/lib/linux-x64/libphonon.so' -d "$TMP/x"
mkdir -p "$OUT/lib/windows-x64" "$OUT/lib/linux-x64"
cp "$TMP/x/steamaudio/lib/windows-x64/phonon.dll" "$OUT/lib/windows-x64/"
cp "$TMP/x/steamaudio/lib/linux-x64/libphonon.so" "$OUT/lib/linux-x64/"
echo "$SA_VERSION" > "$OUT/VERSION"
echo "steam audio $SA_VERSION in $OUT (re-run cmake to copy phonon.dll next to game.exe)"
