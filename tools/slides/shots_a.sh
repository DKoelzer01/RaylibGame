#!/bin/bash
# Slideshow shots, batch A (#1071): the shots the existing --shot-* flags and bench options can take. Each line runs the
# game once (HUD off, the native window size) and keeps its screenshot under OUT as <name>.png; at the end a labelled
# contact sheet (sheet.png) is made for picking. Make slides from the chosen ones with
#   game.exe --to-slide OUT/<name>.png resources/menu/slides/<name>.jpg
# Usage: tools/slides/shots_a.sh [OUT] [EXE]   (default OUT ./slides-a, EXE build/check/bin/game.exe); run from the repo root.
# Needs the PC's GPU and a quiet machine (no other game running): screenshots of a half-loaded world are useless.
set -u
OUT="${1:-slides-a}"
EXE="${2:-build/check/bin/game.exe}"
mkdir -p "$OUT"
COMMON="--no-hud --native-window"
# name | environment | arguments | the screenshot it writes
SHOTS=(
  "valley||--seed 42 --bench 12 --day|bench.png"
  "mountain||--seed 42 --bench 20 --shot-mountain --day|bench_mountain.png"
  "biomes||--seed 42 --bench 90 --shot-biomes --day|bench_biome_*.png"
  "sunset||--seed 42 --bench 12 --sun-elev 2|bench.png"
  "sunrise||--seed 42 --bench 12 --sun-elev -1|bench.png"
  "night||--seed 42 --bench 12 --sun-elev -20|bench.png"
  "cave-mouth||--seed 42 --bench 40 --shot-cave-mouth --day|bench_mouth_*.png"
  # --test-caves takes no shot of its own: the bench's last frame (bench.png) shows where the fly-in ends.
  "lava-tube|RAYLIBGAME_CAVES_TUBE=1 RAYLIBGAME_CAMERA_LIGHT=1|--seed 42 --bench 40 --test-caves --day|bench.png"
  "karst|RAYLIBGAME_CAVES_KARST=1 RAYLIBGAME_CAMERA_LIGHT=1|--seed 42 --bench 40 --test-caves --day|bench.png"
  "room||--seed 42 --bench 20 --shot-room|bench_room.png"
  "machines||--seed 42 --bench 20 --shot-machines --day|bench_machines.png"
  "campfire||--seed 42 --bench 20 --shot-campfire --sun-elev 1|bench_campfire.png"
  "animals||--seed 42 --bench 30 --shot-animals --day|bench_animals_models.png"
)
for line in "${SHOTS[@]}"; do
  IFS='|' read -r name envs args pattern <<< "$line"
  rm -f $pattern
  t0=$SECONDS
  env $envs timeout 240 "$EXE" $COMMON $args > "$OUT/$name.log" 2>&1
  rc=$?
  n=0
  found=( $(ls $pattern 2>/dev/null) ) # counted once, before any moves: one shot keeps the plain name, several get suffixes
  for f in "${found[@]}"; do
    base="${f%.png}"; suffix="${base#bench}"; suffix="${suffix#_}"
    if [ ${#found[@]} -gt 1 ]; then mv "$f" "$OUT/$name-$suffix.png"; else mv "$f" "$OUT/$name.png"; fi
    n=$((n + 1))
  done
  echo "shots: $name: exit $rc, $n shot(s), $((SECONDS - t0)) s"
done
# The contact sheet: every shot, 4 across, 480 px wide each, its file name under it.
powershell -NoProfile -ExecutionPolicy Bypass -File "$(dirname "$0")/sheet.ps1" -Dir "$OUT" -Out "$OUT/sheet.png"
