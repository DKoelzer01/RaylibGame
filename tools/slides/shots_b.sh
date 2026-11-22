#!/bin/bash
# Slideshow shots, batch B (#1071): more --shot-* views (lit, dusk, night, karst and lava-tube mouths), then the scenarios that frame a subject well and take a picture with c.grab when RAYLIBGAME_SLIDE_SHOTS is set. Each line runs the
# game once (HUD off, the native window size) and keeps its screenshot under OUT as <name>.png; at the end a labelled
# contact sheet (sheet.png) is made for picking. Make slides from the chosen ones with
#   game.exe --to-slide OUT/<name>.png resources/menu/slides/<name>.jpg
# Usage: tools/slides/shots_a.sh [OUT] [EXE]   (default OUT ./slides-b, EXE build/check/bin/game.exe); run from the repo root.
# Needs the PC's GPU and a quiet machine (no other game running): screenshots of a half-loaded world are useless.
set -u
OUT="${1:-slides-b}"
EXE="${2:-build/check/bin/game.exe}"
mkdir -p "$OUT"
COMMON="--no-hud --native-window"
# name | environment | arguments | the screenshot it writes
SHOTS=(
  "room-day||--seed 42 --bench 20 --shot-room --day|bench_room.png"
  "shelter-dusk||--seed 42 --bench 20 --shot-shelter --sun-elev -2|bench_shelter*.png"
  "campfire-night||--seed 42 --bench 20 --shot-campfire --sun-elev -8|bench_campfire.png"
  "glow-night||--seed 42 --bench 20 --shot-glow --sun-elev -8|bench_glow*.png"
  "terrain||--seed 42 --bench 30 --shot-terrain --day|bench_terrain*.png"
  "far-forest||--seed 42 --bench 30 --shot-far-forest --day|bench_far_forest.png"
  "clouds||--seed 42 --bench 30 --shot-clouds --day|bench_clouds_*.png"
  "karst-mouth|RAYLIBGAME_CAVES_KARST=1|--seed 42 --bench 40 --shot-cave-mouth --day|bench_mouth_*.png"
  "tube-mouth|RAYLIBGAME_CAVES_TUBE=1|--seed 42 --bench 40 --shot-cave-mouth --day|bench_mouth_*.png"
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
# Scenarios that frame a subject well take their picture themselves (c.grab, only with RAYLIBGAME_SLIDE_SHOTS set):
# a cave chamber, a karst sinkhole from above, the rotation line turning, the ladder climb, a thornback charging.
for sc in caves/chamber-view caves/sinkhole factory/rotation world/ladder creatures/arena-thornback; do
  t0=$SECONDS
  before=$(ls "$OUT"/*.png 2>/dev/null | wc -l)
  RAYLIBGAME_SLIDE_SHOTS="$OUT" timeout 400 "$EXE" $COMMON --scenario "$sc" --scenario-no-record --scenario-out "$OUT/scenario-runs" \
    > "$OUT/${sc//\//_}.log" 2>&1
  rc=$?
  echo "shots: $sc: exit $rc, $(( $(ls "$OUT"/*.png 2>/dev/null | wc -l) - before )) shot(s), $((SECONDS - t0)) s"
done
# The contact sheet: every shot, 4 across, 480 px wide each, its file name under it.
powershell -NoProfile -ExecutionPolicy Bypass -File "$(dirname "$0")/sheet.ps1" -Dir "$OUT" -Out "$OUT/sheet.png"
