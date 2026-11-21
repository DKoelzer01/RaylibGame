#!/bin/bash
# Slideshow shots, batch C (#1071): the hand-placed views of tests/scenarios/slides.cpp (a river valley, a lit cave chamber,
# a lava tube, the planet, its ring and the moon from space), each taking its pictures through c.grab. Then a contact sheet.
# Usage: tools/slides/shots_c.sh [OUT] [EXE]   (default OUT ./slides-c, EXE build/check/bin/game.exe); run from the repo root.
# Needs the PC's GPU and a quiet machine.
set -u
OUT="${1:-slides-c}"
EXE="${2:-build/check/bin/game.exe}"
mkdir -p "$OUT"
for sc in ${SLIDES:-slides/river-valley slides/cave-chamber slides/lava-tube slides/space}; do
  t0=$SECONDS
  before=$(ls "$OUT"/*.png 2>/dev/null | wc -l)
  RAYLIBGAME_SLIDE_SHOTS="$OUT" timeout 420 "$EXE" --no-hud --native-window --scenario "$sc" --scenario-no-record \
    --scenario-out "$OUT/scenario-runs" > "$OUT/${sc//\//_}.log" 2>&1
  rc=$?
  echo "shots: $sc: exit $rc, $(( $(ls "$OUT"/*.png 2>/dev/null | wc -l) - before )) shot(s), $((SECONDS - t0)) s"
done
powershell -NoProfile -ExecutionPolicy Bypass -File "$(dirname "$0")/sheet.ps1" -Dir "$OUT" -Out "$OUT/sheet.png"
