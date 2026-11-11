#!/bin/bash
# Regression: the scripted --test-* gates and benchmarks, run one after another.
#   tools/regression.sh [planet|space|all]        (default: all)
# Environment:
#   EXE=path      the game to run (default build/check/bin/game.exe)
#   SNAPSHOT=1    copy the exe and resources/ to a temp folder first, so you can keep editing and
#                 building while it runs (resources are read at run time)
# Prints each test's "bench:" lines and ends with a PASS/FAIL count. Exit code: 1 if anything failed.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
set_name="${1:-all}"
exe="${EXE:-$here/build/check/bin/game.exe}"
dir="$here"
if [ "${SNAPSHOT:-0}" = "1" ]; then
    dir="$(mktemp -d)"
    cp -r "$here/resources" "$dir/"
    cp "$exe" "$dir/game.exe"
    exe="$dir/game.exe"
    echo "snapshot: $dir"
fi
cd "$dir"
A="--bodies resources/data/bodies_asteroids.txt"
log="$(mktemp)"
run() { local t0=$(date +%s); echo "=== $*" | tee -a "$log"; timeout 1200 "$exe" "$@" 2>&1 | grep "^bench" | grep -v "world ready" | tee -a "$log"; echo "--- $(( $(date +%s) - t0 )) s"; }

planet() {
    run --seed 42 --bench 20 --bench-walk
    run --seed 42 --bench 10
    run --seed 42 --bench 5 --test-gravity
    run --seed 42 --bench 5 --test-mantle
    run --seed 42 --bench 5 --test-moon
    run --seed 42 --bench 5 --test-precision
    run --seed 42 --bench 5 --test-airjump
    run --seed 42 --bench 5 --test-stand 10
    run --seed 42 --bench 5 --test-mine 20
    run --seed 42 --bench 5 --test-suffocate
    run --seed 42 --bench 3 --look-down --test-dig 10 --test-saveload
    run --seed 42 --bench 5 --test-stand 3 --look-down --crack-check
    run --seed 9 --bench 8 --test-stand 6 --telescope 16
    run --seed 42 --bench 5 --test-planet-walk 1000
    run --seed 42 --bench 3 --test-hole
    run --seed 9 --bench 3 --test-trees
    run --seed 9 --bench 3 --test-chop
    run --seed 42 --bench 3 --test-persist 100000
    run --seed 9 --bench 6 --test-build 10000
    run --seed 9 --bench 3 --test-factory 10000
    run --seed 42 --bench 3 --test-progression
    run --seed 42 --bench 3 --test-progression full
    run --seed 42 --bench 3 --test-wall
    run --seed 9 --bench 3 --test-wall
    run --seed 42 --bench 3 --test-menus
    run --seed 42 --bench 3 --test-water
    run --seed 9 --bench 3 --test-water
    run --seed 42 --bench 3 --test-flyfollow
    run --seed 42 --bench 3 --test-hotbar
    run --seed 9 --bench 20 --bench-fly 200
    run --seed 9 --bench 10 --bench-fly 200 --no-foliage
    run --seed 42 --bench 3 --test-smallbody
    run --seed 42 --bench 2 --test-inventory
    run --seed 42 --bench 3 --test-map
    run --seed 42 --bench 3 --test-viewmodel
    run --seed 42 --bench 2 --test-quests
    run --seed 42 --bench 3 --test-hunger
    run --seed 42 --bench 3 --test-farm
    run --seed 42 --bench 3 --test-animals
    run --seed 42 --bench 3 --test-suit
    run --seed 42 --bench 3 --test-tools
    run --seed 42 --bench 3 --test-rivers
    run --seed 9 --bench 3 --test-rivers
    run --seed 42 --bench 3 --test-threats
    run --seed 9 --bench 3 --test-threats
    run --seed 42 --bench 2 --test-recipes
    run --seed 42 --bench 3 --test-crafting
}

space() {
    run --seed 42 --bench 3 --test-ring
    run --seed 42 --bench 3 --test-blueprint
    run --seed 42 --bench 3 --time-scale 8 --test-rocket
    run --seed 9 --bench 3 --time-scale 8 --test-rocket
    run $A --bench 5 --test-orbit 10
    run $A --bench 5 --test-ride
}

case "$set_name" in
    planet) planet ;;
    space) space ;;
    all) planet; space ;;
    *) echo "usage: $0 [planet|space|all]"; exit 2 ;;
esac
pass=$(grep -c "PASS" "$log"); fail=$(grep -c "FAIL" "$log")
echo "=== regression ($set_name): $pass PASS, $fail FAIL"
grep "FAIL" "$log"
rm -f "$log"
[ "$fail" -eq 0 ]
