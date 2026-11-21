#!/bin/bash
# Regression: the scripted --test-* gates and benchmarks, run one after another.
#   tools/regression.sh [planet|eroded|space|multiplayer|all] [--only headless|local] [--list]
#   (default set: all; erosion is the default terrain since #362: "eroded" keeps its own gates)
# Every test is tagged (#1028): H = headless, runs on the Linux headless/EGL build (the deweyserver runner's nightly job,
# .github/workflows/selfhosted.yml): simulation, fluids, building, saves, inventory and crafting, multiplayer, terrain
# logic. L = local, needs the Windows client on the PC: rendering and judged screenshots, camera and movement feel,
# perf benches, verdicts that read frame times or time the machine (ms gates), menus. When unsure: L.
#   --only headless   only the H tests (the runner's nightly)
#   --only local      only the L tests (the PC's midnight run)
#   --list            print the tag and arguments of every test in the set, run nothing
#   --shard I/N       only every Nth of the selected tests, starting with the Ith (1..N): the runner splits the headless
#                     half into jobs short enough to let other runs in between (#1028)
#   (no --only: every test, as before the split)
# Environment:
#   EXE=path      the game to run (default build/check/bin/game.exe)
#   SNAPSHOT=1    copy the exe and resources/ to a temp folder first, so you can keep editing and
#                 building while it runs (resources are read at run time)
#   RUNDIR=path   run in this folder instead of the checkout (it needs resources/ and assets/; saves/ and the erosion
#                 bake cache go there, so a persistent RUNDIR keeps the bake between runs)
#   EXE_GNS=path  multiplayer: the GameNetworkingSockets build (default build/gns/bin/game.exe, else EXE)
#   GNS=1         multiplayer: EXE itself has GNS (the Linux server build), so the tools/test_*.sh join tests run too
# Prints each test's "bench:" lines and ends with a PASS/FAIL count. Exit code: 1 if anything failed.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
set_name="all"; only=""; list=""; shard_i=1; shard_n=1; picked=0
while [ $# -gt 0 ]; do
    case "$1" in
        --only) only="${2:-}"; shift ;;
        --only=*) only="${1#--only=}" ;;
        --list) list=1 ;;
        --shard) shard_i="${2%/*}"; shard_n="${2#*/}"; shift ;;
        -*) echo "unknown option $1"; set_name="?" ;;
        *) set_name="$1" ;;
    esac
    shift
done
case "$only" in ""|headless|local) ;; *) echo "--only takes headless or local"; exit 2 ;; esac
case "$shard_i/$shard_n" in *[!0-9/]*|/*|*/) echo "--shard takes I/N"; exit 2 ;; esac
[ "$shard_n" -ge 1 ] && [ "$shard_i" -ge 1 ] && [ "$shard_i" -le "$shard_n" ] || { echo "--shard takes I/N with 1 <= I <= N"; exit 2; }
exe="${EXE:-$here/build/check/bin/game.exe}"
dir="${RUNDIR:-$here}"
if [ "${SNAPSHOT:-0}" = "1" ] && [ -z "$list" ]; then
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
# want TAG: does this run include tests tagged TAG? (with --list: print the line instead, the caller skips the test)
selected() { # TAG: in --only's half and in --shard's share (counts the tests of the half, so shards interleave)
    [ -z "$only" ] || [ "$only" = "$1" ] || return 1
    picked=$((picked + 1))
    [ $(( (picked - 1) % shard_n + 1 )) -eq "$shard_i" ]
}
want() {
    if [ -n "$list" ]; then return 1; fi
    selected "$1"
}
listed() { # TAG ARGS...: print one --list line, with the RAYLIBGAME_* variables the test runs under
    local tag="$1"; shift
    selected "$tag" || return 0
    printf '%-9s %s%s\n' "$tag" "$(env | grep -o '^RAYLIBGAME_[A-Z0-9_]*=[^ ]*' | sort | tr '\n' ' ')" "$*"
}
H() { if [ -n "$list" ]; then listed headless "$@"; elif want headless; then run "$@"; fi; }
L() { if [ -n "$list" ]; then listed local "$@"; elif want local; then run "$@"; fi; }
# A tools/test_*.sh script in the multiplayer set (tag, name for --list); the caller runs it when this returns 0.
script() { if [ -n "$list" ]; then listed "$1" "$2"; return 1; fi; want "$1"; }
# (#1367) A script that joins over the network: without a GameNetworkingSockets build (no EXE_GNS / build/gns, no GNS=1) it
# is reported as SKIP with the reason, never left out of the log (the net tests silently did not run for weeks).
netscript() {
    script "$@" || return 1
    if [ "$gexe" = "$exe" ] && [ "${GNS:-0}" != 1 ]; then
        echo "=== $2" | tee -a "$log"
        echo "bench: $2 -> SKIP (needs a GameNetworkingSockets build: EXE_GNS, build/gns/bin/game.exe or GNS=1)" | tee -a "$log"
        return 1
    fi
    return 0
}

# Comments after a test: why it is local (the H tests read game state, not pixels, frames or the clock's speed).
# "wall clock" marks a headless test whose phases or timeouts run on GetTime(): a much slower machine could change it.
# The flag-mode fluid tests run the fluid job at the machine's speed (only scenarios set RAYLIBGAME_FLUID_LOCKSTEP=1): the
# headless ones gave the PC's verdicts on the runner (#1028), --test-leak did not and is local until its scenario replaces it.
planet() {
    L --seed 42 --bench 20 --bench-walk # perf bench
    L --seed 42 --bench 10 # perf bench
    H --seed 42 --bench 5 --test-gravity
    H --seed 42 --bench 5 --test-mantle
    H --seed 42 --bench 5 --test-moon
    H --seed 42 --bench 5 --test-precision
    H --seed 42 --bench 5 --test-airjump
    L --seed 42 --bench 5 --test-stand 10 # camera jitter per frame, shadow texels slid per second
    H --seed 42 --bench 5 --test-mine 20
    L --seed 42 --bench 5 --test-suffocate # wall clock: must die within 15 s (failed on the slower runner, #1028)
    H --seed 42 --bench 3 --look-down --test-dig 10 --test-saveload
    L --seed 42 --bench 5 --test-stand 3 --look-down --crack-check # reads the screen for magenta holes
    L --seed 9 --bench 8 --test-stand 6 --telescope 16 # the telescope view, camera jitter
    H --seed 42 --bench 5 --test-planet-walk 1000 # wall clock: 45 s per leg (a lost leg is no failure)
    H --seed 42 --bench 3 --test-hole
    H --seed 9 --bench 3 --test-trees
    H --seed 9 --bench 3 --test-chop
    L --seed 42 --bench 3 --test-persist 100000 # gates on save/load/chunk-build milliseconds
    H --seed 9 --bench 6 --test-build 10000
    L --seed 9 --bench 3 --test-factory 10000 # gates on the factory tick (< 5 ms)
    H --seed 42 --bench 3 --test-machines
    H --test-machine-stack
    L --seed 42 --bench 3 --test-torch # ground luminance from screenshots
    L --seed 42 --bench 5 --test-lamps # ground luminance from screenshots
    H --seed 42 --bench 5 --test-fill --day
    L --seed 9 --bench 20 --test-structures --day # a screenshot tour that waits for streaming at each site (no verdict on the runner, #1028)
    H --seed 42 --bench 5 --test-fishing
    H --seed 42 --bench 5 --test-trapping
    H --seed 42 --bench 3 --test-combat
    H --seed 42 --bench 3 --test-bed
    H --seed 42 --bench 3 --test-save-split # (#1110: records never silently replaced; the MP <-> SP round trip)
    H --seed 42 --bench 5 --test-meteor
    H --seed 42 --bench 3 --test-animal-path
    H --seed 42 --bench 3 --test-progression # wall clock: 120 s per step (flaky: #372)
    H --seed 42 --bench 3 --test-progression full # wall clock: 120 s per step
    H --seed 42 --bench 3 --test-wall
    H --seed 9 --bench 3 --test-wall
    H --seed 42 --bench 3 --test-doorway # 2 m openings pass, 1.5 m ones stop the 1.7 m player (#255)
    H --seed 42 --bench 2 --test-facing # parts face the player until R/V (#196)
    L --seed 42 --bench 3 --test-ledge # movement feel: the fall off the edge timed on the wall clock
    H --seed 42 --bench 3 --test-fall
    H --seed 42 --bench 5 --test-hover
    L --seed 42 --bench 3 --test-telescope-close # worst frame after closing (< 25 ms)
    L --seed 1 --bench 3 --test-tor # movement feel: ground flicker under 0.15 s on the wall clock
    L --seed 42 --bench 3 --test-menus # the title screen and menus
    H --seed 42 --bench 3 --test-water
    H --seed 9 --bench 3 --test-water
    L --seed 42 --bench 3 --test-flyfollow # the fly camera following a spinning body
    H --seed 42 --bench 3 --test-hotbar
    H --seed 42 --bench 3 --test-netids # stable network ids over a save and load (#317)
    H --seed 42 --bench 3 --test-net # loopback connect via GameNetworkingSockets; SKIP without ENABLE_GNS (#310; the Linux build has it)
    H --seed 42 --bench 3 --test-two-players # a second Player 2 km away in its own frame (#314); wall clock: 15 m in 5 s
    H --seed 42 --bench 3 --test-fingerprint # world fingerprint, join info, gameplay RNG (#316)
    H --seed 42 --bench 3 --test-interest # tree tiles and animals around a second centre 3 km away (#318)
    L --seed 9 --bench 20 --bench-fly 200 # perf bench
    L --seed 9 --bench 10 --bench-fly 200 --no-foliage # perf bench
    L --seed 42 --bench 3 --test-smallbody # wall clock: landing on the pebble timed in seconds (failed on the runner, #1028)
    H --seed 42 --bench 2 --test-inventory
    H --seed 42 --bench 3 --test-map
    L --seed 42 --bench 3 --test-viewmodel # the hands' bob and screenshots of them
    H --seed 42 --bench 2 --test-quests
    H --seed 42 --bench 3 --test-hunger
    H --seed 42 --bench 3 --test-farm
    L --seed 42 --bench 3 --test-animals # gates on the update time (200 animals < 1 ms; flee distance flaky: #373)
    H --seed 42 --bench 3 --test-suit
    H --seed 42 --bench 3 --test-tools
    H --seed 42 --bench 3 --test-rivers
    H --seed 9 --bench 3 --test-rivers
    H --seed 42 --bench 3 --test-threats
    H --seed 9 --bench 3 --test-threats
    H --seed 42 --bench 2 --test-recipes
    H --seed 42 --bench 3 --test-crafting
    H --seed 42 --bench 3 --test-materials
    H --seed 42 --bench 3 --test-decor
    H --seed 42 --bench 3 --test-weather
    L --seed 42 --bench 3 --test-flow # gates on the fluid step time (worst 1 s mean < 30 ms)
    RAYLIBGAME_FLUID_CELL=1 L --seed 42 --bench 3 --test-flow # the 1 m fluid grid (opt-in since #714) still runs; step-time gate
    RAYLIBGAME_EROSION=1 L --seed 42 --bench 5 --test-flow-sea # the erosion bake's rivers (#190; the noise terrain's river mouths stand over the sea); --test-flow's step-time gate
    RAYLIBGAME_EROSION=1 H --seed 42 --bench 3 --test-dam # #361/#432: a low dam backs a released river reach up and it spills over (eroded valleys; up to 7 min); wall clock
    H --seed 42 --day --bench 3 --test-dam-divert # #656: a dam that holds a river dries it 200 m downstream, a second dam upstream cuts the first reach's inflow, the state survives save/load, dug out the river comes back: flow past 0.5 within 60 s and still rising at 90 s (~3 min; --dam-divert-full waits for the full return); wall clock
    H --seed 42 --bench 3 --test-pond
    H --seed 42 --bench 3 --test-part-wall # #391: a wall of build blocks holds water; taken down, the water runs on
    L --sandbox --seed 42 --day --bench 3 --test-leak # the fluid job runs at the machine's speed: the runner reached 0.55 of the orifice rate, the PC 0.79 (#1028); #680/#691: a 1 m hole in a dug wall drains like an orifice, then a weir; no water on the grass
    H --sandbox --seed 42 --day --bench 3 --test-leak-parts # #692/#694: the same through a dam of placed blocks (at 1 m cells, opt-in, a one-block hole off the fluid grid stays shut)
    H --sandbox --seed 42 --day --bench 3 --test-dam-break # #680: Ritter front, depth and discharge (reported; the pass is the model work on #680)
    H --seed 42 --bench 3 --test-flow-meter # #392: the flow readout reads a 1 m^3/s stream
    H --seed 42 --bench 3 --test-pressure # #393: hydrostatic pressure in a filled shaft and a pond
    H --seed 42 --day --bench 3 --test-bucket
    H --seed 42 --bench 3 --test-pump
    H --seed 42 --bench 5 --test-pump-hand
    H --seed 42 --bench 3 --test-pump-boiler # a full boiler takes no world water (#229)
    H --seed 42 --bench 3 --test-oil
    H --seed 42 --bench 3 --test-spring
}

# The erosion-baked world (Plan 7 gate, #362): the gates the #355 review asked for, with RAYLIBGAME_EROSION=1. The first
# run bakes the world (~1-2 min, saved as saves/bench.sav.erosion in the run folder); the rest load it (~0.1 s).
eroded() {
    export RAYLIBGAME_EROSION=1
    L --seed 42 --bench 10 # perf bench
    L --seed 42 --bench 5 --test-stand 3 --look-down --crack-check # reads the screen
    H --seed 42 --bench 5 --test-planet-walk 1000 # wall clock: 45 s per leg
    H --seed 42 --bench 3 --test-pond
    H --seed 42 --bench 3 --test-rivers
    L --seed 42 --bench 5 --test-flow-sea # step-time gate
    H --seed 42 --bench 3 --test-build 2000
    H --seed 42 --bench 3 --test-placement
    L --seed 42 --bench 20 --bench-walk # perf bench
    L --seed 42 --bench 20 --bench-fly 60 # perf bench
    unset RAYLIBGAME_EROSION
}

# Asteroid and space tests: opt-in (the user leaves them out of normal runs); kept on the PC as before the split.
space() {
    L --seed 42 --bench 3 --test-ring
    L --seed 42 --bench 3 --test-blueprint
    L --seed 42 --bench 3 --time-scale 8 --test-rocket
    L --seed 9 --bench 3 --time-scale 8 --test-rocket
    L --seed 42 --bench 3 --test-damage
    L --seed 42 --bench 3 --time-scale 8 --test-station
    L $A --bench 5 --test-orbit 10
    L $A --bench 5 --test-ride
}

# Multiplayer (#347): a host and in-process clients over loopback. The network tests need a GNS build (EXE_GNS, default
# build/gns/bin/game.exe; configure with -DENABLE_GNS=ON after tools/build_gns.sh), else they SKIP on the plain exe.
# The Linux server build has GNS: there EXE_GNS = EXE.
multiplayer() {
    local gexe="${EXE_GNS:-$here/build/gns/bin/game.exe}"; [ -x "$gexe" ] || gexe="$exe"
    local sv; sv="$(mktemp -d)/mp.sav"
    local keep="$exe"; exe="$gexe"
    H --host "$sv" --seed 42 --bench 3 --port 27941 --test-host
    H --host "$sv" --bench 3 --port 27941 --test-net-move
    H --host "$sv" --bench 3 --port 27941 --test-net-dig
    H --host "$sv" --bench 3 --port 27941 --test-net-build
    H --host "$sv" --bench 3 --port 27941 --test-net-sync # (#776) joining a world with many edits: indexed, hashed, rate-bounded
    RAYLIBGAME_NET_LAG_MS=50 RAYLIBGAME_NET_JITTER_MS=10 RAYLIBGAME_NET_LOSS=2 RAYLIBGAME_NET_SEED=1 H --host "$sv" --bench 3 --port 27941 --test-mp-bots 6
    exe="$keep"
    H --seed 42 --bench 3 --test-time-follow
    H --test-net-interest
    H --test-net-guard # (#349) per-client message limits: a normal client never dropped; floods and garbage kicked
    exe="$gexe"; RAYLIBGAME_MAX_PER_ADDRESS=6 H --host "$sv" --sandbox --bench 3 --port 27941 --test-net-abuse; exe="$keep" # (#349, #780, #781; names, chat, the per-address cap) needs GNS (else SKIP)
    H --test-server-admin
    H --seed 42 --bench 3 --test-player-cost 8 # (#768) player cost per update (reported), and determinism: a warm vs a poisoned collision cache stay bit-identical
    # The dedicated server (#326): config file, console commands, log file, "stop" with a final save (sandbox world).
    if script headless tools/test_server.sh; then
        echo "=== tools/test_server.sh" | tee -a "$log"
        EXE="$gexe" bash "$here/tools/test_server.sh" 2>&1 | grep "=== test_server" | sed 's/=== test_server: \([A-Z]*\).*/bench: server (#326): config, console, log, stop -> \1/' | tee -a "$log"
    fi
    # Admin end to end (#326 part 2, #345): chat, motd, an op's ban from chat kicks, refusals, the whitelist, the console's kick.
    if netscript headless tools/test_server_admin.sh; then
        echo "=== tools/test_server_admin.sh" | tee -a "$log"
        EXE="$gexe" bash "$here/tools/test_server_admin.sh" 2>&1 | grep "=== test_server_admin" | sed 's/=== test_server_admin: \([A-Z]*\).*/bench: server admin end to end (#345): -> \1/' | tee -a "$log"
    fi
    # --join end to end: a host and two clients as separate processes (#341 core); only with a GNS build.
    if netscript headless tools/test_join.sh; then
        echo "=== tools/test_join.sh" | tee -a "$log"
        EXE="$gexe" bash "$here/tools/test_join.sh" 2>&1 | grep "^bench: join (" | tee -a "$log"
    fi
    if netscript headless "tools/test_join.sh on a non-default world (#757)"; then
        echo "=== tools/test_join.sh on a non-default world (#757)" | tee -a "$log"
        WORLDGEN="ocean 0.4 mountains 1.6 cliffs 0 karst 0 erosion 0" PORT=27952 EXE="$gexe" bash "$here/tools/test_join.sh" 2>&1 | grep "^bench: join (" | tee -a "$log"
    fi
    # Join-in-progress sync (#776): a client joins a host world with 3000 edits (fingerprints must match), then rejoins from its cache (fetching nothing).
    if netscript headless tools/test_join_sync.sh; then
        echo "=== tools/test_join_sync.sh" | tee -a "$log"
        EXE="$gexe" bash "$here/tools/test_join_sync.sh" 2>&1 | grep "=== test_join_sync" | sed 's/=== test_join_sync: \([A-Z]*\).*/bench: join sync (#776): fresh join, then a rejoin from the cache -> \1/' | tee -a "$log"
    fi
    # Plants and farming on the server (#335), on the quick old terrain and on the default terrain (#1107: plants near the
    # spawn of a menu-hosted game couldn't be felled)
    if netscript headless "tools/test_net_plants.sh (#335)"; then
        echo "=== tools/test_net_plants.sh (#335)" | tee -a "$log"
        EXE="$gexe" bash "$here/tools/test_net_plants.sh" 2>&1 | grep "=== test_net_plants" | sed 's/=== test_net_plants: \([A-Z]*\).*/bench: plants and farming on the server (#335): -> \1/' | tee -a "$log"
    fi
    if netscript headless "tools/test_net_plants.sh on the default terrain (#1107)"; then
        echo "=== tools/test_net_plants.sh on the default terrain (#1107)" | tee -a "$log"
        TERRAIN=default PORT=27969 EXE="$gexe" bash "$here/tools/test_net_plants.sh" 2>&1 | grep "=== test_net_plants" | sed 's/=== test_net_plants: \([A-Z]*\).*/bench: plants near the spawn on the default terrain (#1107): -> \1/' | tee -a "$log"
    fi
    # A menu-hosted save loaded in single player (#901): the server's record holds where the player left, and single player
    # puts them there (as the same player)
    if netscript headless "tools/test_hosted_load.sh (#901)"; then
        echo "=== tools/test_hosted_load.sh (#901)" | tee -a "$log"
        PORT=27976 EXE="$gexe" bash "$here/tools/test_hosted_load.sh" 2>&1 | grep "=== test_hosted_load\|  FAIL" | sed 's/=== test_hosted_load: \([A-Z]*\).*/bench: a hosted save loaded in single player keeps the position (#901): -> \1/' | tee -a "$log"
    fi
    if netscript headless "tools/test_animals_net.sh (#334)"; then
        echo "=== tools/test_animals_net.sh (#334)" | tee -a "$log"
        EXE="$gexe" bash "$here/tools/test_animals_net.sh" 2>&1 | grep "=== test_animals_net" | sed 's/=== test_animals_net: \([A-Z]*\).*/bench: animals on the server (#334): same herd, a hit applied for both -> \1/' | tee -a "$log"
    fi
    if netscript headless "tools/test_animals_bite.sh (#334)"; then
        echo "=== tools/test_animals_bite.sh (#334)" | tee -a "$log"
        EXE="$gexe" bash "$here/tools/test_animals_bite.sh" 2>&1 | grep "=== test_animals_bite" | sed 's/=== test_animals_bite: \([A-Z]*\).*/bench: an animal bites a remote player through the server (#334): -> \1/' | tee -a "$log"
    fi
}

case "$set_name" in
    multiplayer) multiplayer ;;
    planet) planet; eroded ;;
    eroded) eroded ;;
    space) space ;;
    all) planet; eroded; space; multiplayer ;;
    *) echo "usage: $0 [planet|eroded|space|multiplayer|all] [--only headless|local] [--list]"; exit 2 ;;
esac
if [ -n "$list" ]; then rm -f "$log"; exit 0; fi
# (#1367) A tools/ script that ended without a PASS, FAIL or SKIP line (it crashed, timed out, or its result line changed
# shape) is a FAIL, not a test that silently did not happen.
awk '/^=== tools\// { if (cur != "" && !got) print "bench: " cur " -> FAIL (no result line)"; cur = substr($0, 5); got = 0; next }
     /PASS|FAIL|SKIP/ { got = 1 }
     END { if (cur != "" && !got) print "bench: " cur " -> FAIL (no result line)" }' "$log" > "$log.missing"
if [ -s "$log.missing" ]; then tee -a "$log" < "$log.missing"; fi
rm -f "$log.missing"
pass=$(grep -c "PASS" "$log"); fail=$(grep -c "FAIL" "$log"); skip=$(grep -c "SKIP" "$log")
what="$set_name${only:+ --only $only}"; [ "$shard_n" -gt 1 ] && what="$what --shard $shard_i/$shard_n"
echo "=== regression ($what): $pass PASS, $fail FAIL, $skip SKIP"
grep "FAIL" "$log"
rm -f "$log"
[ "$fail" -eq 0 ]
