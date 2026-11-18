#!/bin/bash
# Animals on the server (#334), end to end: a --host with a test herd beside the spawn, and two --join clients as separate
# processes. Both must see the same herd (the same animal ids, from the server); the hunter walks to one and swings at it -
# requests the server checks and applies - until it falls; the watcher must see that animal's health drop and then the animal
# go; the host must credit the drops to the hunter. EXE: a GNS build (default build/gns/bin/game.exe). PORT (27961).
# SHOT=<file name>: the watcher takes a screenshot of the herd at the end (in the working directory). Prints the checks and
# "=== test_animals_net: PASS|FAIL". Capped at 10 minutes.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
exe="${EXE:-$here/build/gns/bin/game.exe}"
port="${PORT:-27961}"
tmp="$(mktemp -d)"
cd "$here"
SECONDS=0
printf 'erosion 0\n' > "$tmp/host.sav.worldgen" # (a quick world: no bake)
. "$here/tools/lib_netmem.sh" # (#837) light clients, and enough memory for every instance first
net_wait_for_memory 12 test_animals_net || exit 0
RAYLIBGAME_HOST_HERD="woolback 4" "$exe" --host "$tmp/host.sav" --seed 42 --bench 600 --port "$port" > "$tmp/host.log" 2>&1 &
hp=$!
# #795: on Windows a signal to the msys PID can miss the native game.exe; stop it by its Windows PID as well.
stop_host() {
    local w; w=$(cat "/proc/$hp/winpid" 2>/dev/null)
    kill "$hp" 2>/dev/null
    for _ in 1 2 3 4 5; do kill -0 "$hp" 2>/dev/null || break; sleep 1; done
    wait "$hp" 2>/dev/null
    # only our own host: the PID still names a game.exe moments after we started it
    if [ -n "$w" ] && tasklist //FI "PID eq $w" //FI "IMAGENAME eq game.exe" 2>/dev/null | grep -q game.exe; then taskkill //F //PID "$w" >/dev/null 2>&1; fi
    true
}
trap stop_host EXIT
until grep -q "test herd" "$tmp/host.log" 2>/dev/null || ! kill -0 "$hp" 2>/dev/null || [ $SECONDS -gt 240 ]; do sleep 1; done
grep -h "test herd" "$tmp/host.log" || { echo "test_animals_net: no herd on the host -> FAIL"; stop_host; exit 1; }
cpids=()
RAYLIBGAME_ANIMALS_SHOT="${SHOT:-}" RAYLIBGAME_NET_ANIMALS_ROLE=watcher RAYLIBGAME_JOIN_LINGER=50 timeout 300 "$exe" --join "127.0.0.1:$port" --name watcher --bench 3 $NET_LIGHT_CLIENT --test-net-animals > "$tmp/watcher.log" 2>&1 &
cpids+=($!)
sleep 3
RAYLIBGAME_NET_ANIMALS_ROLE=hunter timeout 300 "$exe" --join "127.0.0.1:$port" --name hunter --bench 3 $NET_LIGHT_CLIENT --test-net-animals > "$tmp/hunter.log" 2>&1 &
cpids+=($!)
for c in "${cpids[@]}"; do wait "$c"; done
sleep 1; stop_host
grep -h "^bench: net animals" "$tmp/hunter.log" "$tmp/watcher.log"
grep -h "felled animal" "$tmp/host.log" | sed 's/^/host: /'
ok=1
target=$(grep -o "target [0-9a-f]* fell" "$tmp/hunter.log" | awk '{print $2}')
hunterIds=$(grep -o "saw [0-9]* animals:.*" "$tmp/hunter.log" | cut -d: -f2)
watcherIds=$(grep -o "saw [0-9]* animals:.*" "$tmp/watcher.log" | cut -d: -f2)
chk() { if [ "$1" = 1 ]; then echo "  ok: $2"; else echo "  MISSING: $2"; ok=0; fi; }
chk "$([ -n "$hunterIds" ] && [ -n "$watcherIds" ] && echo 1)" "both clients saw animals"
same=1; for id in $hunterIds; do echo " $watcherIds " | grep -q " $id " || same=0; done
chk "$([ -n "$hunterIds" ] && [ $same = 1 ] && echo 1)" "the same herd: every animal the hunter saw, the watcher saw too (ids from the server)"
chk "$([ -n "$target" ] && echo 1)" "the hunter's target fell (the server applied its swings)"
wl=$(grep "(watcher): $target health" "$tmp/watcher.log")
first=$(echo "$wl" | awk '{print $7}'); low=$(echo "$wl" | awk '{print $9}' | tr -d ,); gone=$(echo "$wl" | awk '{print $11}')
chk "$([ -n "$wl" ] && awk -v a="$first" -v b="$low" 'BEGIN{exit !(b < a)}' && echo 1)" "the watcher saw its health drop ($first -> $low)"
chk "$([ "$gone" = 1 ] && echo 1)" "and saw it go"
chk "$(grep -q "felled animal $target" "$tmp/host.log" && echo 1)" "the host credited the drops to the hunter"
echo "=== test_animals_net: $([ $ok = 1 ] && echo PASS || echo FAIL) after ${SECONDS} s ($tmp)"
[ $ok = 1 ]
