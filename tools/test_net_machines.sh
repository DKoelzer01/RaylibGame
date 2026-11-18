#!/bin/bash
# Machines on the server (#332). A peaceful --host gives every joiner a campfire, 2 roots and 2 oak logs. "feeder" places the
# campfire through the server and puts the roots in and the logs in its fuel slot, by request; the server smelts them (a
# root a 10 s); "taker" takes the outputs (O) after 40 s. Checks: the feeder's items left by deltas; both clients heard the
# campfire work (MachineStates) and the feeder got its view (MachineView); the taker's 2 roasted roots arrived by deltas;
# the feeder got none of them.
# EXE: a GNS build (default build/gns/bin/game.exe). PORT (27967). Prints "=== test_net_machines: PASS|FAIL". Git Bash.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
exe="${EXE:-$here/build/gns/bin/game.exe}"
port="${PORT:-27967}"
tmp="$(mktemp -d)"
cd "$here"
SECONDS=0
printf 'erosion 0\n' > "$tmp/host.sav.worldgen"
. "$here/tools/lib_netmem.sh"
net_wait_for_memory 12 test_net_machines || exit 0
RAYLIBGAME_HOST_DIFFICULTY=0 RAYLIBGAME_HOST_GIVE="campfire 1,roots 2,oak_log 2" "$exe" --host "$tmp/host.sav" --seed 42 --bench 600 --port "$port" > "$tmp/host.log" 2>&1 &
hp=$!
stop_host() { # #795: by the Windows PID too, only while it is still our game.exe
    local w; w=$(cat "/proc/$hp/winpid" 2>/dev/null)
    kill "$hp" 2>/dev/null
    for _ in 1 2 3 4 5; do kill -0 "$hp" 2>/dev/null || break; sleep 1; done
    wait "$hp" 2>/dev/null
    if [ -n "$w" ] && tasklist //FI "PID eq $w" //FI "IMAGENAME eq game.exe" 2>/dev/null | grep -q game.exe; then taskkill //F //PID "$w" >/dev/null 2>&1; fi
    true
}
trap stop_host EXIT
until grep -q "listening on port" "$tmp/host.log" 2>/dev/null || ! kill -0 "$hp" 2>/dev/null || [ $SECONDS -gt 240 ]; do sleep 1; done
until grep -q "world ready" "$tmp/host.log" 2>/dev/null || [ $SECONDS -gt 280 ]; do sleep 1; done
join() { RAYLIBGAME_MACHINES_SECONDS=55 RAYLIBGAME_MACHINES_TAKE_AT=42 RAYLIBGAME_NET_MACHINES_ROLE=$1 timeout 300 "$exe" --join "127.0.0.1:$port" --name "$1" --bench 3 $NET_LIGHT_CLIENT --test-net-machines > "$tmp/$1.log" 2>&1; }
join taker & tp=$!
sleep 2
join feeder
wait $tp
stop_host
ok=1
chk() { if [ "$1" = 1 ]; then echo "  ok: $2"; else echo "  MISSING: $2"; ok=0; fi; }
f=$(grep -h "^bench: net machines (feeder" "$tmp/feeder.log"); t=$(grep -h "^bench: net machines (taker" "$tmp/taker.log"); echo "$f"; echo "$t"
chk "$(echo "$f" | grep -q "campfire seen 1, placed 1, fed 1" && echo 1)" "the feeder placed the campfire through the server and fed it"
chk "$(echo "$f" | grep -q "roots 0, logs 0" && echo 1)" "the roots and logs left its inventory (the server's deltas)"
chk "$(echo "$f" | grep -q "saw it working 1" && echo "$t" | grep -q "saw it working 1" && echo 1)" "both clients heard it working (MachineStates)"
chk "$(echo "$f" | grep -q "views [1-9]" && echo 1)" "the feeder got the campfire's view (MachineView, its watch)"
chk "$(echo "$t" | grep -q "took 1: roasted roots 0 -> 2" && echo 1)" "the taker took the 2 roasted roots out (O): they arrived by deltas"
chk "$(echo "$f" | grep -q "roasted roots -1 -> 0" && echo 1)" "and the feeder got none of them"
echo "=== test_net_machines: $([ $ok = 1 ] && echo PASS || echo FAIL) after ${SECONDS} s ($tmp)"
[ $ok = 1 ]
