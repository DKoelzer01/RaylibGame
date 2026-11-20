#!/bin/bash
# A joined player rides the planet through gaps in its inputs (#844). A peaceful --host at --time-scale 60 (terra's day in
# a minute) and one light client that joins and stands still for 100 s. The server steps a player only by its inputs; it
# used to leave the player's world-frame velocity behind while the client loaded, and the first step "landed" it at
# 28-36 m/s: a death on a peaceful world. Checks: the server logged no death ("host: <name> died"), and the client joined.
# RAYLIBGAME_REST_GAP=0 turns the fix off (A/B: it dies again). EXE: a GNS build (default build/gns/bin/game.exe).
# PORT (27991). Prints "=== test_net_timescale: PASS|FAIL". Git Bash.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
exe="${EXE:-$here/build/gns/bin/game.exe}"
port="${PORT:-27991}"
tmp="$(mktemp -d)"
cd "$here"
SECONDS=0
printf 'erosion 0\n' > "$tmp/host.sav.worldgen"
. "$here/tools/lib_netmem.sh"
net_wait_for_memory 8 test_net_timescale || exit 0
RAYLIBGAME_HOST_DIFFICULTY=0 "$exe" --host "$tmp/host.sav" --seed 42 --bench 600 --time-scale 60 --port "$port" > "$tmp/host.log" 2>&1 &
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
timeout 300 "$exe" --join "127.0.0.1:$port" --name walker --bench 100 $NET_LIGHT_CLIENT > "$tmp/walker.log" 2>&1
stop_host
ok=1
chk() { if [ "$1" = 1 ]; then echo "  ok: $2"; else echo "  MISSING: $2"; ok=0; fi; }
grep -h "^host: .* died" "$tmp/host.log"
chk "$(grep -q "^join: joined" "$tmp/walker.log" && echo 1)" "the client joined"
chk "$(grep -q "^host: walker joined" "$tmp/host.log" && ! grep -q "^host: .* died" "$tmp/host.log" && echo 1)" "no death on the peaceful host at 60x"
echo "=== test_net_timescale: $([ $ok = 1 ] && echo PASS || echo FAIL) after ${SECONDS} s ($tmp)"
