#!/bin/bash
# Nights survived count for the team on the server (#338). A peaceful --host (no night threats) declares night from 12 to
# 24 s after the first join (RAYLIBGAME_HOST_NIGHT: a real night is half an hour); a --join "sleeper" stays 40 s. The server sees dark fall and day come back for that player without a death, and counts a night for their team;
# the client's book (sent by the server) shows it.
# EXE: a GNS build (default build/gns/bin/game.exe). PORT (27965). Prints "=== test_team_nights: PASS|FAIL". Git Bash.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
exe="${EXE:-$here/build/gns/bin/game.exe}"
port="${PORT:-27965}"
tmp="$(mktemp -d)"
cd "$here"
SECONDS=0
printf 'erosion 0\n' > "$tmp/host.sav.worldgen"
. "$here/tools/lib_netmem.sh"
net_wait_for_memory 9 test_team_nights || exit 0
RAYLIBGAME_HOST_DIFFICULTY=0 RAYLIBGAME_HOST_NIGHT="12,24" "$exe" --host "$tmp/host.sav" --seed 42 --bench 600 --port "$port" > "$tmp/host.log" 2>&1 &
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
RAYLIBGAME_QUESTS_SECONDS=40 RAYLIBGAME_NET_QUESTS_ROLE=sleeper timeout 240 "$exe" --join "127.0.0.1:$port" --name sleeper --bench 3 $NET_LIGHT_CLIENT --test-net-quests > "$tmp/sleeper.log" 2>&1
stop_host
ok=1
chk() { if [ "$1" = 1 ]; then echo "  ok: $2"; else echo "  MISSING: $2"; ok=0; fi; }
line=$(grep -h "^bench: net nights" "$tmp/sleeper.log"); echo "$line"
n=$(echo "$line" | grep -o "): [0-9]* nights" | grep -o "[0-9]*")
chk "$(echo "$line" | grep -q "deaths 0" && echo 1)" "the sleeper didn't die (peaceful)"
chk "$([ "$n" = 1 ] && echo 1)" "the server counted the night survived (dark from 12 to 24 s) for its team, once, and the book it sent shows it ($n)"
echo "=== test_team_nights: $([ $ok = 1 ] && echo PASS || echo FAIL) after ${SECONDS} s ($tmp)"
[ $ok = 1 ]
