#!/bin/bash
# Team exploring on the server (#338 slice 3). A peaceful --host puts "walker" and "stroller" in one team; the walker
# sprints ~600 m east, the stroller ~100 m west, both from the spawn. The server reveals the map tiles around its own view
# of each player (the local map's tiles and radius) into one set per team per body and raises the team's "explored" to its
# size. Checks: the team's total (in each client's book) exceeds what either player revealed alone (its own map), and is
# less than the two together (the tiles round the spawn they share count once).
# EXE: a GNS build (default build/gns/bin/game.exe). PORT (27966). Prints "=== test_team_explore: PASS|FAIL". Git Bash.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
exe="${EXE:-$here/build/gns/bin/game.exe}"
port="${PORT:-27966}"
tmp="$(mktemp -d)"
cd "$here"
SECONDS=0
printf 'erosion 0\n' > "$tmp/host.sav.worldgen"
. "$here/tools/lib_netmem.sh"
net_wait_for_memory 12 test_team_explore || exit 0
RAYLIBGAME_HOST_DIFFICULTY=0 RAYLIBGAME_HOST_TEAM="walker,stroller" "$exe" --host "$tmp/host.sav" --seed 42 --bench 600 --port "$port" > "$tmp/host.log" 2>&1 &
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
join() { RAYLIBGAME_QUESTS_SECONDS=90 RAYLIBGAME_NET_QUESTS_ROLE=$1 timeout 300 "$exe" --join "127.0.0.1:$port" --name "$1" --bench 3 $NET_LIGHT_CLIENT --test-net-quests > "$tmp/$1.log" 2>&1; }
join stroller & sp=$!
join walker
wait $sp
stop_host
ok=1
chk() { if [ "$1" = 1 ]; then echo "  ok: $2"; else echo "  MISSING: $2"; ok=0; fi; }
w=$(grep -h "^bench: net explore (walker" "$tmp/walker.log"); s=$(grep -h "^bench: net explore (stroller" "$tmp/stroller.log"); echo "$w"; echo "$s"
grep "revealed .* map tiles itself" "$tmp/host.log"
# The server's own numbers: what each player revealed alone (its tiles), and the team's total when each left.
sess() { grep -o "$1 joined as session [0-9]*" "$tmp/host.log" | head -1 | awk '{print $NF}'; }
alone() { grep "host: session $1 revealed" "$tmp/host.log" | head -1 | grep -o "revealed [0-9]*" | awk '{print $2}'; }
ws=$(sess walker); ss=$(sess stroller); aw=$(alone "$ws"); as=$(alone "$ss")
team=$(grep -o "its team explored [0-9]*" "$tmp/host.log" | awk '{print $NF}' | sort -n | tail -1)
wt=$(echo "$w" | grep -o "team explored [0-9]*" | awk '{print $NF}'); st=$(echo "$s" | grep -o "team explored [0-9]*" | awk '{print $NF}')
chk "$([ -n "$aw" ] && [ -n "$as" ] && [ "$aw" -gt "$as" ] && echo 1)" "the walker revealed more ground than the stroller on the server ($aw vs $as tiles)"
chk "$([ -n "$team" ] && [ "$team" -gt "$aw" ] && [ "$team" -gt "$as" ] && echo 1)" "the team's total ($team) exceeds what either revealed alone ($aw, $as): their tiles are merged"
chk "$([ -n "$team" ] && [ "$team" -lt $((aw + as)) ] && echo 1)" "and is less than the two together ($((aw + as))): the tiles round the spawn they share count once"
chk "$([ -n "$wt" ] && [ -n "$st" ] && [ "$wt" -gt "$aw" ] && [ "$st" -gt "$as" ] && echo 1)" "both clients' books show the team's total, more than their own share ($wt, $st)"
echo "=== test_team_explore: $([ $ok = 1 ] && echo PASS || echo FAIL) after ${SECONDS} s ($tmp)"
[ $ok = 1 ]
