#!/bin/bash
# Team quest books on the server (#338, #816). A --host puts alice and bob in one team (RAYLIBGAME_HOST_TEAM) and gives
# every joiner berries, hungry; carol is a team of her own. alice eats once through the server: the server counts it for
# her team, finishes the "eat something" node (e0f), pays its 4 berries to alice by deltas and tells the team. Checks:
# alice's book has e0f done and her berries went 3 -> 6 (one eaten, four paid); bob's book (the same team's) has it done
# and his berries are untouched; carol's book doesn't; the team heard the news.
# EXE: a GNS build (default build/gns/bin/game.exe). PORT (27964). Prints "=== test_team_quests: PASS|FAIL". Git Bash.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
exe="${EXE:-$here/build/gns/bin/game.exe}"
port="${PORT:-27964}"
tmp="$(mktemp -d)"
cd "$here"
SECONDS=0
printf 'erosion 0\n' > "$tmp/host.sav.worldgen"
. "$here/tools/lib_netmem.sh" # (#837) light clients, and enough memory for every instance first
net_wait_for_memory 15 test_team_quests || exit 0
RAYLIBGAME_HOST_TEAM="alice,bob" RAYLIBGAME_HOST_HUNGER=50 RAYLIBGAME_HOST_GIVE="berries 3,oak_log 2" "$exe" --host "$tmp/host.sav" --seed 42 --day --bench 600 --port "$port" > "$tmp/host.log" 2>&1 &
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
[ -d /proc ] || echo "test_team_quests: warning: not Git Bash (no /proc): the host may be left running" >&2
until grep -q "listening on port" "$tmp/host.log" 2>/dev/null || ! kill -0 "$hp" 2>/dev/null || [ $SECONDS -gt 240 ]; do sleep 1; done
until grep -q "world ready" "$tmp/host.log" 2>/dev/null || [ $SECONDS -gt 280 ]; do sleep 1; done
sleep 2
join() { RAYLIBGAME_QUESTS_SECONDS=30 RAYLIBGAME_NET_QUESTS_ROLE=$2 timeout 200 "$exe" --join "127.0.0.1:$port" --name "$1" --bench 3 $NET_LIGHT_CLIENT --test-net-quests > "$tmp/$1.log" 2>&1; }
join bob mate & bp=$!
join carol other & cp=$!
sleep 4
join alice eater
wait $bp $cp
stop_host
ok=1
chk() { if [ "$1" = 1 ]; then echo "  ok: $2"; else echo "  MISSING: $2"; ok=0; fi; }
for p in alice bob carol; do grep -h "^bench: net quests" "$tmp/$p.log"; done
a=$(grep -h "^bench: net quests" "$tmp/alice.log"); b=$(grep -h "^bench: net quests" "$tmp/bob.log"); c=$(grep -h "^bench: net quests" "$tmp/carol.log")
team() { echo "$1" | grep -o "team [0-9]*" | awk '{print $2}'; }
ber() { echo "$1" | grep -o "berries [0-9-]* -> [0-9-]*"; }
chk "$([ -n "$(team "$a")" ] && [ "$(team "$a")" = "$(team "$b")" ] && [ "$(team "$a")" != "$(team "$c")" ] && echo 1)" "alice and bob share a team, carol has her own ($(team "$a") $(team "$b") $(team "$c"))"
chk "$(echo "$a" | grep -q "e0f done 1" && echo 1)" "alice ate once: her team's 'eat something' node is done on the server"
chk "$(ber "$a" | awk '{exit !($4 == $2 + 3)}' && echo 1)" "and its 4 berries were paid to her by deltas ($(ber "$a"): one eaten, four paid)"
chk "$(echo "$b" | grep -q "e0f done 1" && echo 1)" "bob (the same team) sees the node done in the book the server sends"
chk "$(ber "$b" | awk '{exit !($4 == $2)}' && echo 1)" "the reward went to the completing player only (bob's $(ber "$b"))"
chk "$(echo "$c" | grep -q "e0f done 0" && echo 1)" "carol (another team) doesn't have it"
chk "$(echo "$b" | grep -q "quest news [1-9]" && echo 1)" "the team heard the news"
# slice 2: crafting on the server
cr=$(grep -h "^bench: net crafting (eater" "$tmp/alice.log"); echo "$cr"
chk "$(echo "$cr" | grep -q "logs 2 -> 0, planks 0 -> 8" && echo 1)" "alice crafted planks twice on the server: 2 logs -> 8 planks, by deltas"
chk "$(grep -q "chat: Your team hasn.t unlocked ore_scanner" "$tmp/alice.log" && echo 1)" "a recipe her team hasn't unlocked (ore scanner) was refused"
echo "=== test_team_quests: $([ $ok = 1 ] && echo PASS || echo FAIL) after ${SECONDS} s ($tmp)"
[ $ok = 1 ]
