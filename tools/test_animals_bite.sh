#!/bin/bash
# Animals on the server, slice 2 (#334): a hostile one bites a remote player through the server. A --host places a thornback
# (a daytime charger) 8 m from the spawn with threats on; a --join "victim" stands there. The thornback must hunt the victim
# (the server's nearest player, not the host's idle stand-in), its bites must reach the victim's client as Damage (its health
# drops), and the victim's feed request must be answered by the server (refused: it has nothing a thornback eats).
# EXE: a GNS build (default build/gns/bin/game.exe). PORT (27962). Prints "=== test_animals_bite: PASS|FAIL". Capped at 6 min.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
exe="${EXE:-$here/build/gns/bin/game.exe}"
port="${PORT:-27962}"
tmp="$(mktemp -d)"
cd "$here"
SECONDS=0
printf 'erosion 0\n' > "$tmp/host.sav.worldgen"
RAYLIBGAME_HOST_HUNGER=50 RAYLIBGAME_HOST_GIVE="berries 3,bed 1" RAYLIBGAME_HOST_HERD="thornback 1" RAYLIBGAME_HOST_HERD_THREATS=1 "$exe" --host "$tmp/host.sav" --seed 42 --day --bench 600 --port "$port" > "$tmp/host.log" 2>&1 &
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
[ -d /proc ] || echo "test_animals_bite: warning: not Git Bash (no /proc): the host may be left running; run it with Git Bash" >&2
until grep -q "test herd" "$tmp/host.log" 2>/dev/null || ! kill -0 "$hp" 2>/dev/null || [ $SECONDS -gt 240 ]; do sleep 1; done
grep -h "test herd" "$tmp/host.log" || { echo "test_animals_bite: no thornback on the host -> FAIL"; stop_host; exit 1; }
RAYLIBGAME_VICTIM_SECONDS=40 RAYLIBGAME_NET_ANIMALS_ROLE=victim timeout 200 "$exe" --join "127.0.0.1:$port" --name victim --bench 3 --test-net-animals > "$tmp/victim.log" 2>&1
stop_host
line=$(grep "^bench: net animals (victim)" "$tmp/victim.log")
echo "$line"
ok=1
bites=$(echo "$line" | grep -o "bitten [0-9]*" | awk '{print $2}')
h0=$(echo "$line" | grep -o "health [0-9]* -> [0-9]*" | awk '{print $2}'); h1=$(echo "$line" | grep -o "health [0-9]* -> [0-9]*" | awk '{print $4}')
chk() { if [ "$1" = 1 ]; then echo "  ok: $2"; else echo "  MISSING: $2"; ok=0; fi; }
chk "$([ -n "$bites" ] && [ "$bites" -gt 0 ] && echo 1)" "the thornback hunted the remote player and its bites reached that client ($bites)"
chk "$(grep -q "chat: It won't eat" "$tmp/victim.log" && echo 1)" "the feed request was answered by the server (refused: it has nothing a thornback eats)"
# #337: survival on the server
sv=$(grep "^bench: net survival (victim" "$tmp/victim.log"); echo "$sv"
num() { echo "$sv" | grep -o "$1 [-0-9.]* -> [-0-9.]*"; }
f0=$(num food | awk '{print $2}'); f1=$(num food | awk '{print $4}'); b0=$(num berries | awk '{print $2}'); b1=$(num berries | awk '{print $4}')
chk "$(echo "$sv" | grep -q "[1-9][0-9]* vitals heard" && echo 1)" "the server's Vitals reached the client"
chk "$([ -n "$f1" ] && awk -v a="$f0" -v b="$f1" 'BEGIN{exit !(b > a)}' && echo 1)" "eating: the server's food went up ($f0 -> $f1)"
chk "$([ -n "$b1" ] && [ "$b1" = "$((b0 - 1))" ] && echo 1)" "eating: one berry gone, by the server's delta ($b0 -> $b1)"
chk "$(echo "$sv" | grep -q "died 1" && echo 1)" "the bites killed it on the server (death in Vitals)"
chk "$(echo "$sv" | grep -q "respawned 1" && echo 1)" "and it respawned with full health"
chk "$(echo "$sv" | grep -q "bitten while dead 0$" && echo 1)" "nothing bit it while it was dead"
# #800: a bed placed through the server, set as the spawn (BedUse); the respawn wakes on it
bd=$(grep "^bench: net bed (victim" "$tmp/victim.log"); echo "$bd"
wake=$(echo "$bd" | grep -o "woke [-0-9.]*" | awk '{print $2}')
chk "$(echo "$bd" | grep -q "placed 1, set 1" && echo 1)" "the bed was placed through the server and set as the spawn"
chk "$(grep -q "chat: Spawn set" "$tmp/victim.log" && echo 1)" "the server answered: spawn set"
chk "$([ -n "$wake" ] && awk -v w="$wake" 'BEGIN{exit !(w >= 0 && w < 1.5)}' && echo 1)" "the respawn woke it on its bed ($wake m)"
echo "=== test_animals_bite: $([ $ok = 1 ] && echo PASS || echo FAIL) after ${SECONDS} s ($tmp)"
[ $ok = 1 ]
