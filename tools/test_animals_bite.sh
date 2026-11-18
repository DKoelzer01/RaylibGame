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
RAYLIBGAME_HOST_HERD="thornback 1" RAYLIBGAME_HOST_HERD_THREATS=1 "$exe" --host "$tmp/host.sav" --seed 42 --day --bench 600 --port "$port" > "$tmp/host.log" 2>&1 &
hp=$!
until grep -q "test herd" "$tmp/host.log" 2>/dev/null || ! kill -0 "$hp" 2>/dev/null || [ $SECONDS -gt 240 ]; do sleep 1; done
grep -h "test herd" "$tmp/host.log" || { echo "test_animals_bite: no thornback on the host -> FAIL"; kill "$hp" 2>/dev/null; exit 1; }
RAYLIBGAME_NET_ANIMALS_ROLE=victim timeout 200 "$exe" --join "127.0.0.1:$port" --name victim --bench 3 --test-net-animals > "$tmp/victim.log" 2>&1
kill "$hp" 2>/dev/null; wait "$hp" 2>/dev/null
line=$(grep "^bench: net animals (victim)" "$tmp/victim.log")
echo "$line"
ok=1
bites=$(echo "$line" | grep -o "bitten [0-9]*" | awk '{print $2}')
h0=$(echo "$line" | grep -o "health [0-9]* -> [0-9]*" | awk '{print $2}'); h1=$(echo "$line" | grep -o "health [0-9]* -> [0-9]*" | awk '{print $4}')
chk() { if [ "$1" = 1 ]; then echo "  ok: $2"; else echo "  MISSING: $2"; ok=0; fi; }
chk "$([ -n "$bites" ] && [ "$bites" -gt 0 ] && echo 1)" "the thornback hunted the remote player and its bites reached that client ($bites)"
chk "$([ -n "$h1" ] && [ "$h1" -lt "$h0" ] && echo 1)" "the victim's health went down ($h0 -> $h1)"
chk "$(grep -q "chat: It won't eat" "$tmp/victim.log" && echo 1)" "the feed request was answered by the server (refused: it has nothing a thornback eats)"
echo "=== test_animals_bite: $([ $ok = 1 ] && echo PASS || echo FAIL) after ${SECONDS} s ($tmp)"
[ $ok = 1 ]
