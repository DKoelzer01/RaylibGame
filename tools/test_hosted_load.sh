#!/bin/bash
# A hosted save loaded in single player (#901 part 2, the user: "a multiplayer save kept inventory when loaded back in single
# player, but the player's position was different"). As the menu's Host does: a --host server, one player joins (fixed
# RAYLIBGAME_PLAYER_UUID), walks and leaves (the server writes her record), the server is stopped with "stop" (its final
# save); then single player loads that save as the same player (--test-load-hold). Compares where she left (the client's
# body-local position), what the server's record says (player_body) and where single player puts her.
# EXE: a GNS build (default build/gns/bin/game.exe). PORT (27975). Uses this worktree's saves/world.sav for the single-player
# load (any one there is moved aside and put back). Prints "=== test_hosted_load: PASS|FAIL". Capped at 15 minutes.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
exe="${EXE:-$here/build/gns/bin/game.exe}"
port="${PORT:-27975}"
uuid="0a9b1c2d-3e4f-4a5b-8c6d-7e8f9a0b1c2d"
tmp="$(mktemp -d)"
cd "$here"
SECONDS=0
# The host's world: seed 42's erosion bake from the benches when there is one (keyed on its parameters: else it bakes)
[ -f saves/bench.sav.erosion ] && cp saves/bench.sav.erosion "$tmp/host.sav.erosion"
# (the record and the world file are written as she leaves: the host is then stopped)
timeout 700 "$exe" --host "$tmp/host.sav" --seed 42 --bench 900 --port "$port" > "$tmp/host.log" 2>&1 &
hp=$!
until grep -q "world ready" "$tmp/host.log" 2>/dev/null || ! kill -0 "$hp" 2>/dev/null || [ $SECONDS -gt 300 ]; do sleep 1; done
mkdir -p saves/mp_cache
[ -f "$tmp/host.sav.erosion" ] && cp "$tmp/host.sav.erosion" saves/mp_cache/client_alice.sav.erosion
RAYLIBGAME_JOIN_LINGER=0 timeout 400 "$exe" --join "127.0.0.1:$port" --name alice --bench 3 --test-join > "$tmp/alice.log" 2>&1
sleep 3; kill "$hp" 2>/dev/null; wait "$hp" 2>/dev/null
left=$(grep "^join: alice leaves" "$tmp/alice.log" | head -1)
rec=$(ls "$tmp"/host.sav.players/*.player 2>/dev/null | head -1) # (a bench client's id is made from its name)
uuid=$(basename "$rec" .player)
recLine=$(grep "^player_body" "$rec" 2>/dev/null | head -1)
echo "client: ${left:-none}"
echo "record: ${recLine:-none}"
grep -h "host: saved on shutdown\|server: stopped" "$tmp/host.log" | head -2
# Single player: the hosted save as saves/world.sav (this worktree's), loaded as the same player
aside=""
if [ -e saves/world.sav ]; then aside="$tmp/aside"; mkdir -p "$aside"; for f in saves/world.sav*; do mv "$f" "$aside/"; done; fi
for f in "$tmp"/host.sav*; do cp -r "$f" "saves/world.sav${f#$tmp/host.sav}"; done
RAYLIBGAME_PLAYER_UUID="$uuid" timeout 300 "$exe" --seed 42 --bench 60 --test-load-hold > "$tmp/sp.log" 2>&1
grep -h "^bench: load-hold" "$tmp/sp.log" | grep -v "waiting\| t=" | tail -2
rm -rf saves/world.sav*
[ -n "$aside" ] && for f in "$aside"/*; do mv "$f" saves/; done
# The verdict: the record holds where she left (within 2 m: her last steps against the server's acknowledged state), and
# single player puts her there (the load-hold test's own PASS)
ok=1
cx=$(echo "$left" | awk '{print $8, $9, $10}'); rx=$(echo "$recLine" | awk '{print $3, $4, $5}') # ("join: alice leaves at body N local x y z")
d=$(awk -v a="$cx" -v b="$rx" 'BEGIN { split(a, p, " "); split(b, q, " "); if (length(p) < 3 || length(q) < 3) { print -1; exit }
     print sqrt((p[1]-q[1])^2 + (p[2]-q[2])^2 + (p[3]-q[3])^2) }')
echo "the record is $d m from where she left"
awk -v d="$d" 'BEGIN { exit !(d >= 0 && d < 2.0) }' || { echo "  FAIL: the server's record isn't where she left"; ok=0; }
grep -q "^bench: load-hold.*PASS" "$tmp/sp.log" || { echo "  FAIL: single player didn't keep her at the record's position"; ok=0; }
echo "=== test_hosted_load: $([ $ok = 1 ] && echo PASS || echo FAIL) after $SECONDS s ($tmp)"
