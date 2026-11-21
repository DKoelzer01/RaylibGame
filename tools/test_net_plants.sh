#!/bin/bash
# Plants and farming on the server (#335). A peaceful --host (crops sown ripe: RAYLIBGAME_FARM_INSTANT) gives every joiner a
# stone shovel and 2 grain seeds. "feller" fells the nearest bush or tree by request, then tills ground ahead with the
# shovel (farmland placed by the server, free), sows grain there and harvests it once the server's fields show it ripe;
# "watcher" only looks. Checks: the drops reached the feller by deltas; the plant is gone on both clients (one removal,
# the same id); the watcher saw the crop come and go (FarmFields); the seeds left and the grain arrived by deltas.
# EXE: a GNS build (default build/gns/bin/game.exe). PORT (27968). TERRAIN=default: the default (erosion) terrain (#1107: plants
# near the spawn of a menu-hosted game). Prints "=== test_net_plants: PASS|FAIL". Git Bash.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
exe="${EXE:-$here/build/gns/bin/game.exe}"
port="${PORT:-27968}"
tmp="$(mktemp -d)"
cd "$here"
SECONDS=0
# TERRAIN=default (#1107): the default (erosion) terrain, as a menu-hosted game has; else the quick old terrain
if [ "${TERRAIN:-}" != default ]; then printf 'erosion 0\n' > "$tmp/host.sav.worldgen"; fi
. "$here/tools/lib_netmem.sh"
net_wait_for_memory 12 test_net_plants || exit 0
RAYLIBGAME_FARM_INSTANT=1 RAYLIBGAME_HOST_DIFFICULTY=0 RAYLIBGAME_HOST_GIVE="stone_shovel 1,grain_seeds 2" "$exe" --host "$tmp/host.sav" --seed 42 --bench $([ "${TERRAIN:-}" = default ] && echo 1200 || echo 600) --port "$port" > "$tmp/host.log" 2>&1 &
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
# (#1367) The clients regenerate the host's world: hand them its erosion bake, as test_join does (the cache is keyed on the
# bake's parameters, so a mismatch just bakes again). Without it each client baked for ~80 s, one after the other.
shared=0
if [ -f "$tmp/host.sav.erosion" ]; then
    mkdir -p saves/mp_cache
    for n in watcher feller; do cp "$tmp/host.sav.erosion" "saves/mp_cache/client_$n.sav.erosion"; done
    shared=1
fi
join() { RAYLIBGAME_PLANTS_SECONDS=${PLANTS_SECONDS:-40} RAYLIBGAME_NET_PLANTS_ROLE=$1 timeout ${CLIENT_TIMEOUT:-300} "$exe" --join "127.0.0.1:$port" --name "$1" --bench 3 $NET_LIGHT_CLIENT --test-net-plants > "$tmp/$1.log" 2>&1; }
if [ "${TERRAIN:-}" = default ] && [ "$shared" = 0 ]; then PLANTS_SECONDS=240 CLIENT_TIMEOUT=600 join watcher & wp=$!; else join watcher & wp=$!; fi
if [ "${TERRAIN:-}" = default ] && [ "$shared" = 0 ]; then
    # (#1107) without a shared bake each client bakes the erosion terrain itself: two at once with the host ran out of memory
    # (bad_alloc). The feller joins once the watcher is in the world ("world ready" prints only at the end; the watcher
    # listens longer).
    until grep -q "loading screen after a spawn" "$tmp/watcher.log" 2>/dev/null || ! kill -0 "$wp" 2>/dev/null || [ $SECONDS -gt 600 ]; do sleep 2; done
else
    sleep 2
fi
join feller
wait $wp
stop_host
ok=1
chk() { if [ "$1" = 1 ]; then echo "  ok: $2"; else echo "  MISSING: $2"; ok=0; fi; }
f=$(grep -h "^bench: net plants (feller" "$tmp/feller.log"); w=$(grep -h "^bench: net plants (watcher" "$tmp/watcher.log"); echo "$f"; echo "$w"
id=$(echo "$f" | sed -n 's/.*: plant \([0-9a-f]*\) standing.*/\1/p')
chk "$(echo "$f" | grep -q "felled 1: plant [0-9a-f]* standing 0" && [ "$id" != 0 ] && echo 1)" "the feller's plant came down (the server's removal reached it)"
drop=$(echo "$f" | grep -o "; [a-z_]* [0-9-]* -> [0-9-]* (+[0-9]* expected)")
chk "$(echo "$drop" | awk '{gsub(/[(+]/, "", $6); exit !($5 + 0 == $3 + $6 && $6 + 0 > 0)}' && echo 1)" "its drops arrived by deltas ($drop)"
chk "$(echo "$w" | grep -q "removals heard [1-9].*: $id\b" && echo "$w" | grep -q " $id" && echo 1)" "the watcher heard the same removal ($id)"
f2=$(grep -h "^bench: net farm (feller" "$tmp/feller.log"); w2=$(grep -h "^bench: net farm (watcher" "$tmp/watcher.log"); echo "$f2"; echo "$w2"
chk "$(echo "$f2" | grep -q "farmland 1; sown 1: seeds 2 -> [23];" && echo 1)" "the feller tilled (server farmland) and sowed (a seed out, 1-2 back with the harvest)"
chk "$(echo "$f2" | grep -q "harvested 1: grain 0 -> [1-9]" && echo 1)" "it harvested the ripe crop: grain arrived by deltas"
chk "$(echo "$w2" | grep -q "crop seen 1, gone 1" && echo 1)" "the watcher saw the crop grow and go (FarmFields)"
# (#920) XP counted on the server and sent as SkillState: the feller's felling and farming; the watcher, who did nothing, none
f3=$(grep -h "^bench: net skills (feller" "$tmp/feller.log"); w3=$(grep -h "^bench: net skills (watcher" "$tmp/watcher.log"); echo "$f3"; echo "$w3"
chk "$(echo "$f3" | awk -F'[;,]' '{ split($3, a, " "); split($4, b, " "); split($5, c, " "); exit !((a[2] + b[2]) > 0 && c[2] > 0) }' && echo 1)" "the feller earned felling and farming XP on the server"
chk "$(echo "$w3" | grep -q "woodcutting 0.0, foraging 0.0, farming 0.0" && echo 1)" "the watcher earned none"
# (#1107) why the server refused a plant, if it did (the host's own log of each refusal)
grep -h "refused for session" "$tmp/host.log" | head -5 | sed 's/^/  host: /'
echo "=== test_net_plants: $([ $ok = 1 ] && echo PASS || echo FAIL) after ${SECONDS} s ($tmp)"
[ $ok = 1 ]
