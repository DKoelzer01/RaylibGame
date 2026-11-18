#!/bin/bash
# --join end to end (#341 core, #339): a --host server and two --join clients, three separate processes on loopback.
# Each client regenerates the server's world from the Welcome (fingerprint must match), starts at the server's state, walks
# 8 s predicting itself, follows the server's world clock, and sees the other client walking (--test-join in main.cpp).
# Needs a GNS build: EXE (default build/gns/bin/game.exe). PORT (default 27951). LAG_MS: simulated lag on the clients (30).
# Prints the clients' "bench: join" lines and "=== test_join: N PASS, M FAIL". Whole run capped at 15 minutes.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
exe="${EXE:-$here/build/gns/bin/game.exe}"
port="${PORT:-27951}"
lag="${LAG_MS:-30}"
tmp="$(mktemp -d)"
cd "$here"
SECONDS=0
# A non-default world (#757): WORLDGEN="ocean 0.4 mountains 1.6 erosion 0" (key value pairs) becomes the host save's sidecar
# before it generates; SANDBOX=1 hosts the flat test plane. The clients get both in the Welcome (the fingerprint must match).
if [ -n "${WORLDGEN:-}" ]; then
    printf '%s %s\n' ${WORLDGEN} > "$tmp/host.sav.worldgen"
    echo "test_join: host world-gen settings: $(tr '\n' ' ' < "$tmp/host.sav.worldgen")"
fi
sandbox=""; [ "${SANDBOX:-0}" = "1" ] && sandbox="--sandbox"
"$exe" --host "$tmp/host.sav" --seed 42 $sandbox --bench 600 --port "$port" > "$tmp/host.log" 2>&1 &
hp=$!
until grep -q "listening on port" "$tmp/host.log" 2>/dev/null || ! kill -0 "$hp" 2>/dev/null || [ $SECONDS -gt 300 ]; do sleep 1; done
if ! grep -q "listening on port" "$tmp/host.log"; then
    echo "test_join: the host never listened ($(tail -1 "$tmp/host.log")) -> FAIL"; kill "$hp" 2>/dev/null; exit 1
fi
until grep -q "world ready" "$tmp/host.log" 2>/dev/null || [ $SECONDS -gt 300 ]; do sleep 1; done
echo "test_join: host up after ${SECONDS} s ($tmp)"
cpids=()
for n in alice bob; do
    shot=""; [ "$n" = bob ] && shot="bench_join_bob.png" # bob turns to alice and takes a picture (in the working directory)
    RAYLIBGAME_JOIN_SHOT="$shot" RAYLIBGAME_NET_LAG_MS="$lag" RAYLIBGAME_JOIN_OTHERS=1 RAYLIBGAME_JOIN_LINGER=20 \
        timeout 600 "$exe" --join "127.0.0.1:$port" --name "$n" --bench 3 --test-join > "$tmp/$n.log" 2>&1 &
    cpids+=($!)
done
for c in "${cpids[@]}"; do wait "$c"; done
sleep 1
kill "$hp" 2>/dev/null; wait "$hp" 2>/dev/null
grep -h "^join:\|^bench: join" "$tmp/alice.log" "$tmp/bob.log"
grep -h "joined as session\|left" "$tmp/host.log" | sed 's/^/host: /'
pass=$(cat "$tmp/alice.log" "$tmp/bob.log" | grep "^bench: join" | grep -c "PASS"); fail=$((2 - pass))
echo "=== test_join: $pass PASS, $fail FAIL (logs in $tmp)"
[ "$fail" -eq 0 ]
