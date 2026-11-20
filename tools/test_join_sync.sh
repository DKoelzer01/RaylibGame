#!/bin/bash
# Join-in-progress sync end to end (#776): a --host server whose world has many edits (--host-edits, 3000 dig strokes
# over 600 m around the start) and a --join client in a separate process, twice:
#   1. with nothing cached: it fetches every region, installs them, and its terrain fingerprint must equal the host's;
#   2. again under the same name: its cache (the client save's regions) matches, so it fetches only what changed since (at
#      most a tenth of the entries: the host's world moves on between the joins), and the
#      fingerprint still matches.
# Needs a GNS build: EXE (default build/gns/bin/game.exe). PORT (default 27961). EDITS (default 3000). LAG_MS (30).
# Prints the "join: synced" lines, the host's fingerprints and "=== test_join_sync: PASS|FAIL". Capped at 15 minutes.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
exe="${EXE:-$here/build/gns/bin/game.exe}"
port="${PORT:-27961}"
edits="${EDITS:-3000}"
lag="${LAG_MS:-30}"
name="syncer"
tmp="$(mktemp -d)"
cd "$here"
SECONDS=0
rm -rf "saves/mp_cache/client_$name.sav" "saves/mp_cache/client_$name.sav."* # (no cache to begin with)
RAYLIBGAME_SYNC_FINGERPRINT=1 "$exe" --host "$tmp/host.sav" --seed 42 --host-edits "$edits" --bench 900 --port "$port" > "$tmp/host.log" 2>&1 &
hp=$!
until grep -q "world ready" "$tmp/host.log" 2>/dev/null || ! kill -0 "$hp" 2>/dev/null || [ $SECONDS -gt 300 ]; do sleep 1; done
if ! grep -q "listening on port" "$tmp/host.log"; then echo "test_join_sync: the host never listened -> FAIL"; kill "$hp" 2>/dev/null; exit 1; fi
echo "test_join_sync: host up after ${SECONDS} s with $(grep -o '[0-9]* dig strokes' "$tmp/host.log") ($tmp)"
ok=1
for run in 1 2; do
    RAYLIBGAME_NET_LAG_MS="$lag" RAYLIBGAME_JOIN_OTHERS=0 timeout 400 "$exe" --join "127.0.0.1:$port" --name "$name" --bench 3 --test-join > "$tmp/client$run.log" 2>&1
    line="$(grep "^join: synced" "$tmp/client$run.log")"
    mine="$(echo "$line" | grep -o 'fingerprint [0-9a-f]*' | cut -d' ' -f2)"
    host="$(grep "terrain fingerprint" "$tmp/host.log" | tail -1 | grep -o '[0-9a-f]*$')"
    fetched="$(echo "$line" | grep -o 'synced: [0-9]* of [0-9]*' | cut -d' ' -f2)"
    echo "run $run: $(grep '^join: cache' "$tmp/client$run.log")"
    echo "run $run: $line"
    echo "run $run: host fingerprint $host"
    [ -n "$mine" ] && [ "$mine" = "$host" ] || { echo "run $run: fingerprints differ (or missing)"; ok=0; }
    # A rejoin fetches only what changed since: the host's world moves on between the joins (structures build near
    # players, animals), so a few entries may differ; with the cache it must be a small part, never everything.
    total="$(echo "$line" | grep -o 'synced: [0-9]* of [0-9]*' | cut -d' ' -f4)"
    if [ "$run" -ge 2 ] && [ $(( ${fetched:-999} * 10 )) -gt "${total:-0}" ]; then echo "run $run: fetched $fetched of $total entries with a cache of the same world (want at most a tenth)"; ok=0; fi
done
kill "$hp" 2>/dev/null; wait "$hp" 2>/dev/null
grep -h "join sync" "$tmp/host.log" | sed 's/^/host: /'
echo "=== test_join_sync: $([ $ok = 1 ] && echo PASS || echo FAIL) (${SECONDS} s; logs in $tmp)"
[ $ok = 1 ]
