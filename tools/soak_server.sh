#!/bin/bash
# Server soak (#348): a --host process and a bot-swarm process (--join --bots N --soak S) on loopback, both with simulated
# lag and loss, so each round trip is ~2 x LAG_MS. The host measures while all N bots (+ the swarm's standing player) are
# online: frame work time (mean/p99/max), the per-system CPU breakdown, traffic per client; this script adds the host's
# peak memory and counts timeouts.
#   tools/soak_server.sh N [SECONDS]      (default 120 s)   env: EXE (GNS build), PORT, LAG_MS (50), LOSS (2), SAVE_DIR,
#   INTEREST (1; 0 = every state to every client, #340 A/B)
# SAVE_DIR (default /tmp/raylibgame_soak) keeps a pristine host save (made once: later runs skip the erosion bake); each run
# starts from a fresh copy of it, so no bot resumes where an earlier run left it.
# Timing numbers are only meaningful in a quiet window (no other benches or builds).
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
n="${1:-8}"; secs="${2:-120}"
exe="${EXE:-$here/build/gns/bin/game.exe}"
port="${PORT:-27981}"; lag="${LAG_MS:-50}"; loss="${LOSS:-2}"
dir="${SAVE_DIR:-/tmp/raylibgame_soak}"; mkdir -p "$dir"
win="$(cygpath -m "$dir" 2>/dev/null || echo "$dir")"
cd "$here"
SECONDS=0
# Fresh sessions every run (#768: a reused save kept each bot's position, so later runs walked into unbuilt terrain and their
# streaming cost wasn't comparable). The pristine copy is the host save as created, before anyone joined (with its erosion
# bake cache), made once; each run starts from a copy of it.
if [ ! -f "$dir/pristine/host.sav" ]; then
    rm -rf "$dir/pristine" "$dir/run"; mkdir -p "$dir/pristine"
    pw="$(cygpath -m "$dir/pristine" 2>/dev/null || echo "$dir/pristine")"
    "$exe" --host "$pw/host.sav" --seed 42 --port "$port" > "$dir/pristine/make.log" 2>&1 &
    mp=$!
    until grep -qi "world ready" "$dir/pristine/make.log" 2>/dev/null || ! kill -0 "$mp" 2>/dev/null || [ $SECONDS -gt 300 ]; do sleep 2; done
    kill "$mp" 2>/dev/null; wait "$mp" 2>/dev/null
    echo "soak: made the pristine host save in ${SECONDS} s"
fi
rm -rf "$dir/run"; cp -r "$dir/pristine" "$dir/run"
win="$(cygpath -m "$dir/run" 2>/dev/null || echo "$dir/run")"
RAYLIBGAME_NET_INTEREST="${INTEREST:-1}" RAYLIBGAME_SOAK_PLAYERS=$((n + 1)) RAYLIBGAME_SOAK_SECONDS="$secs" RAYLIBGAME_NET_LAG_MS="$lag" RAYLIBGAME_NET_LOSS="$loss" RAYLIBGAME_NET_SEED=1 \
    "$exe" --host "$win/host.sav" --seed 42 --port "$port" --max-players $((n + 1)) > "$dir/host_${n}_i${INTEREST:-1}.log" 2>&1 &
hp=$!
# Up: listening and the world streamed in ("World ready in ...").
until grep -qi "world ready" "$dir/host_${n}_i${INTEREST:-1}.log" 2>/dev/null || ! kill -0 "$hp" 2>/dev/null || [ $SECONDS -gt 300 ]; do sleep 2; done
echo "soak: host up after ${SECONDS} s; $n bots for ${secs} s, ${lag} ms lag each way, ${loss}% loss, interest ${INTEREST:-1}"
RAYLIBGAME_NET_LAG_MS="$lag" RAYLIBGAME_NET_LOSS="$loss" RAYLIBGAME_NET_SEED=2 \
    timeout $((secs + 400)) "$exe" --join "127.0.0.1:$port" --name observer --bots "$n" --soak "$secs" > "$dir/swarm_${n}_i${INTEREST:-1}.log" 2>&1 &
sp=$!
# The host's peak memory, sampled while the soak runs (Windows pid behind the MSYS one).
wp=$(cat /proc/$hp/winpid 2>/dev/null || echo "")
until grep -q "bench: server soak" "$dir/host_${n}_i${INTEREST:-1}.log" 2>/dev/null || ! kill -0 "$sp" 2>/dev/null || [ $SECONDS -gt $((secs + 600)) ]; do sleep 3; done
mem=""; [ -n "$wp" ] && mem=$(powershell -NoProfile -Command "\$p = Get-Process -Id $wp -ErrorAction SilentlyContinue; if (\$p) { '{0:N0} MB working set, {1:N0} MB peak, {2:N0} MB private' -f (\$p.WorkingSet64/1MB), (\$p.PeakWorkingSet64/1MB), (\$p.PrivateMemorySize64/1MB) }" 2>/dev/null)
wait "$sp" 2>/dev/null
sleep 2; kill "$hp" 2>/dev/null; wait "$hp" 2>/dev/null
grep -h "^bench: server soak\|^bench: cpu per frame" "$dir/host_${n}_i${INTEREST:-1}.log"
grep -h "^bench: mp-soak" "$dir/swarm_${n}_i${INTEREST:-1}.log"
echo "bench: server memory with $n bots: ${mem:-unknown}"
echo "bench: timeouts: $(grep -c "timed out" "$dir/host_${n}_i${INTEREST:-1}.log") on the host (logs: $dir)"
