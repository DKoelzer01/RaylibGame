#!/bin/bash
# The #348 timing windows (run in a quiet machine window: no other benches or builds). Each run: a host and a bot swarm,
# 50 ms lag each way (~100 ms round trip), 2% loss, 100 s of measurement; interest management on and off (#340 A/B).
#   tools/soak_windows.sh A   -> 8 and 16 bots (4 runs, ~10 min)
#   tools/soak_windows.sh B   -> 32 bots (2 runs, ~6 min)
# Results go to stdout and SAVE_DIR (default /tmp/raylibgame_soak)/window_<A|B>.txt. Hard cap: 12 minutes.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
w="${1:-A}"
dir="${SAVE_DIR:-/tmp/raylibgame_soak}"; mkdir -p "$dir"
out="$dir/window_$w.txt"; : > "$out"
case "$w" in A) runs="8:1 8:0 16:1 16:0" ;; B) runs="32:1 32:0" ;; *) echo "usage: $0 A|B"; exit 2 ;; esac
timeout 720 bash -c '
for r in '"$runs"'; do
    n=${r%%:*}; i=${r##*:}
    echo "=== $n bots, interest $i" | tee -a "'"$out"'"
    INTEREST=$i bash "'"$here"'/tools/soak_server.sh" "$n" 100 2>&1 | grep -E "^bench:" | tee -a "'"$out"'"
done'
echo "=== window $w done (results in $out)"
