#!/bin/bash
# (#1447) One GPU (RTX 2080 SUPER) serves every runner's job: with three GL contexts time-slicing it, a trivial
# full-screen pass measured 20 ms and the fluid readback 56 ms a tick, at full clocks and ~50% utilisation. Every GPU
# step runs under this lock shared (they still share the GPU with each other); a timing step takes it exclusive and runs
# alone. A waiting exclusive step closes the gate (gpu.want) so new shared steps queue behind it instead of starving it.
#
#   tools/ci/gpu_lock.sh shared|exclusive [WAIT_S] -- command args...
#
# WAIT_S (default 3600): if the lock isn't free by then, the command runs without it, and the step summary says so.
# The command inherits the lock's descriptor: it is held until the command and its children end.
set -u
mode=${1:-}; shift
wait=3600
if [ "${1:-}" != "--" ]; then wait=$1; shift; fi
[ "${1:-}" = "--" ] && shift
case "$mode" in shared) flag=-s ;; exclusive) flag=-x ;; *) echo "usage: gpu_lock.sh shared|exclusive [WAIT_S] -- cmd..." >&2; exit 2 ;; esac
dir="${RUNNER_WORKSPACE:-/tmp}/ci-persist"
mkdir -p "$dir"
exec 8>> "$dir/gpu.want" 9>> "$dir/gpu.lock"
t0=$SECONDS
taken=1
if flock $flag -w "$wait" 8; then # the gate: shared passes unless an exclusive step is waiting or running
    left=$((wait - (SECONDS - t0))); [ $left -lt 1 ] && left=1
    flock $flag -w "$left" 9 || taken=0
else
    taken=0
fi
[ "$mode" = shared ] && flock -u 8 # (a shared step only passes the gate; an exclusive one keeps it closed while it runs)
if [ $taken = 1 ]; then msg="gpu lock: $mode, after waiting $((SECONDS - t0)) s"; else msg="gpu lock: $mode NOT taken in $wait s: running without it"; fi
echo "$msg ($(date -u +%H:%M:%S))" >&2
[ -n "${GITHUB_STEP_SUMMARY:-}" ] && echo "$msg ($(date -u +%H:%M:%S))" >> "$GITHUB_STEP_SUMMARY"
"$@"
