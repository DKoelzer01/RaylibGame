#!/bin/bash
# The self-hosted jobs' paths (#1016, #1028, #1410): appends OUT, RUN, ccache and the version to $GITHUB_ENV.
# Outside the checkout (actions/checkout cleans it): GNS, FetchContent, ccache and the run folder (the erosion bake
# cache in saves/) stay warm between runs. Needs REF (the branch) in the environment.
#
# (#1410) Several runners take jobs at once and may share this volume: each job claims a slot, ci-persist/slotN with its
# own out/ (the build and the binary) and run/ (saves, the bake, scenario and join scratch), for the whole job - a
# background flock holds slotN.lock until the job ends (the runner ends the job's leftover processes; an ephemeral
# runner's container goes with them). Two jobs never share a slot; a slot stays warm for the next job that takes it.
# ccache stays shared (it is safe under concurrency).
set -eu
root="$RUNNER_WORKSPACE/ci-persist"
mkdir -p "$root/ccache"
echo "runner: ${RUNNER_NAME:-?}; ci-persist: $(stat -c 'device %d inode %i' "$root" 2>/dev/null || echo '?'); $(df -h "$root" 2>/dev/null | tail -1)"
slot="" holder=""
for i in 1 2 3 4 5 6; do
    mkdir -p "$root/slot$i"
    # A holder that keeps the slot's lock for the job (6 h at most); -n: a slot another job holds fails at once
    nohup flock -n "$root/slot$i.lock" sleep 21600 > /dev/null 2>&1 &
    holder=$!
    sleep 1
    if kill -0 "$holder" 2> /dev/null; then slot=$i; break; fi
done
if [ -z "$slot" ]; then echo "::error::no free ci-persist slot (6 in use): more concurrent jobs than slots"; exit 1; fi
dir="$root/slot$slot"
mkdir -p "$dir/out" "$dir/run"
echo "ci-persist slot $slot ($dir), held by pid $holder"
branch="${REF#refs/heads/}"
{
    echo "OUT=$dir/out"
    echo "RUN=$dir/run"
    echo "CI_SLOT=$slot"
    echo "CCACHE_DIR=$root/ccache"
    echo "CCACHE_MAXSIZE=2G"
    echo "RAYLIBGAME_VERSION=$branch $(git rev-parse --short HEAD)"
} >> "$GITHUB_ENV"
if [ -n "${GITHUB_STEP_SUMMARY:-}" ]; then echo "runner ${RUNNER_NAME:-?}, ci-persist slot $slot" >> "$GITHUB_STEP_SUMMARY"; fi
if command -v ccache > /dev/null; then
    # CMake (3.17+) takes the launcher from the environment: tools/build_linux.sh needs no change.
    { echo "CMAKE_C_COMPILER_LAUNCHER=ccache"; echo "CMAKE_CXX_COMPILER_LAUNCHER=ccache"; } >> "$GITHUB_ENV"
fi
