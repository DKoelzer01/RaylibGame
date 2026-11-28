#!/bin/bash
# Does a task's test catch a known-bad diff? (#2095)  mutant.sh <task-id> <diff file>  -> "caught" when the test fails.
set -uo pipefail
B=$HOME/llm/bench; tid=$1; diff=$2; T=$B/tasks/$tid
exec 7> "$B/slot1.lock"; flock -w 7200 7 || { echo "bench busy"; exit 3; }
WORK=$B/work1; OUT=$B/out1; base=$(jq -r .base "$T/task.json")
git -C "$WORK" fetch -q "$B/repo" "$base"
(cd "$WORK" && git checkout -q -f --detach "$base" && git clean -qfdx && git apply --whitespace=nowarn "$T/test.patch" "$diff") || { echo "setup failed"; exit 1; }
RUN=$(mktemp -d); : > $RUN/rounds.tsv
docker run --rm --cpus 4 --cpu-shares 512 --memory 12g --security-opt no-new-privileges --cap-drop ALL \
  --gpus device=GPU-26c7bdf7-5620-6610-af3f-b4b131814584 -e NVIDIA_DRIVER_CAPABILITIES=graphics,compute,utility \
  -v "$WORK:/work" -v "$OUT:/cache/out" -v "$HOME/rg-sandbox/cache/ccache:/cache/ccache" -v "$RUN:/run" -v "$B/bin:/bench:ro" \
  -v "$HOME/runner-stack/cache-1:/gpulock" -e CCACHE_DIR=/cache/ccache -e CMAKE_C_COMPILER_LAUNCHER=ccache -e CMAKE_CXX_COMPILER_LAUNCHER=ccache \
  -e RAYLIBGAME_FAST_BUILD=1 -e JOBS=4 -e KIND="$(jq -r .kind $T/task.json)" -e TEST_CMD="$(jq -r .test_cmd $T/task.json)" \
  -w /work llm-bench:1 bash /bench/testcmd.sh mutant > $RUN/out.txt 2>&1
rc=$?
[ $rc = 0 ] && echo "$tid $(basename $diff): NOT caught (test passes)" || echo "$tid $(basename $diff): caught (exit $rc)"
grep -E "Failure|FAILED|Expected|BUILD FAILED" $RUN/out.txt | head -4
rm -rf "$RUN"
