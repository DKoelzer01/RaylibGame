#!/bin/bash
# Build the work tree and run the task's test (#2095). Runs inside the llm-bench container in /work. Prints a short,
# model-friendly result (errors or the test's tail) and appends one timing line per call to /run/rounds.tsv:
#   <label> <start> <end> <build_s> <exit>
# Usage: testcmd.sh [label]   (TEST_CMD and KIND come from the environment)
label=${1:-round}
OUT=/cache/out; JOBS=${JOBS:-4}
t0=$(date +%s.%N)
bash tools/build_linux.sh --inside --no-gns --out "$OUT" --jobs "$JOBS" > /run/build-last.log 2>&1
rc=$?
if [ $rc = 0 ] && [ "$KIND" = unit ]; then
  cmake --build "$OUT/build" -j "$JOBS" --target game_tests >> /run/build-last.log 2>&1; rc=$?
fi
t1=$(date +%s.%N)
if [ $rc != 0 ]; then
  echo "BUILD FAILED"
  grep -E "error|Error" /run/build-last.log | grep -v -E "^make|Error [0-9]" | head -25
  printf '%s\t%s\t%s\t%s\t%s\n' "$label" "$t0" "$(date +%s.%N)" "$(awk -v a="$t0" -v b="$t1" 'BEGIN{print b-a}')" build_failed >> /run/rounds.tsv
  exit 2
fi
export GAME="$OUT/bin/game" TESTS="$OUT/bin/game_tests"
if [ "$KIND" = scenario ]; then exec 9>> /gpulock/gpu.lock; flock -s -w 1200 9 || true; rm -rf /run/scenario-out; fi
timeout 600 bash -c "$TEST_CMD" > /run/test-last.log 2>&1
trc=$?
if [ "$KIND" = scenario ] && [ $trc = 0 ]; then
  # (#2137) a single-name scenario run exits 0 even when its verdict is FAIL: the verdict in results.json decides
  n=0; bad=0
  for f in $(find /run/scenario-out -name results.json 2>/dev/null); do
    n=$((n + 1)); [ "$(jq -r .verdict "$f")" = PASS ] || bad=$((bad + 1))
  done
  [ $n = 0 ] && { echo "no scenario results.json"; trc=3; }
  [ $bad -gt 0 ] && trc=1
fi
grep -E "^scenario: check .* FAIL|verdict" /run/test-last.log | head -10
tail -30 /run/test-last.log
echo "exit $trc"
printf '%s\t%s\t%s\t%s\t%s\n' "$label" "$t0" "$(date +%s.%N)" "$(awk -v a="$t0" -v b="$t1" 'BEGIN{print b-a}')" "$trc" >> /run/rounds.tsv
exit $trc
