#!/bin/bash
# Hidden-suite baseline for a task (#2095): base + tests + the reference, the full game_tests; the names of the tests that
# fail there go to tasks/<id>/hidden_baseline.txt. A model run's hidden check counts only failures not in that list.
#   baseline.sh <task-id> ...
set -uo pipefail
B=$HOME/llm/bench; REPO=$B/repo; WORK=$B/work1; OUT=$B/out1
exec 7> "$B/slot1.lock"; flock -w 7200 7 || { echo "bench busy"; exit 3; }
for tid in "$@"; do
  T=$B/tasks/$tid; base=$(jq -r .base "$T/task.json")
  [ -d "$WORK/.git" ] || git clone -q --shared "$REPO" "$WORK"; git -C "$WORK" fetch -q "$REPO" "$base"
  (cd "$WORK" && git checkout -q -f --detach "$base" && git clean -qfdx && git apply --whitespace=nowarn "$T/test.patch" "$T/ref.patch") \
    || { echo "$tid: setup failed"; continue; }
  docker run --rm --name "llm-bench-base-$$" --cpus 4 --cpu-shares 512 --memory 12g --security-opt no-new-privileges --cap-drop ALL \
    -v "$WORK:/work" -v "$OUT:/cache/out" -v "$HOME/rg-sandbox/cache/ccache:/cache/ccache" \
    -e CCACHE_DIR=/cache/ccache -e CMAKE_C_COMPILER_LAUNCHER=ccache -e CMAKE_CXX_COMPILER_LAUNCHER=ccache -e RAYLIBGAME_FAST_BUILD=1 \
    -w /work llm-bench:1 bash -c 'bash tools/build_linux.sh --inside --no-gns --out /cache/out --jobs 4 > /dev/null 2>&1 \
      && cmake --build /cache/out/build -j 4 --target game_tests > /dev/null 2>&1 \
      && timeout 300 /cache/out/bin/game_tests --gtest_brief=1 2>&1 | grep -E "^\[  FAILED  \] [A-Za-z_0-9]+\.[A-Za-z_0-9]+" | awk "{print \$4}" | sort -u' \
    > "$T/hidden_baseline.txt"
  echo "$tid: $(wc -l < "$T/hidden_baseline.txt") failing at the reference: $(tr '\n' ' ' < "$T/hidden_baseline.txt")"
done
