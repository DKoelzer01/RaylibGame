#!/bin/bash
# Local-model benchmark: one task, one model (#2095). Runs on deweyserver as claude.
#   run.sh <task-id> <label>             the agent implements the task; label names the model/config (e.g. qwen3.5-9b)
#   run.sh <task-id> validate            no model: base + tests must FAIL, + the reference must PASS (checks the task)
# One job at a time (flock); builds at low priority next to CI (4 cpus, cpu-shares 512 = the runners').
# Result: one JSON row appended to results.jsonl; the run's files in runs/<time>-<task>-<label>/.
set -uo pipefail
B=$HOME/llm/bench; tid=${1:?task}; label=${2:?label}
T=$B/tasks/$tid; [ -f "$T/task.json" ] || { echo "no task $tid"; exit 1; }
API=${LLM_API:-http://host.docker.internal:8095/v1}
exec 7> "$B/slot1.lock"; flock -w 3600 7 || { echo "bench busy"; exit 3; }
REPO=$B/repo; WORK=$B/work1; OUT=$B/out1
[ -d "$REPO/.git" ] || { git clone -q "$HOME/rg-sandbox/repo" "$REPO" && git -C "$REPO" remote set-url origin https://github.com/DKoelzer01/RaylibGame.git; }
base=$(jq -r .base "$T/task.json"); kind=$(jq -r .kind "$T/task.json"); test_cmd=$(jq -r .test_cmd "$T/task.json")
git -C "$REPO" cat-file -e "$base^{commit}" 2>/dev/null || git -C "$REPO" fetch -q origin
# a standalone clone with its own objects (hardlinked): a worktree's .git and a --shared clone's alternates point
# outside the container's mount, and Aider reads the repo through GitPython
[ -d "$WORK/.git" ] && [ -f "$WORK/.git/objects/info/alternates" ] && rm -rf "$WORK"
[ -d "$WORK/.git" ] || git clone -q "$REPO" "$WORK"; git -C "$WORK" fetch -q "$REPO" "$base"
mkdir -p "$OUT"
cd "$WORK" && git checkout -q -f --detach "$base" && git clean -qfdx && git apply --whitespace=nowarn "$T/test.patch" \
  && git add -A && git -c user.name=bench -c user.email=bench@local commit -qm "bench: task tests" || { echo "setup failed"; exit 1; }
start=$(git rev-parse HEAD)
RUN=$B/runs/$(date +%Y%m%d-%H%M%S)-$tid-$label; mkdir -p "$RUN"; : > "$RUN/rounds.tsv"
metric() { curl -s localhost:8095/metrics | awk -v k="llamacpp:$1" '$1==k{print $2}'; }
m0p=$(metric prompt_tokens_total); m0g=$(metric tokens_predicted_total); m0ps=$(metric prompt_seconds_total); m0gs=$(metric tokens_predicted_seconds_total)
t0=$(date +%s)
if [ "$label" = validate ]; then
  inner='bash /bench/testcmd.sh base > /run/base.out; echo "base exit $?"; git apply --whitespace=nowarn /task/ref.patch && bash /bench/testcmd.sh ref > /run/ref.out; echo "ref exit $?"'
else
  # AGENT=mini (default: bin/mini_agent.py) or aider (bin/aider_agent.sh); CTX_CHARS = the mini agent's conversation budget
  if [ "${AGENT:-mini}" = aider ]; then agentcmd='timeout 2700 bash /bench/aider_agent.sh > /run/agent.out 2>&1'
  else agentcmd='timeout 2700 python3 /bench/mini_agent.py --api "$API" --model "$LABEL" --ctx-chars "${CTX_CHARS:-80000}" > /run/agent.out 2>&1'; fi
  inner="$agentcmd"'; echo "agent exit $?"
         bash /bench/testcmd.sh final > /run/final.out; echo "final exit $?"
         [ "$KIND" = unit ] || cmake --build /cache/out/build -j "$JOBS" --target game_tests > /run/hidden-build.log 2>&1
         timeout 300 /cache/out/bin/game_tests --gtest_brief=1 > /run/hidden.log 2>&1; echo "hidden exit $?"'
fi
docker run --rm --name "llm-bench-$$" --add-host host.docker.internal:host-gateway \
  --gpus device=GPU-26c7bdf7-5620-6610-af3f-b4b131814584 -e NVIDIA_DRIVER_CAPABILITIES=graphics,compute,utility \
  --cpus 4 --cpu-shares 512 --memory 12g --pids-limit 4096 --security-opt no-new-privileges --cap-drop ALL \
  -v "$WORK:/work" -v "$OUT:/cache/out" -v "$HOME/rg-sandbox/cache/ccache:/cache/ccache" -v "$T:/task:ro" -v "$RUN:/run" \
  -v "$B/bin:/bench:ro" -v "$HOME/runner-stack/cache-1:/gpulock" \
  -e CCACHE_DIR=/cache/ccache -e CCACHE_MAXSIZE=3G -e CMAKE_C_COMPILER_LAUNCHER=ccache -e CMAKE_CXX_COMPILER_LAUNCHER=ccache \
  -e RAYLIBGAME_FAST_BUILD=1 -e JOBS=4 -e KIND="$kind" -e TEST_CMD="$test_cmd" -e API="$API" -e LABEL="$label" -e CTX_CHARS="${CTX_CHARS:-80000}" -e AIDER_EDIT_FORMAT="${AIDER_EDIT_FORMAT:-diff}" \
  -w /work llm-bench:1 bash -c "$inner" > "$RUN/inner.log" 2>&1
wall=$(( $(date +%s) - t0 ))
cat "$RUN/inner.log"
# ---- measures (host side)
d() { awk -v a="$1" -v b="$2" 'BEGIN{printf "%.1f", b-a}'; }
m1p=$(metric prompt_tokens_total); m1g=$(metric tokens_predicted_total); m1ps=$(metric prompt_seconds_total); m1gs=$(metric tokens_predicted_seconds_total)
code_of() { grep -E "^$1 exit" "$RUN/inner.log" | awk '{print $3}'; }
allowed=$(jq -r '.files[]' "$T/task.json")
touched=$(git -C "$WORK" diff --name-only "$start")
scope_bad=$(comm -13 <(echo "$allowed" | sort) <(echo "$touched" | sed '/^$/d' | sort) | tr '\n' ' ')
lines=$(git -C "$WORK" diff --numstat "$start" | awk '{s+=$1+$2} END{print s+0}')
# added lines whose line ending differs from their file's (the repo is LF on Linux checkouts, CRLF on Windows ones);
# added debug output
crlf_bad=0
for f in $touched; do
  [ -f "$WORK/$f" ] || continue
  if git -C "$WORK" show "$start:$f" 2>/dev/null | head -50 | grep -q $'\r$'; then
    n=$(git -C "$WORK" diff "$start" -- "$f" | grep -E '^\+[^+]|^\+$' | grep -vc $'\r$')
  else
    n=$(git -C "$WORK" diff "$start" -- "$f" | grep -E '^\+[^+]|^\+$' | grep -c $'\r$')
  fi
  crlf_bad=$((crlf_bad + ${n:-0}))
done
debug=$(git -C "$WORK" diff "$start" | grep -cE '^\+.*(printf\(|std::cout|std::cerr|fprintf\(stderr)' || true)
warn=0; for f in $touched; do n=$(grep -c "$f.*warning:" "$RUN/build-last.log" 2>/dev/null); warn=$((warn + ${n:-0})); done
git -C "$WORK" diff "$start" > "$RUN/model.diff"
build_s=$(awk -F'\t' '{s+=$4} END{printf "%.0f", s}' "$RUN/rounds.tsv")
agent='{}'; [ -f "$RUN/agent.json" ] && agent=$(cat "$RUN/agent.json")
jq -nc --arg task "$tid" --arg label "$label" --arg run "$(basename "$RUN")" --argjson wall "$wall" \
  --arg final "$(code_of final)" --arg hidden "$(code_of hidden)" --arg basex "$(code_of base)" --arg refx "$(code_of ref)" \
  --arg scope "$scope_bad" --argjson lines "$lines" --argjson crlf "$crlf_bad" --argjson debug "$debug" --argjson warn "$warn" \
  --argjson build_s "${build_s:-0}" --argjson agent "$agent" \
  --arg pt "$(d "${m0p:-0}" "${m1p:-0}")" --arg gt "$(d "${m0g:-0}" "${m1g:-0}")" --arg ps "$(d "${m0ps:-0}" "${m1ps:-0}")" --arg gs "$(d "${m0gs:-0}" "${m1gs:-0}")" \
  --slurpfile tj "$T/task.json" \
  '{time: (now|todate), task: $task, label: $label, run: $run, kind: $tj[0].kind, difficulty: $tj[0].difficulty,
    pass: ($final == "0"), hidden_pass: ($hidden == "0"), base_exit: $basex, ref_exit: $refx,
    wall_s: $wall, build_test_s: $build_s, rounds: ($agent.tests // null), steps: ($agent.steps // null), finished: ($agent.finished // null),
    edit_errors: ($agent.edit_errors // null), out_of_scope_attempts: ($agent.out_of_scope // null),
    scope_violations: $scope, lines_changed: $lines, ref_lines: $tj[0].ref_lines_changed, crlf_bad_lines: $crlf, debug_lines: $debug, new_warnings: $warn,
    prompt_tokens: ($pt|tonumber), gen_tokens: ($gt|tonumber), prompt_s: ($ps|tonumber), gen_s: ($gs|tonumber),
    gen_tok_s: (if ($gs|tonumber) > 0 then (($gt|tonumber)/($gs|tonumber)*10|round/10) else null end),
    prompt_tok_s: (if ($ps|tonumber) > 0 then (($pt|tonumber)/($ps|tonumber)|round) else null end)}' | tee -a "$B/results.jsonl"
