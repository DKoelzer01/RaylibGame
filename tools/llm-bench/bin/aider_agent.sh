#!/bin/bash
# Aider as the benchmark's agent (#2095): headless, one message (the spec), its own test loop (--auto-test runs
# testcmd.sh after each edit and feeds failures back, up to Aider's 3 reflections). Runs in /work inside llm-bench.
# Env: API (the llama.cpp /v1 endpoint), AIDER_EDIT_FORMAT (diff | whole | udiff; default diff).
export HOME=/tmp/aider-home; mkdir -p $HOME
git config --global user.name bench; git config --global user.email bench@local
files=$(jq -r '.files[]' /task/task.json)
{ cat /task/spec.md
  echo
  echo "Files you may edit: $(echo $files | tr ' ' ', ')."
  echo "The test is already in the tree; it is built and run automatically after your edits. Keep the existing style;"
  echo "no debug prints, no unrelated changes, don't change the test."
} > /run/aider-message.md
t0=$(date +%s)
aider --model openai/local --openai-api-base "$API" --openai-api-key x \
  --model-metadata-file /bench/aider-models.json --model-settings-file /bench/aider-settings.yml \
  --edit-format "${AIDER_EDIT_FORMAT:-diff}" --map-tokens 2048 \
  --message-file /run/aider-message.md --test-cmd "bash /bench/testcmd.sh aider" --auto-test \
  --yes-always --no-auto-commits --no-dirty-commits --no-gitignore --no-show-model-warnings --no-check-update \
  --analytics-disable --no-pretty --no-stream --no-suggest-shell-commands --no-detect-urls --timeout 600 \
  --chat-history-file /run/aider.history.md --input-history-file /run/aider.input --llm-history-file /run/aider.llm.log \
  $files
rc=$?
tests=$(grep -c '^aider' /run/rounds.tsv 2>/dev/null)
sent=$(grep -oE 'Tokens: [0-9.]+k? sent' /run/aider.history.md 2>/dev/null | wc -l)
jq -n --argjson tests "${tests:-0}" --argjson rc $rc --argjson s $(( $(date +%s) - t0 )) --argjson turns "${sent:-0}" \
  '{tests: $tests, steps: $turns, finished: ($rc == 0), stop: "aider exit \($rc)", agent_s: $s, edit_errors: null, out_of_scope: null}' > /run/agent.json
exit $rc
