#!/bin/bash
# Print the GPU job's runner and time window for each run id given.
for r in "$@"; do
  gh api "repos/DKoelzer01/RaylibGame/actions/runs/$r/jobs" --jq ".jobs[]|select(.name|test(\"GPU\"))|\"$r \(.runner_name) \(.started_at[11:19])-\(.completed_at[11:19]) \(.conclusion)\""
done
