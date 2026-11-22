#!/bin/bash
# Wait (max 60 min) for a workflow run to complete; print its conclusion and failed jobs.
R=repos/DKoelzer01/RaylibGame; id=$1; end=$((SECONDS+3600))
until [ "$(gh api $R/actions/runs/$id --jq .status)" = completed ] || [ $SECONDS -gt $end ]; do sleep 60; done
gh api $R/actions/runs/$id --jq '"run \(.id): \(.status)/\(.conclusion)"'
gh api "$R/actions/runs/$id/jobs" --jq '.jobs[] | select(.conclusion!="success") | "  job \(.name): \(.conclusion)"'
# Exit non-zero unless the run succeeded, so `wait-run.sh X && merge` stops on a red run (it merged #1453 on a red run once).
[ "$(gh api $R/actions/runs/$id --jq .conclusion)" = success ]
