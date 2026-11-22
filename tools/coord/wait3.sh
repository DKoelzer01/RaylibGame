#!/bin/bash
# Wait for CI ok on PRs given as args (max 25 min), then print each result.
end=$((SECONDS+1500))
for p in "$@"; do
  sha=$(gh api repos/DKoelzer01/RaylibGame/pulls/$p --jq .head.sha)
  while [ $SECONDS -lt $end ]; do
    c=$(gh api "repos/DKoelzer01/RaylibGame/commits/$sha/check-runs?check_name=CI%20ok" --jq '[.check_runs[]|"\(.status)/\(.conclusion)"]|first // "none"')
    case "$c" in completed/*) break;; esac; sleep 60
  done
  echo "#$p ${sha:0:7}: $c"
done
