#!/bin/bash
# Merge ready PRs into m84-integration in order; retarget stacked ones first. Stops at the first conflict
# or at a PR whose "CI ok" check (#1016) hasn't passed on its head commit (NOCI=1 skips that gate).
R=repos/DKoelzer01/RaylibGame
for p in "$@"; do
  base=$(gh api $R/pulls/$p --jq .base.ref)
  if [ "${KEEPBASE:-0}" != "1" ] && [ "$base" != "m84-integration" ]; then
    gh api -X PATCH $R/pulls/$p -f base=m84-integration --silent && echo "#$p retargeted from $base"
    sleep 4
  fi
  if [ "${NOCI:-0}" != "1" ]; then
    sha=$(gh api $R/pulls/$p --jq .head.sha)
    ci=$(gh api "$R/commits/$sha/check-runs?check_name=CI%20ok" --jq '[.check_runs[] | "\(.status)/\(.conclusion)"] | first // "none"')
    case "$ci" in
      completed/success) ;;
      *) echo "#$p STOP: CI ok is $ci on ${sha:0:7} (wait for it, or NOCI=1 if the PR predates CI)"; exit 1;;
    esac
  fi
  st=""
  for i in $(seq 1 12); do
    st=$(gh api $R/pulls/$p --jq '"\(.mergeable)/\(.mergeable_state)"')
    case "$st" in true/*|false/*) break;; esac
    sleep 5
  done
  case "$st" in
    true/clean|true/unstable|true/blocked|true/behind)
      gh pr merge $p --merge --subject "Merge #$p ($(date +%Y-%m-%d))" >/dev/null 2>&1
      echo "#$p $(gh api $R/pulls/$p --jq '"merged=\(.merged) \(.merge_commit_sha[0:7])"')";;
    *) echo "#$p STOP: $st"; exit 1;;
  esac
done
echo "ALL DONE"
