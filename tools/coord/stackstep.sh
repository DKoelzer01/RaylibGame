#!/bin/bash
# Merge one stacked PR into staging: retarget, refresh its branch (starts CI), wait for CI ok, merge.
# Usage: bash tools/coord/stackstep.sh <PR>
p="$1"
base=$(gh api repos/DKoelzer01/RaylibGame/pulls/$p --jq .base.ref)
if [ "$base" != "staging" ]; then
  gh api -X PATCH repos/DKoelzer01/RaylibGame/pulls/$p -f base=staging --jq '"retargeted #\(.number) to \(.base.ref)"'
  sleep 5
fi
gh api -X PUT repos/DKoelzer01/RaylibGame/pulls/$p/update-branch --jq .message 2>&1 | tail -1
sleep 40
bash tools/coord/wait3.sh "$p" && bash tools/coord/merge-ready.sh "$p" 2>&1 | tail -2
