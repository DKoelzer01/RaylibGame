#!/bin/bash
# Delete remote branches whose every commit is already in origin/staging.
# Usage: bash tools/coord/delmerged.sh <branch>...   (DRY=1 to only report)
git fetch -q origin staging
for b in "$@"; do
  if ! git fetch -q origin "$b" 2>/dev/null; then echo "$b: not on origin"; continue; fi
  n=$(git rev-list --count "origin/staging..origin/$b")
  if [ "$n" != "0" ]; then echo "$b: KEEP ($n commits not in staging)"; continue; fi
  if [ -n "$DRY" ]; then echo "$b: would delete"; continue; fi
  gh api -X DELETE "repos/DKoelzer01/RaylibGame/git/refs/heads/$b" >/dev/null 2>&1 && echo "$b: deleted" || echo "$b: delete failed"
done
