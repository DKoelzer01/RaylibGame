#!/bin/bash
# busybox (w64devkit on PATH) is not bash: re-run under Git's bash (#1477)
if [ -z "${BASH_VERSINFO:-}" ]; then for b in "C:/Program Files/Git/bin/bash.exe" "/c/Program Files/Git/bin/bash.exe" "/usr/bin/bash"; do [ -x "$b" ] && exec "$b" "$0" "$@"; done; echo "needs bash" >&2; exit 2; fi
# Runs (or finds) the self-hosted workflow for a branch, waits for it with a deadline, and lists its failures,
# each marked "known (#N)" from tools/ci/known_failures.txt or "NEW" (#1482).
#   runner_check.sh <branch> [--dispatch] [--regression]     or     runner_check.sh --run <run id>
R=repos/DKoelzer01/RaylibGame
root=$(git rev-parse --show-toplevel)
known="$root/tools/ci/known_failures.txt"
if [ "$1" = "--run" ]; then id=$2
else
  br=$1; shift
  if [ "$1" = "--dispatch" ] || [ "$2" = "--dispatch" ]; then
    extra=(); [[ " $* " == *" --regression "* ]] && extra=(-f regression=headless)
    gh workflow run selfhosted.yml -r "$br" "${extra[@]}" >/dev/null || exit 1
    sleep 8
  fi
  id=$(gh api "$R/actions/workflows/selfhosted.yml/runs?branch=$br&per_page=1" --jq '.workflow_runs[0].id')
  [ -n "$id" ] || { echo "runner-check: no self-hosted run for $br (use --dispatch)"; exit 1; }
fi
echo "run $id: https://github.com/DKoelzer01/RaylibGame/actions/runs/$id"
end=$((SECONDS+3600))
until [ "$(gh api $R/actions/runs/$id --jq .status)" = completed ] || [ $SECONDS -gt $end ]; do sleep 60; done
concl=$(gh api $R/actions/runs/$id --jq .conclusion)
echo "conclusion: $concl"
gh api "$R/actions/runs/$id/jobs" --jq '.jobs[] | "  job \(.name): \(.conclusion)"'
[ "$concl" = success ] && exit 0
new=0
while IFS= read -r line; do
  tag="NEW"
  while IFS='|' read -r pat issue note; do
    pat=$(echo "$pat" | sed 's/[[:space:]]*$//'); [[ -z "$pat" || "$pat" == \#* ]] && continue
    if echo "$line" | grep -Eq "$pat"; then tag="known $(echo $issue)"; break; fi
  done < "$known"
  [ "$tag" = NEW ] && new=$((new+1))
  echo "  [$tag] $line"
done < <(gh run view "$id" --log 2>/dev/null | grep -E -- '-> FAIL' | grep -v '^.*=== ' | sed 's/.*Z //' | cut -c1-220 | sort -u)
echo "new failures: $new"
[ $new -eq 0 ] && exit 3 || exit 1   # 3 = red, but only known failures
