#!/bin/bash
# The failures of a self-hosted run, each marked "known (#N)" from tools/ci/known_failures.txt or "NEW" (#1479), into the
# job summary (or stdout). The matching is the runner-check skill's (.claude/skills/runner-check/scripts/runner_check.sh):
# a "-> FAIL" line is known when it matches a line's extended regex (the text before the first '|').
#   tools/ci/mark_failures.sh <log or folder>...      (folders are searched for *.log)
# Exit 0 always: the steps that ran the tests already carry the verdict; this only labels it.
root=$(cd "$(dirname "$0")/../.." && pwd)
# (#1946) this checkout's list plus staging's, so entries added there since the branch was cut count too
knownText=$(bash "$root/tools/ci/known_list.sh")
out="${GITHUB_STEP_SUMMARY:-/dev/stdout}"
files=()
for a in "$@"; do
  if [ -d "$a" ]; then while IFS= read -r f; do files+=("$f"); done < <(find "$a" -name '*.log' -type f 2>/dev/null)
  elif [ -f "$a" ]; then files+=("$a"); fi
done
new=0; old=0; body=""
if [ ${#files[@]} -gt 0 ]; then
  while IFS= read -r line; do
    tag="NEW"
    while IFS='|' read -r pat issue note; do
      pat=$(echo "$pat" | sed 's/[[:space:]]*$//'); [[ -z "$pat" || "$pat" == \#* ]] && continue
      if echo "$line" | grep -Eq "$pat"; then tag="known $(echo $issue)"; break; fi
    done <<< "$knownText"
    if [ "$tag" = NEW ]; then new=$((new + 1)); else old=$((old + 1)); fi
    body+="- **[$tag]** \`$line\`"$'\n'
  done < <(cat "${files[@]}" | tr -d '\r' | grep -E -- '-> FAIL' | grep -v '=== ' | cut -c1-220 | sort -u)
fi
{
  echo ""
  echo "### Failures: $new NEW, $old known (#1479)"
  if [ $((new + old)) -eq 0 ]; then echo "_None in the logs._"; else printf '%s' "$body"; fi
  [ $new -gt 0 ] && echo "" && echo "NEW failures need a diagnosis, or an issue and a line in \`tools/ci/known_failures.txt\`."
} >> "$out"
echo "mark_failures: $new NEW, $old known"
exit 0
