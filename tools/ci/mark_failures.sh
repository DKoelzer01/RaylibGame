#!/bin/bash
# The failures of a self-hosted run, each marked "known (#N)" from tools/ci/known_failures.txt or "NEW" (#1479), into the
# job summary (or stdout). The matching is the runner-check skill's (.claude/skills/runner-check/scripts/runner_check.sh):
# a "-> FAIL" line is known when it matches a line's extended regex (the text before the first '|').
#   tools/ci/mark_failures.sh [--verdict <rc dir>] <log or folder>...      (folders are searched for *.log)
# Without --verdict: exit 0 always (labels only).
# (#2120) With --verdict: the scenario steps only report (continue-on-error), and this decides the job. It exits 1 when
# a scenario failed with no known line matching its result line, including a scenario that died before printing one
# (the runner's "scenario: === NAME FAIL" with no "scenario: NAME -> " line), or when a scenario step's exit code
# (<rc dir>/rc-*.txt, "<log name> <code>") is non-zero while its log lists no failed scenario (the run itself broke).
root=$(cd "$(dirname "$0")/../.." && pwd)
# (#1946) this checkout's list plus staging's, so entries added there since the branch was cut count too
knownText=$(bash "$root/tools/ci/known_list.sh")
out="${GITHUB_STEP_SUMMARY:-/dev/stdout}"
rcDir=""
if [ "${1:-}" = --verdict ]; then rcDir="$2"; shift 2; fi
files=()
for a in "$@"; do
  if [ -d "$a" ]; then while IFS= read -r f; do files+=("$f"); done < <(find "$a" -name '*.log' -type f 2>/dev/null)
  elif [ -f "$a" ]; then files+=("$a"); fi
done
# The known line's issue ("#N") a line matches, or nothing.
knownIssue() {
  local pat issue note
  while IFS='|' read -r pat issue note; do
    pat=$(echo "$pat" | sed 's/[[:space:]]*$//'); [[ -z "$pat" || "$pat" == \#* ]] && continue
    if echo "$1" | grep -Eq "$pat"; then echo $issue; return; fi
  done <<< "$knownText"
}
new=0; old=0; body=""
if [ ${#files[@]} -gt 0 ]; then
  while IFS= read -r line; do
    issue=$(knownIssue "$line")
    if [ -z "$issue" ]; then tag="NEW"; new=$((new + 1)); else tag="known $issue"; old=$((old + 1)); fi
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
[ -z "$rcDir" ] && exit 0

# (#2120) The verdict, per failed scenario (the runner's "=== NAME FAIL" line), so a crash or a kill counts too.
bad=0; why=""
for f in "${files[@]}"; do
  log=$(tr -d '\r' < "$f")
  while IFS= read -r name; do
    [ -z "$name" ] && continue
    known=""
    while IFS= read -r l; do
      [ -z "$l" ] && continue
      known=$(knownIssue "$l"); [ -n "$known" ] && break
    done < <(grep -F -- "scenario: $name -> " <<< "$log" | grep -E -- '-> (FAIL|TIMEOUT|ERROR)')
    if [ -z "$known" ]; then bad=$((bad + 1)); why+="- \`$name\` ($(basename "$f"))"$'\n'; fi
  done < <(sed -nE 's/^scenario: === ([^ ]+) FAIL in .*/\1/p' <<< "$log" | sort -u)
done
for r in "$rcDir"/rc-*.txt; do
  [ -f "$r" ] || continue
  read -r logName rc < "$r"
  [ "${rc:-0}" = 0 ] && continue
  if [ "$rc" = started ]; then # the step's own timeout ended it: the scenarios after the last one listed never ran
    bad=$((bad + 1)); why+="- \`$logName\`: the step did not finish (its timeout ended it)"$'\n'; continue
  fi
  if ! tr -d '\r' < "$rcDir/$logName" 2>/dev/null | grep -qE '^scenario: === [^ ]+ FAIL in '; then
    bad=$((bad + 1)); why+="- \`$logName\` exited $rc with no failed scenario listed (the run broke)"$'\n'
  fi
done
{
  echo ""
  if [ $bad -eq 0 ]; then echo "**Verdict (#2120): green** - every scenario failure is known."
  else echo "**Verdict (#2120): red** - $bad NEW:"; printf '%s' "$why"; fi
} >> "$out"
echo "mark_failures verdict: $bad NEW"
[ $bad -eq 0 ]
