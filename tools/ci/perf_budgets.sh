#!/bin/bash
# (#1690) Judge the budget scenes against tools/ci/perf_budgets.txt.
#
#   tools/ci/perf_budgets.sh OUT_DIR...     one scenario-out dir per run (3 runs: the median of each metric is judged)
#
# Reads every results.json under the dirs (the scene from its "scenario" field; needs jq), prints a markdown table (scene,
# metric, value, limit, margin) to $GITHUB_STEP_SUMMARY and stdout, and decides per row:
#   limit 0          "not judged"  (a placeholder: reported only; so is a metric the scene did not report)
#   value <= limit   "ok"
#   up to 20% over   "flag"        (a warning, exit 0)
#   more than 20%    "over"        (exit 1), unless PERF_BUDGETS_ENFORCE is not 1 (the report-only week): then "would block", exit 0
#
# Rerun before failing (not done here): the workflow can rerun only the scene that is "over"
# (`game --scenario perf/<name> --scenario-no-record --scenario-out DIR4`) and pass DIR4 as one more dir; the median is then
# taken over all the dirs given.
# Env: PERF_BUDGETS_FILE (default perf_budgets.txt next to this script), PERF_BUDGETS_NOTE (a line printed above the table).
set -u
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
file="${PERF_BUDGETS_FILE:-$here/perf_budgets.txt}"
sum="${GITHUB_STEP_SUMMARY:-/dev/null}"
enforce="${PERF_BUDGETS_ENFORCE:-0}"
[ $# -ge 1 ] || { echo "usage: perf_budgets.sh OUT_DIR..." >&2; exit 2; }
dirs=("$@")
command -v jq > /dev/null || { echo "Perf budgets: jq missing, not judged" | tee -a "$sum"; exit 0; }

median() { sort -g | awk '{ v[NR] = $1 } END { if (NR) print (NR % 2) ? v[(NR + 1) / 2] : (v[NR / 2] + v[NR / 2 + 1]) / 2 }'; }

out="| scene | metric | value (median) | limit (ms) | margin | verdict |\n|---|---|---|---|---|---|\n"
flags=0; overs=0
while IFS= read -r line; do
    line="${line%%#*}"
    read -r -a w <<< "$line"
    n=${#w[@]}
    [ "$n" -ge 3 ] || continue
    scene="${w[0]}"; limit="${w[$((n - 1))]}"; metric="${w[*]:1:$((n - 2))}"
    val="$(for d in "${dirs[@]}"; do find "$d" -name results.json 2> /dev/null; done | while IFS= read -r f; do
        jq -r --arg s "$scene" --arg m "$metric" 'select(.scenario == $s) | .metrics[]? | select(.name == $m) | .value' "$f" 2> /dev/null
    done | median)"
    if [ -z "$val" ]; then out+="| $scene | $metric | - | $limit | - | not judged (no value) |\n"; continue; fi
    verdict="$(awk -v v="$val" -v l="$limit" 'BEGIN { if (l + 0 <= 0) print "nj"; else if (v <= l) print "ok"; else if (v <= l * 1.2) print "flag"; else print "over" }')"
    margin="$(awk -v v="$val" -v l="$limit" 'BEGIN { if (l + 0 <= 0) print "-"; else printf "%+.0f%%", (v / l - 1) * 100 }')"
    case "$verdict" in
        nj)   text="not judged" ;;
        ok)   text="ok" ;;
        flag) text="flag (up to 20% over)"; flags=$((flags + 1)); echo "::warning::perf budget: $scene $metric = $val ms, $margin over $limit ms" ;;
        over) overs=$((overs + 1))
              if [ "$enforce" = 1 ]; then text="**over (blocks)**"; else text="over (would block)"; fi
              echo "::warning::perf budget: $scene $metric = $val ms, $margin over $limit ms (more than 20%)" ;;
    esac
    out+="| $scene | $metric | $(printf '%.2f' "$val") | $limit | $margin | $text |\n"
done < "$file"

{
    echo ""
    echo "### Perf budgets (#1690): ${#dirs[@]} run dir(s), medians; $([ "$enforce" = 1 ] && echo enforcing || echo report-only)"
    [ -n "${PERF_BUDGETS_NOTE:-}" ] && echo "$PERF_BUDGETS_NOTE"
    printf '%b' "$out"
    echo ""
    echo "$overs over, $flags flagged. A limit of 0 is a placeholder (not judged yet)."
} | tee -a "$sum"
[ "$enforce" = 1 ] && [ "$overs" -gt 0 ] && exit 1
exit 0
