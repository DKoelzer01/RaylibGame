#!/bin/bash
# Saves a compact per-run summary of a self-hosted job (#1476) into ci-persist/summaries/<run id>/<job>/, outside the
# slot folders, so the deweyserver MCP tools (ci_summary, scenario_history) can read results without downloading logs.
#   save_summary.sh <job name>      (runs from the checkout; RUN, REF, GITHUB_RUN_ID, GITHUB_SHA come from the job)
# Keeps 30 days. Never fails the job.
set -u
job=${1:-job}
root="${RUNNER_WORKSPACE:-}/ci-persist/summaries"
[ -n "${RUNNER_WORKSPACE:-}" ] && [ -n "${GITHUB_RUN_ID:-}" ] || exit 0
dest="$root/$GITHUB_RUN_ID/$job"
mkdir -p "$dest" || exit 0
run="${RUN:-}"
# Only files written during this job count: newer than the slot file, which the job touches when it claims its slot. (#1786:
# a whole-minute age with 2 minutes' slack took another run's regression summary, written 39 s before this job began.)
# Without the slot file: the last 180 minutes.
fresh=(-mmin -180)
lock="${RUNNER_WORKSPACE}/ci-persist/slot${CI_SLOT:-}.lock"
[ -n "${CI_SLOT:-}" ] && [ -f "$lock" ] && fresh=(-newer "$lock")
{
    echo "{\"run\": $GITHUB_RUN_ID, \"attempt\": ${GITHUB_RUN_ATTEMPT:-1}, \"job\": \"$job\", \"branch\": \"${REF#refs/heads/}\","
    echo " \"sha\": \"${GITHUB_SHA:-}\", \"event\": \"${GITHUB_EVENT_NAME:-}\", \"saved\": \"$(date -u +%FT%TZ)\"}"
} > "$dest/meta.json"
cp tools/ci/known_failures.txt "$dest/" 2> /dev/null
if [ -n "$run" ] && [ -d "$run" ]; then
    # One line per scenario: name, verdict, reason, times (the full results.json files stay in the job's artifact).
    : > "$dest/scenarios.tsv"
    while IFS= read -r f; do
        # results.json is written by the game with one key per line: pull the few fields without jq
        name=$(grep -m1 '"scenario"' "$f" | sed 's/.*: *"\(.*\)",*/\1/')
        verdict=$(grep -m1 '"verdict"' "$f" | sed 's/.*: *"\(.*\)",*/\1/')
        reason=$(grep -m1 '"reason"' "$f" | sed 's/.*: *"\(.*\)",*/\1/' | tr '\t' ' ' | cut -c1-300)
        wall=$(grep -m1 '"wall_seconds"' "$f" | sed 's/.*: *\([0-9.]*\).*/\1/')
        started=$(grep -m1 '"started"' "$f" | sed 's/.*: *"\(.*\)",*/\1/')
        printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$verdict" "$wall" "$started" "$reason" >> "$dest/scenarios.tsv"
    done < <(find "$run/scenarios" -name results.json "${fresh[@]}" 2> /dev/null)
    # Every FAIL line this job printed (scenario logs, the GPU check, the join tests, the regression). The slot's run
    # folder outlives jobs, so only files written during this job count ("${fresh[@]}", above).
    find "$run" -maxdepth 3 \( -name 'scenarios-*.log' -o -name 'terrain-gpu.log' -o -name 'regression-*.log' -o -path '*/join-logs/*.log' \) \
        "${fresh[@]}" -print0 2> /dev/null | xargs -0 -r cat | grep -E -- '-> FAIL|^=== ' | cut -c1-400 | sort -u > "$dest/fail_lines.txt"
    find "$run" -maxdepth 1 -name 'regression-summary-*.md' "${fresh[@]}" -exec cp {} "$dest/" \; 2> /dev/null
fi
# Retention: 30 days
find "$root" -mindepth 1 -maxdepth 1 -type d -mtime +30 -exec rm -rf {} + 2> /dev/null
echo "summary saved: $dest ($(wc -l < "$dest/scenarios.tsv" 2> /dev/null || echo 0) scenarios)"
exit 0
