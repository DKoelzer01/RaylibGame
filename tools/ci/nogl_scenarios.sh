#!/bin/bash
# (#1960) The sandbox-world scenarios that pass without a GPU, for GitHub's hosted Linux runners (pr-checks.yml): the
# null GL (RAYLIBGAME_NO_GL=1) draws nothing and turns the GPU paths off, so this covers the sim, input, UI-state, save
# and net logic on every PR without self-hosted runner time. The GPU runner still runs everything (selfhosted.yml).
#
#   tools/ci/nogl_scenarios.sh GAME SHARD COUNT [OUTDIR]   run shard SHARD (1-based) of COUNT
#   tools/ci/nogl_scenarios.sh GAME list                   print the selected scenarios
#
# Selected: world "sandbox" (no planet: no erosion bake), minus tags that need pixels or timing (perf, visual, light,
# fx, local) and minus tools/ci/nogl_exclude.txt (scenarios that need the GPU path). Shards take every COUNT-th name
# of the sorted list, so a new scenario lands in one shard without reshuffling the rest much.
set -uo pipefail
game="$1"; shift
here="$(cd "$(dirname "$0")" && pwd)"
export RAYLIBGAME_NO_GL=1
all="$("$game" --scenario-list 2>/dev/null)" || { echo "nogl: --scenario-list failed" >&2; exit 1; }
exclude="$(sed -e 's/#.*//' -e 's/[[:space:]]*$//' -e 's/\r$//' "$here/nogl_exclude.txt" | grep -v '^$' || true)"
names="$(echo "$all" | awk '$2 == "sandbox" && $0 !~ /tags [^ ]*(perf|visual|light|fx|local)/ {print $1}' | sort |
    grep -vxF -f <(echo "${exclude:-#none}") || true)"
[ -n "$names" ] || { echo "nogl: no scenarios selected" >&2; exit 1; }
if [ "${1:-}" = list ]; then echo "$names"; exit 0; fi
shard="$1"; count="$2"; out="${3:-build/nogl-scenarios}"
pick="$(echo "$names" | awk -v s="$shard" -v c="$count" '(NR - 1) % c == s - 1' | paste -sd, -)"
[ -n "$pick" ] || { echo "nogl: shard $shard of $count is empty"; exit 0; }
echo "nogl: shard $shard of $count: $(echo "$pick" | tr ',' '\n' | wc -l) of $(echo "$names" | wc -l) scenarios"
mkdir -p "$out"
"$game" --scenario "$pick" --scenario-no-record --scenario-out "$out" > "$out/run.log" 2>&1
rc=$?
grep -E -- '^scenario: [^ ]+ -> ' "$out/run.log"
passes="$(grep -cE -- '^scenario: [^ ]+ -> PASS' "$out/run.log" || true)"
# (#2123) Each failure, "known (#N)" when it matches a line of known_failures_nogl.txt (this checkout's plus staging's),
# else NEW. Only NEW failures fail the shard: "CI ok" needs the shards, and a tracked flake must not block every merge.
known="$(bash "$here/known_list.sh" known_failures_nogl.txt)"
new=0 old=0 marked=""
while IFS= read -r line; do
    [ -z "$line" ] && continue
    tag=NEW
    while IFS='|' read -r pat issue _; do
        pat="$(echo "$pat" | sed 's/[[:space:]]*$//')"; [[ -z "$pat" || "$pat" == \#* ]] && continue
        if echo "$line" | grep -Eq -- "$pat"; then tag="known $(echo $issue)"; break; fi
    done <<< "$known"
    if [ "$tag" = NEW ]; then new=$((new + 1)); else old=$((old + 1)); fi
    marked+="- **[$tag]** \`$(echo "$line" | cut -c1-200)\`"$'\n'
done < <(grep -E -- '^scenario: [^ ]+ -> (FAIL|TIMEOUT|ERROR)' "$out/run.log" | tr -d '\r')
echo "nogl: $passes passed, $new NEW and $old known failures (exit $rc)"
printf '%s' "$marked"
if [ -n "${GITHUB_STEP_SUMMARY:-}" ]; then
    { echo "### Null-GL scenarios, shard $shard of $count: $passes passed, $new NEW and $old known failures"
      printf '%s' "$marked"
      [ $new -gt 0 ] && echo "" && echo "NEW failures need a fix, or an issue and a line in \`tools/ci/known_failures_nogl.txt\`."; } >> "$GITHUB_STEP_SUMMARY"
fi
# Red unless every picked scenario printed a verdict (a crash cuts the list short) and the exit code is 0, or 1 (the
# game's "some scenario failed") with only known failures.
picked=$(echo "$pick" | tr ',' '\n' | grep -c .)
if [ $((passes + new + old)) -ne "$picked" ]; then echo "nogl: $((passes + new + old)) verdicts for $picked scenarios (exit $rc): a crash?"; exit 1; fi
[ $new = 0 ] && { [ "$rc" = 0 ] || { [ "$rc" = 1 ] && [ $old -gt 0 ]; }; }
