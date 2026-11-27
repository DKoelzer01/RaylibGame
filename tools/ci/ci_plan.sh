#!/bin/bash
# (#1871) What a self-hosted run does: its tier, and for a PR's run the scenarios its changes need. Writes to $GITHUB_OUTPUT
# (or prints): tier=affected|batch|full, skip=true|false, scenarios=<--scenario list, "ALL", or "">, server=true|false,
# regression=true|false, and says why in the job summary.
#
#   affected  a dispatch on a PR's branch (its runner evidence): the build, the unit tests, the GPU check when its paths
#             changed (tools/ci/gpu_changes.sh), and only the scenarios tools/ci/scenario_map.txt maps its changes to
#             (diffed against the merge base with the integration branch or ladder, whichever is nearer).
#   batch     the integration branch every few hours (the schedule), or a dispatch on it: all headless scenarios, the
#             server and join tests, the nightly-only steps. Skipped when the branch hasn't moved since the last batch.
#   full      the nightly, a push to develop, or a dispatch with tier=full (before staging -> develop and develop -> main):
#             the batch plus the headless regression (the nightly, or when asked).
#
# Env: EVENT, SCHEDULE (github.event.schedule), REF, TIER_INPUT (auto|affected|batch|full), REGRESSION_INPUT (none|headless),
# INTEGRATION (the integration branch: vars.INTEGRATION_BRANCH, default staging),
# PERSIST (ci-persist on the runners: the last batched commit per branch), HEAD_SHA (tests: another commit than HEAD).
# --dry: print, don't write markers.
set -u
dry=0; [ "${1:-}" = --dry ] && dry=1
out="${GITHUB_OUTPUT:-/dev/stdout}"
summary="${GITHUB_STEP_SUMMARY:-/dev/null}"
integration="${INTEGRATION:-staging}"
branch="${REF#refs/heads/}"
nightlyCron="0 5 * * *"
tier="" skip=false regression=false why="" kind=""
case "${EVENT:-workflow_dispatch}" in
    schedule)
        if [ "${SCHEDULE:-}" = "$nightlyCron" ]; then tier=full; regression=true; kind=nightly; else tier=batch; kind=batch; fi ;;
    push)
        tier=full; why="a push to $branch" ;;
    *)
        tier="${TIER_INPUT:-auto}"
        # scheduled: a batch started from outside (a cron on deweyserver until this file's schedule is on main, which is the
        # only copy GitHub's schedules read): the same skip as the schedule's
        if [ "$tier" = scheduled ]; then tier=batch; kind=batch
        elif [ "$tier" = auto ] || [ -z "$tier" ]; then
            case "$branch" in "$integration"|develop|main) tier=batch; why="a dispatch on $branch" ;; *) tier=affected; why="a dispatch on $branch (a PR's evidence)" ;; esac
        else why="a dispatch asking for tier $tier"; fi ;;
esac
if [ -n "$kind" ]; then # a timed run: skipped when the branch hasn't moved since the last one of its kind
    head=$(git rev-parse "${HEAD_SHA:-HEAD}" 2>/dev/null || echo "")
    marker="${PERSIST:-/tmp}/ci-last-$kind-${branch//\//_}.sha"
    if [ -n "$head" ] && [ -f "$marker" ] && [ "$(cat "$marker")" = "$head" ]; then
        skip=true; why="$branch hasn't moved since the last $kind run (${head:0:8})"
    else
        last=$(cut -c1-8 "$marker" 2>/dev/null || true)
        lastFull=$(tr -d '[:space:]' 2>/dev/null < "$marker" || true) # (#1947) the batch's server gate diffs from here
        why="the $kind run of $branch at ${head:0:8}$([ -n "$last" ] && echo ", covering $last..${head:0:8}: $(git rev-list --count "$last..$head" 2>/dev/null || echo '?') commits")"
        [ $dry = 1 ] || { mkdir -p "$(dirname "$marker")" && echo "$head" > "$marker"; }
    fi
fi
[ "${REGRESSION_INPUT:-none}" = headless ] && regression=true

scenarios="" server=false
if [ "$tier" = affected ]; then
    # The changed paths against the nearer base (the integration branch or ladder: a PR targets one of them).
    files="" base=""
    for b in "$integration" ladder; do
        git fetch --quiet --depth=300 origin "$b" 2>/dev/null || continue
        mb=$(git merge-base "${HEAD_SHA:-HEAD}" "origin/$b" 2>/dev/null) || continue
        f=$(git diff --name-only "$mb" "${HEAD_SHA:-HEAD}" 2>/dev/null) || continue
        if [ -z "$base" ] || [ "$(printf '%s\n' "$f" | grep -c .)" -lt "$(printf '%s\n' "$files" | grep -c .)" ]; then files="$f"; base="$b@${mb:0:8}"; baseSha="$mb"; fi
    done
    if [ -z "$base" ]; then
        scenarios=ALL; server=true; why="$why; no merge base with $integration or ladder: everything"
    else
        map="$(dirname "$0")/scenario_map.txt"
        declare -A want=()
        all=0 unmatched=0
        add() { local t; IFS=',' read -ra ts <<< "$1"; for t in "${ts[@]}"; do t="${t// /}"; [ -z "$t" ] && continue
                case "$t" in ALL) all=1 ;; SERVER) server=true ;; *) want["$t"]=1 ;; esac; done; }
        mapText=$(tr -d '\015' < "$map") # a Windows checkout has CRLF
        default=$(grep -E '^DEFAULT \| ' <<< "$mapText" | head -1); default="${default##* | }"
        while read -r f; do
            [ -z "$f" ] && continue
            hit=0
            # the scenarios in a changed scenario file, unless they are local-only (static check 5's list: they need the PC)
            if [[ "$f" == tests/scenarios/*.cpp ]] && [ -f "$f" ] && ! grep -o 's.tags = {[^}]*}' "$f" | grep -qE '"(perf|visual|fx|light|local)"'; then
                for n in $(grep -o 'SCENARIO("[^"]*")' "$f" | sed 's/SCENARIO("\(.*\)")/\1/'); do add "$n"; hit=1; done
            fi
            while IFS= read -r line; do
                case "$line" in ""|"#"*|"DEFAULT | "*) continue ;; esac
                re="${line% | *}" tg="${line##* | }"
                printf '%s\n' "$f" | grep -Eq "$re" && { add "$tg"; hit=1; }
            done <<< "$mapText"
            [ $hit = 0 ] && case "$f" in src/*|resources/*) unmatched=1 ;; esac
        done <<< "$files"
        [ $unmatched = 1 ] && add "$default"
        # A mapped name this branch doesn't have (the map serves staging and ladder) is dropped: the game exits 2 when a
        # --scenario list matches nothing, and a list with one unknown name among others would hide it.
        for n in "${!want[@]}"; do
            case "$n" in tag:*) continue ;; esac
            grep -rqF "SCENARIO(\"$n\"" tests/scenarios 2>/dev/null || { unset "want[$n]"; why="$why; $n is not on this branch"; }
        done
        if [ $all = 1 ]; then scenarios=ALL; server=true
        else scenarios=$(printf '%s\n' "${!want[@]}" | sort | paste -sd, -); fi
        why="$why; $(printf '%s\n' "$files" | grep -c .) files changed against $base"
    fi
elif [ "$tier" = batch ] && [ -n "${lastFull:-}" ] && changed=$(git diff --name-only "$lastFull" "${HEAD_SHA:-HEAD}" 2>/dev/null); then
    # (#1947) The staging batch runs every scenario, but the dedicated server and join tests only when the commits since
    # the last batch touch a path the map sends to SERVER or ALL (the full tier always runs them). No marker, or a marker
    # outside the checkout's history: they run.
    scenarios=ALL; server=false
    mapText=$(tr -d '\015' < "$(dirname "$0")/scenario_map.txt")
    hitBy=""
    while read -r f; do
        [ -z "$f" ] && continue
        while IFS= read -r line; do
            case "$line" in ""|"#"*|"DEFAULT | "*) continue ;; esac
            re="${line% | *}" tg="${line##* | }"
            printf '%s\n' "$tg" | grep -Eq '(^|,) *(SERVER|ALL) *(,|$)' || continue
            printf '%s\n' "$f" | grep -Eq "$re" && { server=true; hitBy="$f"; break 2; }
        done <<< "$mapText"
    done <<< "$changed"
    why="$why; server and join tests: $([ $server = true ] && echo "$hitBy changed" || echo "no network or host path changed since ${lastFull:0:8}")"
else
    scenarios=ALL; server=true
fi
{
    echo "tier=$tier"
    echo "skip=$skip"
    echo "scenarios=$scenarios"
    echo "server=$server"
    echo "regression=$regression"
    echo "base=${baseSha:-}" # (the affected tier: the merge base, so the GPU check is gated on the PR's own changes)
} >> "$out"
line="CI tier: **$tier**$([ $skip = true ] && echo ' (skipped)') - $why. Scenarios: ${scenarios:-none}; server and join tests: $server; regression: $regression"
echo "$line"
echo "$line" >> "$summary"
