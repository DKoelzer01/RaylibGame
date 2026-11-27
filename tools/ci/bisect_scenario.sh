#!/bin/bash
# busybox (w64devkit on PATH) is not bash: re-run under Git's bash (#1477)
if [ -z "${BASH_VERSINFO:-}" ]; then for b in "C:/Program Files/Git/bin/bash.exe" "/c/Program Files/Git/bin/bash.exe" "/usr/bin/bash"; do [ -x "$b" ] && exec "$b" "$0" "$@"; done; echo "needs bash" >&2; exit 2; fi
# (#1871) Which merge broke a scenario: a binary search over the first-parent commits in <good>..<bad> (on staging these
# are the PR merges), each one built and run in the deweyserver build sandbox (~/rg-sandbox/rb.sh takes a commit as well
# as a branch). Nothing runs on this PC and the work tree is not touched. For the batched staging run: its summary names
# the commit range it covers; give the last green batch as <good> and the red one as <bad>.
#
#   tools/ci/bisect_scenario.sh <scenario> <good> <bad> [extra game args]
#   tools/ci/bisect_scenario.sh fluid/flow-meter 1a2b3c4d 5e6f7a8b
#   tools/ci/bisect_scenario.sh --dry fluid/flow-meter origin/develop origin/staging   (prints the commits and the steps)
#
# A commit whose build fails is skipped (its neighbour is tried). Each step is one sandbox job, ~5-15 min with a warm
# cache, so a range of N merges takes about log2(N) of them. Commits must be pushed (GitHub serves a commit by its id).
set -u
dry=0; [ "${1:-}" = --dry ] && { dry=1; shift; }
[ $# -ge 3 ] || { sed -n '4,13p' "$0"; exit 2; }
scen="$1" goodName="$2" badName="$3"
good=$(git rev-parse --verify "$2^{commit}") && bad=$(git rev-parse --verify "$3^{commit}") || exit 2
shift 3; extra="$*"
git merge-base --is-ancestor "$good" "$bad" || { echo "bisect: $goodName is not an ancestor of $badName"; exit 2; }
# oldest first; the last entry is <bad>
mapfile -t revs < <(git rev-list --first-parent --reverse "$good..$bad")
n=${#revs[@]}
[ "$n" -gt 0 ] || { echo "bisect: no commits in $goodName..$badName"; exit 2; }
steps=0; x=$n; while [ "$x" -gt 1 ]; do x=$(( (x + 1) / 2 )); steps=$((steps + 1)); done
echo "bisect: $scen over $n first-parent commits in ${good:0:8}..${bad:0:8}, about $steps sandbox runs"
ssh="C:/Windows/System32/OpenSSH/ssh.exe"; [ -x "$ssh" ] || ssh=ssh
log="build/bisect-${scen//\//_}-$(date +%Y%m%d-%H%M%S).log"; mkdir -p build
if [ $dry = 1 ]; then for r in "${revs[@]}"; do git log -1 --format='  %h %s' "$r"; done; exit 0; fi
# run <sha>: 0 pass, 1 fail, 125 skip (build failed or the scenario isn't there), exits on a busy sandbox
run() {
    local rc
    if [ -n "${BISECT_FAKE:-}" ]; then # tests of the search itself: a command that gets the commit and exits 0/1/125
        $BISECT_FAKE "$1" > "$log.$2" 2>&1; rc=$?
    else
        local q="export RAYLIBGAME_SHADER_CACHE=/cache/shaders; \$GAME --scenario '$scen' --scenario-no-record $extra"
        q=${q//\'/\'\\\'\'}
        "$ssh" deweyserver "~/rg-sandbox/rb.sh '$1' -- '$q'" > "$log.$2" 2>&1; rc=$?
    fi
    local step="$log.$2" verdict=1
    cat "$step" >> "$log"
    if [ $rc = 0 ]; then verdict=0
    elif grep -q "^sandbox: busy" "$step"; then rm -f "$step"; echo "bisect: the sandbox stayed busy for 30 min; stopped (log: $log)"; exit 3
    elif [ $rc = 125 ] || grep -q "BUILD FAILED\|no scenario matches" "$step"; then verdict=125; fi
    rm -f "$step"
    return $verdict
}
# Check the ends first: a scenario that fails at <good> too (flaky, or the environment) makes a bisect meaningless.
lo=-1 hi=$((n - 1)) # revs[lo] passes (lo=-1 is <good>), revs[hi] fails
echo "bisect: checking ${bad:0:8} fails ..."
run "$bad" bad; rc=$?
[ $rc = 0 ] && { echo "bisect: $scen passes at ${bad:0:8}: nothing to find (flaky? see $log)"; exit 1; }
[ $rc = 125 ] && { echo "bisect: ${bad:0:8} doesn't build or lacks $scen (see $log)"; exit 2; }
echo "bisect: checking ${good:0:8} passes ..."
run "$good" good; rc=$?
[ $rc = 0 ] || { echo "bisect: $scen fails at ${good:0:8} too: pick an older <good> (see $log)"; exit 1; }
declare -A skipped=()
while [ $((hi - lo)) -gt 1 ]; do
    mid=$(( (lo + hi) / 2 ))
    # a skipped commit: try its untested neighbours inside (lo, hi)
    probe=""
    for ((d = 0; mid + d < hi || mid - d > lo; d++)); do
        for p in $((mid + d)) $((mid - d)); do
            if [ $p -gt $lo ] && [ $p -lt $hi ] && [ -z "${skipped[$p]:-}" ]; then probe=$p; break 2; fi
        done
    done
    [ -z "$probe" ] && break # everything between lo and hi failed to build
    r=${revs[$probe]}
    echo "bisect: $(git log -1 --format='%h %s' "$r") ..."
    run "$r" "s$probe"; rc=$?
    case $rc in
        0) echo "  pass"; lo=$probe ;;
        1) echo "  FAIL"; hi=$probe ;;
        *) echo "  skipped (build failed)"; skipped[$probe]=1 ;;
    esac
done
echo
if [ $((hi - lo)) -eq 1 ]; then
    echo "bisect: first failing commit for $scen:"
else
    echo "bisect: the first failing commit is one of these (the others in between didn't build):"
    for ((i = lo + 1; i < hi; i++)); do git log -1 --format='  %h %s' "${revs[$i]}"; done
fi
git log -1 --format='  %h %s (%an, %ad)' --date=short "${revs[$hi]}"
echo "(log: $log)"
