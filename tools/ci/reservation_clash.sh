#!/bin/bash
# (#2009) Static check 7: a kNetProtocol or kSaveVersion number this branch reserves in tools/ci/protocol_reservations.txt
# must not be reserved for another issue on the other integration base or in an open PR. Check 6 reads only this
# branch's own copy, so #2006 (into ladder) and #2008 (into staging) both took kNetProtocol 56.
#
#   tools/ci/reservation_clash.sh [BASE]   BASE = the commit to compare against: only lines added since BASE count as
#                                          this branch's (a clash among them is an error); without BASE every line is
#                                          checked and a clash is a warning
#
# Sources: origin/staging and origin/ladder (fetched, depth 1), and the lines that open PRs into staging and ladder add to
# the file (GitHub API; needs gh with a token, else that part is skipped with a note). A number reserved twice for the
# same issue is no clash (a sync between the bases repeats lines). Prints ::error / ::warning lines; exit 1 on an error.
set -uo pipefail
export LC_ALL=C
REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_DIR"
RES=tools/ci/protocol_reservations.txt
BASE="${1:-}"
GH_REPO="${GH_REPO:-DKoelzer01/RaylibGame}"
[ -f "$RES" ] || { echo "reservation_clash: no $RES"; exit 0; }

# "<key>\t<issue>\t<line>" for each reservation line: key = the number, or "save <n>"
parse() {
    tr -d '\r' | sed -nE 's/^\+?[[:space:]]*((save[[:space:]]+)?[0-9]+)[[:space:]]*\|[^|]*\|[[:space:]]*(#[0-9]+)[[:space:]]*\|.*$/\1\t\3\t&/p' |
        grep -v '(unused)' | sed -E 's/^save[[:space:]]+/save /; s/\t\+/\t/'
}

mine_all="$(parse < "$RES")"
if [ -n "$BASE" ]; then
    base_lines="$(git show "$BASE:$RES" 2>/dev/null | parse)"
    mine="$(grep -vxF -f <(printf '%s\n' "${base_lines:-#none}") <<< "$mine_all" || true)"
    level=error
else
    mine="$mine_all"; level=warning
fi
[ -n "$mine" ] || { echo "reservation_clash: this branch adds no reservation"; exit 0; }

others=""
for b in staging ladder; do
    git fetch -q --no-tags --depth=1 origin "$b" 2>/dev/null
    lines="$(git show "origin/$b:$RES" 2>/dev/null | parse)"
    [ -n "$lines" ] && others+="$(sed "s/\$/\torigin\/$b/" <<< "$lines")"$'\n'
done
notes=""
if command -v gh > /dev/null && gh auth status > /dev/null 2>&1 || [ -n "${GH_TOKEN:-}" ]; then
    self="${PR_NUMBER:-}"
    for b in staging ladder; do
        for n in $(gh api "repos/$GH_REPO/pulls?state=open&base=$b&per_page=100" -q '.[].number' 2>/dev/null); do
            [ "$n" = "$self" ] && continue
            patch="$(gh api "repos/$GH_REPO/pulls/$n/files?per_page=100" \
                -q ".[] | select(.filename == \"$RES\") | .patch" 2>/dev/null)"
            [ -n "$patch" ] || continue
            lines="$(grep -E '^\+[^+]' <<< "$patch" | parse)"
            [ -n "$lines" ] && others+="$(sed "s/\$/\tPR #$n/" <<< "$lines")"$'\n'
        done
    done
else
    notes="open PRs not checked (no gh token)"
fi

clashes=0
while IFS=$'\t' read -r key issue line; do
    [ -n "$key" ] || continue
    hits="$(awk -F'\t' -v k="$key" -v i="$issue" '$1 == k && $2 != i {print $4 " (" $2 ")"}' <<< "$others" | sort -u | paste -sd ';' -)"
    if [ -n "$hits" ]; then
        what="kNetProtocol $key"; [[ "$key" == save* ]] && what="kSaveVersion ${key#save }"
        echo "::$level file=$RES::static check 7: $what is reserved here for $issue and also in $hits: take the next free number (#2009)"
        clashes=$((clashes + 1))
    fi
done <<< "$mine"
echo "reservation_clash: $(wc -l <<< "$mine") reservation(s) checked, $clashes clash(es)${notes:+; $notes}"
[ "$clashes" -eq 0 ] || [ "$level" != error ]
