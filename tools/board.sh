#!/usr/bin/env bash
# tools/board.sh - cheap project-board and issue updates for RaylibGame (#582).
#
# GitHub gives the user 5000 GraphQL points an hour, shared by every session. `gh project item-edit --url` spends one
# request (plus lookups) per field, and `gh project item-list` reads the whole board. This tool batches every board write
# into one GraphQL request and does everything else through REST (a separate 5000/hour budget).
#
#   tools/board.sh get  N [N...]                        Status/Priority/Size/Kind per issue (one query per 50 issues)
#   tools/board.sh set  [-s S] [-p P] [-z Z] [-k K] N [N...]   set those fields on every issue (adds it to the board)
#   tools/board.sh fill [-s S] [-p P] [-z Z] [-k K] N [N...]   like set, but only fields that are still empty
#   tools/board.sh new  -t TITLE (-b BODY | -F FILE) [-l a,b] [--parent N] [-s S] [-p P] [-z Z] [-k K]
#                                                       create via REST, link the sub-issue via REST, set the board
# Values: Status Backlog|Ready|"In progress"|"In review"|Done; Priority P0-P3; Size XS|S|M|L|XL;
#         Kind Feature|Bug|Performance|Tracking.  Every query prints "graphql: cost C, remaining R" to stderr.
set -euo pipefail
REPO_OWNER=DKoelzer01; REPO_NAME=RaylibGame; REPO="$REPO_OWNER/$REPO_NAME"
PROJECT=PVT_kwHOAe7Hf84BrZzd
declare -A FIELD=( [Status]=PVTSSF_lAHOAe7Hf84BrZzdzhpssRE [Priority]=PVTSSF_lAHOAe7Hf84BrZzdzhpssSs
                   [Size]=PVTSSF_lAHOAe7Hf84BrZzdzhpssSw [Kind]=PVTSSF_lAHOAe7Hf84BrZzdzhpssZg )
declare -A OPT=(
  ["Status:Backlog"]=bf767b79 ["Status:Ready"]=2c8b23de ["Status:In progress"]=4d99a4b6 ["Status:In review"]=ecd2f2aa
  ["Status:Done"]=e280656f
  ["Priority:P0"]=e93b5bb8 ["Priority:P1"]=026eb06a ["Priority:P2"]=2906a46d ["Priority:P3"]=c36f5ea9
  ["Size:XS"]=39a448a0 ["Size:S"]=032ae8bd ["Size:M"]=6bd53109 ["Size:L"]=d7a23c47 ["Size:XL"]=ee6d351f
  ["Kind:Feature"]=dcf904b0 ["Kind:Bug"]=fae89327 ["Kind:Performance"]=8172895f ["Kind:Tracking"]=d260109d )
CHUNK=40   # mutations / issues per GraphQL request

die() { echo "board.sh: $*" >&2; exit 1; }

# Parsed options: WANT[Field]=Value for the fields given.
declare -A WANT=(); TITLE=""; BODY=""; BODYFILE=""; LABELS=""; PARENT=""; NUMS=()
parse() {
  while [ $# -gt 0 ]; do
    case "$1" in
      -s) WANT[Status]="$2"; shift 2;;   -p) WANT[Priority]="$2"; shift 2;;
      -z) WANT[Size]="$2"; shift 2;;     -k) WANT[Kind]="$2"; shift 2;;
      -t) TITLE="$2"; shift 2;;          -b) BODY="$2"; shift 2;;
      -F) BODYFILE="$2"; shift 2;;       -l) LABELS="$2"; shift 2;;
      --parent) PARENT="$2"; shift 2;;
      -*) die "unknown option $1";;
      *) NUMS+=("${1#\#}"); shift;;
    esac
  done
  for f in "${!WANT[@]}"; do [ -n "${OPT[$f:${WANT[$f]}]:-}" ] || die "bad $f value '${WANT[$f]}'"; done
}

# One GraphQL request, filtered with gh's built-in jq (there is no jq binary on this machine).
gql() { gh api graphql -f query="$1" --jq "$2"; }

# get: "N<TAB>itemId<TAB>Status<TAB>Priority<TAB>Size<TAB>Kind" per issue (item id empty when not on the board).
get_rows() {
  local nums=("$@") i
  for ((i = 0; i < ${#nums[@]}; i += CHUNK)); do
    local q="query { repository(owner:\"$REPO_OWNER\", name:\"$REPO_NAME\") {"
    for n in "${nums[@]:i:CHUNK}"; do
      q+=" i$n: issue(number:$n) { number projectItems(first:5) { nodes { id project { id } fieldValues(first:12) {"
      q+=" nodes { ... on ProjectV2ItemFieldSingleSelectValue { name field { ... on ProjectV2SingleSelectField { name } } } } } } } }"
    done
    q+=" } rateLimit { cost remaining resetAt } }"
    # One filter: a "COST ..." line (sent to stderr below) and one TSV row per issue.
    gql "$q" '("COST graphql: cost \(.data.rateLimit.cost), remaining \(.data.rateLimit.remaining), resets \(.data.rateLimit.resetAt)"),
      (.data.repository | to_entries[] | .value as $v
      | ([$v.projectItems.nodes[] | select(.project.id == "'"$PROJECT"'")] | first) as $it
      | (($it.fieldValues.nodes // []) | map(select(.field.name != null)) | map({(.field.name): .name}) | add // {}) as $f
      | [($v.number|tostring), ($it.id // ""), ($f.Status // ""), ($f.Priority // ""), ($f.Size // ""), ($f.Kind // "")] | @tsv)' \
    | while IFS= read -r line; do case "$line" in COST*) echo "${line#COST }" >&2;; *) printf '%s\n' "$line";; esac; done
  done
}

# add: ensure every issue is on the board; prints "N<TAB>itemId". One REST call per issue for the node id, then ONE
# GraphQL mutation for all of them (addProjectV2ItemById returns the existing item when it is already there).
add_items() {
  local nums=("$@") i
  for ((i = 0; i < ${#nums[@]}; i += CHUNK)); do
    local m="mutation {"
    for n in "${nums[@]:i:CHUNK}"; do
      local nid; nid=$(gh api "repos/$REPO/issues/$n" --jq .node_id)
      m+=" a$n: addProjectV2ItemById(input:{projectId:\"$PROJECT\", contentId:\"$nid\"}) { item { id } }"
    done
    m+=" }"
    gql "$m" '.data | to_entries[] | [(.key[1:]), .value.item.id] | @tsv'
  done
}

# set_fields ROWS: lines "itemId<TAB>Field<TAB>Value"; one batched mutation per CHUNK of updates.
set_fields() {
  local lines=() m="" k=0
  mapfile -t lines
  [ ${#lines[@]} -eq 0 ] && { echo "nothing to set"; return; }
  for line in "${lines[@]}"; do
    IFS=$'\t' read -r item field value <<< "$line"
    m+=" u$k: updateProjectV2ItemFieldValue(input:{projectId:\"$PROJECT\", itemId:\"$item\", fieldId:\"${FIELD[$field]}\","
    m+=" value:{singleSelectOptionId:\"${OPT[$field:$value]}\"}}) { projectV2Item { id } }"
    k=$((k + 1))
    if [ $((k % CHUNK)) -eq 0 ]; then gql "mutation {$m }" '.data | length' >/dev/null; m=""; fi
  done
  [ -n "$m" ] && gql "mutation {$m }" '.data | length' >/dev/null
  echo "set $k field values"
}

cmd="${1:-}"; shift || true
case "$cmd" in
  get)
    parse "$@"; [ ${#NUMS[@]} -gt 0 ] || die "get N [N...]"
    get_rows "${NUMS[@]}" | while IFS=$'\t' read -r n item s p z k; do
      printf '#%s  %-12s %-3s %-3s %-12s%s\n' "$n" "${s:--}" "${p:--}" "${z:--}" "${k:--}" "$([ -z "$item" ] && echo '  (not on board)')"; done;;
  set|fill)
    parse "$@"; [ ${#NUMS[@]} -gt 0 ] || die "$cmd [-s S] [-p P] [-z Z] [-k K] N [N...]"
    [ ${#WANT[@]} -gt 0 ] || die "give at least one of -s -p -z -k"
    declare -A ITEM=() HAVE=()
    if [ "$cmd" = fill ]; then
      while IFS=$'\t' read -r n item s p z k; do ITEM[$n]="$item"; HAVE[$n:Status]="$s"; HAVE[$n:Priority]="$p"; HAVE[$n:Size]="$z"; HAVE[$n:Kind]="$k"
      done < <(get_rows "${NUMS[@]}")
    fi
    missing=(); for n in "${NUMS[@]}"; do [ -n "${ITEM[$n]:-}" ] || missing+=("$n"); done
    if [ ${#missing[@]} -gt 0 ]; then while IFS=$'\t' read -r n item; do ITEM[$n]="$item"; done < <(add_items "${missing[@]}"); fi
    for n in "${NUMS[@]}"; do for f in "${!WANT[@]}"; do
      if [ "$cmd" = set ] || [ -z "${HAVE[$n:$f]:-}" ]; then printf '%s\t%s\t%s\n' "${ITEM[$n]}" "$f" "${WANT[$f]}"; fi
    done; done | set_fields;;
  new)
    parse "$@"; [ -n "$TITLE" ] || die "new -t TITLE (-b BODY | -F FILE) ..."
    [ -n "$BODYFILE" ] && BODY=$(cat "$BODYFILE")
    args=(-f "title=$TITLE" -f "body=$BODY"); IFS=',' read -ra ls <<< "$LABELS"; for l in "${ls[@]}"; do [ -n "$l" ] && args+=(-f "labels[]=$l"); done
    read -r num dbid < <(gh api "repos/$REPO/issues" "${args[@]}" --jq '"\(.number) \(.id)"')
    echo "created #$num https://github.com/$REPO/issues/$num"
    [ -n "$PARENT" ] && gh api "repos/$REPO/issues/$PARENT/sub_issues" -F "sub_issue_id=$dbid" --silent && echo "linked under #$PARENT"
    if [ ${#WANT[@]} -gt 0 ]; then
      item=$(add_items "$num" | cut -f2)
      for f in "${!WANT[@]}"; do printf '%s\t%s\t%s\n' "$item" "$f" "${WANT[$f]}"; done | set_fields
    fi;;
  *) sed -n '2,16p' "$0"; exit 1;;
esac
