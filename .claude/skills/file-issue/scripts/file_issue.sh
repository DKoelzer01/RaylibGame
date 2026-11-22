#!/bin/bash
# busybox (w64devkit on PATH) is not bash: re-run under Git's bash (#1477)
if [ -z "${BASH_VERSINFO:-}" ]; then for b in "C:/Program Files/Git/bin/bash.exe" "/c/Program Files/Git/bin/bash.exe" "/usr/bin/bash"; do [ -x "$b" ] && exec "$b" "$0" "$@"; done; echo "needs bash" >&2; exit 2; fi
# Search for an existing issue, then (with --create) file a new one (#1476). REST only, except gh issue create.
#   file_issue.sh search "<keywords>"
#   file_issue.sh create "<title>" <body-file> [--labels a,b] [--parent N]
R=DKoelzer01/RaylibGame
case "$1" in
  search)
    gh api -X GET search/issues -f q="repo:$R is:issue $2" -f per_page=10 \
      --jq '.items[] | "#\(.number) [\(.state)] \(.title)"'
    ;;
  create)
    title=$2; body=$3; shift 3
    args=(--title "$title" --body-file "$body")
    while [ $# -gt 0 ]; do
      case "$1" in
        --labels) IFS=, read -ra ls <<< "$2"; for l in "${ls[@]}"; do args+=(--label "$l"); done; shift 2 ;;
        --parent) args+=(--parent "$2"); shift 2 ;;
        *) echo "unknown option $1" >&2; exit 2 ;;
      esac
    done
    url=$(gh issue create "${args[@]}" 2>&1 | tail -1)
    echo "$url"
    case "$url" in https://*) echo "issue #${url##*/}";; *) exit 1;; esac
    ;;
  *) sed -n '2,5p' "$0"; exit 2 ;;
esac
