#!/bin/bash
# busybox (w64devkit on PATH) is not bash: re-run under Git's bash (#1477)
if [ -z "${BASH_VERSINFO:-}" ]; then for b in "C:/Program Files/Git/bin/bash.exe" "/c/Program Files/Git/bin/bash.exe" "/usr/bin/bash"; do [ -x "$b" ] && exec "$b" "$0" "$@"; done; echo "needs bash" >&2; exit 2; fi
# Mechanical pre-review of a diff (#1877): what a script can catch before a Claude session spends review turns on it.
# Used by sessions before they ask for review or label a PR ready, on subagent diffs, and as the guardrail of the planned
# local coding tool (local_code: --files and --max-lines are its chunk scope).
#
#   tools/ci/prereview.sh [BASE] [--files a,b,dir/,*.glob] [--max-lines N] [--committed]
#
#   BASE         what to compare against (default origin/staging); the diff starts at merge-base(BASE, HEAD)
#   --files      the chunk's scope: changed files outside it are blockers (comma list of paths, dir/ prefixes, globs)
#   --max-lines  changed lines (added + removed) above N are a blocker; without it, a warning above 600
#   --committed  only committed changes (default: committed + uncommitted in the work tree; untracked files are listed)
#
# Checks: scope; diff size; CR in .sh files (they must stay LF: bash fails on CR) and mixed line endings in other changed
# files; debug leftovers in added lines (std::cout, #if 0, DEBUG/dbg prints, TODO remove, XXX, HACK); changes Claude must
# look at (kNetProtocol / NetMsg, save format, new headers, CMake, workflows, resources/data); tools/ci/static_checks.sh.
# Exit 0 clean, 1 warnings only, 2 blockers. Output is one short report.
set -uo pipefail
export LC_ALL=C
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"

base="" files="" maxl="" committed=0
while [ $# -gt 0 ]; do
    case "$1" in
        --files) files="$2"; shift 2 ;;
        --max-lines) maxl="$2"; shift 2 ;;
        --committed) committed=1; shift ;;
        -h|--help) sed -n '4,19p' "$0"; exit 0 ;;
        -*) echo "prereview: unknown option $1" >&2; exit 2 ;;
        *) base="$1"; shift ;;
    esac
done
base=${base:-origin/staging}
mb=$(git merge-base "$base" HEAD 2>/dev/null) || { echo "prereview: no merge base with $base" >&2; exit 2; }
if [ $committed = 1 ]; then range=("$mb" HEAD); else range=("$mb"); fi

blockers=() warnings=() notes=()
block() { blockers+=("$*"); }
warn() { warnings+=("$*"); }
note() { notes+=("$*"); }

mapfile -t changed < <(git diff --name-only --diff-filter=ACMRT "${range[@]}")
mapfile -t deleted < <(git diff --name-only --diff-filter=D "${range[@]}")
mapfile -t added_files < <(git diff --name-only --diff-filter=A "${range[@]}")
untracked=()
[ $committed = 0 ] && mapfile -t untracked < <(git ls-files --others --exclude-standard)
read -r nadd ndel < <(git diff --numstat "${range[@]}" | awk '$1 != "-" { a += $1; d += $2 } END { print a + 0, d + 0 }')
total=$((nadd + ndel))

# --- scope
in_scope() { # path -> 0 if it matches one of the --files entries
    local p="$1" f IFS=,
    for f in $files; do
        f=${f# }
        case "$f" in
            */) [[ "$p" == "$f"* ]] && return 0 ;;
            *\**|*\?*) [[ "$p" == $f ]] && return 0 ;;
            *) [ "$p" = "$f" ] && return 0 ;;
        esac
    done
    return 1
}
if [ -n "$files" ]; then
    for p in "${changed[@]}" "${deleted[@]}" "${untracked[@]}"; do
        [ -n "$p" ] || continue
        in_scope "$p" || block "out of scope: $p"
    done
fi

# --- size
if [ -n "$maxl" ]; then
    [ "$total" -gt "$maxl" ] && block "diff too large: $total changed lines (max $maxl)"
elif [ "$total" -gt 600 ]; then
    warn "large diff: $total changed lines; consider splitting it"
fi

# --- line endings (work-tree files: git stores LF and converts on checkout, so the diff can't show them)
for p in "${changed[@]}" "${untracked[@]}"; do
    [ -f "$p" ] || continue
    case "$p" in *.png|*.jpg|*.wav|*.ogg|*.ttf|*.bin|*.gz|*.zip|*.exe|*.dll|tests/saves/*) continue ;; esac
    # count CRs with tr: a native Windows awk reads in text mode and never sees them
    crlf=$(tr -cd '\r' < "$p" | wc -c); lf=$(( $(wc -l < "$p") - crlf ))
    case "$p" in
        *.sh) [ "$crlf" -gt 0 ] && block "$p: $crlf CRLF lines (shell scripts must be LF: .gitattributes eol=lf)" ;;
        *) [ "$crlf" -gt 0 ] && [ "$lf" -gt 0 ] && warn "$p: mixed line endings ($crlf CRLF, $lf LF): Edit-tool matching may fail" ;;
    esac
done

# --- debug leftovers in added lines of code
leftovers=$(git diff -U0 "${range[@]}" -- '*.cpp' '*.h' '*.hpp' '*.c' '*.glsl' '*.vs' '*.fs' '*.sh' ':!tools/ci/prereview.sh' |
    awk '/^\+\+\+ / { f = substr($0, 7); next }
         /^@@/ { split($3, a, ","); n = substr(a[1], 2) + 0; next }
         /^\+/ { s = substr($0, 2)
                 if (s ~ /std::cout|std::cerr|^#if 0|printf\(" *(DEBUG|DBG|dbg|XXX)|fprintf\(stderr, *" *(DEBUG|DBG|dbg)|TODO:? *(remove|delete|revert)|FIXME:? *(remove|temp)|[^A-Za-z]XXX[^A-Za-z]|[^A-Za-z]HACK[^A-Za-z]|\/\/ *(debug|DEBUG)[: ]/)
                     { gsub(/^[ \t]+/, "", s); printf "%s:%d: %s\n", f, n, substr(s, 1, 120) }
                 n++ }')
if [ -n "$leftovers" ]; then
    while IFS= read -r l; do warn "leftover? $l"; done <<< "$leftovers"
fi

# --- changes a Claude session must look at
for p in "${changed[@]}" "${deleted[@]}"; do
    case "$p" in
        src/engine/netproto.h) note "netproto.h changed: NetMsg / kNetProtocol (protocol changes are Claude's; reserve the number)" ;;
        src/engine/saveversion.h|src/engine/savemigrate.*|*/save*.cpp|*/save*.h) note "save code changed: $p (save format?)" ;;
        CMakeLists.txt|*/CMakeLists.txt|cmake/*) note "build files changed: $p" ;;
        .github/*) note "workflow/CI change: $p (needs the owner's review)" ;;
        resources/data/*) note "data changed: $p" ;;
        CLAUDE.md|.claude/*) note "agent config changed: $p" ;;
    esac
done
for p in "${added_files[@]}" "${untracked[@]}"; do
    case "$p" in src/*.h|src/*.hpp) note "new header: $p (new interface?)" ;; esac
done
[ ${#untracked[@]} -gt 0 ] && note "untracked files (not in the diff yet): ${untracked[*]}"

# --- static checks (NetMsg ids, bake version, protocol reservations, scenario steps...)
sc_out=$(GITHUB_STEP_SUMMARY=/dev/null bash tools/ci/static_checks.sh "$mb" 2>&1)
sc_rc=$?
# warnings about files this diff doesn't touch are already on the base branch: not this diff's problem
changed_set=" ${changed[*]:-} ${untracked[*]:-} "
old_warn=0
while IFS= read -r l; do
    case "$l" in
        ::error*) block "static check: ${l#*::}" ;;
        ::warning*)
            f=$(printf '%s' "$l" | sed -n 's/^::warning file=\([^,]*\),.*/\1/p')
            if [ -z "$f" ] || [[ "$changed_set" == *" $f "* ]]; then warn "static check: ${l#*::}"; else old_warn=$((old_warn + 1)); fi ;;
    esac
done <<< "$sc_out"
[ $old_warn -gt 0 ] && note "$old_warn static-check warnings on files this diff doesn't touch (already on the base)"
[ $sc_rc -ne 0 ] && ! printf '%s\n' "${blockers[@]:-}" | grep -q 'static check' && block "static checks failed (exit $sc_rc)"

# --- report
verdict=clean; rc=0
[ ${#warnings[@]} -gt 0 ] && { verdict=warnings; rc=1; }
[ ${#blockers[@]} -gt 0 ] && { verdict=BLOCKED; rc=2; }
echo "prereview: $verdict - ${#changed[@]} files changed, ${#deleted[@]} deleted, +$nadd -$ndel vs $(git rev-parse --short "$mb") ($base)"
for x in "${blockers[@]:-}"; do [ -n "$x" ] && echo "  BLOCK  $x"; done
n=0; for x in "${warnings[@]:-}"; do [ -n "$x" ] || continue; n=$((n + 1)); [ $n -le 20 ] && echo "  WARN   $x"; done
[ $n -gt 20 ] && echo "  WARN   ... $((n - 20)) more"
for x in "${notes[@]:-}"; do [ -n "$x" ] && echo "  NOTE   $x"; done
exit $rc
