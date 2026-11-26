#!/bin/bash
# busybox (w64devkit on PATH) is not bash: re-run under Git's bash (#1477)
if [ -z "${BASH_VERSINFO:-}" ]; then for b in "C:/Program Files/Git/bin/bash.exe" "/c/Program Files/Git/bin/bash.exe" "/usr/bin/bash"; do [ -x "$b" ] && exec "$b" "$0" "$@"; done; echo "needs bash" >&2; exit 2; fi
# Compact output for builds, tests, benches and sandbox runs (#1876). Everything a session runs lands in its context and
# is re-read on every later turn, so print only what decides the next step and keep the full log in a file.
#
#   tools/quiet.sh [--label L] [--tail N] [--max N] [--all] -- <command...>
#
#   tools/quiet.sh -- cmake --build build/check -j 4
#   tools/quiet.sh --label pump -- build/check/bin/game.exe --seed 42 --bench 3 --test-pump
#   tools/quiet.sh -- tools/sandbox.sh '$TESTS --gtest_brief=1'
#   tools/quiet.sh -- build/check/bin/game.exe --scenario tag:fluid --scenario-no-record
#
# Prints: the exit code, the time and the log path; compiler and linker errors (the first few with their context lines)
# and warnings (deduplicated, capped); every FAIL verdict (bench "-> FAIL", "FAILED", scenario FAIL, gtest [  FAILED  ]);
# PASS verdicts as a count (all of them with --all, or when there are few); bench and scenario summary lines; and, when
# the command failed without any of those, the last --tail lines. Exits with the command's exit code, so gating
# (`tools/quiet.sh -- cmake --build ... && ...`) still works. Logs: build/logs/<label>-<time>.log, pruned after 7 days.
# Read the log (Read with an offset, or grep it) only when the summary isn't enough.
set -u
label="" tail_n=15 max_n=30 all=0
while [ $# -gt 0 ]; do
    case "$1" in
        --label) label="$2"; shift 2 ;;
        --tail) tail_n="$2"; shift 2 ;;
        --max) max_n="$2"; shift 2 ;;
        --all) all=1; shift ;;
        --) shift; break ;;
        -h|--help) sed -n '4,20p' "$0"; exit 0 ;;
        *) break ;;
    esac
done
[ $# -gt 0 ] || { echo "quiet: no command (tools/quiet.sh [--label L] -- <command...>)" >&2; exit 2; }

repo="$(cd "$(dirname "$0")/.." && pwd)"
logdir="$repo/build/logs"
mkdir -p "$logdir"
find "$logdir" -name '*.log' -mtime +7 -delete 2>/dev/null
if [ -z "$label" ]; then
    label=$(basename "$1"); label=${label%.exe}; label=${label%.sh}
    # a game run is named by its first --test-*/--scenario argument, a build by its target dir
    for a in "$@"; do case "$a" in --test-*|--bench-*) label="$label-${a#--}"; break ;; esac; done
fi
label=$(printf '%s' "$label" | tr -c 'A-Za-z0-9._-' '_')
log="$logdir/$label-$(date +%Y%m%d-%H%M%S).log"

start=$(date +%s)
"$@" > "$log" 2>&1
rc=$?
secs=$(( $(date +%s) - start ))
lines=$(wc -l < "$log" | tr -d ' ')

# Windows paths for the log, so the Read tool can open it as printed
shown="$log"
command -v cygpath > /dev/null 2>&1 && shown=$(cygpath -m "$log")
echo "quiet: $label exit=$rc in ${secs}s, $lines lines -> $shown"

LC_ALL=C awk -v max="$max_n" -v all="$all" -v rc="$rc" -v tailn="$tail_n" '
function trim(s) { sub(/\r$/, "", s); return s }
{ line = trim($0); L[NR] = line }
# compiler and linker problems (gcc/clang "file:line:col: error:", ld "undefined reference", make/ninja stops)
line ~ /(: |^)(fatal )?error: |undefined reference to|multiple definition of|collect2: error|ld returned/ {
    nerr++; if (nerr <= 3) { errat[nerr] = NR } ; next_ctx = 0 }
line ~ /: warning: / {
    w = line; sub(/^.*[\/\\]/, "", w)                  # shorten the path to the file name
    if (!(w in seenw)) { seenw[w] = 1; nw++; if (nw <= 10) warn[nw] = w } }
# verdicts
line ~ /-> FAIL|[^A-Z]FAILED|^FAILED|\[  FAILED  \]|scenario: .* FAIL|: FAIL( |$|:)|^FAIL/ && line !~ /0 failed/ {
    nfail++; if (nfail <= max) fail[nfail] = line; next }
line ~ /-> PASS|\[  PASSED  \]|scenario: === .* PASS|: PASS( |$|:)|[^A-Z]PASSED/ {
    npass++; if (npass <= max) pass[npass] = line; next }
line ~ /^scenario: [0-9]+ passed|^bench: .*(ms|fps|FPS|p50|p95|p99)|^\[==========\]|tests? passed|^summary:|^SUMMARY/ {
    nsum++; if (nsum <= 12) summ[nsum] = line }
END {
    if (nerr) {
        printf "errors: %d (first %d):\n", nerr, (nerr < 3 ? nerr : 3)
        for (i = 1; i <= 3 && i <= nerr; i++) for (j = errat[i]; j <= errat[i] + 2 && j <= NR; j++) {
            if (j > errat[i] && L[j] ~ /(error|warning|note): /) break
            print "  " substr(L[j], 1, 300) }
    }
    if (nw) { printf "warnings: %d distinct%s\n", nw, (nw > 10 ? " (first 10)" : ""); for (i = 1; i <= nw && i <= 10; i++) print "  " substr(warn[i], 1, 240) }
    if (nfail) { printf "FAIL: %d\n", nfail; for (i = 1; i <= nfail && i <= max; i++) print "  " substr(fail[i], 1, 300); if (nfail > max) printf "  ... %d more in the log\n", nfail - max }
    if (npass) {
        if (all || npass <= 5) { printf "PASS: %d\n", npass; for (i = 1; i <= npass && i <= max; i++) print "  " substr(pass[i], 1, 200) }
        else printf "PASS: %d verdict lines (--all lists them)\n", npass }
    for (i = 1; i <= nsum; i++) print "  " substr(summ[i], 1, 240)
    if (rc != 0 && !nerr && !nfail) {
        printf "last %d lines:\n", tailn
        for (i = (NR > tailn ? NR - tailn + 1 : 1); i <= NR; i++) print "  " substr(L[i], 1, 300) }
}' "$log"
exit $rc
