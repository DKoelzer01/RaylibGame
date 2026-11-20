#!/bin/bash
# Markdown summary of a tools/regression.sh log (stdout of the run, with its "=== <test>" and "--- N s" lines) (#1028):
# the totals line(s), then one row per test with its verdict and duration: new failures first, then known ones, then
# the rest. KNOWN (optional, tools/regression_known.txt) lists known failures as "<args> | <issue>".
#   tools/ci/regression_summary.sh regression.log [KNOWN] >> "$GITHUB_STEP_SUMMARY"
# Exit code: 1 when a test failed that KNOWN doesn't list, else 0.
set -u
log="${1:?usage: regression_summary.sh <log> [known-failures file]}"
known="${2:-/dev/null}"
grep "^=== regression (" "$log" | sed 's/^=== /**/; s/$/**/' || echo "**regression: no totals line (the run did not finish)**"
echo
awk -v knownfile="$known" '
BEGIN {
    n = 0
    while ((getline line < knownfile) > 0) {
        sub(/\r$/, "", line) # (a Windows checkout has CRLF)
        if (line ~ /^#/ || line !~ / \| /) continue
        k = index(line, " | ")
        issue[substr(line, 1, k - 1)] = substr(line, k + 3)
    }
}
function flush() {
    if (name == "") return
    v = fail ? "FAIL" : pass ? "PASS" : skip ? "SKIP" : "-"
    cls = !fail ? 2 : (name in issue) ? 1 : 0 # 0 new failure, 1 known failure, 2 the rest
    if (cls == 1) v = "FAIL (known " issue[name] ")"
    if (cls == 0) { v = "**FAIL (new)**"; ++newFails }
    if (name in issue) { seen[name] = 1; if (fail) failedOnce[name] = 1 }
    shown = name; gsub(/\|/, "\\|", shown)
    rows[n] = sprintf("| %s | %s | %s |", v, secs, shown); cl[n++] = cls
    total += secs
    name = ""
}
/^=== regression \(/ { flush(); done = 1; next } # the totals of a set, then its FAIL lines again: not a test
/^=== / { flush(); done = 0; name = substr($0, 5); pass = fail = skip = 0; secs = "?"; next }
done { next }
/^--- [0-9]+ s$/ { secs = $2; next }
/FAIL/ { fail = 1 } /PASS/ { pass = 1 } /SKIP/ { skip = 1 }
END {
    flush()
    known = 0; for (i = 0; i < n; ++i) if (cl[i] == 1) ++known
    # listed, run, never failed (every run of those arguments), and not marked flaky: maybe fixed
    for (t in seen) if (!(t in failedOnce) && issue[t] !~ /flaky/) fixed = fixed "\n- `" t "` (" issue[t] ")"
    printf "%d tests, %.0f min in the tests; **%d new failures**, %d known failures (tools/regression_known.txt)\n", n, total / 60, newFails, known
    if (fixed != "") printf "\nKnown failures that passed (fixed? remove them from tools/regression_known.txt):%s\n", fixed
    printf "\n| verdict | s | test |\n| --- | ---: | --- |\n"
    for (c = 0; c <= 2; ++c) for (i = 0; i < n; ++i) if (cl[i] == c) print rows[i]
    exit newFails > 0
}' "$log"
rc=$?
if grep -q "FAIL" "$log"; then
    echo; echo "<details><summary>FAIL lines</summary>"; echo; echo '```'
    grep "FAIL" "$log" | grep -v "^=== regression" | sort -u | cut -c1-400
    echo '```'; echo "</details>"
fi
exit $rc
