#!/bin/bash
# Markdown summary of a tools/regression.sh log (stdout of the run, with its "=== <test>" and "--- N s" lines) (#1028):
# the totals line(s), then one row per test with its verdict and duration, failures first.
#   tools/ci/regression_summary.sh regression.log >> "$GITHUB_STEP_SUMMARY"
set -u
log="${1:?usage: regression_summary.sh <log>}"
grep "^=== regression (" "$log" | sed 's/^=== /**/; s/$/**/' || echo "**regression: no totals line (the run did not finish)**"
echo
awk '
BEGIN { n = 0 }
function flush() {
    if (name == "") return
    v = fail ? "FAIL" : pass ? "PASS" : skip ? "SKIP" : "-"
    rows[n] = sprintf("| %s | %s | %s |", v, secs, name); failed[n++] = fail
    total += secs
    name = ""
}
/^=== regression \(/ { flush(); done = 1; next } # the totals of a set, then its FAIL lines again: not a test
/^=== / { flush(); done = 0; name = substr($0, 5); gsub(/\|/, "\\|", name); pass = fail = skip = 0; secs = "?"; next }
done { next }
/^--- [0-9]+ s$/ { secs = $2; next }
/FAIL/ { fail = 1 } /PASS/ { pass = 1 } /SKIP/ { skip = 1 }
END {
    flush()
    printf "%d tests, %.0f min in the tests\n\n| verdict | s | test |\n| --- | ---: | --- |\n", n, total / 60
    for (i = 0; i < n; ++i) if (failed[i]) print rows[i]
    for (i = 0; i < n; ++i) if (!failed[i]) print rows[i]
}' "$log"
if grep -q "FAIL" "$log"; then
    echo; echo "<details><summary>FAIL lines</summary>"; echo; echo '```'
    grep "FAIL" "$log" | grep -v "^=== regression" | sort -u | cut -c1-400
    echo '```'; echo "</details>"
fi
