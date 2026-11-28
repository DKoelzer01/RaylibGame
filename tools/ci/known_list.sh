#!/bin/bash
# (#1946) The known-failures list to match a run's failures against: this checkout's tools/ci/known_failures.txt plus the
# integration branch's (origin/$INTEGRATION, default staging). A branch cut before an entry was added would otherwise tag
# that failure NEW, and a PR's own new entry still counts. Prints the merged lines (CR stripped, duplicates dropped).
# The fetch failing (no network, no remote) is silent: then only the checkout's list. KNOWN_FAILURES_LOCAL_ONLY=1 skips it.
root=$(cd "$(dirname "$0")/../.." && pwd)
b="${INTEGRATION:-staging}"
{
  tr -d '\r' < "$root/tools/ci/known_failures.txt" 2>/dev/null
  if [ "${KNOWN_FAILURES_LOCAL_ONLY:-0}" != 1 ]; then
    # A shallow clone (the runner's checkout) fetches one commit; a full clone must not become shallow.
    depth=(); [ "$(git -C "$root" rev-parse --is-shallow-repository 2>/dev/null)" = true ] && depth=(--depth=1)
    git -C "$root" fetch --quiet "${depth[@]}" origin "$b" > /dev/null 2>&1
    git -C "$root" show "origin/$b:tools/ci/known_failures.txt" 2>/dev/null | tr -d '\r'
  fi
} | awk '!seen[$0]++'
