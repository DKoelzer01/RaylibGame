#!/bin/bash
# busybox (w64devkit on PATH) is not bash: re-run under Git's bash (#1477)
if [ -z "${BASH_VERSINFO:-}" ]; then for b in "C:/Program Files/Git/bin/bash.exe" "/c/Program Files/Git/bin/bash.exe" "/usr/bin/bash"; do [ -x "$b" ] && exec "$b" "$0" "$@"; done; echo "needs bash" >&2; exit 2; fi
# Build and test on deweyserver instead of this PC (user 2026-11-22): for compiling and headless checks while the user
# playtests, or whenever this PC should stay free. Pushes the current commit (HEAD, not the work tree's uncommitted
# changes) to a scratch ref, then deweyserver builds it (Linux headless,
# no GNS, its own ccache) and runs the command there with $GAME / $TESTS set. Nothing runs on this PC.
#
#   tools/sandbox.sh                          build only (compile errors + warning count)
#   tools/sandbox.sh '$TESTS --gtest_brief=1'                     unit tests
#   tools/sandbox.sh '$GAME --seed 42 --bench 1 --test-terrain-gpu'  CPU = GPU terrain (bit-exact)
#   tools/sandbox.sh '$GAME --seed 42 --cave-stats 800'           cave counts + the density hash
#   tools/sandbox.sh '$GAME --scenario tag:caves --scenario-no-record'
# Single-quote the command so $GAME stays for the server. One sandbox job at a time; ~5-15 min with a warm cache.
# Visual checks (captures, screenshots) still need this PC and a machine slot.
set -u
br=$(git rev-parse --abbrev-ref HEAD)
[ "$br" = HEAD ] && { echo "sandbox: check out a branch first (detached HEAD)"; exit 2; }
# (#2127) The commit goes to a scratch ref, not a branch: pushing the branch itself left stray branches on GitHub (a local
# name nobody meant to publish). One ref per commit, so two sessions never overwrite each other's (and nothing is
# force-pushed); it is deleted when the job ends. Push the branch yourself when you want it on GitHub.
sha=$(git rev-parse HEAD)
ref="refs/sandbox/$br-${sha:0:8}"
git push -q origin "HEAD:$ref" 2>&1 | grep -v '^remote:'
git ls-remote --exit-code origin "$ref" > /dev/null || { echo "sandbox: push of ${sha:0:8} to $ref failed"; exit 2; }
trap 'git push -q origin ":$ref" 2>&1 | grep -v "^remote:"' EXIT
trap 'exit 130' INT TERM # (a timeout's TERM too: the EXIT trap still deletes the ref)
ssh="C:/Windows/System32/OpenSSH/ssh.exe"; [ -x "$ssh" ] || ssh=ssh
# (#1769) Linked shader programs in the sandbox's persistent /cache (the job's work tree is cleaned each time): ~16 s of
# compiles per game process otherwise. A command can still set RAYLIBGAME_SHADER_CACHE itself (0: a cold start).
q=${1:-}; [ -n "$q" ] && q="export RAYLIBGAME_SHADER_CACHE=\${RAYLIBGAME_SHADER_CACHE:-/cache/shaders}; $q"
# the command goes to the server as one argument, single-quoted for its shell
q=${q//\'/\'\\\'\'}
"$ssh" deweyserver "~/rg-sandbox/rb.sh '$ref' -- '$q'"
