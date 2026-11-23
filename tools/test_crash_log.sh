#!/bin/bash
# Crash entries in log.txt (#1526, #1527): each --test-crash kind runs in its own scratch folder (it ends before any
# window opens) and its log.txt must say what happened. Then the stop from outside: a run sleeping in --test-crash wait,
# a second run starting in the same folder (log.txt is now the second's), and the first stopped by SIGINT: its line must
# land in log.prev.txt, never inside the second run's log.txt.
# Usage: tools/test_crash_log.sh [EXE]   (default build/check/bin/game.exe). Prints "=== test_crash_log: PASS|FAIL".
set -u
exe="$(cd "$(dirname "${1:-build/check/bin/game.exe}")" && pwd)/$(basename "${1:-build/check/bin/game.exe}")"
tmp="$(mktemp -d)"
win=0; case "$(uname -s)" in MINGW*|MSYS*|CYGWIN*) win=1 ;; esac
export RAYLIBGAME_FATAL_BOX=0 # no error box, whatever the console
ok=1
chk() { if grep -qE "$2" "$tmp/$1/log.txt" 2>/dev/null; then echo "  ok: $1: $3"; else echo "  MISSING: $1: $3"; ok=0; fi; }
run() { mkdir -p "$tmp/$1"; (cd "$tmp/$1" && timeout 60 "$exe" --test-crash "$1" > out.txt 2> err.txt); echo "  $1: exit $?"; }
for k in segv thread-segv assert terminate fatal; do run $k; done
if [ $win = 1 ]; then
    chk segv "FATAL crash: exception 0xc0000005 \(access violation\) at .* = game.exe\+0x" "the exception, its module and offset"
    chk segv "FATAL crash: writing address 0+\b" "the address it wrote"
    chk segv "FATAL stack: game.exe\+0x" "a stack"
    chk segv "FATAL minidump: saves/logs/crash-.*\.dmp" "a minidump"
    ls "$tmp"/segv/saves/logs/crash-*.dmp > /dev/null 2>&1 && echo "  ok: segv: the .dmp exists" || { echo "  MISSING: segv: the .dmp file"; ok=0; }
    chk thread-segv "FATAL crash: exception 0xc0000005 .*thread [0-9]+ 'crasher'" "the crashing thread's name"
    chk terminate "FATAL uncaught exception: std::terminate with no exception" "std::terminate (a joinable std::thread destroyed)"
else
    chk segv "FATAL crash: signal 11" "the signal"
    chk thread-segv "FATAL crash: signal 11 .*thread 'crasher'" "the crashing thread's name"
    chk terminate "FATAL crash: signal 6" "std::terminate's abort"
fi
chk segv "build [^ ]+ [0-9a-f]{7}" "the build"
chk assert "FATAL assertion failed: --test-crash assert, file .*logger\.cpp, line [0-9]+" "the assertion's text, file and line"
chk assert "FATAL abort\(\) after the failed assertion above" "abort() names the assertion"
chk fatal "FATAL .*\[Log\] --test-crash fatal" "the LOG_FATAL line"
chk fatal "FATAL abort\(\) after the FATAL line above" "abort() names the FATAL line"

# (#1527) stopped by Ctrl+C with a thread still running (the factory's was: std::terminate, the error box, a hang), while
# a newer run owns log.txt. (Git bash's and busybox's timeout end a native process hard, so the run sends Ctrl+C to itself.)
mkdir -p "$tmp/stop"
(cd "$tmp/stop" && RAYLIBGAME_TEST_CRASH_DELAY=4 exec "$exe" --test-crash ctrl-c > a.txt 2>&1) & a=$!
for _ in $(seq 1 50); do grep -q "test-crash ctrl-c" "$tmp/stop/log.txt" 2>/dev/null && break; sleep 0.2; done
(cd "$tmp/stop" && timeout 60 "$exe" --test-crash fatal > b.txt 2>&1)
end=$((SECONDS + 20)); while kill -0 $a 2>/dev/null && [ $SECONDS -lt $end ]; do sleep 0.2; done
if kill -0 $a 2>/dev/null; then echo "  MISSING: stop: the run ended on its Ctrl+C (still alive after 20 s: killed)"; ok=0; kill -9 $a 2>/dev/null; else echo "  ok: stop: the run ended on its Ctrl+C"; fi
wait $a 2>/dev/null
grep -q "FATAL" "$tmp/stop/log.prev.txt" && { echo "  MISSING: stop: no FATAL line (std::terminate) from the stopped run"; ok=0; } || echo "  ok: stop: no FATAL line (std::terminate) from the stopped run"
grep -q "stopped from outside" "$tmp/stop/log.prev.txt" && echo "  ok: stop: its line went to its own log (log.prev.txt)" || { echo "  MISSING: stop: 'stopped from outside' in log.prev.txt"; ok=0; }
grep -q "stopped from outside\|test-crash wait" "$tmp/stop/log.txt" && { echo "  MISSING: stop: the newer run's log.txt has none of the old run's lines"; ok=0; } || echo "  ok: stop: the newer run's log.txt has none of the old run's lines"
grep -q "after the FATAL line above" "$tmp/stop/log.txt" && echo "  ok: stop: the newer run's log is whole" || { echo "  MISSING: stop: the newer run's own lines"; ok=0; }
echo "=== test_crash_log: $([ $ok = 1 ] && echo PASS || echo FAIL) ($tmp)"
[ $ok = 1 ]
