#!/bin/bash
# busybox (w64devkit on PATH) is not bash: re-run under Git's bash (#1477)
if [ -z "${BASH_VERSINFO:-}" ]; then for b in "C:/Program Files/Git/bin/bash.exe" "/c/Program Files/Git/bin/bash.exe" "/usr/bin/bash"; do [ -x "$b" ] && exec "$b" "$0" "$@"; done; echo "needs bash" >&2; exit 2; fi
# Sets up this worktree's agent build (#1482): copies the gitignored CMakeUserPresets.json from the main checkout,
# configures build/check with its own GAME_OUTPUT_DIR (never the user's root game.exe), retries once when the first
# configure fails on the antivirus scan of a.exe, then builds game and game_tests at low priority with -j 4.
#   setup.sh [--no-build]
export PATH="/c/Program Files/CMake/bin:/c/raylib/w64devkit/bin:$PATH"
wt=$(git rev-parse --show-toplevel) || exit 1
main=$(dirname "$(git rev-parse --path-format=absolute --git-common-dir)")
[ "$wt" = "$main" ] && { echo "setup: this is the main checkout, not a worktree; refusing (its game.exe is the user's)"; exit 1; }
[ -f "$wt/CMakeUserPresets.json" ] || cp "$main/CMakeUserPresets.json" "$wt/" 2>/dev/null && echo "setup: presets copied"
out="$wt/build/check/bin"
conf() { cmake -S "$wt" -B "$wt/build/check" -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DGAME_OUTPUT_DIR="$out" > "$wt/build/configure.log" 2>&1; }
mkdir -p "$wt/build"
conf || { echo "setup: configure failed once (often the a.exe antivirus scan); retrying"; sleep 5; conf; } || { tail -20 "$wt/build/configure.log"; exit 1; }
grep -q "GAME_OUTPUT_DIR:.*=$out" "$wt/build/check/CMakeCache.txt" || { echo "setup: GAME_OUTPUT_DIR is not $out - stop"; exit 1; }
echo "setup: configured build/check -> $out"
[ "$1" = "--no-build" ] && exit 0
bash "$wt/tools/slot.sh" wait-quiet || exit 1
powershell -NoProfile -Command "\$p = Start-Process cmake -ArgumentList '--build','$wt/build/check','-j','4','--target','game','game_tests' -NoNewWindow -PassThru -RedirectStandardOutput '$wt/build/build.log'; \$p.PriorityClass = 'BelowNormal'; \$p.WaitForExit(); exit \$p.ExitCode"
rc=$?
grep -E "error:|warning:" "$wt/build/build.log" | head -20
[ $rc -eq 0 ] && echo "setup: built $out/game.exe" || { echo "setup: build failed (build/build.log)"; exit 1; }
