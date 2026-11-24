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
# (#1758) -O1: an agent build checks behaviour; timing on the PC needs a -O3 build (RAYLIBGAME_FAST_BUILD=OFF)
# (#1804) The preset's compilers (CMakeUserPresets.json "w64devkit"): w64devkit's own gcc/g++/make, not whatever cc/c++
# come first on PATH (MSYS2's do on this PC).
wk="C:/raylib/w64devkit/bin"
conf() { cmake -S "$wt" -B "$wt/build/check" -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DRAYLIBGAME_FAST_BUILD=ON -DGAME_OUTPUT_DIR="$out" \
    -DCMAKE_C_COMPILER="$wk/gcc.exe" -DCMAKE_CXX_COMPILER="$wk/g++.exe" -DCMAKE_MAKE_PROGRAM="$wk/mingw32-make.exe" > "$wt/build/configure.log" 2>&1; }
# (#1804) The antivirus scan of CMake's test programs can break the compiler ABI detection without failing the configure
# (CMAKE_<LANG>_ABI_COMPILED stays 0, no implicit link libraries): the link then puts the C libraries on the line and
# logger.cpp's __imp__assert (#1526) collides with msvcrt's. Detected here; the compiler files are deleted and redone.
abi_ok() {
    local f
    for f in "$wt"/build/check/CMakeFiles/*/CMakeCCompiler.cmake "$wt"/build/check/CMakeFiles/*/CMakeCXXCompiler.cmake; do
        [ -f "$f" ] || return 1
        grep -Eq 'set\(CMAKE_(C|CXX)_ABI_COMPILED (TRUE|1)\)' "$f" || return 1
    done
}
mkdir -p "$wt/build"
conf || { echo "setup: configure failed once (often the a.exe antivirus scan); retrying"; sleep 5; conf; } || { tail -20 "$wt/build/configure.log"; exit 1; }
for try in 1 2 3; do
    abi_ok && break
    echo "setup: compiler ABI detection failed (antivirus scan, #1804); redoing it ($try/3)"
    rm -rf "$wt"/build/check/CMakeFiles/[0-9]*.[0-9]*
    sleep 5
    conf || { sleep 5; conf; } || { tail -20 "$wt/build/configure.log"; exit 1; }
done
# Some worktrees lose the C++ test program every time (its .bin never appears): take the detection from another worktree
# with the same CMake and g++ versions whose detection worked (it names only the compiler, rewritten to ours), and keep it.
if ! abi_ok; then
    ver=$(basename "$(ls -d "$wt"/build/check/CMakeFiles/[0-9]*.[0-9]* | head -1)")
    gv=$("$wk/g++.exe" -dumpfullversion)
    for lang in C CXX; do
        mine="$wt/build/check/CMakeFiles/$ver/CMake${lang}Compiler.cmake"
        grep -Eq "set\(CMAKE_${lang}_ABI_COMPILED (TRUE|1)\)" "$mine" && continue
        donor=""
        while read -r o; do
            f="$o/build/check/CMakeFiles/$ver/CMake${lang}Compiler.cmake"
            [ "$o" != "$wt" ] && [ -f "$f" ] && grep -Eq "set\(CMAKE_${lang}_ABI_COMPILED (TRUE|1)\)" "$f" && grep -q "CMAKE_${lang}_COMPILER_VERSION \"$gv\"" "$f" && { donor="$f"; break; }
        done < <(git -C "$wt" worktree list --porcelain | sed -n 's/^worktree //p')
        [ -n "$donor" ] || { echo "setup: the $lang ABI detection fails and no other worktree has a good one (#1804)"; exit 1; }
        cc=$([ $lang = C ] && echo gcc.exe || echo g++.exe)
        sed "1s#.*#set(CMAKE_${lang}_COMPILER \"$wk/$cc\")#" "$donor" > "$mine"
        cp "$(dirname "$donor")/CMakeDetermineCompilerABI_${lang}.bin" "$(dirname "$mine")/" 2>/dev/null
        echo "setup: $lang compiler detection taken from $(dirname "$donor") (#1804)"
    done
    conf || { tail -20 "$wt/build/configure.log"; exit 1; }
    abi_ok || { echo "setup: the compiler ABI detection still fails (#1804)"; exit 1; }
fi
grep -q "GAME_OUTPUT_DIR:.*=$out" "$wt/build/check/CMakeCache.txt" || { echo "setup: GAME_OUTPUT_DIR is not $out - stop"; exit 1; }
echo "setup: configured build/check -> $out"
[ "$1" = "--no-build" ] && exit 0
bash "$wt/tools/slot.sh" wait-quiet || exit 1
powershell -NoProfile -Command "\$p = Start-Process cmake -ArgumentList '--build','$wt/build/check','-j','4','--target','game','game_tests' -NoNewWindow -PassThru -RedirectStandardOutput '$wt/build/build.log'; \$p.PriorityClass = 'BelowNormal'; \$p.WaitForExit(); exit \$p.ExitCode"
rc=$?
grep -E "error:|warning:" "$wt/build/build.log" | head -20
[ $rc -eq 0 ] && echo "setup: built $out/game.exe" || { echo "setup: build failed (build/build.log)"; exit 1; }
