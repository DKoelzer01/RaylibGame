# Tests for machine-guard-classify.ps1 (#1864). Run: powershell -NoProfile -File .claude/hooks/machine-guard.tests.ps1
# Exit code = number of failures.
. (Join-Path $PSScriptRoot 'machine-guard-classify.ps1')

$nl = "`n"
$cases = @(
    # --- runs (must stay caught)
    @('run',   'build/check/bin/game.exe --seed 42 --bench 3 --test-pond'),
    @('run',   './game.exe --seed 1'),
    @('run',   'cd build/check/bin && ./game --seed 1'),
    @('run',   'export PATH="/c/Program Files/CMake/bin:$PATH" && cmake --build build/check -j 4 > b.log && build/check/bin/game.exe --seed 42'),
    @('run',   "timeout 900 bash -c 'cd x; build/check/bin/game_tests.exe --gtest_brief=1'"),
    @('run',   '& "C:\Users\dewey\RaylibGame-lod\build\check\bin\game.exe" --seed 42'),
    @('run',   'Start-Process build\check\bin\game.exe -ArgumentList "--seed 42"'),
    @('run',   'nice -n 10 build/check/bin/game_tests --gtest_filter=Belt*'),
    @('run',   'RAYLIBGAME_FRAMES=1 build/check/bin/game.exe --seed 42'),
    @('run',   'ctest --test-dir build/check'),
    @('run',   "bash <<'EOF'${nl}build/check/bin/game.exe --seed 42${nl}EOF"),
    @('run',   'cat list.txt | build/check/bin/game.exe --stdin'),
    # --- builds (must stay caught)
    @('build', 'cmake --build build/check -j 4'),
    @('build', '"C:\Program Files\CMake\bin\cmake.exe" --build build/check'),
    @('build', 'cd src && g++ -c foo.cpp'),
    @('build', "bash -c 'mingw32-make -C build/check'"),
    @('build', 'build/check/bin/shader_tool.exe x.glsl'),
    # --- false positives (#1864): text, read-only and remote commands
    @('',      'cd /c/x && gh issue create --title "game crashes" --body "the game.exe crashes in game_tests; make sure"'),
    @('',      'export PATH="/c/Program Files/CMake/bin:$PATH"; gh pr view 1794'),
    @('',      "git commit -m `"fix: game.exe hang (#12)`""),
    @('',      "cat > build/checklist.md <<'EOF'${nl}Run game.exe and check the map.${nl}make sure to try game_tests${nl}EOF"),
    @('',      "git commit -m @'${nl}game.exe no longer hangs${nl}'@"),
    @('',      'sha256sum build/playtest/bin/game.exe'),
    @('',      'ls -la build/playtest/bin/game.exe && stat build/playtest/bin/game.exe'),
    @('',      "ssh deweyserver 'cd /cache/out && ./game_tests --gtest_brief=1'"),
    @('',      "tools/sandbox.sh '`$TESTS --gtest_brief=1'"),
    @('',      "bash tools/sandbox.sh '`$GAME --seed 42 --bench 1 --test-terrain-gpu'"),
    @('',      'gh api repos/o/r/issues/5/comments -f body="the game.exe build"'),
    @('',      'Get-Process game -ErrorAction SilentlyContinue'),
    @('',      'taskkill /IM game.exe /F'),
    @('',      'cmake --build build/check --target game ; echo done') # a build, checked below
)
$cases[-1][0] = 'build'

$fail = 0
foreach ($c in $cases) {
    $got = Get-MachineGuardKind $c[1]
    if ($got -ne $c[0]) {
        $fail++
        Write-Output ("FAIL want '{0}' got '{1}': {2}" -f $c[0], $got, ($c[1] -replace "`n", '\n'))
    }
}
Write-Output ("{0} cases, {1} failed" -f $cases.Count, $fail)
exit $fail
