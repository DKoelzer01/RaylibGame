# PostToolUse hook (Edit|Write): after a C++ source edit, rebuild the agent build and report compiler
# errors back to Claude (exit 2 + stderr). Quiet when the build is fine or the file isn't C++.
$ErrorActionPreference = 'SilentlyContinue'
$raw = [Console]::In.ReadToEnd()
try { $evt = $raw | ConvertFrom-Json } catch { exit 0 }
$file = $evt.tool_input.file_path
if (-not $file -or $file -notmatch '\.(cpp|cc|c|h|hpp)$') { exit 0 }

# Build the checkout that holds the edited file (#1197), not the session's start folder: walk up from the file to the
# first folder with both CMakeLists.txt and a .git entry (a worktree's .git is a file). No such folder: skip.
$repo = $null
$dir = Split-Path -Parent ([IO.Path]::GetFullPath($file))
while ($dir) {
    if ((Test-Path (Join-Path $dir 'CMakeLists.txt')) -and (Test-Path (Join-Path $dir '.git'))) { $repo = $dir; break }
    $parent = Split-Path -Parent $dir
    if ($parent -eq $dir) { break }
    $dir = $parent
}
if (-not $repo) { exit 0 }
$build = Join-Path $repo 'build/check'
if (-not (Test-Path (Join-Path $build 'CMakeCache.txt'))) { exit 0 } # not configured here: skip (the capture skill configures it)

# Someone else's quiet window (#1477, tools/slot.sh): don't build now; the next edit (or a manual build) catches up.
$slot = Join-Path $repo 'tools/slot.sh'
$bash = 'C:\Program Files\Git\bin\bash.exe'
if ((Test-Path $slot) -and (Test-Path $bash)) {
    Push-Location $repo
    $q = (& $bash $slot quiet-other 2>$null) -join ''
    $quiet = ($LASTEXITCODE -eq 0)
    Pop-Location
    if ($quiet) { exit 0 }
}

# Below-normal priority (#1197): cmake and the compilers inherit it, so a build never stutters a running game.
try { (Get-Process -Id $PID).PriorityClass = 'BelowNormal' } catch {}
$env:PATH = "C:\Program Files\CMake\bin;C:\raylib\w64devkit\bin;$env:PATH"
$out = & cmake --build $build -j 4 2>&1 # -j 4 (#1477): several sessions build on this PC
if ($LASTEXITCODE -eq 0) {
    # Warnings in the file just edited (a missed -Wformat once crashed a test, #1378): report them, capped.
    $leaf = [IO.Path]::GetFileName($file)
    $warns = $out | Where-Object { $_ -match 'warning:' -and $_ -match [regex]::Escape($leaf) } | Select-Object -First 10
    if ($warns) {
        [Console]::Error.WriteLine("build-check: build ok, but warnings in ${leaf}:")
        $warns | ForEach-Object { [Console]::Error.WriteLine($_) }
        exit 2
    }
    exit 0
}

# Errors only (warnings are noise here), capped so the context stays small.
$errs = $out | Where-Object { $_ -match 'error:|undefined reference|Error [0-9]' } | Select-Object -First 25
if (-not $errs) { $errs = $out | Select-Object -Last 25 }
[Console]::Error.WriteLine("build-check: build/check failed after editing $([IO.Path]::GetFileName($file)):")
$errs | ForEach-Object { [Console]::Error.WriteLine($_) }
exit 2
