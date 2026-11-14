# PostToolUse hook (Edit|Write): after a C++ source edit, rebuild the agent build and report compiler
# errors back to Claude (exit 2 + stderr). Quiet when the build is fine or the file isn't C++.
$ErrorActionPreference = 'SilentlyContinue'
$raw = [Console]::In.ReadToEnd()
try { $evt = $raw | ConvertFrom-Json } catch { exit 0 }
$file = $evt.tool_input.file_path
if (-not $file -or $file -notmatch '\.(cpp|cc|c|h|hpp)$') { exit 0 }

$repo = $env:CLAUDE_PROJECT_DIR
if (-not $repo) { $repo = (Get-Location).Path }
$build = Join-Path $repo 'build/check'
if (-not (Test-Path (Join-Path $build 'CMakeCache.txt'))) { exit 0 } # not configured here: skip (the capture skill configures it)

$env:PATH = "C:\Program Files\CMake\bin;C:\raylib\w64devkit\bin;$env:PATH"
$out = & cmake --build $build -j 8 2>&1
if ($LASTEXITCODE -eq 0) { exit 0 }

# Errors only (warnings are noise here), capped so the context stays small.
$errs = $out | Where-Object { $_ -match 'error:|undefined reference|Error [0-9]' } | Select-Object -First 25
if (-not $errs) { $errs = $out | Select-Object -Last 25 }
[Console]::Error.WriteLine("build-check: build/check failed after editing $([IO.Path]::GetFileName($file)):")
$errs | ForEach-Object { [Console]::Error.WriteLine($_) }
exit 2
