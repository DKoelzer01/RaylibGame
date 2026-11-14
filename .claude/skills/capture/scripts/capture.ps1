# Capture skill: build the agent build, run the game with in-game recording (#375), report the outputs.
param(
    [Parameter(Mandatory = $true)][string]$Name,
    [string]$GameArgs = "--seed 42 --bench 10",
    [int]$Frames = 60,
    [int]$Every = 5,
    [switch]$NoHud,
    [string]$BuildDir = "build/check",
    [int]$TimeoutSec = 300
)
$ErrorActionPreference = 'Stop'

# Repo root: the skill lives in <repo>/.claude/skills/capture/scripts.
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..\..')).Path
Set-Location $repo
$env:PATH = "C:\Program Files\CMake\bin;C:\raylib\w64devkit\bin;$env:PATH"

# 1. Build (never test a stale exe).
if (-not (Test-Path "$BuildDir/CMakeCache.txt")) {
    Write-Host "capture: configuring $BuildDir"
    cmake -S . -B $BuildDir -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release "-DGAME_OUTPUT_DIR=$repo/$BuildDir/bin" | Out-Host
    if ($LASTEXITCODE -ne 0) { Write-Error "capture: configure failed"; exit 1 }
}
cmake --build $BuildDir -j 8 2>&1 | Select-Object -Last 15 | Out-Host
if ($LASTEXITCODE -ne 0) { Write-Error "capture: build failed - fix it before capturing"; exit 1 }
$exe = Join-Path $repo "$BuildDir/bin/game.exe"
if (-not (Test-Path $exe)) { Write-Error "capture: $exe not found"; exit 1 }

# 2. Output dir.
$stamp = Get-Date -Format 'yyyyMMdd-HHmmss'
$out = Join-Path $repo "captures/$Name-$stamp"
New-Item -ItemType Directory -Force $out | Out-Null

# 3. Does this build have the recording flags (#375)?
$hasRecord = (Select-String -Path $exe -Pattern '--record' -SimpleMatch -Quiet)
$argList = @($GameArgs -split ' ' | Where-Object { $_ })
if ($hasRecord) {
    $argList += @('--record', $out, '--frames', "$Frames", '--record-every', "$Every", '--gif', '--sheet')
    if ($NoHud) { $argList += @('--record-hud', '0') }
} else {
    Write-Warning "capture: this build has no --record yet (#375); falling back to the bench screenshot."
}

# 4. Run from the repo root (resources/ and assets/ resolve there), log to the capture dir.
Write-Host "capture: $exe $($argList -join ' ')"
$log = Join-Path $out 'game.log'
$p = Start-Process -FilePath $exe -ArgumentList $argList -WorkingDirectory $repo -RedirectStandardOutput $log -RedirectStandardError "$log.err" -PassThru
if (-not $p.WaitForExit($TimeoutSec * 1000)) {
    Stop-Process -Id $p.Id -Force
    Write-Error "capture: timed out after $TimeoutSec s (see $log)"; exit 1
}
Get-Content $log -Tail 10 | Out-Host

# 5. Report.
if (-not $hasRecord) {
    Get-ChildItem $repo -Filter 'bench*.png' | Where-Object { $_.LastWriteTime -gt $p.StartTime } | ForEach-Object { Copy-Item $_.FullName $out }
}
$frames = @(Get-ChildItem $out -Filter 'frame_*.png' -ErrorAction SilentlyContinue)
Write-Host "capture: dir     $out"
Write-Host "capture: frames  $($frames.Count)"
foreach ($f in 'sheet.png', 'capture.gif') { if (Test-Path "$out/$f") { Write-Host "capture: $($f.PadRight(7)) $out/$f" } }
Get-ChildItem $out -Filter 'bench*.png' -ErrorAction SilentlyContinue | ForEach-Object { Write-Host "capture: shot    $($_.FullName)" }
exit $p.ExitCode
