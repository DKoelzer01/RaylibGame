# Capture skill: build the agent build, run the game with in-game recording (#375, src/space/recorder.*),
# and report the outputs from the game's "record:" summary lines.
param(
    [Parameter(Mandatory = $true)][string]$Name,      # capture name (the game files it under <root>/<date>/<name>-<time>/)
    [string]$GameArgs = "--seed 42 --bench 10",       # the scenario
    [int]$Frames = 120,
    [int]$Every = 1,
    [double]$Start = 0,                               # extra seconds after loading + WorldIdle before frame 0
    [switch]$NoHud,
    [switch]$NoGif,
    [string]$Compare = "",                            # an earlier capture's name or path: per-frame pixel diff
    [switch]$Continue,                                # keep running after the last frame (default: the game exits)
    [string]$BuildDir = "build/check",
    [int]$TimeoutSec = 600
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
    if ($LASTEXITCODE -ne 0) { Write-Error "capture: configure failed (rerun once: the first configure can fail on an AV scan)"; exit 1 }
}
# Merge stderr inside cmd, not PowerShell: PS 5.1 turns each stderr line (any compiler warning) into a
# NativeCommandError, which $ErrorActionPreference = 'Stop' makes fatal (#492). Success = the exit code only.
cmd /c "cmake --build `"$BuildDir`" -j 8 2>&1" | Select-Object -Last 15 | Out-Host
if ($LASTEXITCODE -ne 0) { Write-Error "capture: build failed - fix it before capturing"; exit 1 }
$exe = Join-Path $repo "$BuildDir/bin/game.exe"
if (-not (Test-Path $exe)) { Write-Error "capture: $exe not found"; exit 1 }

# 2. Does this build have the recorder (#375)?
$hasRecord = (Select-String -Path $exe -Pattern '--record-every' -SimpleMatch -Quiet)
$argList = @($GameArgs -split ' ' | Where-Object { $_ })
$capName = "$Name-$(Get-Date -Format 'HHmmss')"
if ($hasRecord) {
    $argList += @('--record', $capName, '--frames', "$Frames", '--record-every', "$Every", '--sheet')
    if (-not $NoGif) { $argList += '--gif' }
    if ($Start -gt 0) { $argList += @('--record-start', "$Start") }
    if ($NoHud) { $argList += @('--record-hud', '0') }
    if ($Continue) { $argList += '--record-continue' }
    if ($Compare) { $argList += @('--record-compare', $Compare) }
} else {
    Write-Warning "capture: this build has no recorder yet (#375 not merged into this branch); falling back to the bench screenshots."
}

# 3. Run from the repo root (resources/ and assets/ resolve there); log to a temp file.
$logDir = Join-Path $env:TEMP "raylib-capture"
New-Item -ItemType Directory -Force $logDir | Out-Null
$log = Join-Path $logDir "$capName.log"
Write-Host "capture: $exe $($argList -join ' ')"
$p = Start-Process -FilePath $exe -ArgumentList $argList -WorkingDirectory $repo -RedirectStandardOutput $log -RedirectStandardError "$log.err" -PassThru
if (-not $p.WaitForExit($TimeoutSec * 1000)) {
    Stop-Process -Id $p.Id -Force
    Write-Error "capture: timed out after $TimeoutSec s (log: $log)"; exit 1
}

# 4. Report: the recorder's summary lines carry the resolved paths (UNC as //DEWEYSERVER/...).
$rec = @(Select-String -Path $log -Pattern '^record:' | ForEach-Object { $_.Line })
if ($rec.Count -gt 0) {
    $rec | ForEach-Object { Write-Host $_ }
    $m = [regex]::Match($rec[0], '^record:\s+(\S+)\s+frames=(\d+)\s+gif=(\S+)\s+sheet=(\S+)')
    if ($m.Success) {
        $dir = $m.Groups[1].Value -replace '/', '\'
        Write-Host "capture: dir    $dir"
        Write-Host "capture: frames $($m.Groups[2].Value)  (frame_0000.png ...)"
        if ($m.Groups[4].Value -ne '-') { Write-Host "capture: sheet  $($m.Groups[4].Value -replace '/', '\')" }
        if ($m.Groups[3].Value -ne '-') { Write-Host "capture: gif    $($m.Groups[3].Value -replace '/', '\')" }
        # A Read of a UNC path always asks for permission: on deweyserver, look with media_get instead (SKILL.md).
        $mm = [regex]::Match($m.Groups[1].Value, '(?i)claude-media/(.+?)/?$')
        if ($mm.Success) { Write-Host "capture: look   mcp__deweyserver__media_get(path=""$($mm.Groups[1].Value)"")  (frame=N for one frame; not Read)" }
    }
} else {
    Get-Content $log -Tail 12 | Out-Host
    if ($hasRecord) { Write-Warning "capture: no 'record:' line - the run ended before the capture finished (log: $log)" }
    $shots = Get-ChildItem $repo -Filter 'bench*.png' | Where-Object { $_.LastWriteTime -gt $p.StartTime }
    $shots | ForEach-Object { Write-Host "capture: shot   $($_.FullName)" }
}
Write-Host "capture: log    $log"
exit $p.ExitCode
