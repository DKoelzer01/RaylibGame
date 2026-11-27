# PreToolUse hook (Bash|PowerShell), #1477: the one-at-a-time machine rule, enforced.
# - Running the game or the tests (game.exe, game_tests.exe, ctest) needs the machine: this session's worktree must
#   hold it (tools/slot.sh request|take). Otherwise the command is refused with how to ask.
# - Building (cmake --build) during someone else's quiet window is refused (tools/slot.sh wait-quiet, then build).
# Commands that only mention those names (issue text, commit messages, heredocs, gh, git, grep, sha256sum, ssh,
# tools/sandbox.sh...) pass; the classifier is in machine-guard-classify.ps1, its tests in machine-guard.tests.ps1.
# Fails open: any error here lets the command through, so a broken hook never blocks work.
$ErrorActionPreference = 'Stop'

# (#1919) The hold is named after the worktree slot.sh ran in, but the event's cwd is the session's start folder, so a
# session (or its subagent) working in another worktree failed `slot.sh holds`. Its hold still counts when it carries
# this session's id (slot.sh records CLAUDE_CODE_SESSION_ID; subagents share their session's id), or, for a hold
# without one, when the command works in the holder's worktree (a path inside it, or a cd into it).
function Test-HoldIsThisSession([string]$bash, [string]$slot, [string]$session, [string]$cmd) {
    $info = (& $bash $slot holder-info 2>$null) -join ''
    if ($LASTEXITCODE -ne 0 -or -not $info) { return $false }
    $name, $sid = $info -split '\|', 2
    if ($sid -and $session) { return $sid -eq $session }
    if (-not $name) { return $false }
    $lname = $name.ToLowerInvariant()
    # one spelling for paths: forward slashes, lower case, /c/... -> c:/...
    $c = (($cmd -replace '\\', '/').ToLowerInvariant()) -replace '(^|[\s"''=])/([a-z])/', '$1$2:/'
    foreach ($line in (git worktree list --porcelain 2>$null)) {
        if ($line -notlike 'worktree *') { continue }
        $wt = ($line.Substring(9) -replace '\\', '/').ToLowerInvariant().TrimEnd('/')
        if ((Split-Path -Leaf $wt) -eq $lname -and $c -match ([regex]::Escape($wt) + '([/\s"'';&|]|$)')) { return $true }
    }
    # a relative cd into it: cd .claude/worktrees/<name>
    return $c -match ('(^|[;&|]\s*)cd\s+["'']?(\S*/)?' + [regex]::Escape($lname) + '/?["'']?(\s|$|;|&)')
}

try {
    $raw = [Console]::In.ReadToEnd()
    $evt = $raw | ConvertFrom-Json
    $cmd = [string]$evt.tool_input.command
    if (-not $cmd) { exit 0 }
    $cwd = [string]$evt.cwd
    if (-not $cwd) { $cwd = (Get-Location).Path }

    # What the command executes, not the words in its text (#1864): see machine-guard-classify.ps1. Direct compiler
    # calls and programs run out of a build folder count as builds (2026-11-23 incident), so a refused cmake build
    # can't be worked around during someone else's quiet window.
    . (Join-Path $PSScriptRoot 'machine-guard-classify.ps1')
    $kind = Get-MachineGuardKind $cmd
    $runs = $kind -eq 'run'
    $builds = $kind -eq 'build'
    if (-not ($runs -or $builds)) { exit 0 }

    $slot = Join-Path $PSScriptRoot '../../tools/slot.sh' # this checkout's slot tool (CLAUDE_PROJECT_DIR may not be the worktree)
    if (-not (Test-Path $slot)) { exit 0 }  # an older checkout without the tool: don't enforce
    $bash = 'C:\Program Files\Git\bin\bash.exe'
    if (-not (Test-Path $bash)) { exit 0 }
    Push-Location $cwd
    try {
        if ($runs) {
            & $bash $slot holds *> $null
            $held = $LASTEXITCODE -eq 0
            if (-not $held) { $held = Test-HoldIsThisSession $bash $slot ([string]$evt.session_id) $cmd }
            if (-not $held) {
                $st = (& $bash $slot status 2>&1) -join "`n"
                [Console]::Error.WriteLine("machine-guard: running the game/tests needs the machine slot (one session at a time, #1477).`n$st`nAsk with: tools/slot.sh request <minutes> [quiet] <what>   (run it in the background; it returns when the machine is yours), then tools/slot.sh free when done.")
                exit 2
            }
            # One run at a time even inside the holding session: its subagents share its worktree and so count as the
            # holder, and two runs at once spoil timing and captures (user, 2026-11-23). A run that needs two processes
            # (a host and a client) opts out with SLOT_PARALLEL=1 in the command.
            if ($cmd -notmatch 'SLOT_PARALLEL=1') {
                $busy = Get-Process -Name 'game', 'game_tests' -ErrorAction SilentlyContinue |
                    Where-Object { $_.Path -and $_.Path -notmatch '(?i)\\RaylibGame\\game\.exe$' } # the user's root game.exe aside
                if ($busy) {
                    $list = ($busy | ForEach-Object { "$($_.Id) $($_.Path)" }) -join "`n  "
                    [Console]::Error.WriteLine("machine-guard: a game/test run is already going (one run at a time, subagents included):`n  $list`nWait for it to finish (or stop your own stale one). A run that needs a host and a client together: put SLOT_PARALLEL=1 in the command.")
                    exit 2
                }
            }
        } elseif ($builds) {
            $q = (& $bash $slot quiet-other 2>$null) -join ''
            if ($LASTEXITCODE -eq 0) {
                [Console]::Error.WriteLine("machine-guard: quiet window held by $q - no builds until it ends. Wait with: tools/slot.sh wait-quiet (then build at low priority, -j 4).")
                exit 2
            }
        }
    } finally { Pop-Location }
    exit 0
} catch { exit 0 }
