# RaylibGame LOOKUP GATE (#1958, the user 2026-11-27: "the lookup should require tool usage unless tool usage can't find
# it"; 2026-11-28: the repo hook blocks too). The user-level ~/.claude/hooks/rg-lookup-gate.ps1 stands down where this
# file carries the words "LOOKUP GATE", so a PC with both runs the gate once. (This file was the "not a block" hint.)
#
# PreToolUse (Bash, PowerShell, Grep): a command that searches or prints source files (grep/rg/cat/sed/head/
# Select-String/Get-Content/git grep|show over src/, tests/ or resources/) is refused unless
#   - a deweyserver code_* lookup in this session found nothing or failed in the last 10 minutes (the "can't find it"
#     case: the index then can't help), or
#   - every source path the command names is a file this branch changed (the index holds staging/ladder/develop only).
# A reader after a pipe filters another command's output and is not gated. Read stays open: a lookup gives the file and
# line, then Read with offset/limit takes only those lines (CLAUDE.md).
# PostToolUse (mcp__deweyserver__code_*): an empty or failed lookup stamps the session's miss file.
# Off switch: RAYLIBGAME_LOOKUP_GATE=0. Decisions are logged to ~/.claude/lookup-gate.log (time, session, decision).
$ErrorActionPreference = 'SilentlyContinue'
$p = [string]$env:CLAUDE_PROJECT_DIR
$off = [string]$env:RAYLIBGAME_LOOKUP_GATE
if ($off -in @('0', 'false', 'off')) { exit 0 }
try { $evt = [Console]::In.ReadToEnd() | ConvertFrom-Json } catch { exit 0 }
$sid = ([string]$evt.session_id) -replace '[^A-Za-z0-9_-]', ''
if (-not $sid) { $sid = 'none' }
$miss = Join-Path $env:TEMP "claude-lookup-miss-$sid.txt"
$log = Join-Path $env:USERPROFILE '.claude\lookup-gate.log'
$now = [DateTimeOffset]::UtcNow.ToUnixTimeSeconds()
function Write-GateLog([string]$what) { Add-Content -Path $log -Value "$([DateTime]::UtcNow.ToString('s')) $sid $what" -Encoding ascii }
$tool = [string]$evt.tool_name

# --- PostToolUse: remember a lookup that came back empty or failed
if ([string]$evt.hook_event_name -eq 'PostToolUseFailure') {
    # a lookup that failed outright (deweyserver down, MCP error) fires this event, not PostToolUse: always a miss
    if ($tool -match '^mcp__deweyserver__code_') { Set-Content -Path $miss -Value $now -Encoding ascii; Write-GateLog "miss-failed $tool" }
    exit 0
}
if ([string]$evt.hook_event_name -eq 'PostToolUse') {
    if ($tool -notmatch '^mcp__deweyserver__code_') { exit 0 }
    $r = $evt.tool_response | ConvertTo-Json -Compress -Depth 8
    $empty = (-not $r) -or $r -match 'total_hits\W{1,6}0\b' -or $r -match '(results|matches|refs|symbols|definitions)\W{1,6}\[\s*\]' -or
             $r -match '(?i)not found|no match|not indexed|unknown branch' -or $r -match '"?error\\?"\s*:'
    if ($empty) { Set-Content -Path $miss -Value $now -Encoding ascii; Write-GateLog "miss $tool" }
    exit 0
}

# --- PreToolUse: is this a read of source?
$srcRe = '\b(src|tests|resources)([/\\]|\s|$|")'
if ($tool -eq 'Grep') {
    $where = "$([string]$evt.tool_input.path) $([string]$evt.tool_input.glob)"
    $isRead = ($where -match $srcRe) -or ([string]$evt.tool_input.type -in @('cpp', 'c', 'glsl'))
    $cmd = $where
} else {
    $cmd = [string]$evt.tool_input.command
    if (-not $cmd) { exit 0 }
    # quoted text is a pattern, a sed script or a message, not a path to read (a sed whose replacement says "src/" over a
    # build/ file was refused)
    $cmd = $cmd -replace "'[^']*'", "''" -replace '"[^"]*"', '""'
    $reader = '(grep|rg|egrep|cat|head|tail|sed|awk|less|Select-String|Get-Content|sls|gc|git\s+(grep|show))'
    # a reader as a command word (start, or after ; & ( or a newline - not after a pipe), then a source path before the
    # next separator
    $isRead = $cmd -match "(^|[;&(\n]\s*)$reader\b[^;&|\n]*$srcRe"
}
if (-not $isRead) { exit 0 }

# a recent empty lookup: the index can't help, read away
$t = 0
if ((Test-Path $miss) -and [long]::TryParse((Get-Content $miss -TotalCount 1), [ref]$t) -and ($now - $t) -le 600) { Write-GateLog "allow miss"; exit 0 }

# deweyserver unreachable (outage: its code_* tools may vanish rather than fail, so no miss is ever stamped): reads are
# open. Probed at most once a minute, 3 s timeout; the answer is cached in the temp dir.
$probe = Join-Path $env:TEMP 'claude-lookup-probe.txt'
$up = $null
if (Test-Path $probe) { $l = (Get-Content $probe -TotalCount 1) -split ' '; if ($l.Count -eq 2 -and ($now - [long]$l[0]) -le 60) { $up = $l[1] -eq 'up' } }
if ($null -eq $up) {
    $up = $false
    try { $null = Invoke-WebRequest -Uri 'https://mcp.tail68a71f.ts.net/mcp' -Method Head -TimeoutSec 3 -UseBasicParsing; $up = $true }
    catch { if ($_.Exception.Response) { $up = $true } } # any HTTP answer (401, 405...) means the server is there
    Set-Content -Path $probe -Value "$now $(if ($up) { 'up' } else { 'down' })" -Encoding ascii
}
if (-not $up) { Write-GateLog "allow server-down"; exit 0 }

# every named source path is one of this branch's own changes (committed since staging, staged, unstaged or untracked)
$paths = @([regex]::Matches(($cmd -replace '\\', '/'), '(?<![\w:])(?:[\w.-]+:)?((?:src|tests|resources)/[\w./-]*)') | ForEach-Object {
    if ($_.Value -match '^[\w.-]+:(src|tests|resources)/') { '#indexed-ref' } else { $_.Groups[1].Value.TrimEnd('/', '.') } })
if ($paths.Count -gt 0 -and $paths -notcontains '#indexed-ref') {
    $cwd = [string]$evt.cwd
    if (-not $cwd) { $cwd = $p }
    # the branch's own base: of staging and ladder, the one it differs from least (a ladder branch against staging
    # would count all of ladder as its own)
    $changed = $null
    foreach ($b in @('origin/staging', 'origin/ladder')) {
        $mb = git -C $cwd merge-base HEAD $b 2>$null
        if (-not $mb) { continue }
        $d = @(git -C $cwd diff --name-only $mb 2>$null)
        if ($null -eq $changed -or $d.Count -lt $changed.Count) { $changed = $d }
    }
    if ($null -eq $changed) { $changed = @() }
    $changed += git -C $cwd ls-files --others --exclude-standard -- src tests resources 2>$null
    $own = $true
    foreach ($q in $paths) { if ($changed -notcontains $q) { $own = $false; break } }
    if ($own) { Write-GateLog "allow own"; exit 0 }
}

Write-GateLog "deny $tool"
$why = "LOOKUP REQUIRED (#1958, the user's rule): this reads source files with a shell or Grep search. Use the deweyserver lookups first: mcp__deweyserver__code_find(name), code_refs(name), code_outline(file), code_grep(pattern, path_glob, branch) (staging, ladder, develop; load once with ToolSearch ""select:mcp__deweyserver__code_find,mcp__deweyserver__code_refs,mcp__deweyserver__code_outline,mcp__deweyserver__code_grep""), then Read the lines they point at with offset/limit. Allowed without a lookup: files your branch changed (name them exactly), and any source read for 10 minutes after a lookup that found nothing or failed."
@{ hookSpecificOutput = @{ hookEventName = 'PreToolUse'; permissionDecision = 'deny'; permissionDecisionReason = $why } } | ConvertTo-Json -Compress -Depth 4
exit 0
