# (#1897) SessionStart hook, matcher "compact": after an auto-compaction, put the session's build/session-handoff.md into the
# new context next to the summary. The safety net for topic sessions launched with --autocompact (tools/coord/cycle_session.sh):
# a compaction summary loses detail, while the handoff is what the session itself wrote down to continue from.
# A SessionStart hook's stdout becomes context. Prints nothing when there is no handoff.
$ErrorActionPreference = 'SilentlyContinue'
$cwd = $null
try { $cwd = ([Console]::In.ReadToEnd() | ConvertFrom-Json).cwd } catch {}
if (-not $cwd) { $cwd = (Get-Location).Path }
$f = Join-Path $cwd 'build/session-handoff.md'
if (-not (Test-Path $f)) {
    $r = git -C $cwd rev-parse --show-toplevel 2>$null
    if ($r) { $f = Join-Path $r 'build/session-handoff.md' }
}
if (-not (Test-Path $f)) { exit 0 }
$item = Get-Item $f
$text = Get-Content -Raw -Encoding utf8 $f
$max = 20000 # characters: a handoff is a page, not a log
if ($text.Length -gt $max) { $text = $text.Substring(0, $max) + "`n[... cut at $max characters; read the file for the rest]" }
$age = [int]((Get-Date) - $item.LastWriteTime).TotalMinutes
Write-Output "Context was compacted. Your handoff file $f (written $age min ago; if it is older than the work in the summary above, the summary wins):"
Write-Output ''
Write-Output $text
exit 0
