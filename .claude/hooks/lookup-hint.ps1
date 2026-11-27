# (#1958) PreToolUse hint (never a block) for Bash/PowerShell commands that search or print source files (grep, rg, cat,
# sed, head, Select-String, Get-Content, git grep/show over src/, tests/ or resources/): the deweyserver code_grep /
# code_find / code_outline answer in a few lines instead of pages of output that stay in context. Rate-limited per session
# (the 1st match, then every 10th) so it doesn't become noise for work on a branch's own changes, which the index can't see.
$ErrorActionPreference = 'SilentlyContinue'
try { $evt = [Console]::In.ReadToEnd() | ConvertFrom-Json } catch { exit 0 }
$cmd = [string]$evt.tool_input.command
if (-not $cmd) { exit 0 }
$reader = '(grep|rg|egrep|cat|head|tail|sed|awk|less|Select-String|Get-Content|sls|gc|git\s+(grep|show))'
$source = '\b(src|tests|resources)([/\\]|\s|$)'
# a reader as a command word (start, or after ; & ( or a newline; not after a pipe, where it filters another command's
# output), then a source path before the next separator
if ($cmd -notmatch "(^|[;&(\n]\s*)$reader\b[^;&|\n]*$source") { exit 0 }
$sid = ([string]$evt.session_id) -replace '[^A-Za-z0-9_-]', ''
if (-not $sid) { $sid = 'none' }
$f = Join-Path $env:TEMP "claude-lookup-hint-$sid.txt"
$n = 0
if (Test-Path $f) { [int]::TryParse((Get-Content $f -TotalCount 1), [ref]$n) | Out-Null }
$n++
Set-Content -Path $f -Value $n -Encoding ascii
if ((($n - 1) % 10) -ne 0) { exit 0 }
$msg = "HINT (#1958, not a block; match $n this session): this command reads source files into context. For code on staging/ladder/develop, mcp__deweyserver__code_grep(pattern, path_glob, branch), code_find(name), code_outline(file) or code_refs(name) answer in a few lines (load them once: ToolSearch ""select:mcp__deweyserver__code_grep,mcp__deweyserver__code_find,mcp__deweyserver__code_outline,mcp__deweyserver__code_refs""). Keep using grep/Read for code your branch changed."
@{ hookSpecificOutput = @{ hookEventName = 'PreToolUse'; additionalContext = $msg } } | ConvertTo-Json -Compress -Depth 4
exit 0
