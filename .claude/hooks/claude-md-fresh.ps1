# (#1911) SessionStart hook: a session started in a checkout on an old branch loads that branch's CLAUDE.md and misses the
# rules added on staging since. Compare the cwd's CLAUDE.md with origin/staging's (the local ref, no fetch: it must be fast,
# and the refs are shared by every worktree, which other sessions fetch all day).
#   same as staging                     -> print nothing
#   an older staging version            -> print staging's CLAUDE.md; a SessionStart hook's stdout becomes context
#   neither (the branch edited it)      -> one line: it differs, and how to see the diff
$ErrorActionPreference = 'SilentlyContinue'
$cwd = $null
try { $cwd = ([Console]::In.ReadToEnd() | ConvertFrom-Json).cwd } catch {}
if (-not $cwd) { $cwd = (Get-Location).Path }
$root = git -C $cwd rev-parse --show-toplevel 2>$null
if (-not $root) { exit 0 }
$f = Join-Path $root 'CLAUDE.md'
if (-not (Test-Path $f)) { exit 0 }
$ref = 'origin/staging'
$want = git -C $root rev-parse --verify -q "${ref}:CLAUDE.md" 2>$null
if (-not $want) { exit 0 }
# hash-object applies the checkout filters (CRLF), so this is the blob the file would commit as
$have = git -C $root hash-object -- $f 2>$null
if (-not $have -or $have -eq $want) { exit 0 }
# Blobs CLAUDE.md had on staging: "--raw" lines read ":100644 100644 <old> <new> M<TAB>CLAUDE.md" (full hashes with --no-abbrev).
$old = @{}
foreach ($line in (git -C $root log -200 --format= --raw --no-abbrev $ref -- CLAUDE.md 2>$null)) {
    $p = $line -split '\s+'
    if ($p.Count -ge 4) { $old[$p[2]] = 1; $old[$p[3]] = 1 }
}
$branch = git -C $root rev-parse --abbrev-ref HEAD 2>$null
if ($old.ContainsKey($have)) {
    $text = (git -C $root show "${ref}:CLAUDE.md" 2>$null) -join "`n"
    Write-Output "NOTE (#1911): the CLAUDE.md loaded above is an OLD version (this checkout is on '$branch'). $ref's CLAUDE.md below is current and supersedes it where they differ:"
    Write-Output ''
    Write-Output $text
} else {
    Write-Output "NOTE (#1911): this checkout's CLAUDE.md ('$branch') differs from ${ref}'s (edited on this branch, or $ref is older). See: git diff $ref -- CLAUDE.md"
}
exit 0
