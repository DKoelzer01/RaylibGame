# Command classifier for machine-guard.ps1 (#1477, false positives fixed in #1864).
# Get-MachineGuardKind returns 'run' (starts the game, game_tests or ctest on this PC), 'build' (a compile) or ''.
# It looks at what a command EXECUTES, not at words in its text:
# - message text is dropped first: heredoc bodies (unless they feed a shell), PowerShell here-strings, and the quoted
#   values of --body/--title/-m and similar flags, so issue text or a commit message naming game.exe passes;
# - the command is split into segments (; && || | & newline $( and backtick), each segment's wrappers are peeled
#   (VAR=x, timeout, nice, bash -c, Start-Process...), and a segment whose command word is read-only or remote
#   (gh, git, grep, sha256sum, ssh, tools/sandbox.sh...) is skipped;
# - a run needs game(_tests).exe or a path to game/game_tests (bash finds ./game as game.exe), not the bare word;
#   a compiler counts only as the command word, cmake --build anywhere in a non-inert segment.
# Splitting inside quotes is deliberate: `bash -c 'cd x; build/check/bin/game.exe'` still yields the game.exe segment.

$script:MgInert = '^(gh|git|grep|egrep|fgrep|rg|ls|cat|head|tail|less|wc|sed|awk|find|echo|printf|cp|mv|rm|mkdir|rmdir|touch|stat|file|du|df|sha256sum|sha1sum|md5sum|cmp|diff|sort|uniq|cut|tr|jq|xxd|strings|objdump|nm|ldd|which|type|pwd|cd|pushd|popd|export|unset|test|\[|\[\[|true|false|sleep|date|basename|dirname|realpath|readlink|tee|ssh|scp|curl|sandbox\.sh|slot\.sh|copy|Copy-Item|Move-Item|Remove-Item|New-Item|Rename-Item|Set-Location|Push-Location|Pop-Location|Get-\w+|Select-String|Select-Object|Where-Object|Sort-Object|Measure-Object|Format-\w+|Out-File|Out-String|Set-Content|Add-Content|Write-Output|Write-Host|Test-Path|Resolve-Path|Split-Path|Join-Path|ConvertFrom-Json|ConvertTo-Json|Stop-Process|taskkill|tasklist)$'

function Remove-MgMessageText([string]$s) {
    # Heredocs: keep the header line, drop the body unless the line's command is a shell (then the body runs).
    $s = [regex]::Replace($s, '(?ms)^([^\n]*?)<<-?[ \t]*([''"]?)(\w+)\2([^\n]*)\n.*?^[ \t]*\3[ \t]*\r?$', {
            param($m)
            if ($m.Groups[1].Value -match '(?i)(^|[\s;&|(])(bash|sh|zsh|pwsh|powershell|cmd)(\.exe)?(\s|$)') { return $m.Value }
            return $m.Groups[1].Value + '<<' + $m.Groups[3].Value + $m.Groups[4].Value
        })
    # PowerShell here-strings are data (commit messages, file text).
    $s = [regex]::Replace($s, '(?s)@''\r?\n.*?\r?\n''@', "''")
    $s = [regex]::Replace($s, '(?s)@"\r?\n.*?\r?\n"@', "''")
    # Quoted values of text flags (gh --body/--title, git -m, gh api -f body=...).
    $q = '("(?:[^"\\]|\\.)*"|''[^'']*'')'
    $s = [regex]::Replace($s, '(?i)(\s(?:--body|--title|--notes|--message|--comment|--subject|--description|-m|-b|-t))(?:\s+|=)' + $q, '$1 ""')
    $s = [regex]::Replace($s, '(\s(?:-f|-F|--field|--raw-field)\s+[\w.\[\]]+=)' + $q, '$1""')
    # The quoted command of ssh / tools/sandbox.sh runs on deweyserver, not on this PC.
    $s = [regex]::Replace($s, '(?i)((?:^|[\s;&|(])(?:ssh(\.exe)?|[\w./\\-]*sandbox\.sh)\s[^''"\r\n;&|]*?)' + $q, '$1""')
    return $s
}

function Get-MgSegmentCore([string]$seg) {
    # Peel wrappers until the real command word leads the segment.
    $prev = $null
    while ($seg -ne $prev) {
        $prev = $seg
        $seg = $seg -replace '^[\s({!]+', ''
        $seg = $seg -replace '^["'']+', ''
        $seg = $seg -replace '^\w+=("[^"]*"|''[^'']*''|\S*)\s+', ''                               # VAR=x cmd
        $seg = $seg -replace '(?i)^(timeout|gtimeout)(\s+-\S+)*\s+\S+\s+', ''
        $seg = $seg -replace '(?i)^nice(\s+-n\s*-?\d+|\s+-+\d+)?\s+', ''
        $seg = $seg -replace '(?i)^(time|env|exec|nohup|command|builtin)\s+', ''
        $seg = $seg -replace '(?i)^(bash|sh|zsh)(\.exe)?(\s+-\w+)*\s+-c\s+', ''
        $seg = $seg -replace '(?i)^(bash|sh|zsh)(\.exe)?\s+(?=\S+\.sh(\s|$|["'']))', ''            # bash tools/x.sh
        $seg = $seg -replace '(?i)^(powershell|pwsh)(\.exe)?(\s+-\w+)*\s+-(Command|c)\s+', ''
        $seg = $seg -replace '(?i)^cmd(\.exe)?\s+/[ck]\s+', ''
        $seg = $seg -replace '(?i)^Start-Process\s+(-FilePath\s+)?', ''
        $seg = $seg -replace '(?i)^Invoke-Expression\s+', ''
    }
    return $seg
}

function Test-MgInertSegment([string]$core) {
    if ($core -match '^\$[\w:]+\s*=') { return $true }                                            # $env:X = 'y'
    $word = (($core -split '\s+')[0]) -replace '["'']', ''
    if (-not $word) { return $true }
    if ($word -match "(?i)$script:MgInert") { return $true }
    $base = ($word -split '[/\\]')[-1]
    return ($base -match "(?i)$script:MgInert")
}

function Get-MachineGuardKind([string]$cmd) {
    if (-not $cmd) { return '' }
    $text = Remove-MgMessageText $cmd
    $segments = [regex]::Split($text, '\r?\n|&&|\|\||[;|&]|\$\(|`')
    $kind = ''
    foreach ($seg in $segments) {
        $core = Get-MgSegmentCore $seg
        if (Test-MgInertSegment $core) { continue }
        if ($core -match '(?i)(^|[\s"''=/\\(])game(_tests)?\.exe(?=$|[\s"''),;])' -or
            $core -match '(?i)(^|[\s"''=(])[\w.:~$/\\-]*[/\\]game(_tests)?(?=$|[\s"''),;])' -or
            $core -match '(?i)^game(_tests)?(\s|$)' -or
            $core -match '(?i)(^|[\s"''(])ctest(\.exe)?(\s|$)') {
            return 'run'
        }
        if ($core -match '(?i)cmake(\.exe)?["'']?\s+--build' -or
            $core -match '(?i)^([\w.:~$/\\-]*[/\\])?(g\+\+|gcc|c\+\+|clang\+\+|clang|cc|ninja|make|mingw32-make)(\.exe)?["'']?(\s|$)' -or
            $core -match '(?i)^[^\s"'']*build[/\\][^\s"'']*\.exe') {
            $kind = 'build'
        }
    }
    return $kind
}
