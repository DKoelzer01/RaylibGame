# (#1958) SessionStart hook: the deweyserver code and issue lookups are deferred MCP tools, so a session has to load them
# with ToolSearch before it can call them, and sessions fell back to grep/cat over the source instead (0 code_* calls in
# four topic sessions, 2026-11-26; shell read/search output is ~40% of the tool-result tokens). A SessionStart hook's stdout
# becomes context: print the exact line to load them and the rule. No network, no git: it must be instant.
$ErrorActionPreference = 'SilentlyContinue'
$tools = 'mcp__deweyserver__code_find,mcp__deweyserver__code_refs,mcp__deweyserver__code_outline,mcp__deweyserver__code_grep,mcp__deweyserver__issue_search,mcp__deweyserver__ci_summary'
Write-Output "CODE LOOKUPS FIRST (#1958): load the deweyserver lookups now with ToolSearch query ""select:$tools"" and use them before grep/rg/cat/sed over src/, tests/ or resources/ (a few lines per answer; they index staging, ladder and develop, refreshed every 15 min). Grep or read the files yourself only for code your branch changed. Run builds and tests through tools/quiet.sh. If deweyserver failed to connect at startup, /mcp reconnects it."
exit 0
