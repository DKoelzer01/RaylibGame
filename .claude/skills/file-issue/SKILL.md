---
name: file-issue
description: File a RaylibGame GitHub issue the project's way - search for a match first (reuse or reopen it), then create one with labels and an optional parent, then set its board fields. Use before starting any task without an issue, and for every bug found along the way.
argument-hint: "<keywords or title>"
---

# File an issue

1. **Search first:** the deweyserver MCP `issue_search("<the problem in a sentence>")` (#1878: ranked, every issue,
   one line each); `bash .claude/skills/file-issue/scripts/file_issue.sh search "<keywords>"` if the MCP is down.
   A match: use it; reopen it (`gh api -X PATCH repos/DKoelzer01/RaylibGame/issues/N -f state=open`) if the problem is back.
2. **Write the body** to a file under build/: the problem or request (the user's words if they gave them), the cause if
   known, the planned fix. Mark guessed values "GUESS".
3. **Create:** `bash .claude/skills/file-issue/scripts/file_issue.sh create "<title naming the problem>" build/issue.md
   --labels bug,water [--parent <tracking issue>]`. Labels: bug or enhancement, plus area labels (rendering, shadows,
   terrain, water, physics, camera, ui, performance, tests, ...). A batch of tasks gets one issue each, under a parent.
4. **Board:** set all four fields with the deweyserver MCP tool `board_set(issue, status, priority, size, kind)`.
   If it says the token is missing, put "#N -> Status / Priority / Size / Kind" in your report instead.
   Priority follows the user's order; never raise one on your own - propose it.
5. Reference it in commits `(#N)` and the PR's "Fixes" list.
