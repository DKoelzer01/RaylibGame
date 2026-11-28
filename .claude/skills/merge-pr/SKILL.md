---
name: merge-pr
description: Coordinator only. Merge reviewed PRs into staging or ladder - auto-merge on "CI ok", or wait for a runner run when the PR's evidence is a runner scenario - then update the board and the playtest board. Use when a session reports a PR ready. Feature sessions never merge.
argument-hint: "<PR> [--runner <run id>] [--ladder]"
---

# Merge a PR (coordinator)

1. Check the PR: base branch (staging, or `ladder` for ladder/combat work), the issue references, the evidence.
   Stacked PR whose base merged: `bash tools/coord/stackstep.sh <PR>` (retargets, updates the branch, waits, merges).
2. Plain PR: `gh pr merge <PR> --auto --merge` - GitHub merges it when "CI ok" passes (both bases are protected).
3. Evidence from a runner scenario: don't auto-merge. Wait for the run, then merge:
   `bash .claude/skills/runner-check/scripts/runner_check.sh --run <id>` (exit 0 or 3 = OK), then
   `KEEPBASE=1 bash tools/coord/merge-ready.sh <PR>` (KEEPBASE keeps a ladder base).
4. After it merges:
   - board: `bash tools/coord/board.sh set <issue> -s "In review"` (Done only when the fix reaches develop) and label
     the issues `fixed in open PR` if they aren't;
   - player-visible change: add a line to the playtest board (build/playtest-board/index.html in the coordinator's
     worktree) and republish it;
   - tell the session (one line), and any session whose branch builds on it.
5. Merges into develop/main, and ladder -> staging, only with the user's explicit OK.
