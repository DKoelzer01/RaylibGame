---
name: morning-report
description: Coordinator only. Build the report the user reads when they come back - merges since they left, each session's progress, decisions waiting for them (as questions with a recommendation), failures and incidents told plainly, and what's next. Use when the user returns after being away, or asks "what happened" / "what's up".
---

# Morning report (coordinator)

Gather, cheaply:
- **Merges, open PRs, issues and runner health in one call:** the deweyserver MCP `morning_data(since)` (#1879): PRs
  merged since then (all bases), open PRs with their "CI ok" state, issues opened and closed, and the runner digest
  (NEW failures, flaky candidates, known failures that now pass). `ci_summary(run)` only for a run you need in detail;
  no log downloads.
- **Decisions waiting:** the "MORNING"/"DECIDE" lines in the coordinator's notes (build/coordinator-notes.md) and the
  sessions' reports since then.
- **Board movement:** issues set to In review / Done since then (the merges' "Fixes" lists).
- **Incidents:** anything that went wrong (a rule broken, a red merge, a near miss) - disclose it plainly with the fix.

Write it in this order, short:
1. One paragraph: what changed for the user (player-visible fixes to playtest, with the playtest board link).
2. Decisions needed - each a question with a recommendation (use AskUserQuestion for up to 4 at once).
3. Per session: done, in progress, next (one line each).
4. Failures and incidents.
5. What happens next without input.

Don't paste logs or long lists; link the PRs and issues.
