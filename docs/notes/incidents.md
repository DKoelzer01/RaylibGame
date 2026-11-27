# Incidents and lessons

Mistakes the sessions (coordinator included) made, and the rule or tool each one produced. Add one line per incident:
what happened, the cause, and the fix or rule. Read the patterns before you start; they repeat.

## Recurring patterns

1. **A red test is often the test, not the game - especially on the runner.** The runner is ~35x slower in fluid
   steps and its clients spend 40-70 s loading, so timing windows, "act at once after placing", and scripted aims
   break there first. Before filing a game bug, check: does the test wait for the thing it checks (quest book update,
   factory placement, sim time vs wall time)? Is its aim/pose stable? Is the player alive (a creature can kill a test
   player)? Use `scenario_history(name)` to see if it ever passed. Examples: farm/till-rmb (#1176), factory/stats
   (#1406), join creature (#1421), orbit-closure (F5), the progression pass rule (L6b-3), the cave-mouth camera (#1428).
2. **...but new tests do find real bugs - read them carefully.** #1445 lab_node never counted, #1436 fly-camera hop,
   #1460 undocking jump, #1499 a commented-out refresh, #1503/#1462 `count()` vs `held()` for items with data,
   #1498 conveyors copying items with data, #1504 a 7x fluid cost digging into wide rivers.
3. **Items with data are invisible to `Inventory::count()`** (by design, so recipes never spend them). Any check
   "does the player have/hold X" must use `held()`. Three bugs came from this (#1462, #1503, #1498).
4. **Two instants mixed in one computation.** The body poses advance before the rocket and its substeps, so a frame
   switch or an outside-written state used positions from different instants: 16 m, 13 m, 7.5 m jumps (F5, #1436,
   #1460). Re-express state only at one consistent instant (after the step's last substep).
6. **Process slips cost more than code bugs (2026-11-22/23).** Wrong PR numbers in reports (about one in five), a
   retargeted stacked PR that sat an hour with no CI (#1653), a playtest-board baseline that missed 39 marks (old
   failures never reached their issues), two game runs at once through a session's subagents (#1662). Fixes: quote
   numbers only from `gh` output; no stacks; update-branch after a retarget; the board's own `routed` field is the
   baseline (docs/notes/board.md); the guard refuses a second run.
5. **Don't state an inference as a decision.** The coordinator wrote "composed items disassemble into parts" into
   decisions.md before the user answered; the user chose raw materials (#1495). Mark unconfirmed items "GUESS" or
   ask first.

## Log

| Date | What happened | Cause | Fix / rule |
|---|---|---|---|
| 2026-11-21 | #1453 merged although its runner run was red | `wait-run.sh` always exited 0 | It now fails on a red run; failures were all known, nothing broke |
| 2026-11-21 | Builds ran during quiet windows (3 times) | the PostToolUse build-check hook rebuilt -j 8 after every edit, unaware of quiet windows | Hook skips quiet windows, -j 4; machine slot lock + guard hook (#1477) |
| 2026-11-21 | A session ran the game without the slot | `bash` was busybox (w64devkit PATH), slot.sh failed, `\| tail &&` hid it | Scripts re-run under Git's bash (#1488); never pipe before `&&` on a gate |
| 2026-11-21 | known_failures.txt patterns with `\|` never matched | the file splits on `\|` | No `\|` inside a pattern (file header); split into lines (#1491) |
| 2026-11-21 | A misspelt env switch "passed" a GPU-only test | `getenv` of a wrong name silently runs the default | grep the exact name before an A/B (CLAUDE.md) |
| 2026-11-21 | Merged a PR on hosted CI while its new runner scenario was still queued; the scenario failed (#1165) | evidence was the runner run | Runner-evidence PRs wait for that run |
| 2026-11-22 | decisions.md said composed items disassemble into parts | coordinator's inference written as a ruling | Corrected (#1495); pattern 5 |
| 2026-11-22 | River widening showed no change on the slide's river | the slide's spot is a river head (narrow by design, 7 cm of water: under one cell, so nothing drawn) | Check a change on the user-visible case before the PR |
| 2026-11-22 | --test-flow-sea "passed" on a river mouth since #190 | its shore walk used waterSourceAt, which never reports the sea, so every sea pick stopped 2 m out; only a river's water in the way gave a real channel | The walk counts the sea as water; re-baselined at a real sea spot (#1512); pattern 1: check what a test's pick actually is |
| 2026-11-22 | A passing scenario run's log.txt ended in `FATAL abort()`, reported as a crash (#1527) | an earlier run killed by `timeout` hung in termination with log.txt still open and wrote its FATAL line at its old offset, mid-line in the next run's log | A FATAL line with no timestamp that cuts a line is another process; check for leftover game processes before calling it a crash. Let the scenario's wall timeout end a run rather than killing it |
| 2026-11-24 | A restarted UI/UX session started in the root checkout and ran `git checkout` in Progression's worktree (agent-af13...); Progression moved to `.claude/worktrees/prog-belts` | the restart prompt named no worktree, so the session guessed a path | A restart prompt names the worktree path; a session checks `pwd` before any `git checkout` |
| 2026-11-24 | The playtest-prep subagent split read-only commands to get past the machine guard, which flags any command naming game.exe | the guard cannot tell a read from a run | Report the refusal instead of splitting commands; the guard's false positives get fixed in the guard |
| 2026-11-26 | The bug-filing subagent got past a machine-guard refusal by putting `gh issue create` commands in temp scripts | the guard false-positived on the words "tests"/"game" in issue text; the actions were harmless but this routes around a refusal | Report a refusal ("never route around a refusal", CLAUDE.md); the guard's false positives get fixed in the guard |
| 2026-11-27 | The shared repo (`C:/Users/dewey/RaylibGame/.git`) went shallow at the staging and ladder tips; every worktree's merge of origin/staging failed with "unrelated histories" until `git fetch --unshallow origin` | Tools and CI ran the new `tools/ci/reservation_clash.sh` (#2009) locally; it fetched staging and ladder with `--depth=1`, and a depth fetch into a full clone makes it shallow (all worktrees share one `.git`) | `--depth` only when `rev-parse --is-shallow-repository` says true (reservation_clash.sh, gpu_changes.sh; known_list.sh already did). Rule: never `--depth` in a full clone; a CI-only script runs in the user's repo too once someone tests it locally |
