---
name: runner-check
description: Run or read the self-hosted deweyserver workflow (GPU bit-exactness, headless scenarios, server and join tests) for a branch or run id, wait with a deadline, and list its failures marked "known (#N)" from tools/ci/known_failures.txt or "NEW". Use instead of grepping runner logs by hand, before reporting a PR ready, and whenever a runner run is red.
argument-hint: "<branch> [--dispatch] [--regression]   or   --run <run id>"
---

# Runner check

Run it in the background (it waits up to 60 min for the run to finish):

```bash
bash .claude/skills/runner-check/scripts/runner_check.sh <branch> --dispatch      # start a run for your branch, then wait
bash .claude/skills/runner-check/scripts/runner_check.sh <branch>                 # the branch's latest run
bash .claude/skills/runner-check/scripts/runner_check.sh --run <id>               # a given run
```

`--regression` adds the headless regression shards (only when the coordinator asks: it costs ~1 h of runner time).
The runner does one GPU job at a time: batch your changes before dispatching.

Exit codes: 0 green; 3 red with only known failures (fine to report, name the issues); 1 red with NEW failures
(diagnose those, or file an issue and add a line to `tools/ci/known_failures.txt` in your PR once it's tracked).

Report the run link, the conclusion and the NEW lines. Don't re-diagnose a known failure: read its issue.
