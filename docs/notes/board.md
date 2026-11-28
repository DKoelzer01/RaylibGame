# Project board: commands and ids (coordinator)

Board: "RaylibGame Roadmap", https://github.com/users/DKoelzer01/projects/5 (project number 5,
owner DKoelzer01, id `PVT_kwHOAe7Hf84BrZzd`). Every open issue belongs on it, with all four fields set. The
board's auto-add workflow adds new issues itself, but sets no fields, so:
- **GraphQL budget (#582):** several sessions running `gh project` exhausted the user's 5000/h GraphQL budget. Feature
  sessions never run `gh project` commands: use REST (`gh api repos/DKoelzer01/RaylibGame/issues/N`, comments, labels,
  close, `search/issues`) and put board changes in your report as "#N -> Status/Priority/Size/Kind"; the coordinator
  applies them in batches (`build/board.sh set N -s/-p/-z/-k`). GraphQL only for `gh pr create`, never retried in a loop.
  The commands below are for the coordinator.
- **New issue:** right after `gh issue create`, set Status, Priority, Size and Kind. `item-add` either adds it or fails
  with "Content already exists" (auto-add was first); then look the item id up instead:
  `gh project item-add 5 --owner DKoelzer01 --url <issue url> --format json --jq .id` gives the item id, then
  `gh project item-edit --project-id PVT_kwHOAe7Hf84BrZzd --id <item id> --field-id <field> --single-select-option-id <option>`
  once per field. Find an existing item's id with
  `gh project item-list 5 --owner DKoelzer01 --limit 500 --format json --jq '.items[] | select(.content.number==N) | .id'`.
- **Status as work moves:** Backlog (not started) -> Ready (requirements clear, unblocked) -> In progress
  (starting work) -> In review (fixed/built in an open PR; also label `fixed in open PR`) -> Done (merged
  into develop, issue closed).
- **Priority:** P0 now (current plan / blocking), P1 next, P2 later, P3 someday. Follow the user's order in
  #158 (features) and #128 (graphics). Don't raise a priority on your own; propose it in the summary instead.
- **Size:** XS < 1 h, S ~ half a session, M ~ 1 session, L 2-3 sessions, XL a milestone of several
  (usually a tracking issue). Re-estimate when the scope turns out different.
- **Kind:** Feature, Bug, Performance, Tracking (an issue with sub-issues).
- Field and option ids (single-select):

  | Field | Field id | Options |
  | --- | --- | --- |
  | Status | `PVTSSF_lAHOAe7Hf84BrZzdzhpssRE` | Backlog `bf767b79`, Ready `2c8b23de`, In progress `4d99a4b6`, In review `ecd2f2aa`, Done `e280656f` |
  | Priority | `PVTSSF_lAHOAe7Hf84BrZzdzhpssSs` | P0 `e93b5bb8`, P1 `026eb06a`, P2 `2906a46d`, P3 `c36f5ea9` |
  | Size | `PVTSSF_lAHOAe7Hf84BrZzdzhpssSw` | XS `39a448a0`, S `032ae8bd`, M `6bd53109`, L `d7a23c47`, XL `ee6d351f` |
  | Kind | `PVTSSF_lAHOAe7Hf84BrZzdzhpssZg` | Feature `dcf904b0`, Bug `fae89327`, Performance `8172895f`, Tracking `d260109d` |

  If an id stops working, re-read them with `gh project field-list 5 --owner DKoelzer01 --format json`.
  ("Type" is reserved by GitHub, hence "Kind".)
- The gh token needs the `project` scope (added 2026-11-12). If a project command fails on scopes, ask the
  user to run `gh auth refresh -h github.com -s project`.

## Playtest board (coordinator)

- Page: https://claude.ai/artifact/MBu4PbuNrgXXP5khnkWjt9, source `build/playtest-board/index.html` (GROUPS by
  category). The user's marks live in its db collection `checks`, one doc per item id: `status` (pass/fail/skip),
  `reason`, `at` (ISO time of the user's last change).
- **The db is the baseline, not a local copy** (a local baseline missed 39 marks on 2026-11-22). After routing a
  mark (fail -> the user's words on the issue + a new issue for a new problem; pass -> `playtest confirmed`; a note ->
  feedback the same way), write `routed` = the doc's `at` into that doc (ArtifactData update, pinned to its version).
  A doc needs handling when `routed` is missing or differs from `at`.
- A new build: retests of items that failed get a fresh id (`1157-2`), passed items leave the page, and the old marks
  are deleted so the page starts clean.
