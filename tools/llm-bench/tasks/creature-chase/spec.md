# A creature keeps chasing its target out to its chase range

**Problem:** in `src/space/creatures.cpp`, `Brain_Step` picks a target among the players: predators and territorial
creatures want a hostile player only while `p.dist < c.aggroRange`, so a target that steps just past the aggro range is
dropped at once (and predators also skip any player at `aggroRange` or more). A creature should take a target within its
aggro range (`range=`) and keep chasing that same target out to its chase range (`chase=`).

**Change (in the player loop of `Brain_Step`):**
- A player counts as in range when `p.dist < c.aggroRange`, or when it is the creature's current target
  (`p.id == s.target && s.target != 0`).
- Predator and territorial `want` use that in-range test instead of `p.dist < c.aggroRange`.
- The final cut-off skips a player when `p.dist >= c.chaseRange`, for every aggro kind.

**Test:** `Creatures.ChaseKeptToChaseRange` (already in the file): not taken at 13 m (range 12), taken at 11 m, kept at
13 m and 29 m, given up at 31 m (chase 30), and not taken again at 20 m.
