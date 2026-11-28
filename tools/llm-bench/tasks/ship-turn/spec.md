# The ship's autopilot turns toward a target dead behind

**Problem:** in `src/space/ship.cpp`, `Ship::update` turns the ship toward `want` with
`QuaternionFromVector3ToVector3(forward(), want)`. raymath returns an all-zero quaternion for opposite vectors, so with a
target dead behind the autopilot never turns, and the slerp decays the orientation's length.

**Change:**
1. Add a free function in `ship.cpp` (above `Ship::renderPos`):
   `Quaternion TurnBetween(Vector3 from, Vector3 to, Vector3 axis)`: the rotation taking unit `from` to unit `to`; when they
   are opposite (dot product below -0.9999) it is a half turn (PI) about `axis`, otherwise raymath's
   `QuaternionFromVector3ToVector3(from, to)`.
2. In `Ship::update`, use `TurnBetween(forward(), want, up())` in place of the raymath call, and normalize the slerped
   orientation (`QuaternionNormalize`) so it stays unit length.

**Test:** `Ship.TurnsTowardATargetDeadBehind` (already in the file).
