# A clearing reaches tree tiles that are still building

**Problem:** in `src/engine/trees.cpp`, `TreeSystem::addClearing` asks again only the tiles that are built
(`kv.second.ready`) and whose `centre` is within `r + kTreeTileSize` of the clearing. A tile still building has no centre
yet (it is all zeros) and its build captured the clearings as they were, so its trees grow inside the new clearing.

**Change:**
1. Add `bool TreeClearingReaches(const TreeTileKey& key, const TreeTile& tile, float planetRadius, Vector3 c, float r)`
   in `trees.cpp` (above `addClearing`) and declare it in `src/engine/trees.h` with the other tile helpers "exposed for
   tests". It uses the tile's `centre` when it has one (not all zero), otherwise the key's centre
   (`Vector3Scale(TreeTileCentreDir(key, planetRadius), planetRadius)`), and returns whether that is within
   `r + kTreeTileSize` of `c`.
2. `addClearing` requests every tile for which `TreeClearingReaches(key, tile, radius, c, r)` is true, built or not.

**Test:** `Trees.ClearingReachesATileStillBuilding` (already in the file).
