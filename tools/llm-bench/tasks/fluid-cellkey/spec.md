# Fluid cell keys with 21-bit fields

**Problem:** in `src/engine/fluid.cpp`, `Key(int i, int j, int k)` packs three coordinates into 11-bit fields (+-1024). It
is meant for brick coordinates, but `LevelSite` and `Fluid::PressureAt` also use it for *cell* coordinates, which reach
past +-1024 at a site's edge (a 256 m site of 0.25 m cells). The fields overflow into each other and different cells get
the same key (`Key(0, 0, 1030) == Key(0, 1, -1018)`), so the levelling pass merges unrelated cells.

**Change:**
1. Add `inline int64_t CellKey(int i, int j, int k)` next to `Key`: three 21-bit fields, each coordinate offset by
   `1 << 20` (i in bits 42+, j in bits 21+, k in bits 0+). Keep `Key` as it is (bricks use it).
2. In `LevelSite` and in `Fluid::PressureAt`, every key built from **cell** coordinates uses `CellKey` instead of `Key`
   (the `seen` sets, the `leafAt` map and its lookups). Brick keys stay on `Key`.

The file is large: use grep to find the `Key(` calls in those two functions.

**Test:** `Fluid.CellKeysStayApartPastTheBrickKeysRange` (already in the file).
