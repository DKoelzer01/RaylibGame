# Site queries scan their whole range

**Problem:** in `src/space/structures.cpp`, `Structures_NearestSite` and `Structures_SitesNear` scan a grid of sample
points around the query position. The samples are half a tile (`kTile * 0.5f`) apart, but the grid's half-width `R` is
counted in whole tiles (`ceilf(range / kTile) + 1`), so a query only sees about half its range: a 3 km query scans
~1.8 km, and /locate, music and animals miss sites.

**Change:** in both functions, count `R` in half-tile steps so the grid covers the whole range. Nothing else changes.

**Test:** the scenario `world/site-query-range` (already in the tree): every site a 6 km query finds within 3 km, a 3 km
query finds too, and `Structures_NearestSite` agrees with the nearest of them.
