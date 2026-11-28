# The Catalog's slots are the Inventory's size

**Problem (the user: "Catalog slots are larger than inventory slots"):** in `src/space/gameui.cpp`, `DrawCreative` (the
Catalog page) draws its item cells at a fixed `48.0f * GameUI_Scale()`, while the Inventory page (`DrawItemsPanel`) uses
the shared slot size `gGridCell` (set as the screen opens, #1148), which shrinks to fit the screen.

**Change:**
1. `DrawCreative` uses `gGridCell` for its cell size when it is set (greater than 0), else the old `48.0f * GameUI_Scale()`.
2. For the test, remember the cell each page last drew with: two file statics (`gDrawnCellInventory`,
   `gDrawnCellCatalog`) set in `DrawItemsPanel` and `DrawCreative`, and a function
   `void GameUI_TestDrawnCells(float& inventoryCell, float& catalogCell)` that returns them, declared in
   `src/space/gameui.h` (next to `GameUI_Scale`).

**Test:** the scenario `ui/catalog-slots` (already in the tree): at two interface sizes the Inventory page and then the
Catalog page are opened; both must draw with the same slot size.
