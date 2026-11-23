# A crypt chamber (#1641): two rows of stone sarcophagi (a block base, a marble slab lid). Room 12 x 10 cells (6 x 5 m).
# Props sit at least 2 cells from the walls; the cell column at the room centre (6, 5) stays free up to y 3.
name Sarcophagus chamber
kind room
theme crypt
weight 1
room 12 10
fill block 2 0 2 5 1 3 stone_brick
fill slab 2 2 2 5 2 3 marble_brick
fill block 6 0 6 9 1 7 stone_brick
fill slab 6 2 6 9 2 7 marble_brick
# Urns between the rows: they break open for ruin_common loot (#1641 slice 2c).
part cube 3 0 5 0 marble_brick loot ruin_common
part cube 8 0 4 0 marble_brick loot ruin_common
part cube 9 0 5 0 marble_brick loot ruin_common
