# A colony bunker's barracks (#1641): rows of bunks (slabs; there is no bed part yet) with an iron panel divider between the rows.
# Room 14 x 10 cells (7 x 5 m); the centre column (7, 5) stays free.
name Barracks
kind room
theme bunker
weight 1
room 14 10
fill slab 2 0 2 5 0 3 white_plaster
fill slab 2 0 6 5 0 7 white_plaster
fill slab 8 0 2 11 0 3 white_plaster
fill slab 8 0 6 11 0 7 white_plaster
fill panel 2 0 4 5 1 4 iron_plate
fill panel 8 0 4 11 1 4 iron_plate
# Footlockers behind the divider: they break open for ruin_common loot (#1641 slice 2c).
part cube 3 0 5 0 oak_planks loot ruin_common
part cube 5 0 5 0 oak_planks loot ruin_common
part cube 10 0 5 0 oak_planks loot ruin_common
