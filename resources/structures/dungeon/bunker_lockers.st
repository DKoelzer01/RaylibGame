# A colony bunker's locker room (#1641): two rows of lockers, iron panels with a glass pane above. Room 12 x 10 cells (6 x 5 m).
# The centre column (6, 5) stays free.
name Locker room
kind room
theme bunker
weight 1
room 12 10
fill panel 2 0 2 9 2 2 iron_plate
fill pane 2 3 2 9 3 2 glass
fill panel 2 0 7 9 2 7 iron_plate
fill pane 2 3 7 9 3 7 glass
# Supply crates between the rows: they break open for ruin_common loot (#1641 slice 2c).
part cube 3 0 4 0 oak_planks loot ruin_common
part cube 8 0 4 0 oak_planks loot ruin_common
part cube 4 0 5 0 oak_planks loot ruin_common
