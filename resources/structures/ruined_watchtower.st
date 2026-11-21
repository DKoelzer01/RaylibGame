# A stone watchtower, half fallen: 4 x 4 m, 7 m high, a doorway, a chest inside (#164).
# Commands: part <type> x y z orient material | fill/walls <type> x0 y0 z0 x1 y1 z1 material | clear x0 y0 z0 x1 y1 z1 | loot x y z table
# Cells are 0.5 m; y is up. Blocks are 2 cells (1 m). Placed: turned by 90-degree steps, sunk `bury` metres,
# `decay` of the parts left out (more near the top), stone swapped for the biome's stone.
name Ruined watchtower
kind ruin
biomes Forest Meadow Savanna Taiga Tundra Alpine Badlands Desert Jungle
weight 1.0
bury 0.5
decay 0.35
dungeon 0.4
walls block 0 0 0 7 13 7 stone_brick
fill slab 2 6 2 5 6 5 oak_planks
clear 2 0 0 5 3 1
loot 3 0 3 ruin_common
