# A broken-open geode: a bowl 9 m across cut into the ground, quartz and amethyst crystals rising from its floor,
# and a cache of gems in the middle (#164: a natural wonder).
# carve x y z r: metres from the site's centre at ground level (x east, y up, z north), cut before the parts go in.
# Parts may go below the ground here (negative y cells): the set keeps their depth.
name Crystal geode
kind wonder
biomes Alpine Tundra Badlands Desert Taiga Forest Meadow Savanna
weight 0.5
bury 0
decay 0
carve 0 -2 0 4.5
fill post 4 -12 0 4 -9 0 quartz_block
fill post -4 -12 2 -4 -10 2 amethyst_block
fill post 0 -12 -5 0 -8 -5 quartz_block
fill post 5 -11 4 5 -9 4 amethyst_block
fill post -5 -11 -4 -5 -8 -4 quartz_block
fill post 2 -11 6 2 -10 6 amethyst_block
fill post -2 -12 -3 -2 -9 -3 amethyst_block
loot 0 -13 0 wonder_gifts
