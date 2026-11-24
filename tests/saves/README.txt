Golden saves (#1693): one directory per save version, tests/saves/v<N>/, each a copy of a save made by that version (world.sav,
world.sav.players/, world.sav.regions/) with facts.txt (what the world holds). The scenario saves/golden-load loads each through the
migration chain and checks the facts. Make a new one with saves/golden-write (see tests/scenarios/saves_golden.cpp).
No v13 yet: produced on the sandbox and committed.
