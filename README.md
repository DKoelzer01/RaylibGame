# RaylibGame
A space open-world survival prototype on a from-scratch raylib voxel planet engine (marching
cubes, Simplex fBm, octree LOD, shadowed lighting). You start on a 60 km home planet with a moon and an
asteroid ring overhead.

## Layout
- `src/engine/` - the reusable engine library (`engine` target): density fields, meshing/LOD, job
  system, player controller, inventory, audio, scenes, debug UI.
- `src/space/` - the space survival game (`game` target): ship, survival, save/load, HUD, `main.cpp`.
- `resources/` - shaders and data (`resources/data/bodies.txt`, `items.txt`, `recipes.txt`).

## Build (CMake)
Needs CMake 3.24+ and a GCC/MinGW toolchain on PATH (e.g. `C:/raylib/w64devkit/bin`).
raylib 5.5 is fetched and built automatically on first configure.

```sh
cmake --preset release          # or: debug, profile (Tracy)
cmake --build --preset release
./game.exe                      # the exe is written to the repo root (assets are loaded relative to it)
```

Options: `-DUSE_SYSTEM_RAYLIB=ON -DRAYLIB_PATH=C:/raylib/raylib` links a prebuilt raylib instead of
fetching; `-DENABLE_TRACY=ON` (the `profile` preset) builds the Tracy v0.11.1 client; connect with
the matching `tracy-profiler` GUI.

VS Code: the CMake Tools extension picks up `CMakePresets.json`; `compile_commands.json` is
generated in the build directory for IntelliSense/clangd.

## Controls
| Key | Action |
| --- | --- |
| Mouse | main menu: New Game, Load Game (saves in `saves/`), Settings (mouse, FOV, volume, foliage, frame cap), Exit |
| Mouse / arrows | look |
| WASD, Space, Ctrl | move (Shift/Alt = faster) |
| LMB / RMB (hold) | dig / fill terrain at the crosshair (filling costs stone, dirt, gravel or sand) |
| Tab | inventory + crafting |
| G | build mode (wheel shape, R / V rotate, RMB place, LMB remove; B saves the plot you look at as a blueprint, N cycles blueprints: RMB places what you can afford) |
| M / , . | aboard a rocket: orbit map / time warp down, up (1x-100x, coasting above the air only) |
| 1-9, 0 | hotbar: the terrain tool, a building material (any shape costs it), a machine, or a consumable (oxygen canister, power cell). In the inventory, hover an item and press 1-0 to bind it |
| L | headlamp (on foot: suit power; piloting: ship nose light) |
| F | toggle free-fly camera (WASD along the view, Space/Ctrl along camera up, mouse turns, Q/E roll, Shift/Alt faster) |
| Space | jump; hold for jetpack |
| E | board / exit the ship; board a rocket (looking at its cockpit) / leave it once landed |
| Mouse wheel | brush size (zoom while the telescope is up) |
| T | raise / lower the telescope |
| F5 / F9 | save / load (the current save file) |
| Esc | pause menu (Resume, Save, Settings, Save & Quit to Menu); Esc again resumes |
| F1 | debug panel (frame graph, streaming stats, sun/shadow/texture tweaks, brush, world regen by seed) |
| F3 | debug overlay (light, shadow box, shadow map, chunk bounds) |

Shaders in `resources/` hot-reload when saved.

## How the terrain works
- Each body has an immutable density field (`src/engine/voxel.*`): asteroids a `PlanetoidField` (seeded fBm density with a radial
  falloff), planets and moons a `PlanetField` (below); each also gives the material at a point.
- Terrain edits (`src/engine/terraindeltas.*`) are spherical density brushes folded straight into a lattice of
  density deltas (1-voxel spacing, 32^3 chunks of 16-bit values, copy-on-write snapshots shared with the workers).
  Density = base field + the trilinear delta, so builds cost the same after a million strokes. Saves keep the
  deltas in region files next to the save (`<save>.regions/`, 16^3 chunks per file, run-length coded), written
  incrementally: only regions changed since the last save. The game autosaves every 5 minutes. Older saves (brush
  lists) are folded into deltas on load.
- Chunks (32^3 cells) are built on a worker pool (`src/engine/jobsystem.*`): sample density, marching cubes,
  normals from the density gradient, material weights in vertex colours. The main thread only
  uploads finished meshes, time-boxed per frame. Edits re-mesh just the touched chunks.
- Each planetoid is an octree of LOD nodes (`UpdateLod` in `src/engine/world.cpp`): a node at level L covers
  32*2^L voxels with 32 cells, so near terrain is fine and far terrain costs a handful of draw calls.
  Nodes overlap their +x/+y/+z neighbours by one cell to hide cracks between levels, and old nodes
  stay on screen until their replacements are built (no holes while moving).
- World positions are doubles (`DVec3`); the scene keeps a floating `renderOrigin` near the camera,
  so rendering stays precise millions of units out.


## Gameplay (space survival prototype)
- **Player** (`src/engine/player.*`): walks on planetoids with gravity toward the nearest body; "up" follows
  the surface so you can walk all the way around. Collision is resolved directly against the density
  field, so tunnels you dig are walkable. Jump (buffered, with coyote time), WASD air control while airborne
  in gravity, jetpack on Space/Ctrl, 6-DOF suit thrusters in open
  space, and you ride along on spinning bodies. **F** toggles the free-fly dev camera.
- **Items & crafting** (`src/engine/inventory.*`, `resources/data/items.txt`, `recipes.txt`): mining yields
  items by voxel material (rock, regolith, ice, ore); building consumes stone/regolith/lunacrete.
  **Tab** opens inventory + crafting. Add items/recipes by editing the text files.
- **Survival** (`src/space/survival.*`): suit oxygen drains outside the ship's life-support field and a planet's breathable air (~5 min a
  tank), suit power pays for recharging the jetpack and the headlamp, health drops when suffocating or from hard landings and
  crashes. Death respawns you at the ship. **1** uses an oxygen canister, **2** a power cell (charges
  the ship instead while piloting). An objectives list guides the loop: dig stone, find iron ore, craft a power cell,
  board the ship, fly to the moon, visit the asteroid ring.
- **Ship** (`src/space/ship.*`): **E** to board/exit. 6-DOF flight (WASD, Space/Ctrl, mouse, Z/C roll,
  Shift boost, X brake) with flight assist, chase camera, terrain collision, crash damage, and an
  energy store. Once at rest on the ground it anchors to that body and rides its spin (position and
  heading) until you thrust off, keeping the surface's velocity. It feels gravity: piloted with energy,
  flight assist hovers against it (costing energy like thrust); unpiloted, it falls, lands and parks.
- **Gravity** (`GravityAt` in `src/engine/world.*`): summed over every body, inverse-square outside a
  body's mean surface and linear to zero at its centre inside. "Up" is against it, and is held where
  gravity is too weak to define it (deep space, the balance point between bodies).
- **Star system** (`src/engine/bodies.*`, `resources/data/bodies.txt`): one body (asteroid, planet, moon)
  or procedural ring per line, with radius, position or orbit, spin or day length, gravity and seed. The
  default is the 60 km home planet, its moon, an asteroid ring 50 km up and three huge rocks in and
  beside it; the old asteroid start is `bodies_asteroids.txt`. Saves embed the system they were
  made in, so editing the file never breaks a save. Errors are reported with line numbers in `log.txt`.
- **Rockets** (`src/space/rocket.*`, `resources/data/ship_parts.txt`): a build plot with a cockpit is a rocket. Build it
  freely from rocket engines, fuel tanks, a guidance computer, landing legs and any structure (parts' mass, thrust,
  engine efficiency and tank sizes are data; structure weighs by material). Looking at it shows the launch check
  (thrust-to-weight above 1.2 here, delta-v for orbit, fuel loaded, what is missing); I loads rocket fuel, E boards.
  Throttle up (Shift; Ctrl down, X cuts) and the plot lifts off as a vehicle: summed gravity, drag against the
  turning air, fuel burnt at the engines' efficiency, in double precision, so coasting above the air is a real
  orbit; W/S, A/D, Z/C steer (the guidance computer holds attitude). The HUD shows altitude, speeds and the orbit
  (apoapsis / periapsis; in orbit once the periapsis clears the air at 35 km). A gentle upright touchdown turns it
  back into a plot on whatever it landed on; hard landings break legs, harder ones engines (damage, not
  destruction). The cockpit is pressurised. Saves keep a rocket in flight. Rocketry is tier 7 of the quest book: parts
  need titanium alloy, aluminium and advanced circuits, and fuel comes from the refinery.
  M opens the orbit map (the orbit drawn from its elements around the body you orbit, apoapsis / periapsis, the air,
  the moon, time to apoapsis); , and . warp time 1x-100x while coasting above the air (any burn, the air or ground
  within 2 km drops it). Fast flight in thick air heats the leading face (re-entry: ~800 C coming down from orbit): the
  part there burns away if it stays past its limit (wood 250 C, glass 900, the rocket parts' alloys 1100-1600, heat
  shield tiles 3000), shielding the parts behind it; the cockpit scorches but holds. The landing autopilot aerobrakes
  above 12 km. Blueprints (resources/blueprints/, your own in blueprints/): a starter rocket ships as a preset.
- **Planetary ring** (`src/engine/ring.*`): 60,000 rocks 5-50 m across (about 1% five times bigger) in a band
  3 km wide and 60 m thick, each on its own circular orbit, derived from the seed and its index (nothing
  stored until you change one). Far away rocks draw as instances of a dozen lumpy rock meshes (skipped when
  under a pixel), and beyond ~40 km the ring is a translucent band in a faint dust halo, dark in the planet's
  shadow, an arc across the sky from the ground. Near a rock (within a reach growing with its size) it becomes
  a full voxel body: dig it, stand on it (low gravity of its own); leave and it turns back into an instance,
  keeping its edits and build plots, which are saved (save v10). Huge rocks orbiting inside the ring clear a
  lane in it (rocks on crossing orbits would collide); the default has goliath and bastion at its edges and
  pebble in a gap in the middle.
- **Audio** (`src/engine/audio.*`): procedural sound effects, no assets.
- **Planets and moons** (`PlanetField` in `src/engine/voxel.*`): a height map over the sphere (continents and
  ocean floor, ridged mountains, hills, river valleys), worm caves, ore veins that get richer with depth, snow
  at high latitudes and on peaks, and an unmineable mantle 1.5% of the radius below sea level (900 m on the
  60 km planet). Coarse LOD chunks leave out detail finer than their sample spacing; nodes over the
  horizon, above the highest terrain or inside the mantle are never built. A sea fills the planet to sea level
  (meshed per chunk where it meets open air; polar seas are frozen into ice shelves), and the planet has an
  atmosphere: sky colour by sun angle and altitude, dusk glow, stars at night and in space, haze with
  distance, and dense blue fog underwater. It is the default start system: you start on temperate
  grassland (picked from the seed), where the air is breathable below ~15 km (the suit refills its oxygen).
- **Build plots** (`src/engine/plots.*`, `src/space/buildtool.*`): G toggles build mode. Parts snap to a 0.5 m grid on
  a plot anchored where the first part goes (up = local gravity, north along the meridian); the hotbar picks the
  material (concrete, stone, wood, iron plate, glass, packed dirt) or a machine, the wheel picks the shape (cube,
  1 m block, slope, corner or a 1/2/4 m beam; any shape from any material), R/V rotate (24 orientations), RMB
  places, LMB removes (refunded). A shape costs its material by size; parts lay their own foundations (terrain
  cleared inside their cells only, small gaps to the ground filled), merge into one mesh per 8 m section with
  hidden faces culled, cast shadows, collide with the player, and are saved with their material.
- **Inland water** (`PlanetField::waterElevation`): rivers in a channel along the river valleys, smaller streams
  on a finer network, and ponds (12-38 m, flat surfaces) on gentle ground; beds are sand and gravel, and digging
  them also turns up clay. Undug grassland facing up is grass all the way through its soil; planets yield dirt
  (moons and asteroids regolith).
- **Fly camera** (F): rides its reference body, the nearest whose sphere of influence holds it: its spin and
  motion within 2 radii, its motion only further out.
- **Factory simulation** (`src/engine/factory.*`, `resources/data/machine_recipes.txt`): machines are build parts
  (chest, combustion generator, electric furnace, macerator, battery) joined by conveyor and cable parts. A
  deterministic simulation runs at 20 ticks a second on its own thread (integer items and energy): conveyor networks
  move outputs to machines that take them, cable networks share generator and battery power (a shortfall slows
  every consumer equally; a machine below the network voltage refuses power). Plots over 1 km away tick once a
  second with a 20x step. Look at a machine for its state; I puts items in, O takes the output. Saved (v9).
- **Tiers and progression** (`resources/data/materials.txt`, `tiers.txt`, `machines.txt`): materials (iron, copper,
  tin, bronze, steel families) generate their ore, dust, ingot and plate items; ore veins come in kinds by 48 m deposit
  (iron, copper, tin, coal). Machines, their power (none, fuel, EU at a voltage, generator, battery) and fuels are data.
  The quest book (top left) walks the tiers, each gated by an item or placed machine: a campfire fires a clay crucible
  (tier 0), a kiln smelts copper, tin and iron and alloys bronze (1), a boiler drives a forge hammer by cable and a coke
  oven and blast furnace make steel (2), glass and vacuum tubes (3), circuits, LV hulls, an LV generator and the electric
  furnace (4). Fuel machines run recipes that don't eat their fuel first.
  From LV, ore processing multiplies ore the GTNH way: a macerator crushes it (x2), an ore washer (water from a water pump)
  purifies it, a centrifuge splits it into dust and a byproduct (iron: nickel, copper: gold, the rest: iron), and the
  electric furnace smelts the dust. MV (5, 128 V): an MV generator, an electrolyzer (water into hydrogen and oxygen;
  bauxite dust into aluminium, the gate), a chemical reactor (polyethylene for good circuits) and a miner that works the
  deposit under it (bauxite and ilmenite deposits join the ores). HV (6, 512 V): an HV generator, an arc furnace
  (titanium from ilmenite and coke, then titanium alloy: the gate), a refinery (rocket fuel) and a precision assembler
  (alloy plates, advanced circuits). A machine on a network above its voltage refuses the power, so each tier runs its
  own cables and generator.
- **Distance sphere** (`src/engine/farsphere.*`): beyond 3 radii (telescope zoom included) planets and moons draw as a
  sphere with colour and normal maps baked on the workers from the height map and biome rules (sea by depth, forest
  canopy), lit by the sun with the night side and a sea glint, under the same haze as the terrain.
- **Foliage** (`src/engine/foliage.*`, `trees.*`): instanced grass tufts on grass terrain (from each finest chunk's mesh,
  to 70 m), and trees and bushes built from prisms, placed on a 64 m cube-sphere tile grid from the seed (a
  forest field, grass or dirt, gentle slopes, solid ground, below the tree line; broadleaf in the temperate band,
  conifers further north and up high), drawn instanced to 1 km (simpler meshes past 200 m). Trunks collide.
  LMB on a tree chops it (wood), on a bush harvests it (fibre, sometimes berries); plants regrow from a sapling
  after 10 minutes and are full size 20 minutes later; removals are saved. Bench runs save to `saves/bench.sav`.
  Beyond the trees' 1 km range, the terrain paints in forest canopy from the same forest field (forests to the
  horizon and from orbit); forests are densest in the tropics and thin out toward the tundra. The start has a
  14 m clearing for the player and the parked ship.
- **Telescope** (`src/space/telescope.*`): you start with one (craft more from glass and iron plate). T raises it, the
  wheel zooms 4x-64x; terrain in view refines as if that much closer (`SetLodZoom`). Bodies more than 20 km away
  draw in a separate far pass with their own depth range, and haze only counts the air the view ray passes
  through (an 8 km layer), so the moon and the ring are visible from the ground.
- Asteroid density is a displaced sphere (mean surface 0.72 R) with a uniform 0.25/voxel gradient: crisp surfaces and
  predictable dig depth.
## Benchmark / test mode
```sh
./game.exe --bench 10                       # 1920x1080 uncapped; waits for the world to load, prints frame stats
./game.exe --bench 5 --seed 42 --size 200 --cam 400 200 -400
./game.exe --bench 2 --seed 1553044913 --cam 60 45 -50 --test-dig 10 --test-saveload
```
Writes `bench.txt` and `bench.png`. `--cam` = camera offset from the first planetoid, `--count` =
planetoid count, `--system R` = belt radius, `--lod F` = LOD factor, `--first x y z` = first body position
(`--size`, `--count`, `--system` and `--first` override the first body and first ring of `bodies.txt`),
`--bodies file` uses another system file,
`--dump-system` prints every body, `--load` loads `saves/world.sav` at start, `--fps N` caps the frame rate,
`--bench-walk` walks the player through the measured window (frame time while terrain streams), `--crack-check`
clears to magenta without a skybox and reports the share of pixels showing through the terrain (look down at it),
`--ui` opens the debug panel, `--overlay` the F3 overlay, `--nocull` disables
frustum culling, `--test-dig N` carves a tunnel at the crosshair, `--test-saveload` saves and reloads,
`--test-walk S` / `--test-stand S` spawn the player and walk/stand (reports grounding + penetration),
`--test-planet-walk M` walks M metres on body 0 and prints PASS/FAIL (milestone 1 gate: 1000),
`--test-hole` (walk into a fresh pit, build a pillar underfoot), `--test-trees` (walk into a trunk), `--test-chop` (fell a tree, save/reload, regrow), `--test-persist N`
(milestone 2 gate: N strokes, flat builds, save and load under a second), `--test-build N` (milestone 3 gate:
N parts, standing on them, save/load, frame time with the factory in view), `--test-factory N` (milestone 4 gate:
N working machines, tick under 5 ms, plates reach the chests, contents survive a reload), `--test-progression`
(milestone 5 gate: a scripted playthrough from an empty inventory to a placed, powered, working LV electric furnace;
`--test-progression full`: the milestone 8 gate, on through ore processing, MV and HV to titanium alloy, a miner and the
refinery),
`--test-wall` (walk into shaft and tunnel walls: no launches), `--test-menus` (main menu, new game, pause, save & quit,
load), `--test-water` (inland water near the start, its bed and clay), `--test-flyfollow` (the fly camera rides the
planet and the moon), `--test-hotbar` (materials, consumables, foundations, materials through a save), `--test-ring`
(milestone 6 gate: 60k rocks, frame time from the ground and inside the ring, a rock promoted, dug, left, revisited
and reloaded with its hole, standing on a huge rock), `--test-blueprint` (place the starter rocket preset, fill in, save and load), `--test-rocket` (milestone 7 gate, with `--time-scale 8`: a rocket built,
fuelled and flown by the autopilot to orbit, saved and loaded in orbit, brought down and landed), `--bench-fly V`
/ `--bench-fly-alt A` (camera flight over body 0, LOD streaming; negative A: underwater), `--no-foliage`,
`--test-mine N` digs N strokes at the feet, crafts everything possible and builds 3 strokes,
`--test-fly S` autopilots the ship to the nearest other asteroid, `--test-suffocate` checks death/respawn,
`--test-moon` (moon on its analytic orbit; standing and jumping on it), `--time-scale k` (bodies spin and orbit
k times faster), `--test-mantle` (planet terrain: height band, mantle refuses digging, ore by depth, caves; on the
default system), `--test-ride` (standing/jumping on a spinning body; `--ride-spin D` sets its rate), `--test-orbit S` (jetpack orbit stability),
`--test-park` (parked ship rides a spinning body), `--test-gravity` (field values, continuity between
bodies, an unpiloted ship dropping and parking, a piloted ship hovering), `--test-precision` (float vs double
placement error standing on body 0; run it with a large body via `--bodies`), `--spin-test` (writes `spin_*.png` frames),
`--sun x y z` sets the sun direction, `--headlamp` / `--inventory` start with the lamp on / inventory
open, `--look-down` tilts the view to the ground (lighting checks).

## Third-party
See `third_party/README.md` (Dear ImGui v1.92.7, rlImGui tag Raylib_5_5).
