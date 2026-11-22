# Decisions and rulings (user)

The user's decisions by date and area, moved out of CLAUDE.md (#1478). Check here before changing a behaviour the user decided.

Player and building (2026-11-14)
- Player collision capsule 1.7 m tall (5 spheres at most 0.23 m apart; `RAYLIBGAME_PLAYER_HEIGHT` for A/B); 2 m
  doorways and tunnels must pass.
- Parts face automatically (from the hit face and view) until the player presses R/V.
- The Help tab is a Factoriopedia-style database; its "Discovered only" choice is saved per world.

Controls and tools
- Hotbar starts empty; an empty slot is bare hands (no "Hands" item). Ground items (dirt, sand,
  gravel, clay, regolith and raw stone: stone, granite, marble...) go into the hand and RMB places them as
  terrain; a tool in hand never places ground. Parts are built from bricks (2 of a stone -> 1 of its brick, by
  hand), planks, logs, concrete and the rest; loose ground and raw stone are not part materials (user, 2026-11-20, #1098).
- LMB breaks terrain and built parts; parts need a 0.45 s hold (also outside build mode). A part in
  front of the ground blocks the terrain tool (no digging/tilling/placing through structures).
- One swing per cooldown (`1/strokesPerSecond`, chops 0.3 s); the cooldown survives button releases,
  the hand animation lasts exactly the cooldown, and a bar under the crosshair shows it (amber when
  breaking a part).
- Digs: the carve sphere sinks half its radius along the view; the highlight is a depth-tested circle
  on the surface (the sphere's mouth), never a see-through sphere.
- Multi-cell parts placed against a face slide out until no cell is behind that face.
- Tools are held upright; hand bob = one side swing per 1.8 m walked.

Materials and progression
- Placed ground keeps its material: delta chunks carry a material channel (`DeltaChunk::m`), regions
  save as RGN4. Placing uses the filler in hand first.
- Stone tools, their repair and the tier-0 quest use stone; loose stones (tree-system variant 9) are
  picked up with LMB for 1-2 stone.
- Gems are rare: 1 in 3000 rock strokes (quartz 1 in ~1900); veins come later. Tests use
  `terrainTool.gemBoost`.

World and rendering
- Lakes = whole filled basins of the river flood, level held 3 m under the land sampled around them;
  water region = basin cells + 60 m margin; a capped (4 m) bank in the outer margin. Rivers in a lake
  build no banks. Chunks find lakes by BFS distance; dry water-lattice points sit under their ground.
- Ring rocks are cut at 50 km true distance (the mid pass's 60 km far plane cut them by depth).
- Clouds: a shell 1.5 km up in atmosphere.fs, in the planet frame, drifting with the wind; cover
  never below ~0.22 so a few are always out.
- Perf: chunks draw nearest-first (early-Z), CPU copies of chunk meshes are freed after upload.
  Scene::KickGpu (glFlush) at pass boundaries: without it CPU and GPU ran serialized (+25% FPS). GPU timer = elapsed queries.
  Reverse-Z (#135, scenetarget.*): the 3D world draws into an offscreen float-depth target (window sample count, glClipControl
  0..1, infinite far, GL_GREATER) in ONE pass (no far/mid/near depth ranges or clears); RAYLIBGAME_REVZ=0 = old passes.
  Terrain arena (#144, terrainarena.*): one multi-draw-indirect per body per pass; RAYLIBGAME_ARENA=0 = per-chunk draws.
  (Old note "remaining time is present/composition" was wrong: it was CPU/GPU serialization.)
- OpenGL 4.3 is required (user, 2026-11-13, #185): raylib is built for a 4.3 core context, so without it there is no window
  (main.cpp shows an error box). There is no 3.3 renderer; RAYLIBGAME_ARENA=0 / RAYLIBGAME_REVZ=0 / RAYLIBGAME_GPU_TERRAIN=0
  are A/B and debugging switches on a 4.3 context, not fallbacks for older hardware.
- Terrain lattice 0.25 m (kVoxelSize; RAYLIBGAME_VOXEL_1M builds the old 1 m one): digs carve exactly the
  previewed sphere (--test-dig-shape IoU ~0.85, pass bar 0.83: #425). Delta chunks hold sparse 4^3 bricks; regions are RGN5.
- Shadows (2026-11-12): two cascades (near +-32 m, far +-200 m, 2048 each), bilinear PCF, normal-offset
  lookup + slope bias; trees and bushes draw into both (tree vertex shader + depth.fs). Grass casts none
  by choice (blades thinner than a near texel).
- Foliage (user, 2026-11-12): grass = 3 instanced meshes (tuft, 14-blade clump, 24-blade dense clump;
  variant rides in scale + 10*variant), 7 per m^2, darker lower half, drawn to 100 m; trees to 2 km,
  bushes 300 m. The user trades a few FPS for view distance (walk bench ~220 FPS). Draw distances are
  runtime values (--grass-dist, --tree-dist); the terrain canopy tint follows treeDistance.

GL features and streaming (2026-11-13, stack #287-#293)
- The driver hands out a 4.3 core context with newer features as extensions; src/engine/glcaps.* (`GL_Caps()`)
  detects them at runtime and each path falls back when missing. `RAYLIBGAME_GL_LIMIT=43` forces the 4.3 paths.
- A/B switches: RAYLIBGAME_STREAM_OLD=1 (old streaming), RAYLIBGAME_GPU_READBACK_COPY=1 (old readback),
  RAYLIBGAME_UPLOAD_STAGING=0 (direct uploads), RAYLIBGAME_TREE_CULL=0 (no GPU tree culling),
  RAYLIBGAME_INDIRECT_COUNT=1 (opt-in multi-draw-indirect-count with chunk_compact.comp: measured no gain).
- Streaming: <= 96 requests per frame, frees capped at 256/frame, nearest-first uploads, telescope keeps the
  zoom-1 selection resident. GPU readback = persistent mapped ring; uploads go through a 16 MB persistent write
  ring with per-upload fences.
- Preload ring + chunk cache (#262: #294-#296, 2026-11-13): every node the LOD would pick from a camera within
  300 m of the player along the ground (0-8 m over it; nodes touching +-12 m of the height map) is built and kept
  resident, not drawn. A new or loaded game shows "Preparing spawn" until it is built (~8 s, ~38k nodes,
  ~27 M vertices = ~2 GB VRAM on seed 42); afterwards it re-centres every 75 m in spare time only (the view's
  unused request budget, a short queue, 1 ms of uploads, one 8-block GPU batch a frame) and pauses above
  15 m/s. A ball-shaped ring reached 300 m down into caves (415k nodes, 80 M vertices) - keep it on the ground.
  Unwanted built nodes go to an LRU cache (8 M vertices / 40k nodes, ring and view not counted); the arena is
  grown once behind the spawn screen so it never regrows (GB copy) mid-walk. WorldIdle() ignores ring builds.
  The cache also evicts while the arena is over 85% full, and a ring install reserves nodes x 800 + cache + 4 M
  vertices up front (#300: a late regrow was a 30 ms frame closing the telescope).
- Terrain detail (#297): a settings choice Low/Medium/High/Ultra = LOD factor 2 / 2.5 / 3 / 4, default Medium. It
  only sharpens mid-range terrain (the planet is selected to the horizon already). Above 15 m/s against the nearest
  body (back under 10) the factor drops to 2: at 3, flying 60 m/s cost 4.5 ms GPU streaming a frame (115 FPS).
  Seed 42 at Medium: still 267 / walk 250 / fly 262 FPS. Benches: `--lod F` sets the detail, `--lod-fast F` the
  flying factor (99 = no drop); benches do not read settings.txt.
- Tree shadows: caster sets per cascade chosen by distance to the sun axis (1.26 -> 0.09 ms); main view
  culled on the GPU (tree_cull.comp + glDrawArraysIndirect).
- Fly camera follows the body whose sphere of influence it is in, rotating with it within 2 radii (#127).
- Gameplay in the same stack: boilers need fuel + water (pumps fill them), generators are demand driven;
  fall damage safe below 6 m/s (8.5 with servo legs); animals walk on the density surface with ledge
  reach/detours. Ravines are generated (voxel.cpp), not streaming holes. There is no start ship.

Water (Plan 4, user answers)
- Sea, rivers and large lakes (> ~1 km^2) are infinite; rivers can be diverted; small lakes/ponds
  are finite once disturbed. The sea can flood a channel dug inland below sea level.
- Multi-level (cave) water must work. Buckets AND pipes/pumps. All of F1-F6 before the next playtest.
- Decided (2026-11-11): F2 is a sparse 3D cell simulation, Terraria-style, 1 m voxels (fill amount
  per cell; down, then sideways; a slow levelling pass for U-bends / cave pressure), surface meshed
  with marching cubes on the fill amounts. Replaces the per-column span design.
- Fluid cells are 0.25 m by default (user, 2026-11-17, #714); `RAYLIBGAME_FLUID_CELL=1` is the old 1 m grid (opt-in).
  A save records its cell size; sites saved at another size are dropped on load (their simulated water starts over).

Rulings 2026-11-19 to 2026-11-21 (user)
- Order: Plan 11 body frames (#667) is P0; the research/tech tree (#1030) follows it. Combat 5-8 (#1219-#1222) are built in
  the `ladder` branch now; the ladder lands in m84 once, after the tiered machines and the user's playtest.
- The ladder (#1271): tiers 0-7 rewritten; Progression builds the anvil (shared with combat 5); research gates hand
  recipes, machine placement AND machine recipes; machines removed by the rework are refunded into chests; tiers 8+
  follow Plan 13 (8 Orbit+Luna, 9 Belt, 10 Tempest+fusion; Space writes them). The macerator is core at tier 4; tier 4
  adds an LV assembler.
- Tiered machines (#1427, part of the ladder rework): every machine gets a variant for each tier after the one it's
  introduced in, each with its own recipe (not crafted from the previous tier's machine). A tier unlocks its recipes;
  lower-tier recipes run in higher-tier machines, overclocked x2 speed for x3 EU/t per tier step, no cap. Non-electric
  machines serve tiers 1-3; from LV, electric tiered successors replace them.
- Disassembly: an assembler mode, one recipe level down (direct inputs). Yield by the item's tier vs the assembler's:
  same 50%, one above 10%, two above 1%, more 0; one below 75%, two or more below 100%. Whole parts (circuits,
  components) come back as a % chance per unit; raw resources come back exactly, scaled by the yield.
- Multiblocks (#1429, later): built from many subcomponents (newer GregTech / Immersive Engineering style), not "a bigger
  cube" except where that makes sense.
- Alchemy: one generic essence fluid with the aspect tagged per container/network; spills fizzle (no world spill yet,
  #1414). No ley thread (#1249).
- Ladders (#1105): made from planks; W climbs; no stamina cost.
- Accessibility (#1017): colour-vision setting and cues (#1432) -> simulated colour blindness check (#1433) -> controllers,
  which wait for Steam Input (#1434).
- Slideshow (#1071): UI/UX builds hand-placed shots. #1170 horizon culling: closed as measured (<= 3% gain).
- Bench/scenario runs should flag frames where the camera or player is inside solid (#1442).
- Weapons (#1219): per-item instance data first (#1446: parts, quality, sockets, wear on unstackable items), then the
  weapons. **All tools are composed items too** (#1449): one item per class (tool_pickaxe...), with head/handle/binding
  materials, anvil quality, sockets and wear on the copy; mining tier from the head material. Disassembling a composed
  item gives back its parts at the tier yield. Factory slots (#1446 step 2b): chests carry items with data for players;
  belts, inserters and machines refuse them, except an assembler in Disassemble mode.
- Machines with several recipes get a recipe picker ("Makes": Any or one output), saved per machine (L6c, #1459).
- Underground plan (#1465, after the ladder and Plan 11): metals at full GregTech scale (+1-2 magic metals), real-ish
  geology with vein families, ruins from two eras (a lost colony; an older native civilisation), lava as a simulated
  fluid (lava + water -> obsidian), hazards = gas pockets, heat with depth, flooding (no cave-ins). Owners split:
  caves/structures/ore generation Lighting, lava/aquifers Water, metals Progression.
- Coordination (#1476): the machine slot is self-serve through tools/slot.sh with hooks; CLAUDE.md holds core rules
  only; auto-merge on "CI ok"; known runner failures listed; protocol numbers reserved in a file.

Deferred by the user: moon mining launch (#116), motion sickness (#117; comfort options are #168). The history of
the earlier playtest batches is in git history and the closed issues.

## Art pass (branch art-pass-5, 2026-11-12, Plan 5 light items)

- Shaders take `#include "file"` (HotShader expands it, relative to the shader, and watches the file).
  `resources/shadow.glsl` = the two-cascade sun shadow lookup + headlamp for lighting/trees/grass/parts/water.
  Shadow maps use hardware compare (sampler2DShadow, linear): never sample them as sampler2D (debug
  previews call SetShadowCompare(id, false) first). Near cascade blends into the far one over 8% of its box.
- Parts receive shadows + headlamp (ctx.bindSunLighting binds slots 6-7 before plot and water draws).
- Part surface patterns: id in the part mesh's vertex-colour alpha (plots.cpp PartPatternOf), drawn in
  section-local plot metres by parts.fs (stone, brick, planks, metal, glass, tiles, marble, gems, logs,
  lamp glow, machine casing, conveyor, cable, farmland, hull).
- Item icons: src/space/itemicons.* paints every item at startup (one atlas); --shot-icons.
- Particles: src/engine/particles.* (body-local, pooled, rlgl billboards); emitters in terraintool
  (dig/place/chop/bucket), rocket.cpp Exhaust, src/space/ambience.* (rain splashes, fireflies, machine
  smoke, campfire flame); --shot-particles.
- Lighting review leftovers for the user: tonemapping (ACES), hemisphere/SH ambient + SSAO, GGX specular.
  Staggering the far cascade was rejected (planet spin makes a stale map shimmer).
- Don't run screenshot benches while another game instance runs (e.g. a background bisect): chunk
  loading slows and shots show a half-loaded world.

## Plan 5 answers (user, 2026-11-12) - all approved, order B -> C -> A -> machine models -> D
- A: procedural code-built animal models (glTF slot later); predators only chase (no kills yet);
  machine models after the animals.
- B: all 10 biomes (jungle, savanna, desert dunes, badlands/mesas, temperate forest, meadow, swamp,
  taiga, tundra, alpine); regional size (2-5 km, several near spawn); save-breaking is fine.
- C: dramatic terrain tied to biomes, with a good variety of flat building spots as well as
  interesting ones; floating islands rare; ravines and cave mouths in the same pass.
- D: some PBR; stylized look - higher resolution than Minecraft but fairly unrealistic.

## Generic fluids (user, 2026-11-12) - on m25-biomes / PR #30
- One fluid model, not per-fluid copies: kinds come from resources/data/fluids.txt (water = 0, crude_oil;
  room for 8, e.g. lava later: colour, opacity, glow, viscosity, min spread, bucket item). The kind rides in
  cell flag bits 5-7; kinds never mix (a fluid enters only empty cells or its own kind). Meshes per kind,
  drawn by water.fs (fluidColor.a > 0 = a non-water kind). Saves: fill + 256 x kind; sources on a "springs" line.
- Sources (Fluid::AddSource): a rate (m^3/s) and a head; fill the column's surface up to the head, then rest
  (equilibrium, bricks sleep) or overflow into a stream. Slow sources add in drops of >= 8/255 (thinner
  specks soak away: kEvaporate). Place them on the real ground (probe density from ~8 m up: the mesh can
  stand metres above the height map).
- Natural sources (ambience.cpp kSourceTypes): oil seeps in warm biomes (carved basin, head under the rim:
  pools) and water springs in wet hill biomes (head 1.2 m over the ground: streams). --test-oil, --test-spring.
- Follow-up ideas from the user: express rivers/lakes (and later magma) as sources too; springs show as
  thin films on steep scree (consider a carved channel or a lower film limit for springs).
