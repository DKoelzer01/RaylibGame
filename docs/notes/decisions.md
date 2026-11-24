# Decisions and rulings (user)

The user's decisions by date and area, moved out of CLAUDE.md (#1478). Check here before changing a behaviour the user decided.

Player and building (2026-11-14)
- Player collision capsule 1.7 m tall (5 spheres at most 0.23 m apart; `RAYLIBGAME_PLAYER_HEIGHT` for A/B); 2 m
  doorways and tunnels must pass.
- Parts face automatically (from the hit face and view) until the player presses R/V.
- The Help tab is a Factoriopedia-style database; its "Discovered only" choice is saved per world.
- Building (user, 2026-11-24, #1771): a part placed on bare ground starts a new plot at the exact aim point, facing the
  player's heading (or turned with the rotate key), not the body's cardinal grid. Near an existing plot, placement snaps
  to it as before.

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
- Grass under water (user, 2026-11-24, #1774): grass that stays under water turns to dirt. Grass, then tall grass, regrow
  after a delay once the water is gone. Trees never regrow under water.

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
  item gives back the **raw materials** its parts were made of (ingots, planks...), scaled by the tier yield - not the
  parts (user, 2026-11-22). Factory slots (#1446 step 2b): chests carry items with data for players;
  belts, inserters and machines refuse them, except an assembler in Disassemble mode.
- Machines with several recipes get a recipe picker ("Makes": Any or one output), saved per machine (L6c, #1459).
- Underground plan (#1465, after the ladder and Plan 11): metals at full GregTech scale (+1-2 magic metals), real-ish
  geology with vein families, ruins from two eras (a lost colony; an older native civilisation), lava as a simulated
  fluid (lava + water -> obsidian), hazards = gas pockets, heat with depth, flooding (no cave-ins). Owners split:
  ore kinds/vein families/caves/gems Lighting, structures (#1468) and hazards (#1472) UI/UX, lava/aquifers Water,
  metals Progression.
- 2026-11-22: ore patterns per kind for colour vision (#1464) as shown; dry river channels (#1438): fade detail and
  rock features inside channels, after Plan 11 F5, together with wider rivers sized to their valleys + a water-cycle
  slider (#637); the ladder merges into m84 after L7 and composed tools T1 (combat 5-8 continues on m84); stations
  become Fixed(station) frames in F5 (with a free-fall term); the asteroid ship's flight moves to frames in F8.
- Coordination (#1476): the machine slot is self-serve through tools/slot.sh with hooks; CLAUDE.md holds core rules
  only; auto-merge on "CI ok"; known runner failures listed; protocol numbers reserved in a file.

Rulings 2026-11-22 (user)
- Machine components (#1614, sub of #1271): GregTech-style layers (microcraft components -> assemblies -> machines),
  all doable by hand (1.5-2x the machine's metal) or by machine; fair cost, no grind; automation and magic favoured over
  brute force; exploration and quests give shortcuts (loot parts, permanent team-wide blueprints ~25% cheaper per
  assembly, research notes up to half a branch node / 25% of a core node). 7 shapes (plate, rod, wire, gear, screw,
  spring, cable); one alloy per tier (brass T3, solder T4, electrum + invar T5, stainless T6, tungsten steel T7);
  migration one pass per tier. **Budget is report only**: CI prints the per-machine minutes table, never fails on cost
  (structure still fails: hand recipe per step, plannable from raw, no cycles). **Yield order: magic + tech together
  best, tech alone = magic alone, hand worst** (a structure test fixes the order).
- Multiblocks (#1429): the controller projects a see-through ghost, the player places the real parts and forms it with
  a wrench or hammer (no auto-build); placed parts stay visible and animate; fixed shapes, a few (tanks, boilers,
  distillation towers, furnaces) stretchable. All 10 design recommendations accepted: tier = lowest of casings and energy
  hatch; up to 4 parallel; break-to-unform; Craft missing via the planner; no maintenance; first six: large tank, large
  bronze boiler, coke oven battery, EBF, large ore line, distillation tower.
- Conduits (#1620, sub of #1429): 1/8-cell sub-grid, straight/right-angle runs, no crossing; amps by cross-section
  (1x1 = 1 A .. 4x4 = 16 A), voltage by material tier; overload trips a breaker at the source, manual reset; no loss
  inside a plot, poles + transformers between plots; pipes carry one fluid, limited by material temperature and
  diameter; conveyors stay cells (superseded 2026-11-23 by #1735: thin belts). All 9 design recommendations accepted.
- Anvil quality (#1613, sub of #1219): quality = tier-gated materials + smithing skill + a voxel-forging minigame
  (2D 12x8 height grid); 75-130% with grades Crude / Common / Fine (>=110) / Masterwork (>=125); machines make composed
  items at a fixed 100%; re-forging at a better anvil raises quality (costs head material, resets wear). All 8 design
  recommendations accepted.
- Basic circuits stay hand-craftable at tier 4, assembler optional (#1531); the planner hand-plans them from tier 4.
- Volcanoes (#1591): new materials (Ash near the vent, an emissive Lava crust on the crater floor); the crust is
  walkable but hot (heat damage while standing on it); digging through it may later expose simulated lava.
- Cave species are leashed underground (give up the chase at the cave mouth); only the Lost Miner may wander out near
  cave mouths at night and returns underground by day (#1227).
- RMB placing ground from a stack in the off hand stays (#836), even with a tool in the main hand (#1577).
- LOD selection stays distance-based (detail setting scales it), no screen-space-error metric (#1603, #1607).

Rulings 2026-11-23 (user)
- Creature sim (#1650): all 5 accepted: 200 active within 300 m / 30 engaged; 2 ms mean / 4 ms p99 on the user's PC; the
  same budget on the host for 4 players far apart, far creatures sent at 2 Hz; far creatures (past ~60 m) may move
  less accurately; slices 2-3 (tiers + threads) before slice 5 (new AI/combat patterns).
- Creature lifecycle (#1740): cave creatures spawn only with a player underground, 24-64 m away, out of sight, in the
  dark (light grid near zero; before #1739: no torch within ~12 m and no sky); beyond ~96 m or an unloaded chunk they are
  removed and counted back into their home cave, which refills slowly (a few minutes each) up to its cap; only tamed,
  named, guardian and boss creatures persist when far. Light-shy creatures flee lit areas, light-drawn ones come.
- Point-light shadows (#1576, 2026-11-22): cube shadow maps; the held light + the 3 nearest torches/lamps/campfires, and
  **creatures and players cast shadows too** (changed from the recommendation). Anti-aliasing default **MSAA 2x** (#1617);
  the 2.75 detail level is named **Balanced** (#1618). #1739 narrows the cube shadows to the held light.
- Light grid (#1739): static lights through a 1 m RGB-coloured grid, smoothed when sampled, smooth falloff (no bands);
  range candles ~6 m, torches ~12, lamps ~20, big lights ~30 (GUESS, tune), light bends around corners losing strength.
  All extras in: crops and saplings need light (sky light counts), creatures react to light, a light-level readout at
  the aimed spot, sky light flows a short way into caves. Held/dynamic lights stay point lights.
- Terrain cost (#1610): items 3 (shadows) and 4 (far terrain shader) both done behind switches with before/after captures
  for the user to judge; **item 3 then dropped** (far cascade caching saves at most 0.3 ms even body-local; avoid the
  every-other-frame cascade that flickered, #1163). Item 2 (upload byte budget) measured: does not pay off, no PR.
- Plan 11 rendering track (#1703): shadows, particles, creature render/cull/light, wind and network interest move to
  body-local space (one per-body transform at draw time); Lighting does it, Water owns the Plan 11 core.
- Armour (#166): three-part armour and per-material resistance profiles, **cap ~95% (not 60%)**, later (space-age) tech
  stacks mitigation further; early armour must matter: **~50-60% for a full hide set or equivalent**. Load classes cost
  stamina only, **no speed penalty** (medium/heavy cost more). Gems as recommended (armour gems give the matching
  resistance; gem weapons split 80/20); weakness discovery after the first hit; order: armour + stats, load classes,
  gems, UI, tech modules.
- Progression walkthrough (#1694): all 3 accepted (~7 h to the first EBF; time costs measured from one 20-min recording
  of T1-T2 play; report-only).
- Logistics (#1695): all 5 accepted: tubes on the conduit sub-grid; storage without type limits; one logistics step per
  tier; cross-plot storage later at EV; the slowest belt sets the budget.
- Power ladder (#1696): solar at LV, ~4x fuel values per tier, pipe-fed steam turbines at MV, the loss numbers as
  proposed. **Added: fission; liquid and gas fuels in generators/turbines; physical power plants** (heat pipes -> heat
  exchangers -> boiling -> steam -> turbines through the rotation system -> electricity; several heat/steam sources =
  several automation chains); **late-game hydro** (a large-scale dam). Closed-loop condensers with an optional open loop;
  heat pipes short runs, coolant loops for distance; uranium a deep ore (~1.5 km, #1465); dams on infinite rivers capped
  by flow and head. **Fission (changed from the recommendation): meltdowns, radiation and explosions are allowed** for
  poorly managed or high-output reactors, **and safe designs that cannot melt down must exist** (~40% output). Also: 12 m
  crater cap; other players' builds safe by default (server setting); far reactors SCRAM unless steady; contamination
  always decays (debris ~1 week, corium ~1 month); river-carried contamination later; creatures avoid and take damage,
  no mutations for now. (Turbine -> shaft -> alternator was taken as accepted; the liquid-combustion generator stays the
  simple MV option.)
- Perf budgets (#1690): after the report-only week, a miss of up to ~10-20% over target only flags; more blocks the merge
  (a PR may raise a budget in the same diff with a reason).
- POI pointers (#1692) and versioned saves (#1693; saves before 13 are dropped): all accepted as recommended.
- Multiblock sharing (#1429): GTNH-style wall sharing, everything shareable (casings, frames and hatches may serve more
  than one structure; the validator must not assume one block = one structure; breaking a shared block unforms every
  structure using it). Shared energy hatch shortfall: everyone slows evenly (proportional); shared output hatch: reserve
  space at recipe start; shared input hatch: both pull from one buffer, round-robin per recipe start.
- Anvil interaction (#1613): RMB with an item in hand puts it on the anvil, RMB with an empty hand takes it off (the
  workpiece sits in the world); using the anvil needs a crafted hammer and **every stroke costs hammer durability**
  (Quick forge charges a fixed number of strokes' wear).
- Thin belts (#1620, #1735): one lane per belt; a belt is **half a cell wide, 1/8 cell tall with 1/8 open space above**
  (a 1/4-cell placement box); **players never ride them**, not even later. Deterministic client sim + checksums; old
  conveyors auto-convert. **This supersedes "conveyors stay full cells"** (#1620 ruling, 2026-11-22).
- Rotation machines (#1734): all 4 accepted: rotary pump, trip hammer, rotary sieve, mixer now; drill/elevator/fan later;
  the trip hammer does bloom -> ingot but never heads; the rotary pump follows the electric pump's draining rules; mixer
  alloy dust at equal yield.
- Shape-mode digging (#1722): aligns to a nearby plot, else a fixed body-local grid (deterministic; digging never creates
  a plot); **Alt+Scroll** cycles Free -> 1x1 -> 2x2 (larger with better tools); **first-person free look moves to the
  middle mouse button; lock-on is unbound by default** (stays rebindable; today's MMB lock-on #1216 moves).
- Tool assembly window (#1723): a pop-out beside crafting; unusable materials greyed with what is missing; forged heads
  come from the anvil (the window never forges); old per-material recipes removed; mods slot hidden until mods exist; a
  'same again' button. **No taking apart in the window**: only the disassembler, into raw materials.
- Map zoom (#1725): all 6 accepted: a rotatable globe; unexplored = dark sphere with silhouette and day/night; telescope
  'sketching'; bodies discovered by telescope or orbit; team maps shared; the Map item still required, a star chart for
  the system level.
- Static models (#1719): the user wants CC0 packs (e.g. opengameart.org) as placeholders; plan accepted ("your suggestions
  work"): .glb in assets/models/<pack>/, THIRDPARTY entry per asset.
- Sessions (CLAUDE.md, PR #1699): audio moved to UI/UX; magic split by feature.

Deferred by the user: moon mining launch (#116), motion sickness (#117; comfort options are #168). The history of
the earlier playtest batches is in git history and the closed issues.

## Rulings 2026-11-24 (user)

- Lighting (user, 2026-11-24): sky light seeded through cave mouths (#1739 PR 3) is deferred to #1765. Rendering already
  darkens caves; crops use block light plus their own sun check.
- Combat (user, 2026-11-24, #1776): parry rework. An RMB tap parries and holding blocks. A parry succeeds when its window
  overlaps the attack's hit window. A parried projectile or melee hit does no damage, and a parried melee hit staggers.
  The parry window shows as a draining bar in a reusable status-effect widget at the left edge (#1777).
- Process (user, 2026-11-24): the coordinator sets board priorities itself, without "GUESS" marks; the user reassigns when
  they want.
- Merges (user, 2026-11-24, #1744): sessions label finished PRs `ready` (plus `runner-evidence` when the proof is a runner
  run), and the deweyserver merge bot merges them. The coordinator keeps m84 -> develop -> main and ladder -> m84.
- Pace (user, 2026-11-23): until a weekly budget reset, every session still finishes a few items a day, infrastructure first.
- Sandbox (user, 2026-11-24, #1744): slot 2 has no load gate (the load is our own runners while the user develops). Sandbox
  jobs get CPU priority over the runners.
- Runner PR runs (user, 2026-11-24, #1768, "implement all 4 now"): PR runs skip the parity check, perf budgets and pose
  twins; those run nightly or on m84.

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
