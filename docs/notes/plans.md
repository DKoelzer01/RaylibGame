# Plans: links and history

- Plan 2 (early game, done): https://claude.ai/artifact/RL22kZ7VBMugAFgGRPS3cE
- Plan 3 (Terra tiers, space): https://claude.ai/artifact/XVw45AC6vwYitbdhaeP2xL - M21-M24 done (M24 on PR #28).
- Plan 4 (flowing water): https://claude.ai/artifact/27xL5j1zkirSe8pDZgPz59
- **Current plans and design docs (2026-11-21)** - read your area's doc before starting:
  - Plan 11, body-relative physics frames (#667, **P0**, Water): https://claude.ai/artifact/2nJjsd4cSHYP6rj4UAnxSK.
    Behind `RAYLIBGAME_FRAMES=1` while it lands (F3+); code reading positions near a body uses `toLocalOf`. F5 (vehicles
    and stations), F7 (save format) next; Space's star/orbit, stations and lander work waits on them.
  - Factory and Progression (the ladder #1271, research #1030, factory sim #997): https://claude.ai/artifact/F68hF9Yd8aXM4J9KLeGPCn.
  - Plan 13, Space (#1262; tiers 8-15): https://claude.ai/artifact/NaNNSfcfJs6r1an5SFU7m7.
  - Combat (#1215-#1222; decisions in its section 14): https://claude.ai/artifact/EH4w7X5sJ6YS6Hm9vMJ3XJ.
  - Plan 12, creatures (#1223-#1228): https://claude.ai/artifact/BAuGD4A613nu9LZyuDvjPG.
  - Magic (#1001, issues #1248-#1259; decisions in section 14; `RAYLIBGAME_MAGIC=1`): https://claude.ai/artifact/QiR81H8phuCeuERszbm5Nf.
  - Audio (#1003, #1238-#1246, Steam Audio #1020; v2 behind `RAYLIBGAME_AUDIO_V2=1`, fetch the SDK with
    tools/fetch_steamaudio.sh; voice chat is OUT of scope): https://claude.ai/artifact/2DUKL3mr2JMLTfZRFN3ccA.
  - Overhangs and arches (#999, #1231-#1237): https://claude.ai/artifact/8nKfS1yZy41f2i2HXsqthH.
  - Caves, Plan 10 (#615): https://claude.ai/artifact/C3DCqMnRnuc2Nni1BeeZjA. Plan 9 (formations, #586, done):
    https://claude.ai/artifact/DtYGaq5ruSCQ2hcfEBKJ2B. Accessibility (#1017): https://claude.ai/artifact/Pkxw3MRbMvZA1PdA4p3yeQ.
  - Trackers: #158 features, #128 graphics, #303 multiplayer, #719 dungeons, #1429 multiblocks (later).
- **Plans and designs of 2026-11-23** (rulings in decisions.md; design text lives in the issue):
  - Creature sim architecture (#1650); creature lifecycle (#1740); armour and weapons (#166).
  - Progression walkthrough (#1694); logistics beyond belts (#1695); power ladder + fission (#1696); anvil forging (#1613);
    rotation machines (#1734); tool assembly window (#1723); multiblock wall sharing (#1429 addendum).
  - Conduits (#1620) and thin belts (#1735; **supersedes "conveyors stay full cells"** from #1620).
  - Light grid (#1739); point-light shadows (#1576, with #1617 MSAA 2x default, #1618 "Balanced" detail level).
  - Plan 11 rendering track, body-local shadows/particles/creatures/wind/network interest (#1703). Terrain cost (#1610):
    item 3 (far shadow cascade caching) **dropped**.
  - Perf budgets in CI (#1690); versioned saves (#1693); POI pointers (#1692); shape-mode digging (#1722); map zoom levels
    (#1725); static .glb model pipeline (#1719).
- Earlier plans (done or superseded): Plan 7 erosion terrain (#351, https://claude.ai/artifact/EXix8zmVSS2WRUSLAsovJg,
  default since #362, `RAYLIBGAME_EROSION=0` = off), Plan 8 (https://claude.ai/artifact/LFnRaLaCPdgEPsuX7L6dJb).
- Plan 7 details: with erosion (the default) the 3D band is off (RAYLIBGAME_3D=0/1 overrides). Gullies (ErosionGullies,
  RAYLIBGAME_GULLIES=0 A/B); glacial (RAYLIBGAME_GLACIAL); rivers on the bake grid, carve-only water (--water-stats raises
  nothing). Tests: --test-flow-sea, --test-pond, --test-pump-hand, --test-hover. Open questions for the user are on #351.
- Plan 6 (GPU terrain generation + GPU-driven drawing, approved 2026-11-12): https://claude.ai/artifact/VJdKtUZEpBg5uefnqewRyZ -
  issues #143 (steps #144-#150), under #128. Done: steps 1-4 + 7 (PRs #153, #157, #179); step 5 deferred (#148), step 6
  not needed (#149: GPU density is bit-identical, collision keeps the CPU field). One path on GL 4.3 (base-instance records).
- Plan 5 (art, life, world variety): https://claude.ai/artifact/5igt2VrsTyNhZ3BZ94teHX - A animal models/AI, B biomes,
  C 3D terrain, D stylized PBR terrain: all built and merged.
