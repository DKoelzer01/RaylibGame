# Multiplayer architecture (ADR draft)

Status: **draft for review** (2026-11-14, #309 under #303). The decided rows come from the user's answers on #303 (2026-11-13).
Everything marked *Proposed* is a recommendation, not yet agreed.

## Decided

| Topic | Decision |
|---|---|
| Authority | Authoritative server: clients send inputs and requests, the server simulates and validates |
| Topology | Listen server first (a player hosts in-game); the same server code builds headless as a dedicated server later |
| Transport | Steam networking sockets: relay (`ConnectP2P`), lobbies, invites, rich presence, **and** connect by IP (`ConnectByIPAddress`) |
| Scale | Medium servers, 8-32 players; leave room for more later (#350) |
| Sharing | Shared world, bases and tier/quest progression; per-player inventory and survival; PvP and land claims |
| World data | The erosion bake (Plan 7) runs once per world on the server, is stored in the save, sent to joining clients and cached by hash: the world fingerprint (#316) |

## Proposed answers to the open choices

**1. Tick rates.** *Proposed:* a fixed 30 Hz server simulation, and snapshots at 20 Hz, with the rate falling off for far entities (#340).
- Today the player steps physics at 120 Hz substeps inside a variable frame. Keep those substeps inside a 30 Hz server tick, so collision behaves the same.
- The factory and fluids already run on their own clocks (factory ticks; fluid every 25 ms). They become server-side jobs on the same sim clock.

**2. Replication style per system.** *Proposed:*

| System | Style | Why |
|---|---|---|
| Players, animals, vehicles, rockets | Snapshot + delta, interpolated on clients | Continuous motion |
| Terrain edits, plots and parts, blueprints | Reliable ordered events (the existing edit and part command structs) | Discrete, must not be lost; they replay on join |
| Inventories, containers, machines | Server state plus reliable events; clients see machine *views* (#332) | Small; only the owner or viewer needs them |
| Fluids | Active bricks streamed near players; static water never sent (#333) | Generated water is a function of the seed and bake |
| Trees, plants, farming | Events (chopped, planted, grown stage) over the seed-generated base (#335) | The base is deterministic |
| Time, weather | A server clock plus a seed: clients compute the sky (#339) | Nearly free |

**3. Frames.** *Proposed (strongly recommended):*
- Every replicated position is **body-local, in doubles, relative to a body's network id**. Never render space; never world space across the spin.
- The code already moves this way: `toLocalD`, ground frames, and fly-follow #127.
- The render origin stays a client-only concern.

**4. Client prediction.** *Proposed:*
- **Predicted:** own movement (with server reconciliation, #327), the dig and placement preview, and the swing timer.
- **Not predicted:** dig results and yields, crafting, machine state, other players.
- A misprediction snaps the preview, not the world.

**5. Vehicles.** *Proposed:*
- The server simulates rockets and ships.
- The pilot's inputs are predicted locally for camera feel only.
- One global rocket flight state exists today (`rocket.cpp`); it becomes one per vehicle with a network id (#317).

**6. Time warp.** *Proposed:* disabled in multiplayer at first. Later: allowed only when every player is in a vehicle or asleep (a server-wide vote).

**7. Client-only systems.** *Proposed:* these never cross the wire:
- particles, ambience, grass, foliage instancing and clouds;
- LOD, streaming and the preload ring;
- the viewmodel, HUD and screenshots.

## Code survey (2026-11-14, `m49-erosion`)

What stands between today's code and a server. Counts are from grep, a rough guide.

- **Single-player globals** (`extern` in headers):
  - `player`, `inventory`, `survival`, `ship`, `buildTool`, `terrainTool`, `telescope`, `settings`, `worldGen`, `worldDifficulty`, `lodSettings`;
  - the hotbar, suit, quest book, tool wear and navigation state (file statics).
  - Per-player ones move into a `PlayerSession`; world ones into a `World` the server owns.
- **Wall clock in simulation:**
  - `GetTime()`: 273 calls in `main.cpp` (mostly benches), 15 in `animals.cpp`, 8 each in `survival.cpp` and `navigation.cpp`, 7 in `world.cpp`, 6 in `weather.cpp`, 5 each in `terraintool.cpp` and `scene.cpp`.
  - Gameplay ones move to a sim clock (`SimTime()`); benches keep wall time.
- **Input read inside logic:**
  - `IsKey*`/`IsMouse*` calls: `rocket.cpp` 22, `camera.cpp` 14, `ship.cpp` 12, `buildtool.cpp` 12, `player.cpp` 11, `terraintool.cpp` 5, `animals.cpp` 4.
  - These become an input struct per player per tick (the player update already takes `wish`, `jetUp`, `jetDown` that way).
- **Render origin in physics:**
  - `renderOrigin` appears in `player.cpp` 15 times, `world.cpp` 17, `rocket.cpp` 12.
  - Physics must run in body-local doubles; the render origin becomes a draw-time offset only.
- **Drawing inside simulation files:**
  - `animals.cpp` (26) and `survival.cpp` (21) mix HUD and drawing with logic; `fluid.cpp` uploads meshes (4).
  - Split into sim and view halves, so the headless server links no GL.
- **Already in good shape:**
  - The world is generated from the seed (terrain, rivers, trees, ore, herds, weather).
  - Edits, parts, factory messages and fluid commands are small command structs.
  - The factory simulation is deterministic on its own thread.
  - Collision reads the density field, so the server needs no meshes.
  - The bench harness can drive scripted clients (#347).

## Toolchain spike (to do first)

The project builds with MinGW. The Steamworks SDK ships MSVC import libraries, but also a **flat C API** (`steam_api_flat.h`) over `steam_api64.dll`. GNU ld can link a DLL directly, so MinGW should work without switching compilers. **To verify:**
1. Link `steam_api64.dll` from MinGW.
2. `SteamAPI_Init` with `steam_appid.txt` (480, Spacewar, for testing).
3. Open an `ISteamNetworkingSockets` listen socket.
4. Connect by IP between two local instances.

If the flat API proves awkward, Valve's open-source GameNetworkingSockets (same API, BSD) is the fallback for direct connect. It has heavy dependencies (protobuf, OpenSSL).

## Proposed first steps (phase 1, single-player keeps working at each one)

1. **Sim clock.** `SimTime()` / `SimDt()` replace `GetTime()` in gameplay code; benches untouched. Small, mechanical, no behaviour change.
2. **Input struct.** `PlayerInput { wish, look, jump, jetUp, jetDown, use, swing, hotbar... }` filled from raylib once a frame. Player, terrain tool and build tool read only the struct.
3. **PlayerSession.** Move `player`, `inventory`, `survival`, the suit, the hotbar and tool wear into one struct. Single-player has one; the server has N.
4. **Network ids** (#317) for bodies, plots, parts, animals and vehicles: stable across save and load.
5. **Headless build target** (`game_server`, no window). It runs the world, factory, fluids, animals and weather on the sim clock, for profiling (#348) before any networking exists.
6. **Toolchain spike** above, in parallel with steps 1 and 2.

Each step gets its own issue and PR, and is verified by the existing benches (single-player stays green).
