#ifndef FACTORY_H
#define FACTORY_H

// Factory simulation (milestone 4, design doc section 5).
//
// A deterministic simulation on its own thread at a fixed 20 ticks a second, separate from the
// renderer: integer item counts and energy, fixed iteration order, no floating point in the logic.
// The main thread sends each plot's layout when it changes (machines, conveyors, cables) and the
// player's position now and then; the simulation keeps machine state across layout changes, ticks,
// and publishes a snapshot the HUD reads. Item transfers by the player go through commands.
//
// Machines: chest (storage), combustion generator (burns wood into EU), electric furnace and
// macerator (recipes from resources/data/machine_recipes.txt), battery (stores EU).
// Networks: conveyor cells joined face to face form an item network that moves output items to
// machines that take them (a per-network throughput budget); cable cells form a power network:
// generators and batteries supply, consumers share any shortfall equally, and a machine whose tier
// is below the network's voltage refuses power (shows "overvoltage") rather than exploding.
// Plots far from the player tick once a second with a 20x step: the same rates and budgets, a
// twentieth of the cost, never more output than the full simulation.

#include "voxel.h" // Int3
#include <cstdint>
#include <cstdio>
#include <functional>
#include <tuple>
#include <raylib.h>
#include <string>
#include <vector>

class Scene;
class Plot;
enum class PartType : uint8_t;

constexpr int kTicksPerSecond = 20;

enum class MachineStatus : uint8_t { Idle, Working, NoPower, Overvoltage, OutputFull, NoFuel, NoLoad, NoWater };
const char* MachineStatusName(MachineStatus s);

// A campfire's fire (#1103): by itself it lights when it has fuel and something to cook; a player can light it with nothing
// to cook (it burns its fuel: light and warmth) or put it out (no burn, no cooking) until it is lit again or given
// something to cook. Lit, it goes back to by itself when its fuel runs out.
enum FireState : uint8_t { kFireAuto = 0, kFireLit = 1, kFireOut = 2 };

struct MachineView {             // published per machine each tick
    MachineStatus status = MachineStatus::Idle;
    int progress = 0, duration = 0; // recipe ticks done / needed
    int energy = 0, energyMax = 0;  // buffer (generators: fuel ticks left)
    std::string input, output;      // "iron_ore x3"
    std::vector<std::pair<std::string, int>> inSlots, outSlots; // per slot ("", 0 when empty); chests use inSlots
    std::vector<int> inCenti;                                   // (#1102) per input slot, a fraction more (hundredths): a chest's 0.4 dirt
    bool hasFuelSlot = false;       // fuel-burning machines and generators (#194)
    std::pair<std::string, int> fuelSlot;
    int water = 0, waterMax = 0;    // litres in its tank (boilers, #229)
    int heat = 0, heatMax = 0;      // boilers warm up before they make full steam
    bool hotOutputs = false;        // its outputs come out hot (a furnace, kiln...: their icons glow, #580)
    bool canLight = false;          // a campfire: the screen offers to light it or put it out (#1103)
    uint8_t fire = kFireAuto;
    bool rotation = false;          // (#1269) a rotation machine: its shaft's speed and its network's load
    float rpm = 0.0f;
    int speedPercent = 100, stressPercent = 0;
    const char* rotStatus = "";     // "Running", "Slowed", "Too slow", "Stopped", "Jammed"
    bool valid = false;
};

// (#1269) A rotation part (shaft, cog, crank, gauge, rotation machine) on a plot near the player, as of the last tick.
struct RotationView {
    bool valid = false;
    float rpm = 0.0f;               // signed, about the plot axis it lies on
    const char* status = "";        // Rotation::StatusName
    const char* jam = "";           // why its network is jammed ("" if not)
    int stressPercent = 0;          // demand / capacity
    int speedPercent = 100;         // the network's stress ratio: every part turns at this share of its nominal speed
    int capacity = 0, demand = 0;   // SU
};

// A cable's power network, for its look-at panel (#570): what is on it and what flowed in the last step.
struct CableView {
    bool valid = false;             // the cable is known to the simulation (false for a tick after asking)
    int voltage = 0;                // the network's (its highest generator or battery); 0: nothing feeds it
    int generators = 0, consumers = 0, batteries = 0;
    float supply = 0.0f, demand = 0.0f; // EU per tick made (with batteries' share) and wanted
};

struct FactoryStats {
    int machines = 0, activeMachines = 0, itemNetworks = 0, powerNetworks = 0;
    double avgTickMs = 0.0, maxTickMs = 0.0;
    long long ticks = 0;
    long long itemsMade = 0;
    // (#1207) the last minute of ticks (up to 1200), refreshed once a second: whole ticks and their phases (commands in,
    // the plots, the view block)
    int recentTicks = 0;
    double meanTickMs = 0.0, p99TickMs = 0.0, worstRecentMs = 0.0, cmdMs = 0.0, tickOnlyMs = 0.0, viewMs = 0.0;
};

// (#1199) Production statistics (the Statistics screen, and the oracle the summary sim is checked against). Per item: made and
// used by machines (crafts, fuel burnt) and put in and taken out by players, per minute over a window and in total; per
// machine kind: how many there are, their status now, ticks spent in each status, EU made and used per tick.
enum class StatWindow : uint8_t { Minute, Hour, TenHours }; // the last minute (by seconds), hour (by minutes), 10 hours
struct StatLine {
    std::string item;
    bool estimated = false;                                        // (#1261) partly from far plots in the summary sim
    double madePerMin = 0, usedPerMin = 0, handInPerMin = 0, handOutPerMin = 0;
    long long made = 0, used = 0, handIn = 0, handOut = 0;      // totals since the world loaded (finished seconds)
    std::vector<float> madeSpark, usedSpark;                      // the window's buckets, oldest first (per bucket)
};
struct StatKind {
    int type = 0;                                                  // PartType
    int count = 0;
    int now[8] = {};                                               // machines per MachineStatus now
    long long ticks[8] = {};                                       // ticks per MachineStatus since the load
    double euMadePerTick = 0, euUsedPerTick = 0;                   // over the window
};
struct FactoryStatsView {
    StatWindow window = StatWindow::Minute;                        // the one shown (a longer one not filled yet falls back)
    double coveredMinutes = 0;                                     // how much of it has data
    int plots = 0, estimatedPlots = 0;                             // (#1261) plots in the summary sim (their figures are estimates)
    std::vector<StatLine> items;                                   // by item name
    std::vector<StatKind> kinds;                                   // by PartType
    double euMadePerTick = 0, euUsedPerTick = 0, euDemandPerTick = 0;
    long long batteryEnergy = 0, batteryCapacity = 0;              // now
};

namespace Factory {
    void Start();                                    // spawns the simulation thread (idempotent)
    void Stop();
    void Reset();                                    // a new world / a load: forget every plot
    // Main thread, once a frame: send layouts of plots that changed, the player's position, and
    // take the results of item commands (into the inventory).
    void Sync(Scene& world);
    MachineView View(uint64_t bodyId, uint32_t plotId, int part);
    CableView Cable(uint64_t bodyId, uint32_t plotId, Int3 cell); // the network through that cable cell (#570); watched like View
    // Whether a machine on a plot near the player is working (any number of calls a frame: models, smoke).
    // View() follows one machine at a time (the open machine screen) - never call it for many machines.
    bool Working(uint64_t bodyId, uint32_t plotId, int part);
    FactoryStats Stats();
    void ResetTickStats();   // the tick-time ring starts empty (a timed window: mean, p99 and worst then cover only it)
    // (#1199) The production statistics over a window: every simulated plot (onePlot false) or one plot. Waits for the tick.
    FactoryStatsView StatsView(StatWindow window, bool onePlot = false, uint64_t bodyId = 0, uint32_t plotId = 0);
    // Player item transfers (queued for the next tick).
    void Insert(uint64_t bodyId, uint32_t plotId, int part, const std::string& item, int count, uint32_t to = 0, int centi = 0); // to (#332): leftovers to that player; centi (#1102): a fraction more (a chest keeps it)
    // The same, into the fuel slot (Shift+click: an item that is both fuel and a recipe input goes to the inputs by default).
    void InsertFuel(uint64_t bodyId, uint32_t plotId, int part, const std::string& item, int count, uint32_t to = 0);
    // (#1193) Research labs, between ticks: a machine's one output (as a miner's ore; "" any), and its output stacks whose
    // item names start with `prefix` taken out (the game counts them). False / empty: not simulated (yet).
    bool SetVariant(uint64_t bodyId, uint32_t plotId, int part, const std::string& item);
    std::vector<std::pair<std::string, int>> DrainOutputs(uint64_t bodyId, uint32_t plotId, int part, const std::string& prefix);
    // (#1103) Light a campfire (lit = true) or put it out; other machines ignore it.
    void SetFire(uint64_t bodyId, uint32_t plotId, int part, bool lit);
    // (#2072) A sealed ruin cache with a power opener (cache_openers.txt) is a consumer on its plot: a cable or a battery bank
    // against it charges its lock, at most its cap a tick and only from a network at its voltage or more. CacheCharged: how
    // full (0..1; 0 until the sim has it), from the last tick, without waiting for the sim. CacheSpend empties it (opened).
    // Factory::Start hands both to CacheOpeners_SetPower; InstallCachePower does it again (after a test's own hooks).
    float CacheCharged(uint64_t bodyId, uint32_t plotId, int part);
    void CacheSpend(uint64_t bodyId, uint32_t plotId, int part);
    void InstallCachePower();
    // (#1269) A player turns a hand crank: it turns for 0.3 s after the last call (call it each frame RMB is held).
    void Crank(uint64_t bodyId, uint32_t plotId, int part);
    // (#1270) A water wheel's or windmill's speed from its river or wind (signed RPM about its plot axis): rotationsources.*
    void SetSourceRpm(uint64_t bodyId, uint32_t plotId, int part, int rpm);
    // (#1269) A rotation part's speed and its network's load (plots near the player; invalid otherwise).
    RotationView RotationAt(uint64_t bodyId, uint32_t plotId, int part);
    // (#1269) Every rotation part of a near plot: (part, signed RPM about its plot axis). Empty if none or far.
    void RotationSpeeds(uint64_t bodyId, uint32_t plotId, std::vector<std::pair<int, float>>& out);
    // Water piped into a machine's tank (boilers), in litres; what doesn't fit is lost.
    void AddWater(uint64_t bodyId, uint32_t plotId, int part, int litres);
    int TankRoom(uint64_t bodyId, uint32_t plotId, int part); // litres a boiler can still take (less water already sent); 0 if unknown
    bool TakesWater(PartType type);
    // Every machine recipe (machine_recipes.txt) by item id, for the in-game database (#354). ticks/eu: 20 ticks a
    // second; generators list their fuels as one-input recipes with eu = EU per tick produced.
    struct RecipeInfo { PartType machine; std::vector<std::pair<std::string, int>> in, out; int ticks = 0, eu = 0; };
    std::vector<RecipeInfo> RecipeList();
    std::vector<std::string> FuelsOf(PartType type); // item ids a machine burns (machines.txt)
    void Extract(uint64_t bodyId, uint32_t plotId, int part, uint32_t to = 0); // to (#332): the outputs to that player (0: this game)
    // Container screen: one slot's whole stack back to the inventory (next tick).
    // which: 0 inputs, 1 outputs, 2 the fuel slot.
    void TakeSlot(uint64_t bodyId, uint32_t plotId, int part, int slot, int which, uint32_t to = 0);
    // The same onto the cursor (#192: all inventories behave alike): `count` of the slot's stack (0 = all), delivered by
    // Sync to the cursor sink (the inventory screen's held stack) instead of the inventory.
    void TakeSlotToCursor(uint64_t bodyId, uint32_t plotId, int part, int slot, int which, int count, uint32_t to = 0);
    // (#862) `count` of item into one storage slot of a chest (empty, or the same item up to its cap), as a held stack is put
    // down in Minecraft; what doesn't fit goes back as Insert's leftovers do. Not a chest: as Insert (machines route inputs).
    void InsertSlot(uint64_t bodyId, uint32_t plotId, int part, const std::string& item, int count, int slot, uint32_t to = 0, int centi = 0); // centi (#1102): a fraction more
    void SetCursorSink(void (*sink)(const std::string& item, int count, float heat)); // heat: 1 fresh out of a heating machine (#580)
    void (*CursorSink())(const std::string& item, int count, float heat); // (#854) the one set (a joined client puts a server take on the cursor)
    // Water pumps (waterworks): water units in a pump's outputs (as of the last tick), and removing
    // them as they go down its pipes.
    int PumpWater(uint64_t bodyId, uint32_t plotId, int part);
    void Drain(uint64_t bodyId, uint32_t plotId, int part, const std::string& item, int count);
    // Crafting and storing from nearby chests: take `count` of one item out of a chest (into the
    // inventory next tick), and what the chests of plots near the player hold (as of the last tick).
    void Take(uint64_t bodyId, uint32_t plotId, int part, const std::string& item, int count, uint32_t to = 0);
    std::vector<std::pair<std::string, int>> ChestContents(uint64_t bodyId, uint32_t plotId, int part);
    // Tests: put items straight into a machine (fuel, ore) before it runs.
    void Seed(uint64_t bodyId, uint32_t plotId, int part, const std::string& item, int count);
    // Whether a machine type takes an item (recipe input, fuel, or a chest).
    bool Accepts(PartType type, const std::string& item);
    void Advance(float simDt);                       // one simulation step (#315): grants the ticks it owes (20 Hz x speed)
    void WaitTicks(int n);                           // run n more ticks now and block until they have
    void SetSpeed(int multiplier);                   // tests: run the 20 Hz clock this many times faster
    // Saves: machine contents per body index (the plots themselves are saved with the body).
    void Save(FILE* f, uint64_t bodyId);
    void Load(FILE* f, uint64_t bodyId);
    // The same in two steps, for bodies that exist only later (ring rocks loaded dormant, a rocket in
    // flight): read a "machines" block as text (empty if none), then attach it to a body id.
    std::vector<std::string> ReadSaved(FILE* f);
    void LoadSaved(uint64_t bodyId, const std::vector<std::string>& lines);
    // A plot moved to another body / id (a rocket launching or landing): its machines go with it.
    void MovePlot(uint64_t fromBody, uint32_t fromPlot, uint64_t toBody, uint32_t toPlot);

    // ---- Multiplayer (#332): the server simulates, clients see views ----
    // Server: more machines to publish a full view of each tick (the ones its players have open), as View does for one.
    using MachineKey = std::tuple<uint64_t, uint32_t, int>; // body id, plot id, part
    void WatchMany(const std::vector<MachineKey>& keys);
    MachineView Watched(const MachineKey& key);       // the last tick's view of a WatchMany machine (invalid until then)
    // Server: more places whose plots count as near (ticked every tick, their working machines published): its players.
    void SetNearPoints(const std::vector<Vector3>& renderPoints);
    void SetCameraNear(bool counts);                 // (#1302) false on a dedicated server: only its players make plots near
                                                    // (with nobody about, every plot runs in the summary sim)
    // The machines working on one plot now (near plots only): what clients get for animations (MachineStates).
    std::vector<int> WorkingParts(uint64_t bodyId, uint32_t plotId);
    // (#1243) The machines on a near plot stopped for a reason now (no power, overvoltage, output full, no fuel, no water;
    // not plain idle): machine sounds' stall cues. Local simulation only (empty on a joined client).
    std::vector<int> StalledParts(uint64_t bodyId, uint32_t plotId);
    // Client: the server's views stand in for the simulation (which then runs no ticks): View / Working read these.
    void SetRemote(bool remote);
    bool Remote();
    void SetRemoteView(const MachineKey& key, const MachineView& view);
    void SetRemoteWorking(uint64_t bodyId, uint32_t plotId, const std::vector<int>& parts); // that plot's working machines

    // ---- (#1254) Fluid in machine slots: alchemy's essences ----
    // Items whose names start with a registered prefix ("essence_") are fluid held in a machine's slots, not things: no
    // machine or chest accepts them from a conveyor or a hand, hands never take them out, and a removed machine's are lost
    // with it. Only MoveFluids moves them (pipes, phials); recipes make and use them as usual. Register before Start.
    void AddFluidPrefix(const std::string& prefix);
    bool IsFluidItem(const std::string& item);
    // Between ticks: each machine's fluid stacks (by item, summed) in its inputs and outputs, and its empty slots.
    struct FluidSlots { bool valid = false; std::vector<std::pair<std::string, int>> in, out; int freeIn = 0, freeOut = 0; };
    std::vector<FluidSlots> ReadFluids(const std::vector<MachineKey>& keys);
    // Between ticks, in order: up to `count` of a fluid item out of one machine's slots (which: 0 inputs, 1 outputs) into
    // another's. Either end may be absent (hasFrom / hasTo false: a phial poured in, a spill). Moves at most what the source
    // holds and what fits the target, keeping the target's total of the item at most `cap` (0: no cap). Returns each move's
    // amount.
    struct FluidMove { MachineKey from{}, to{}; bool hasFrom = false, hasTo = false; int fromWhich = 1, toWhich = 0; std::string item; int count = 0, cap = 0; };
    std::vector<int> MoveFluids(const std::vector<FluidMove>& moves);
}
// (#1254) More machine recipes, made by the game from its own data (alchemy: the distillery from aspects.txt, the cauldron
// from alchemy.txt): called once when the factory first loads its tables, after machine_recipes.txt. Set it before then.
inline std::function<void(std::vector<Factory::RecipeInfo>& out)> gFactoryExtraRecipes;
// Server: results for a player (Extract / TakeSlot / Take / leftovers with `to`) go here instead of this game's inventory.
inline std::function<void(uint32_t to, const std::string& item, double amount, float heat)> gFactoryDeliver; // (#1102: an amount: a chest's fraction too)
// Client (#332): the machine the screen or look-at panel views changed (Factory::View with SetRemote): ask the server for it.
inline std::function<void(uint64_t bodyId, uint32_t plotId, int part)> gFactoryRemoteWatch;

#endif
