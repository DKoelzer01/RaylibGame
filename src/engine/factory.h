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

struct MachineView {             // published per machine each tick
    MachineStatus status = MachineStatus::Idle;
    int progress = 0, duration = 0; // recipe ticks done / needed
    int energy = 0, energyMax = 0;  // buffer (generators: fuel ticks left)
    std::string input, output;      // "iron_ore x3"
    std::vector<std::pair<std::string, int>> inSlots, outSlots; // per slot ("", 0 when empty); chests use inSlots
    bool hasFuelSlot = false;       // fuel-burning machines and generators (#194)
    std::pair<std::string, int> fuelSlot;
    int water = 0, waterMax = 0;    // litres in its tank (boilers, #229)
    int heat = 0, heatMax = 0;      // boilers warm up before they make full steam
    bool hotOutputs = false;        // its outputs come out hot (a furnace, kiln...: their icons glow, #580)
    bool valid = false;
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
    // Player item transfers (queued for the next tick).
    void Insert(uint64_t bodyId, uint32_t plotId, int part, const std::string& item, int count, uint32_t to = 0); // to (#332): leftovers to that player
    // The same, into the fuel slot (Shift+click: an item that is both fuel and a recipe input goes to the inputs by default).
    void InsertFuel(uint64_t bodyId, uint32_t plotId, int part, const std::string& item, int count, uint32_t to = 0);
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
    void TakeSlotToCursor(uint64_t bodyId, uint32_t plotId, int part, int slot, int which, int count);
    void SetCursorSink(void (*sink)(const std::string& item, int count, float heat)); // heat: 1 fresh out of a heating machine (#580)
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
    // The machines working on one plot now (near plots only): what clients get for animations (MachineStates).
    std::vector<int> WorkingParts(uint64_t bodyId, uint32_t plotId);
    // Client: the server's views stand in for the simulation (which then runs no ticks): View / Working read these.
    void SetRemote(bool remote);
    bool Remote();
    void SetRemoteView(const MachineKey& key, const MachineView& view);
    void SetRemoteWorking(uint64_t bodyId, uint32_t plotId, const std::vector<int>& parts); // that plot's working machines
}
// Server: results for a player (Extract / TakeSlot / Take / leftovers with `to`) go here instead of this game's inventory.
inline std::function<void(uint32_t to, const std::string& item, int count, float heat)> gFactoryDeliver;
// Client (#332): the machine the screen or look-at panel views changed (Factory::View with SetRemote): ask the server for it.
inline std::function<void(uint64_t bodyId, uint32_t plotId, int part)> gFactoryRemoteWatch;

#endif
