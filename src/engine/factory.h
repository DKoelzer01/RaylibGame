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

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

class Scene;
class Plot;
enum class PartType : uint8_t;

constexpr int kTicksPerSecond = 20;

enum class MachineStatus : uint8_t { Idle, Working, NoPower, Overvoltage, OutputFull, NoFuel };
const char* MachineStatusName(MachineStatus s);

struct MachineView {             // published per machine each tick
    MachineStatus status = MachineStatus::Idle;
    int progress = 0, duration = 0; // recipe ticks done / needed
    int energy = 0, energyMax = 0;  // buffer (generators: fuel ticks left)
    std::string input, output;      // "iron_ore x3"
    bool valid = false;
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
    FactoryStats Stats();
    // Player item transfers (queued for the next tick).
    void Insert(uint64_t bodyId, uint32_t plotId, int part, const std::string& item, int count);
    void Extract(uint64_t bodyId, uint32_t plotId, int part);
    // Tests: put items straight into a machine (fuel, ore) before it runs.
    void Seed(uint64_t bodyId, uint32_t plotId, int part, const std::string& item, int count);
    // Whether a machine type takes an item (recipe input, fuel, or a chest).
    bool Accepts(PartType type, const std::string& item);
    void WaitTicks(int n);                           // block until n more ticks have run
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
}

#endif
