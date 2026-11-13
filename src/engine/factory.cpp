#include "factory.h"
#include "inventory.h"
#include "plots.h"
#include "utils/logger.h"
#include "worldtypes.h"
#include "scene.h"
#include <raymath.h>

#include <raylib.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <fstream>
#include <map>
#include <mutex>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <unordered_set>

const char* MachineStatusName(MachineStatus s) {
    switch (s) {
        case MachineStatus::Idle: return "idle";
        case MachineStatus::Working: return "working";
        case MachineStatus::NoPower: return "no power";
        case MachineStatus::Overvoltage: return "overvoltage: refuses power";
        case MachineStatus::OutputFull: return "output full";
        case MachineStatus::NoFuel: return "no fuel";
        case MachineStatus::NoLoad: return "no load (fire banked)";
        case MachineStatus::NoWater: return "no water";
    }
    return "?";
}

namespace {

// ---------------------------------------------------------------------------------------------
// Data
// ---------------------------------------------------------------------------------------------

struct Stack { int item = -1, count = 0; };
constexpr int kSlots = 4;
constexpr int kStackMax = 64;
constexpr int kChestSlots = 16;

struct Recipe {
    PartType machine;
    std::vector<Stack> in, out;
    int ticks = 0, eu = 0;   // generators: ticks of burn per fuel item, EU per tick produced
};

struct Tables {
    std::vector<Recipe> recipes;
    // Item ids are interned from both threads (the main thread's inserts, the simulation's loaded
    // states) and names read on both: a new name grew the vector under a reader. Locked; names are
    // returned by value.
    int item(const std::string& n) {
        std::lock_guard<std::mutex> lock(mutex);
        auto it = itemIds.find(n);
        if (it != itemIds.end()) return it->second;
        itemIds[n] = (int)itemNames.size();
        itemNames.push_back(n);
        return (int)itemNames.size() - 1;
    }
    int find(const std::string& n) { // -1 if never seen
        std::lock_guard<std::mutex> lock(mutex);
        auto it = itemIds.find(n);
        return it == itemIds.end() ? -1 : it->second;
    }
    std::string name(int id) {
        std::lock_guard<std::mutex> lock(mutex);
        return id >= 0 && id < (int)itemNames.size() ? itemNames[id] : std::string();
    }
private:
    std::mutex mutex;
    std::vector<std::string> itemNames;
    std::unordered_map<std::string, int> itemIds;
} tables;

PartType MachineByName(const std::string& n) {
    for (int t = 0; t < (int)PartType::Count; ++t) {
        const char* m = PartMachineName((PartType)t);
        if (m && n == m) return (PartType)t;
    }
    return PartType::Count;
}

// machines.txt: how each machine is powered.
enum class Power : uint8_t { None, Fuel, Eu, Generator, Storage, Battery };
struct MachineDef {
    int volts = 0;
    Power power = Power::None;
    int euOut = 0;
    std::vector<int> fuels;  // item ids it burns
    int waterPerTick = 0;    // litres of water per tick of burning (boilers, #229); 0: no water tank
};
constexpr int kWaterTank = 4000;   // litres a boiler holds (four buckets)
constexpr int kBucketLitres = 1000;
constexpr int kHeatMax = 400;      // ticks of burning from cold to full steam (20 s)
MachineDef machineDefs[(int)PartType::Count];
std::unordered_map<int, int> fuelTicks; // item id -> burn ticks
const MachineDef& Def(PartType t) { return machineDefs[(int)t]; }

void LoadMachineDefs() {
    std::ifstream f("resources/data/machines.txt");
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("fuel ", 0) == 0) {
            std::stringstream ss(line.substr(5));
            std::string item; int ticks = 0;
            if (ss >> item >> ticks) fuelTicks[tables.item(item)] = ticks;
            continue;
        }
        std::vector<std::string> cols;
        std::stringstream ss(line);
        std::string c;
        while (std::getline(ss, c, '|')) {
            const size_t a = c.find_first_not_of(" \t"), b = c.find_last_not_of(" \t\r");
            cols.push_back(a == std::string::npos ? "" : c.substr(a, b - a + 1));
        }
        if (cols.size() < 4) continue;
        const PartType t = MachineByName(cols[0]);
        if (t == PartType::Count) continue;
        MachineDef& d = machineDefs[(int)t];
        d.volts = atoi(cols[1].c_str());
        const std::string& p = cols[2];
        d.power = p == "fuel" ? Power::Fuel : p == "eu" ? Power::Eu : p == "generator" ? Power::Generator
                : p == "storage" ? Power::Storage : p == "battery" ? Power::Battery : Power::None;
        d.euOut = atoi(cols[3].c_str());
        if (cols.size() > 4) { std::stringstream fs(cols[4]); std::string fuel; while (fs >> fuel) d.fuels.push_back(tables.item(fuel)); }
        if (cols.size() > 5) d.waterPerTick = std::max(0, atoi(cols[5].c_str()));
    }
}

void LoadRecipes() {
    if (!tables.recipes.empty()) return;
    LoadMachineDefs();
    std::ifstream f("resources/data/machine_recipes.txt");
    std::string line;
    auto parseStacks = [](const std::string& s, std::vector<Stack>& out) {
        std::stringstream ss(s);
        std::string part;
        while (std::getline(ss, part, ',')) {
            std::stringstream ps(part);
            std::string name; int n = 0;
            if (ps >> name >> n && name != "-") out.push_back({ tables.item(name), n });
        }
    };
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> cols;
        std::stringstream ss(line);
        std::string c;
        while (std::getline(ss, c, '|')) {
            const size_t a = c.find_first_not_of(" \t"), b = c.find_last_not_of(" \t\r");
            cols.push_back(a == std::string::npos ? "" : c.substr(a, b - a + 1));
        }
        if (cols.size() < 5) continue;
        Recipe r;
        r.machine = MachineByName(cols[0]);
        if (r.machine == PartType::Count) continue;
        parseStacks(cols[1], r.out);
        parseStacks(cols[2], r.in);
        r.ticks = std::max(1, atoi(cols[3].c_str()));
        r.eu = std::max(0, atoi(cols[4].c_str()));
        tables.recipes.push_back(r);
    }
    logger.logf("[Factory] %zu machine recipes\n", tables.recipes.size());
}

// ---------------------------------------------------------------------------------------------
// Simulation state (owned by the simulation thread, guarded by simMutex)
// ---------------------------------------------------------------------------------------------

constexpr int kTierLV = 32;                 // volts
constexpr int kBatteryCapacity = 100000;    // EU
constexpr int kBatteryRate = 64;            // EU per tick in or out
constexpr int kConveyorItemsPerSecond = 8;

struct Machine {
    int part = -1;
    PartType type = PartType::Count;
    std::vector<Int3> cells;
    Stack in[kChestSlots], out[kSlots];     // chests use in[] as storage
    Stack fuel;                             // fuel-burning machines and generators: their own fuel slot (#194)
    int water = 0, heat = 0;                // boilers (#229): litres in the tank; heat 0..kHeatMax
    int recipe = -1, progress = 0;          // generators use progress as the banked-fire tick counter
    int energy = 0;                         // battery charge / generator fuel ticks left
    int fuelEU = 0;                         // generator: EU per tick of the current fuel
    int tier = kTierLV;
    int variant = -1;                       // item id: a miner only makes the ore of the deposit under it
    MachineStatus status = MachineStatus::Idle;
    int powerNet = -1, itemNet = -1;
    int received = 0;                       // EU this tick (consumers)
};

struct Net { std::vector<int> machines; int credit = 0, rr = 0, voltage = kTierLV; };

struct PlotSim {
    uint64_t body = 0;
    uint32_t plot = 0, version = 0;
    bool nearPlayer = true;
    std::vector<Machine> machines;
    std::unordered_map<int, int> byPart;
    std::vector<Net> itemNets, powerNets;
};

struct LayoutMsg {
    uint64_t body; uint32_t plot, version;
    std::vector<std::pair<int, PartType>> machines;      // part index, type
    std::vector<std::vector<Int3>> machineCells;
    std::vector<std::string> machineVariant;              // per machine: its one output item ("" any)
    std::vector<Int3> conveyors, cables;
};
struct ItemMsg { int kind; uint64_t body; uint32_t plot; int part; int item; int count; }; // 0 insert, 1 extract, 2 seed, 3 take (one item from a chest), 4 take slot (item = slot, count: 0 input, 1 output, 2 fuel), 5 drain outputs, 6 insert into the fuel slot, 7 water (count = litres)

std::mutex simMutex;                          // state (sim thread holds it for a tick)
std::vector<PlotSim> plots;
std::map<std::tuple<uint64_t, uint32_t, int>, std::vector<std::string>> pendingStates; // from a load

std::mutex queueMutex;                        // main thread -> sim
std::vector<LayoutMsg> layoutQueue;
std::vector<ItemMsg> itemQueue;
std::vector<std::pair<uint64_t, std::pair<uint32_t, bool>>> nearQueue;
struct MoveMsg { uint64_t fromBody; uint32_t fromPlot; uint64_t toBody; uint32_t toPlot; };
std::vector<MoveMsg> moveQueue;               // plots that changed body / id (rocket launch, landing)
bool resetRequested = false;
std::vector<std::pair<int, int>> results;     // sim -> main: (item, count) back to the inventory
std::vector<std::pair<int, int>> refunds;     // sim thread: contents of removed machines (to results)

std::mutex viewMutex;                         // snapshot for the HUD
std::tuple<uint64_t, uint32_t, int> watchKey{ 0, 0, -1 };
MachineView watchView;
std::map<std::tuple<uint64_t, uint32_t, int>, std::vector<std::pair<int, int>>> chestView; // chests of plots near the player: (item, count)
std::map<std::tuple<uint64_t, uint32_t, int>, int> pumpView;     // water pumps (every plot): water units in their outputs
std::map<std::tuple<uint64_t, uint32_t, int>, int> pumpPending;  // pump water drained (queued) but not yet out of pumpView (viewMutex)
inline uint64_t WorkingKey(uint64_t body, uint32_t plot, int part) { return (body * 1315423911ull) ^ ((uint64_t)plot << 24) ^ (uint64_t)(uint32_t)part; }
std::unordered_set<uint64_t> workingView;                          // machines working on near plots (viewMutex)
FactoryStats stats;

std::thread simThread;
std::atomic<bool> running{ false };
std::atomic<long long> tickCount{ 0 };
std::atomic<int> simSpeed{ 1 };
std::condition_variable tickCv;
std::mutex tickCvMutex;

std::unordered_map<std::string, uint32_t> sentVersions; // main thread: layouts already sent

PlotSim* FindPlot(uint64_t body, uint32_t plot) {
    for (auto& p : plots) if (p.body == body && p.plot == plot) return &p;
    return nullptr;
}

// Adds items to a stack array; returns the count that did not fit.
int AddItems(Stack* s, int n, int item, int count) {
    for (int i = 0; i < n && count > 0; ++i)
        if (s[i].item == item && s[i].count < kStackMax) { const int k = std::min(count, kStackMax - s[i].count); s[i].count += k; count -= k; }
    for (int i = 0; i < n && count > 0; ++i)
        if (s[i].count == 0) { s[i].item = item; const int k = std::min(count, kStackMax); s[i].count = k; count -= k; }
    return count;
}
int CountItems(const Stack* s, int n, int item) { int c = 0; for (int i = 0; i < n; ++i) if (s[i].item == item) c += s[i].count; return c; }
void TakeItems(Stack* s, int n, int item, int count) {
    for (int i = 0; i < n && count > 0; ++i)
        if (s[i].item == item) { const int k = std::min(count, s[i].count); s[i].count -= k; count -= k; if (!s[i].count) s[i].item = -1; }
}
int InSlots(const Machine& m) { return m.type == PartType::Chest ? kChestSlots : kSlots; }

bool IsFuel(const MachineDef& def, int item) { for (int f : def.fuels) if (f == item) return true; return false; }
bool RecipeWants(PartType type, int item) {
    for (const Recipe& r : tables.recipes)
        if (r.machine == type) for (const Stack& s : r.in) if (s.item == item) return true;
    return false;
}
bool HasFuelSlot(PartType type) { return !Def(type).fuels.empty(); }
int WaterBucketItem() { static const int id = tables.item("water_bucket"); return id; }
int BucketItem() { static const int id = tables.item("bucket"); return id; }

bool Accepts(PartType type, int item) {
    if (Def(type).power == Power::Storage) return true;
    if (Def(type).waterPerTick > 0 && item == WaterBucketItem()) return true;
    return IsFuel(Def(type), item) || RecipeWants(type, item);
}

// Puts items into a machine: water buckets empty into a boiler's tank (the empty bucket comes back), fuel that no
// recipe of the machine wants goes to its fuel slot, the rest to its inputs (toFuel: the fuel slot whatever, e.g. logs
// for a kiln that also takes them as an input). Returns how many did not fit; `buckets` counts empties given back.
int InsertInto(Machine& m, int item, int count, bool toFuel, int& buckets) {
    const MachineDef& def = Def(m.type);
    buckets = 0;
    if (def.waterPerTick > 0 && item == WaterBucketItem()) {
        while (count > 0 && m.water + kBucketLitres <= kWaterTank) { m.water += kBucketLitres; --count; ++buckets; }
        return count;
    }
    // Fuel goes to the fuel slot. An item that is also a recipe input (logs in a campfire, which chars them) lights an
    // empty fuel slot first and goes to the inputs after that, unless toFuel (Shift+click) asks for the fuel slot.
    const bool dual = RecipeWants(m.type, item);
    if (HasFuelSlot(m.type) && IsFuel(def, item) && (toFuel || !dual || m.fuel.count == 0)) {
        int k = 0;
        if (m.fuel.count == 0 || m.fuel.item == item) { // one kind of fuel at a time
            k = std::min(count, kStackMax - m.fuel.count);
            if (k > 0) { m.fuel.item = item; m.fuel.count += k; }
        }
        count -= k;
        if (count <= 0 || toFuel || !dual) return count;
    }
    if (toFuel || !Accepts(m.type, item)) return count;
    return AddItems(m.in, InSlots(m), item, count);
}

// Burns up to maxTicks of fuel (taking new fuel items from the fuel slot as needed); returns the
// ticks actually burnt. m.energy holds the ticks left on the current item.
int BurnFuel(Machine& m, const MachineDef& def, int maxTicks) {
    int burnt = 0;
    while (burnt < maxTicks) {
        if (m.energy <= 0) {
            if (m.fuel.count <= 0 || !IsFuel(def, m.fuel.item)) break;
            const int f = m.fuel.item;
            if (--m.fuel.count == 0) m.fuel.item = -1;
            auto it = fuelTicks.find(f);
            m.energy = it != fuelTicks.end() ? it->second : 200;
        }
        const int k = std::min(maxTicks - burnt, m.energy);
        m.energy -= k;
        burnt += k;
    }
    return burnt;
}

// A generator's step (#229): it burns only for the EU wanted of it (needEU over these dt ticks: its network's
// consumers and battery room), not every tick whatever the load. A boiler also needs water (litres per tick of
// burning) and heats up: from cold it burns flat out until hot (with or without a cable), making steam in proportion
// to its heat; hot with nothing to power it banks the fire (one tick in eight keeps the heat). Without fuel or water
// it cools. Returns the EU made.
long long GeneratorStep(Machine& m, const MachineDef& def, int dt, long long needEU) {
    const bool steam = def.waterPerTick > 0;
    int want;
    if (steam && m.heat < kHeatMax) want = dt;                                             // warming up
    else want = needEU > 0 ? (int)std::min<long long>(dt, (needEU + def.euOut - 1) / std::max(1, def.euOut)) : 0;
    bool banked = false;
    if (steam && want == 0) { m.progress += dt; want = m.progress / 8; m.progress %= 8; banked = true; }
    if (steam) want = std::min(want, m.water / def.waterPerTick);
    const int burnt = BurnFuel(m, def, want);
    long long eu = 0;
    if (steam) {
        m.water -= burnt * def.waterPerTick;
        if (!banked) eu = (long long)burnt * def.euOut * std::min(m.heat + burnt, kHeatMax) / kHeatMax;
        m.heat = std::min(kHeatMax, m.heat + burnt);
    } else {
        eu = (long long)burnt * def.euOut;
    }
    const bool noWater = steam && m.water < def.waterPerTick;
    const bool noFuel = m.energy <= 0 && m.fuel.count <= 0;
    if (steam && (noWater || noFuel)) m.heat = std::max(0, m.heat - dt);             // the fire is out: it cools
    m.status = noWater ? MachineStatus::NoWater : (noFuel && burnt == 0) ? MachineStatus::NoFuel
             : eu > 0 || (steam && burnt > 0 && !banked) ? MachineStatus::Working : MachineStatus::NoLoad;
    return eu;
}

void BuildNets(PlotSim& ps, const LayoutMsg& L) {
    ps.itemNets.clear();
    ps.powerNets.clear();
    std::unordered_map<Int3, int> machineAt;
    for (size_t i = 0; i < ps.machines.size(); ++i) {
        ps.machines[i].itemNet = ps.machines[i].powerNet = -1;
        for (const Int3& c : ps.machines[i].cells) machineAt[c] = (int)i;
    }
    auto flood = [&](const std::vector<Int3>& cellsIn, std::vector<Net>& nets, bool power) {
        std::unordered_set<Int3> cells(cellsIn.begin(), cellsIn.end()), seen;
        const Int3 dirs[6] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
        for (const Int3& start : cellsIn) {
            if (seen.count(start)) continue;
            Net net;
            std::vector<Int3> stack{ start };
            std::unordered_set<int> members;
            seen.insert(start);
            while (!stack.empty()) {
                const Int3 c = stack.back(); stack.pop_back();
                for (const Int3& d : dirs) {
                    const Int3 n = { c.x + d.x, c.y + d.y, c.z + d.z };
                    if (cells.count(n) && !seen.count(n)) { seen.insert(n); stack.push_back(n); }
                    auto m = machineAt.find(n);
                    if (m != machineAt.end()) members.insert(m->second);
                }
            }
            std::vector<int> sorted(members.begin(), members.end());
            std::sort(sorted.begin(), sorted.end()); // fixed order: determinism
            const int id = (int)nets.size();
            for (int m : sorted) {
                int& slot = power ? ps.machines[m].powerNet : ps.machines[m].itemNet;
                if (slot < 0) { slot = id; net.machines.push_back(m); } // a machine joins one network of each kind
            }
            nets.push_back(std::move(net));
        }
    };
    flood(L.conveyors, ps.itemNets, false);
    flood(L.cables, ps.powerNets, true);
    // A network runs at the highest voltage of its generators and batteries (lower machines refuse it).
    for (Net& net : ps.powerNets) {
        net.voltage = 0;
        for (int m : net.machines) {
            const MachineDef& d = Def(ps.machines[m].type);
            if (d.power == Power::Generator || d.power == Power::Battery) net.voltage = std::max(net.voltage, d.volts);
        }
    }
}

void ApplyLayout(const LayoutMsg& L) {
    PlotSim* ps = FindPlot(L.body, L.plot);
    if (!ps) { plots.push_back({}); ps = &plots.back(); ps->body = L.body; ps->plot = L.plot; }
    ps->version = L.version;
    // Keep the state of machines still there (same part index and type).
    std::vector<Machine> next;
    std::unordered_map<int, int> byPart;
    for (size_t i = 0; i < L.machines.size(); ++i) {
        const int part = L.machines[i].first;
        const PartType type = L.machines[i].second;
        Machine m;
        auto it = ps->byPart.find(part);
        if (it != ps->byPart.end() && ps->machines[it->second].type == type) m = ps->machines[it->second];
        m.part = part; m.type = type; m.cells = L.machineCells[i];
        m.variant = i < L.machineVariant.size() && !L.machineVariant[i].empty() ? tables.item(L.machineVariant[i]) : -1;
        // A loaded save's contents arrive before the layout: apply them once.
        auto pend = pendingStates.find({ L.body, L.plot, part });
        if (pend != pendingStates.end()) {
            const auto& s = pend->second; // recipe progress energy fuelEU, then "in" item count pairs, then "out" pairs
            size_t k = 0;
            auto num = [&]() { return k < s.size() ? atoi(s[k++].c_str()) : 0; };
            m.recipe = num(); m.progress = num(); m.energy = num(); m.fuelEU = num();
            const int nin = num();
            for (int j = 0; j < nin && k + 1 < s.size(); ++j) { const int item = tables.item(s[k]); const int c = atoi(s[k + 1].c_str()); k += 2; AddItems(m.in, InSlots(m), item, c); }
            const int nout = num();
            for (int j = 0; j < nout && k + 1 < s.size(); ++j) { const int item = tables.item(s[k]); const int c = atoi(s[k + 1].c_str()); k += 2; AddItems(m.out, kSlots, item, c); }
            // Then tagged extras: "f item count" (the fuel slot), "w litres", "h heat".
            while (k < s.size()) {
                const std::string tag = s[k++];
                if (tag == "f" && k + 1 < s.size()) { m.fuel.item = tables.item(s[k]); m.fuel.count = atoi(s[k + 1].c_str()); k += 2; }
                else if (tag == "w") m.water = num();
                else if (tag == "h") m.heat = num();
            }
            pendingStates.erase(pend);
        }
        byPart[part] = (int)next.size();
        next.push_back(std::move(m));
    }
    // Machines gone from the layout (the part was removed, or replaced by another kind): what they
    // held goes back to the player instead of vanishing with them.
    for (const Machine& old : ps->machines) {
        auto kept = byPart.find(old.part);
        if (kept != byPart.end() && next[kept->second].type == old.type) continue;
        for (int i = 0; i < InSlots(old); ++i) if (old.in[i].count > 0) refunds.push_back({ old.in[i].item, old.in[i].count });
        for (int i = 0; i < kSlots; ++i) if (old.out[i].count > 0) refunds.push_back({ old.out[i].item, old.out[i].count });
        if (old.fuel.count > 0) refunds.push_back({ old.fuel.item, old.fuel.count });
    }
    ps->machines = std::move(next);
    ps->byPart = std::move(byPart);
    BuildNets(*ps, L);
}

// One step of dt ticks for a plot.
void TickPlot(PlotSim& ps, int dt, long long& made, int& active) {
    // Power: generators burn fuel into EU at their voltage; batteries buffer; consumers on a network
    // above their voltage refuse it.
    for (Net& net : ps.powerNets) {
        long long supply = 0, demand = 0, batteryAvail = 0, batteryRoom = 0;
        for (int mi : net.machines) {
            Machine& m = ps.machines[mi];
            const MachineDef& def = Def(m.type);
            m.received = 0;
            if (def.power == Power::Battery) {
                batteryAvail += std::min(m.energy, kBatteryRate * dt);
                batteryRoom += std::min(kBatteryCapacity - m.energy, kBatteryRate * dt);
            } else if (def.power == Power::Eu) {
                if (def.volts < net.voltage) { m.status = MachineStatus::Overvoltage; continue; }
                if (m.recipe >= 0) demand += (long long)tables.recipes[m.recipe].eu * dt;
            }
        }
        // Generators burn for the consumers and the batteries' room, in a fixed order (#229: not whatever the load).
        for (int mi : net.machines) {
            Machine& m = ps.machines[mi];
            const MachineDef& def = Def(m.type);
            if (def.power == Power::Generator) supply += GeneratorStep(m, def, dt, std::max(0LL, demand + batteryRoom - supply));
        }
        // Share: generation first, then the batteries; a shortfall slows every consumer equally.
        const long long fromBatteries = std::max(0LL, std::min(batteryAvail, demand - supply));
        const long long avail = supply + fromBatteries;
        for (int mi : net.machines) {
            Machine& m = ps.machines[mi];
            const MachineDef& def = Def(m.type);
            if (def.power == Power::Eu && m.recipe >= 0 && def.volts >= net.voltage) {
                const long long need = (long long)tables.recipes[m.recipe].eu * dt;
                m.received = (int)(demand <= avail ? need : need * avail / std::max(1LL, demand));
            }
        }
        long long surplus = supply - std::min(supply, demand);
        long long drain = fromBatteries;
        for (int mi : net.machines) {
            Machine& m = ps.machines[mi];
            if (Def(m.type).power != Power::Battery) continue;
            if (drain > 0) { const long long k = std::min<long long>(drain, std::min(m.energy, kBatteryRate * dt)); m.energy -= (int)k; drain -= k; }
            else if (surplus > 0) { const long long k = std::min<long long>(surplus, std::min(kBatteryCapacity - m.energy, kBatteryRate * dt)); m.energy += (int)k; surplus -= k; }
            m.status = m.energy > 0 ? MachineStatus::Working : MachineStatus::Idle;
        }
    }
    // Machines: start recipes, progress (EU received, fuel burnt, or by themselves), finish into
    // the outputs.
    for (Machine& m : ps.machines) {
        const MachineDef& def = Def(m.type);
        if (def.power == Power::Generator || def.power == Power::Battery || def.power == Power::Storage) {
            // A generator on no network still lights: a boiler warms up and banks its fire (#229).
            if (def.power == Power::Generator && m.powerNet < 0) GeneratorStep(m, def, dt, 0);
            if (def.power == Power::Generator && m.status == MachineStatus::Working) ++active;
            continue;
        }
        // Inputs are only checked when a craft starts and taken when it finishes (#193): take an input out
        // meanwhile and the craft is cancelled, its progress lost. Fuel has its own slot (#194).
        auto inputsHere = [&](const Recipe& r) {
            for (const Stack& s : r.in) if (CountItems(m.in, kSlots, s.item) < s.count) return false;
            return true;
        };
        if (m.recipe >= 0 && !inputsHere(tables.recipes[m.recipe])) { m.recipe = -1; m.progress = 0; }
        if (m.recipe < 0) {
            for (size_t ri = 0; ri < tables.recipes.size(); ++ri) {
                const Recipe& r = tables.recipes[ri];
                if (r.machine != m.type) continue;
                if (m.variant >= 0 && (r.out.empty() || r.out[0].item != m.variant)) continue;
                if (!inputsHere(r)) continue;
                m.recipe = (int)ri; m.progress = 0;
                break;
            }
            if (m.recipe < 0) { m.status = MachineStatus::Idle; continue; } // asleep until items arrive
        }
        const Recipe& r = tables.recipes[m.recipe];
        if (m.progress < r.ticks) {
            int ticks = 0;
            if (def.power == Power::Eu) {
                if (m.powerNet < 0) { m.status = MachineStatus::NoPower; continue; }
                if (def.volts < ps.powerNets[m.powerNet].voltage) { m.status = MachineStatus::Overvoltage; continue; }
                ticks = r.eu > 0 ? m.received / r.eu : dt;
            } else if (def.power == Power::Fuel) {
                ticks = BurnFuel(m, def, std::min(dt, r.ticks - m.progress));
            } else {
                ticks = dt; // runs by itself (coke oven)
            }
            m.progress = std::min(r.ticks, m.progress + ticks);
            m.status = ticks > 0 ? MachineStatus::Working : def.power == Power::Fuel ? MachineStatus::NoFuel : MachineStatus::NoPower;
            if (ticks > 0) ++active;
        }
        if (m.progress >= r.ticks) {
            // Finish only if every output fits.
            Stack trial[kSlots];
            std::copy(m.out, m.out + kSlots, trial);
            bool fits = true;
            for (const Stack& s : r.out) if (AddItems(trial, kSlots, s.item, s.count) > 0) { fits = false; break; }
            if (!fits) { m.status = MachineStatus::OutputFull; continue; }
            for (const Stack& s : r.in) TakeItems(m.in, kSlots, s.item, s.count); // consumed now, at the end (#193)
            std::copy(trial, trial + kSlots, m.out);
            for (const Stack& s : r.out) made += s.count;
            m.recipe = -1; m.progress = 0;
        }
    }
    // Items: each conveyor network moves items from outputs (and chests) to machines that take them.
    for (Net& net : ps.itemNets) {
        net.credit = std::min(net.credit + kConveyorItemsPerSecond * dt, kConveyorItemsPerSecond * kTicksPerSecond * 2);
        const int n = (int)net.machines.size();
        while (net.credit >= kTicksPerSecond && n > 1) {
            bool moved = false;
            for (int s = 0; s < n && !moved; ++s) {
                Machine& src = ps.machines[net.machines[(net.rr + s) % n]];
                const bool chest = src.type == PartType::Chest;
                Stack* from = chest ? src.in : src.out;
                const int fromN = chest ? kChestSlots : kSlots;
                for (int i = 0; i < fromN && !moved; ++i) {
                    if (from[i].count <= 0) continue;
                    const int item = from[i].item;
                    for (int d = 1; d < n && !moved; ++d) {
                        Machine& dst = ps.machines[net.machines[(net.rr + s + d) % n]];
                        if (&dst == &src || (chest && dst.type == PartType::Chest) || !Accepts(dst.type, item)) continue;
                        if (dst.type == PartType::Chest && !chest) { if (AddItems(dst.in, kChestSlots, item, 1) == 0) moved = true; }
                        else if (dst.type != PartType::Chest && CountItems(dst.in, kSlots, item) < kStackMax) {
                            int buckets = 0; // fuel to the fuel slot, water buckets into a boiler's tank (its empty bucket to its outputs)
                            if (InsertInto(dst, item, 1, false, buckets) == 0) { moved = true; if (buckets > 0) AddItems(dst.out, kSlots, BucketItem(), buckets); }
                        }
                        if (moved) { TakeItems(from, fromN, item, 1); }
                    }
                }
            }
            net.rr = (net.rr + 1) % std::max(1, n);
            if (!moved) break;
            net.credit -= kTicksPerSecond;
        }
    }
}

void FillView(const PlotSim& ps, int part, MachineView& v) {
    v = MachineView{};
    auto it = ps.byPart.find(part);
    if (it == ps.byPart.end()) return;
    const Machine& m = ps.machines[it->second];
    v.valid = true;
    v.status = m.status;
    if (m.recipe >= 0) { v.progress = m.progress; v.duration = tables.recipes[m.recipe].ticks; }
    v.energy = m.energy;
    v.energyMax = m.type == PartType::Battery ? kBatteryCapacity : 0;
    v.hasFuelSlot = HasFuelSlot(m.type);
    if (m.fuel.count > 0) v.fuelSlot = { tables.name(m.fuel.item), m.fuel.count };
    if (Def(m.type).waterPerTick > 0) { v.water = m.water; v.waterMax = kWaterTank; v.heat = m.heat; v.heatMax = kHeatMax; }
    auto list = [](const Stack* s, int n) {
        std::string out;
        for (int i = 0; i < n; ++i) if (s[i].count > 0) { if (!out.empty()) out += ", "; out += tables.name(s[i].item) + " x" + std::to_string(s[i].count); }
        return out.empty() ? std::string("-") : out;
    };
    v.input = list(m.in, InSlots(m));
    v.output = list(m.out, kSlots);
    for (int i = 0; i < InSlots(m); ++i) v.inSlots.push_back({ m.in[i].count > 0 ? tables.name(m.in[i].item) : std::string(), m.in[i].count });
    if (m.type != PartType::Chest) for (int i = 0; i < kSlots; ++i) v.outSlots.push_back({ m.out[i].count > 0 ? tables.name(m.out[i].item) : std::string(), m.out[i].count });
}

void SimLoop() {
    using clk = std::chrono::steady_clock;
    auto next = clk::now();
    double emaMs = 0.0, maxMs = 0.0;
    while (running) {
        next += std::chrono::microseconds(1000000 / (kTicksPerSecond * std::max(1, simSpeed.load()))); // tests may run the clock faster
        const auto t0 = clk::now();
        {
            std::lock_guard<std::mutex> lock(simMutex);
            // Commands.
            std::vector<LayoutMsg> layouts;
            std::vector<ItemMsg> items;
            std::vector<std::pair<uint64_t, std::pair<uint32_t, bool>>> nears;
            std::vector<MoveMsg> moves;
            bool reset = false;
            {
                std::lock_guard<std::mutex> q(queueMutex);
                layouts.swap(layoutQueue); items.swap(itemQueue); nears.swap(nearQueue); moves.swap(moveQueue);
                reset = resetRequested; resetRequested = false;
            }
            if (reset) { plots.clear(); emaMs = maxMs = 0.0; }
            // Moved plots first, so the new plot's first layout finds its machines (same part indices).
            for (const MoveMsg& mv : moves)
                if (PlotSim* p = FindPlot(mv.fromBody, mv.fromPlot); p && !FindPlot(mv.toBody, mv.toPlot)) { p->body = mv.toBody; p->plot = mv.toPlot; }
            for (const LayoutMsg& L : layouts) ApplyLayout(L);
            for (const auto& n : nears) if (PlotSim* p = FindPlot(n.first, n.second.first)) p->nearPlayer = n.second.second;
            std::vector<std::pair<int, int>> back;
            back.swap(refunds); // removed machines' contents (from this tick's layouts)
            std::vector<std::pair<std::tuple<uint64_t, uint32_t, int>, int>> drained; // pump drains applied this tick (pumpPending)
            for (const ItemMsg& im : items) {
                PlotSim* p = FindPlot(im.body, im.plot);
                auto it = p ? p->byPart.find(im.part) : decltype(p->byPart.end()){};
                if (!p || it == p->byPart.end()) { if (im.kind == 0) back.push_back({ im.item, im.count }); if (im.kind == 5) drained.push_back({ std::make_tuple(im.body, im.plot, im.part), im.count }); continue; }
                Machine& m = p->machines[it->second];
                if (im.kind == 0 || im.kind == 2 || im.kind == 6) {
                    int buckets = 0;
                    const int left = InsertInto(m, im.item, im.count, im.kind == 6, buckets);
                    if (left > 0 && im.kind != 2) back.push_back({ im.item, left });
                    if (buckets > 0 && im.kind != 2) back.push_back({ BucketItem(), buckets }); // the emptied buckets
                } else if (im.kind == 7) {
                    if (Def(m.type).waterPerTick > 0) m.water = std::min(kWaterTank, m.water + im.count);
                } else if (im.kind == 3) {
                    if (m.type == PartType::Chest) { const int k = std::min(im.count, CountItems(m.in, kChestSlots, im.item)); TakeItems(m.in, kChestSlots, im.item, k); if (k > 0) back.push_back({ im.item, k }); }
                } else if (im.kind == 5) {
                    TakeItems(m.out, kSlots, im.item, im.count); // a pump's water leaves through its pipes (nothing to the inventory)
                    drained.push_back({ std::make_tuple(im.body, im.plot, im.part), im.count });
                } else if (im.kind == 4) {
                    Stack* s = im.count == 2 ? &m.fuel : im.count ? m.out : m.in;
                    const int n = im.count == 2 ? 1 : im.count ? kSlots : InSlots(m);
                    if (im.item >= 0 && im.item < n && s[im.item].count > 0) { back.push_back({ s[im.item].item, s[im.item].count }); s[im.item] = Stack{}; }
                } else {
                    Stack* s = m.type == PartType::Chest ? m.in : m.out;
                    const int n = m.type == PartType::Chest ? kChestSlots : kSlots;
                    for (int i = 0; i < n; ++i) if (s[i].count > 0) { back.push_back({ s[i].item, s[i].count }); s[i] = Stack{}; }
                }
            }
            // Tick: near plots every tick; far plots once a second with a 20x step.
            const long long t = tickCount.load();
            long long made = 0;
            int active = 0, machines = 0, inets = 0, pnets = 0;
            for (size_t i = 0; i < plots.size(); ++i) {
                PlotSim& ps = plots[i];
                machines += (int)ps.machines.size(); inets += (int)ps.itemNets.size(); pnets += (int)ps.powerNets.size();
                if (ps.nearPlayer) TickPlot(ps, 1, made, active);
                else if ((t + (long long)i) % kTicksPerSecond == 0) TickPlot(ps, kTicksPerSecond, made, active);
            }
            const double ms = std::chrono::duration<double, std::milli>(clk::now() - t0).count();
            emaMs = emaMs == 0.0 ? ms : emaMs * 0.95 + ms * 0.05;
            maxMs = std::max(maxMs, ms);
            {
                std::lock_guard<std::mutex> v(viewMutex);
                stats.machines = machines; stats.activeMachines = active; stats.itemNetworks = inets; stats.powerNetworks = pnets;
                stats.avgTickMs = emaMs; stats.maxTickMs = maxMs; stats.ticks = t + 1; stats.itemsMade += made;
                const auto& [b, pl, part] = watchKey;
                if (PlotSim* p = FindPlot(b, pl)) FillView(*p, part, watchView); else watchView = MachineView{};
                for (auto& r : back) results.push_back(r);
                chestView.clear();
                pumpView.clear();
                for (const auto& [k, n] : drained) if (auto pit = pumpPending.find(k); pit != pumpPending.end() && (pit->second -= n) <= 0) pumpPending.erase(pit);
                if (reset) pumpPending.clear();
                workingView.clear(); // machines working on the plots near the player (models, smoke: Factory::Working)
                for (const PlotSim& ps : plots)
                    if (ps.nearPlayer) for (const Machine& m : ps.machines) if (m.status == MachineStatus::Working) workingView.insert(WorkingKey(ps.body, ps.plot, m.part));
                for (const PlotSim& ps : plots)
                    for (const Machine& m : ps.machines)
                        if (m.type == PartType::WaterPump) { int w = 0; for (int i = 0; i < kSlots; ++i) if (m.out[i].count > 0) w += m.out[i].count; pumpView[std::make_tuple(ps.body, ps.plot, m.part)] = w; }
                for (const PlotSim& ps : plots) {
                    if (!ps.nearPlayer) continue;
                    for (const Machine& m : ps.machines) {
                        if (m.type != PartType::Chest) continue;
                        auto& v = chestView[std::make_tuple(ps.body, ps.plot, m.part)];
                        for (int i = 0; i < kChestSlots; ++i) if (m.in[i].count > 0) v.push_back({ m.in[i].item, m.in[i].count });
                    }
                }
            }
        }
        tickCount++;
        tickCv.notify_all();
        std::this_thread::sleep_until(next);
        if (clk::now() > next + std::chrono::milliseconds(250)) next = clk::now(); // don't race to catch up after a stall
    }
}

std::string PlotKey(uint64_t body, uint32_t plot) { return std::to_string(body) + ":" + std::to_string(plot); }

} // namespace

void Factory::Start() {
    if (running) return;
    LoadRecipes();
    running = true;
    simThread = std::thread(SimLoop);
}

void Factory::Stop() {
    if (!running) return;
    running = false;
    if (simThread.joinable()) simThread.join();
}

void Factory::Reset() {
    std::lock_guard<std::mutex> q(queueMutex);
    resetRequested = true;
    layoutQueue.clear();
    itemQueue.clear();
    sentVersions.clear();
}

void Factory::Sync(Scene& world) {
    if (!running) Start();
    static double lastNear = 0.0;
    const bool nearCheck = GetTime() - lastNear > 0.5;
    if (nearCheck) lastNear = GetTime();
    for (auto& obj : world.objects) {
        auto* body = dynamic_cast<Planetoid*>(obj.get());
        if (!body || !body->plots) continue;
        for (auto& pl : body->plots->plots) {
            const std::string key = PlotKey(body->id, pl->id);
            auto sv = sentVersions.find(key);
            if (sv == sentVersions.end() || sv->second != pl->version) {
                LayoutMsg L{ body->id, pl->id, pl->version, {}, {}, {}, {}, {} };
                std::vector<std::pair<Int3, CellShape>> fp;
                for (size_t i = 0; i < pl->parts.size(); ++i) {
                    const Part& part = pl->parts[i];
                    if (part.type == PartType::Count) continue;
                    if (IsMachine(part.type)) {
                        Plot::Footprint(part, fp);
                        std::vector<Int3> cells;
                        for (auto& c : fp) cells.push_back(c.first);
                        L.machines.push_back({ (int)i, part.type });
                        L.machineCells.push_back(std::move(cells));
                        // A miner works the ore deposit 16 m under it.
                        std::string variant;
                        if (part.type == PartType::Miner) {
                            const Vector3 at = pl->toBody({ (part.origin.x + 1.0f) * kCellSize, (part.origin.y - 32.0f) * kCellSize, (part.origin.z + 1.0f) * kCellSize });
                            // The deposit within 40 m under it, or (without one) the scattered common ore there.
                            const auto* pf = dynamic_cast<const PlanetField*>(body->field.get());
                            const int dep = pf ? pf->depositNear(at, 40.0f) : -1;
                            variant = dep >= 0 ? PlanetField::OreKindName(dep) : pf ? (OreKindOf(*body, at) ? OreKindOf(*body, at) : "iron_ore") : OreKindAt(at, body->seed);
                        }
                        L.machineVariant.push_back(variant);
                    } else if (part.type == PartType::Conveyor) L.conveyors.push_back(part.origin);
                    else if (part.type == PartType::Cable) L.cables.push_back(part.origin);
                }
                std::lock_guard<std::mutex> q(queueMutex);
                layoutQueue.push_back(std::move(L));
                sentVersions[key] = pl->version;
            }
            if (nearCheck) {
                const Vector3 o = body->toWorld({ (float)pl->origin.x, (float)pl->origin.y, (float)pl->origin.z });
                const bool near = Vector3Distance(o, world.camera.position) < 1000.0f;
                std::lock_guard<std::mutex> q(queueMutex);
                nearQueue.push_back({ body->id, { pl->id, near } });
            }
        }
    }
    // Items back from the simulation.
    std::vector<std::pair<int, int>> back;
    {
        std::lock_guard<std::mutex> v(viewMutex);
        back.swap(results);
    }
    for (auto& r : back) inventory.add(tables.name(r.first), r.second);
}

MachineView Factory::View(uint64_t body, uint32_t plot, int part) {
    std::lock_guard<std::mutex> v(viewMutex);
    const auto key = std::make_tuple(body, plot, part);
    if (watchKey != key) { watchKey = key; return MachineView{}; } // filled from the next tick on
    return watchView;
}

FactoryStats Factory::Stats() {
    std::lock_guard<std::mutex> v(viewMutex);
    return stats;
}

void Factory::Insert(uint64_t body, uint32_t plot, int part, const std::string& item, int count) {
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back({ 0, body, plot, part, tables.item(item), count });
}
void Factory::Extract(uint64_t body, uint32_t plot, int part) {
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back({ 1, body, plot, part, -1, 0 });
}
void Factory::TakeSlot(uint64_t body, uint32_t plot, int part, int slot, int which) {
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back({ 4, body, plot, part, which == 2 ? 0 : slot, which });
}
void Factory::InsertFuel(uint64_t body, uint32_t plot, int part, const std::string& item, int count) {
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back({ 6, body, plot, part, tables.item(item), count });
}
void Factory::AddWater(uint64_t body, uint32_t plot, int part, int litres) {
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back({ 7, body, plot, part, -1, litres });
}
bool Factory::TakesWater(PartType type) { LoadRecipes(); return Def(type).waterPerTick > 0; }
void Factory::Seed(uint64_t body, uint32_t plot, int part, const std::string& item, int count) {
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back({ 2, body, plot, part, tables.item(item), count });
}

void Factory::WaitTicks(int n) {
    const long long target = tickCount.load() + n;
    std::unique_lock<std::mutex> lk(tickCvMutex);
    tickCv.wait_for(lk, std::chrono::seconds(30), [&] { return tickCount.load() >= target; });
}

void Factory::Save(FILE* f, uint64_t bodyId) {
    std::lock_guard<std::mutex> lock(simMutex);
    int n = 0;
    for (const PlotSim& ps : plots) if (ps.body == bodyId) n += (int)ps.machines.size();
    fprintf(f, "machines %d\n", n);
    for (const PlotSim& ps : plots) {
        if (ps.body != bodyId) continue;
        for (const Machine& m : ps.machines) {
            fprintf(f, "m %u %d %d %d %d %d", ps.plot, m.part, m.recipe, m.progress, m.energy, m.fuelEU);
            int nin = 0; for (int i = 0; i < InSlots(m); ++i) if (m.in[i].count > 0) ++nin;
            fprintf(f, " %d", nin);
            for (int i = 0; i < InSlots(m); ++i) if (m.in[i].count > 0) fprintf(f, " %s %d", tables.name(m.in[i].item).c_str(), m.in[i].count);
            int nout = 0; for (int i = 0; i < kSlots; ++i) if (m.out[i].count > 0) ++nout;
            fprintf(f, " %d", nout);
            for (int i = 0; i < kSlots; ++i) if (m.out[i].count > 0) fprintf(f, " %s %d", tables.name(m.out[i].item).c_str(), m.out[i].count);
            if (m.fuel.count > 0) fprintf(f, " f %s %d", tables.name(m.fuel.item).c_str(), m.fuel.count);
            if (m.water > 0) fprintf(f, " w %d", m.water);
            if (m.heat > 0) fprintf(f, " h %d", m.heat);
            fprintf(f, "\n");
        }
    }
}

std::vector<std::string> Factory::ReadSaved(FILE* f) {
    std::vector<std::string> lines;
    int n = 0;
    if (fscanf(f, " machines %d", &n) != 1) return lines; // none here (older saves)
    char line[4096];
    if (!fgets(line, sizeof line, f)) return lines; // rest of the header line
    for (int i = 0; i < n; ++i) {
        if (!fgets(line, sizeof line, f)) break;
        std::string s = line;
        while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
        lines.push_back(s);
    }
    return lines;
}

void Factory::LoadSaved(uint64_t bodyId, const std::vector<std::string>& lines) {
    LoadRecipes();
    std::lock_guard<std::mutex> lock(simMutex);
    for (const std::string& line : lines) {
        std::stringstream ss(line);
        std::string tag; uint32_t plot; int part;
        if (!(ss >> tag >> plot >> part)) continue;
        std::vector<std::string> rest;
        std::string w;
        while (ss >> w) rest.push_back(w);
        pendingStates[{ bodyId, plot, part }] = rest;
    }
}

void Factory::Load(FILE* f, uint64_t bodyId) { LoadSaved(bodyId, ReadSaved(f)); }

void Factory::MovePlot(uint64_t fromBody, uint32_t fromPlot, uint64_t toBody, uint32_t toPlot) {
    std::lock_guard<std::mutex> q(queueMutex);
    moveQueue.push_back({ fromBody, fromPlot, toBody, toPlot });
}

bool Factory::Accepts(PartType type, const std::string& item) {
    LoadRecipes();
    const int id = tables.find(item);
    return id >= 0 && ::Accepts(type, id);
}

void Factory::SetSpeed(int ticksPerTick) { simSpeed = std::max(1, ticksPerTick); }

void Factory::Take(uint64_t body, uint32_t plot, int part, const std::string& item, int count) {
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back({ 3, body, plot, part, tables.item(item), count });
}

std::vector<std::pair<std::string, int>> Factory::ChestContents(uint64_t body, uint32_t plot, int part) {
    std::vector<std::pair<std::string, int>> out;
    std::lock_guard<std::mutex> v(viewMutex);
    auto it = chestView.find(std::make_tuple(body, plot, part));
    if (it == chestView.end()) return out;
    for (const auto& s : it->second) {
        const std::string name = tables.name(s.first);
        bool merged = false;
        for (auto& o : out) if (o.first == name) { o.second += s.second; merged = true; }
        if (!merged) out.push_back({ name, s.second });
    }
    return out;
}

int Factory::PumpWater(uint64_t body, uint32_t plot, int part) {
    std::lock_guard<std::mutex> v(viewMutex);
    auto it = pumpView.find(std::make_tuple(body, plot, part));
    const int have = it == pumpView.end() ? 0 : it->second;
    auto pit = pumpPending.find(std::make_tuple(body, plot, part));
    return std::max(0, have - (pit == pumpPending.end() ? 0 : pit->second)); // not what is already on its way out
}

void Factory::Drain(uint64_t body, uint32_t plot, int part, const std::string& item, int count) {
    { std::lock_guard<std::mutex> v(viewMutex); pumpPending[std::make_tuple(body, plot, part)] += count; } // PumpWater won't offer it twice
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back({ 5, body, plot, part, tables.item(item), count });
}

bool Factory::Working(uint64_t body, uint32_t plot, int part) {
    std::lock_guard<std::mutex> v(viewMutex);
    return workingView.count(WorkingKey(body, plot, part)) > 0;
}
