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
};
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
    int recipe = -1, progress = 0;
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
struct ItemMsg { int kind; uint64_t body; uint32_t plot; int part; int item; int count; }; // 0 insert, 1 extract, 2 seed

std::mutex simMutex;                          // state (sim thread holds it for a tick)
std::vector<PlotSim> plots;
std::map<std::tuple<uint64_t, uint32_t, int>, std::vector<std::string>> pendingStates; // from a load

std::mutex queueMutex;                        // main thread -> sim
std::vector<LayoutMsg> layoutQueue;
std::vector<ItemMsg> itemQueue;
std::vector<std::pair<uint64_t, std::pair<uint32_t, bool>>> nearQueue;
bool resetRequested = false;
std::vector<std::pair<int, int>> results;     // sim -> main: (item, count) back to the inventory
std::vector<std::pair<int, int>> refunds;     // sim thread: contents of removed machines (to results)

std::mutex viewMutex;                         // snapshot for the HUD
std::tuple<uint64_t, uint32_t, int> watchKey{ 0, 0, -1 };
MachineView watchView;
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

bool Accepts(PartType type, int item) {
    if (Def(type).power == Power::Storage) return true;
    for (int f : Def(type).fuels) if (f == item) return true;
    for (const Recipe& r : tables.recipes)
        if (r.machine == type) for (const Stack& s : r.in) if (s.item == item) return true;
    return false;
}

// Burns up to maxTicks of fuel (taking new fuel items from the inputs as needed); returns the
// ticks actually burnt. m.energy holds the ticks left on the current item.
int BurnFuel(Machine& m, const MachineDef& def, int maxTicks) {
    int burnt = 0;
    while (burnt < maxTicks) {
        if (m.energy <= 0) {
            bool lit = false;
            for (int f : def.fuels)
                if (CountItems(m.in, kSlots, f) > 0) {
                    TakeItems(m.in, kSlots, f, 1);
                    auto it = fuelTicks.find(f);
                    m.energy = it != fuelTicks.end() ? it->second : 200;
                    lit = true;
                    break;
                }
            if (!lit) break;
        }
        const int k = std::min(maxTicks - burnt, m.energy);
        m.energy -= k;
        burnt += k;
    }
    return burnt;
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
        long long supply = 0, demand = 0, batteryAvail = 0;
        for (int mi : net.machines) {
            Machine& m = ps.machines[mi];
            const MachineDef& def = Def(m.type);
            m.received = 0;
            if (def.power == Power::Generator) {
                const int burn = BurnFuel(m, def, dt);
                supply += (long long)burn * def.euOut;
                m.status = burn > 0 ? MachineStatus::Working : MachineStatus::NoFuel;
            } else if (def.power == Power::Battery) {
                batteryAvail += std::min(m.energy, kBatteryRate * dt);
            } else if (def.power == Power::Eu) {
                if (def.volts < net.voltage) { m.status = MachineStatus::Overvoltage; continue; }
                if (m.recipe >= 0) demand += (long long)tables.recipes[m.recipe].eu * dt;
            }
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
            if (def.power == Power::Generator && m.status == MachineStatus::Working) ++active;
            continue;
        }
        if (m.recipe < 0) {
            // Fuel machines: recipes that don't eat fuel come first (a campfire with wood and a
            // crucible mould fires the mould rather than charring its fuel), and a recipe that does
            // eat fuel must leave something burning or to burn.
            auto isFuel = [&](int item) { for (int f : def.fuels) if (f == item) return true; return false; };
            for (int pass = 0; pass < 2 && m.recipe < 0; ++pass)
                for (size_t ri = 0; ri < tables.recipes.size(); ++ri) {
                    const Recipe& r = tables.recipes[ri];
                    if (r.machine != m.type) continue;
                    if (m.variant >= 0 && (r.out.empty() || r.out[0].item != m.variant)) continue;
                    bool eatsFuel = false;
                    for (const Stack& s : r.in) eatsFuel |= isFuel(s.item);
                    if (def.power == Power::Fuel && eatsFuel != (pass == 1)) continue;
                    if (def.power != Power::Fuel && pass == 1) continue;
                    bool ok = true;
                    // Nothing burning: one more of the fuel it eats must be left to light (only the fuel; an
                    // extra of every input stalled recipes holding exactly enough of the rest).
                    for (const Stack& s : r.in) if (CountItems(m.in, kSlots, s.item) < s.count + (isFuel(s.item) && m.energy <= 0 && def.power == Power::Fuel ? 1 : 0)) { ok = false; break; }
                    if (!ok) continue;
                    for (const Stack& s : r.in) TakeItems(m.in, kSlots, s.item, s.count);
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
                        else if (dst.type != PartType::Chest && CountItems(dst.in, kSlots, item) < kStackMax && AddItems(dst.in, kSlots, item, 1) == 0) moved = true;
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
    auto list = [](const Stack* s, int n) {
        std::string out;
        for (int i = 0; i < n; ++i) if (s[i].count > 0) { if (!out.empty()) out += ", "; out += tables.name(s[i].item) + " x" + std::to_string(s[i].count); }
        return out.empty() ? std::string("-") : out;
    };
    v.input = list(m.in, InSlots(m));
    v.output = list(m.out, kSlots);
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
            bool reset = false;
            {
                std::lock_guard<std::mutex> q(queueMutex);
                layouts.swap(layoutQueue); items.swap(itemQueue); nears.swap(nearQueue);
                reset = resetRequested; resetRequested = false;
            }
            if (reset) { plots.clear(); emaMs = maxMs = 0.0; }
            for (const LayoutMsg& L : layouts) ApplyLayout(L);
            for (const auto& n : nears) if (PlotSim* p = FindPlot(n.first, n.second.first)) p->nearPlayer = n.second.second;
            std::vector<std::pair<int, int>> back;
            back.swap(refunds); // removed machines' contents (from this tick's layouts)
            for (const ItemMsg& im : items) {
                PlotSim* p = FindPlot(im.body, im.plot);
                auto it = p ? p->byPart.find(im.part) : decltype(p->byPart.end()){};
                if (!p || it == p->byPart.end()) { if (im.kind == 0) back.push_back({ im.item, im.count }); continue; }
                Machine& m = p->machines[it->second];
                if (im.kind == 0 || im.kind == 2) {
                    const int left = Accepts(m.type, im.item) ? AddItems(m.in, InSlots(m), im.item, im.count) : im.count;
                    if (left > 0 && im.kind == 0) back.push_back({ im.item, left });
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
                            variant = OreKindAt(at, body->seed);
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
            fprintf(f, "\n");
        }
    }
}

void Factory::Load(FILE* f, uint64_t bodyId) {
    LoadRecipes();
    int n = 0;
    if (fscanf(f, " machines %d", &n) != 1) return;
    char line[4096];
    fgets(line, sizeof line, f); // rest of the header line
    std::lock_guard<std::mutex> lock(simMutex);
    for (int i = 0; i < n; ++i) {
        if (!fgets(line, sizeof line, f)) break;
        std::stringstream ss(line);
        std::string tag; uint32_t plot; int part;
        ss >> tag >> plot >> part;
        std::vector<std::string> rest;
        std::string w;
        while (ss >> w) rest.push_back(w);
        pendingStates[{ bodyId, plot, part }] = rest;
    }
}

bool Factory::Accepts(PartType type, const std::string& item) {
    LoadRecipes();
    const int id = tables.find(item);
    return id >= 0 && ::Accepts(type, id);
}

void Factory::SetSpeed(int ticksPerTick) { simSpeed = std::max(1, ticksPerTick); }
