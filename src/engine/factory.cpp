#include "factory.h"
#include "inventory.h"
#include "plots.h"
#include "rotationplot.h"
#include "rotationsources.h"
#include "utils/logger.h"
#include "worldtypes.h"
#include "simclock.h"
#include "scene.h"
#include <raymath.h>

#include <raylib.h>
#include <algorithm>
#include <array>
#include <cstdlib>
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
// #580: results taken from a heating machine's output carry this bit on their item index: they reach the inventory hot.
constexpr int kHotItem = 1 << 30;
static bool HeatsItems(PartType t) {
    switch (t) {
        case PartType::Furnace: case PartType::Kiln: case PartType::CokeOven: case PartType::BlastFurnace:
        case PartType::AlloySmelter: case PartType::ArcFurnace: case PartType::Campfire: return true;
        default: return false;
    }
}
constexpr int kSlots = 4;
// A slot holds as much as a player inventory cell (#566: StackLimit: 8192 resources, 256 blocks, 1 tool), but conveyors and
// inserters keep only kFeedMax of an item in a machine's inputs (a buffer, not a store).
constexpr int kFeedMax = 64;
constexpr int kChestSlots = Inventory::kGrid; // as many cells as the player's main inventory (#495: 50, drawn 10 wide)

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
        const int id = (int)itemNames.size();
        if (id < kMaxCapItems) caps[id] = StackLimit(n); // written before the id is handed out (under the lock)
        itemIds[n] = id;
        itemNames.push_back(n);
        return id;
    }
    // Whole units of the item a slot holds (#566: the player inventory's StackLimit). Lock-free: the sim thread asks
    // per stack move, and an id's cap never changes once it is interned.
    int cap(int id) const { return id >= 0 && id < kMaxCapItems ? caps[id] : 64; }
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
    static constexpr int kMaxCapItems = 4096;
    int caps[kMaxCapItems] = {};
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
enum class Power : uint8_t { None, Fuel, Eu, Generator, Storage, Battery, Rotation }; // Rotation (#1269): at its shaft's speed
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
                : p == "storage" ? Power::Storage : p == "battery" ? Power::Battery : p == "rotation" ? Power::Rotation : Power::None;
        d.euOut = atoi(cols[3].c_str());
        if (cols.size() > 4) { std::stringstream fs(cols[4]); std::string fuel; while (fs >> fuel) d.fuels.push_back(tables.item(fuel)); }
        if (cols.size() > 5) d.waterPerTick = std::max(0, atoi(cols[5].c_str()));
    }
}

// (#1207) built once by LoadRecipes
std::array<std::vector<int>, 256> recipesOf;      // by PartType: indices into tables.recipes
std::array<std::vector<uint8_t>, 256> wantsOf;   // by PartType, by item id: a recipe of it takes the item
// (#1214) sleeping machines and networks: on unless RAYLIBGAME_FACTORY_SLEEP=0 (A/B; the unit test runs both)
bool sleepOn = [] { const char* e = getenv("RAYLIBGAME_FACTORY_SLEEP"); return !(e && e[0] == '0'); }();

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
    // (#1207) Per machine type: its recipes in file order, and which items they take (a recipe start and Accepts looked
    // through every recipe of every machine, each tick).
    for (size_t ri = 0; ri < tables.recipes.size(); ++ri) {
        const Recipe& r = tables.recipes[ri];
        recipesOf[(int)r.machine].push_back((int)ri);
        std::vector<uint8_t>& w = wantsOf[(int)r.machine];
        for (const Stack& s : r.in) { if (s.item >= (int)w.size()) w.resize(s.item + 1, 0); w[s.item] = 1; }
    }
    LOG_INFO("Factory", "%zu machine recipes\n", tables.recipes.size());
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
    // (#1214) asleep: its recipe search found nothing and its contents have not changed since; dirty: they changed (an insert,
    // a take, fuel taken, a craft, a variant, fire or layout change) - it wakes, and its conveyor network wakes
    bool asleep = false, dirty = true;
    PartType type = PartType::Count;
    std::vector<Int3> cells;
    std::vector<Stack> in;                  // InSlots(): kSlots inputs, or a chest's kChestSlots of storage (#495: only chests pay for 50)
    Stack out[kSlots];
    Stack fuel;                             // fuel-burning machines and generators: their own fuel slot (#194)
    int water = 0, heat = 0;                // boilers (#229): litres in the tank; heat 0..kHeatMax
    uint8_t fire = kFireAuto;               // campfires (#1103): lit or put out by a player, else by itself
    int recipe = -1, progress = 0;          // generators use progress as the banked-fire tick counter
    int energy = 0;                         // battery charge / generator fuel ticks left
    int fuelEU = 0;                         // generator: EU per tick of the current fuel
    int tier = kTierLV;
    int variant = -1;                       // item id: a miner only makes the ore of the deposit under it
    MachineStatus status = MachineStatus::Idle;
    int powerNet = -1, itemNet = -1;
    int received = 0;                       // EU this tick (consumers)
    int rotAcc = 0;                         // (#1269) rotation machines: permille of a tick of progress carried over
};

struct Net {
    std::vector<int> machines; int credit = 0, rr = 0, voltage = kTierLV;
    bool asleep = false; // (#1214) its last pass moved nothing and none of its machines changed since
    long long supply = 0, demand = 0; int ticks = 1; // the last step's EU made and wanted, over its ticks (cable tooltips, #570)
};

// (#1199) Production statistics: a count kept in Factorio's windows, rolled once a second of factory time: the last minute by
// seconds, the last hour by minutes, the last 10 hours by 10 minutes, and the total. Sim thread only (under simMutex).
constexpr int kStatBuckets = 60;
struct Series {
    std::array<int32_t, kStatBuckets> sec{}, min{};      // (a minute of a big HV plot is ~6e8 EU: fits)
    std::array<int64_t, kStatBuckets> ten{};             // (10 minutes of it does not)
    int32_t acc = 0, minAcc = 0;
    int64_t tenAcc = 0;
    long long total = 0;
    void add(long long n) { acc += (int32_t)n; }
};
struct ItemSeries { Series made, used, handIn, handOut; }; // machine crafts and fuel burnt; a player's inserts and takes
struct PlotStats {
    std::map<int, ItemSeries> items;                        // by item id (ordered: deterministic)
    std::map<int, Series> euMade, euUsed;                   // EU by machine kind (PartType)
    Series euDemand;                                        // EU the consumers asked for (satisfaction = used / demand)
    std::map<int, std::array<long long, 8>> statusTicks;    // ticks per status, by machine kind (the bottleneck finder)
};

struct PlotSim {
    uint64_t body = 0;
    uint32_t plot = 0, version = 0;
    bool nearPlayer = true;
    std::vector<Machine> machines;
    std::unordered_map<int, int> byPart;
    std::vector<Net> itemNets, powerNets;
    std::unordered_map<Int3, int> cableNet; // cable cell -> its power network (#570)
    PlotStats stats;                        // (#1199)
    RotationPlot rot;                       // (#1269) its shaft networks (per plot: ticks on any worker)
    // (#1207) what this plot last put in the near-plot views (working machines, chests and a hash of their contents)
    struct Published { bool near = false; std::vector<int> working, chests; uint64_t chestHash = 0; } pub;
};
thread_local PlotStats* gCount = nullptr;  // the plot ticking now on this thread (BurnFuel counts the fuel it takes into it; thread_local: plots may tick on several threads later, #997)
long long statSeconds = 0;    // seconds rolled since the reset (how much of each window is filled)
// (#1207) The last minute of factory ticks' times (sim thread), for the mean, p99 and worst and the phase split.
struct TickTimes { float total = 0, cmd = 0, tick = 0, view = 0; };
constexpr size_t kTickRing = 1200;
std::array<TickTimes, kTickRing> tickRing{};
size_t tickRingAt = 0;

// Once a second: each count's second into the rings, minutes and 10 minutes as they complete.
void RollSeries(Series& s) {
    const int slot = (int)(statSeconds % kStatBuckets);
    s.sec[slot] = s.acc; s.minAcc += s.acc; s.total += s.acc; s.acc = 0;
    if ((statSeconds + 1) % 60 == 0) { s.min[(statSeconds / 60) % kStatBuckets] = s.minAcc; s.tenAcc += s.minAcc; s.minAcc = 0; }
    if ((statSeconds + 1) % 600 == 0) { s.ten[(statSeconds / 600) % kStatBuckets] = s.tenAcc; s.tenAcc = 0; }
}
void RollStats(PlotStats& st) {
    for (auto& [id, is] : st.items) { RollSeries(is.made); RollSeries(is.used); RollSeries(is.handIn); RollSeries(is.handOut); }
    for (auto& [k, s] : st.euMade) RollSeries(s);
    for (auto& [k, s] : st.euUsed) RollSeries(s);
    RollSeries(st.euDemand);
}

struct LayoutMsg {
    uint64_t body; uint32_t plot, version;
    std::vector<std::pair<int, PartType>> machines;      // part index, type
    std::vector<std::vector<Int3>> machineCells;
    std::vector<std::string> machineVariant;              // per machine: its one output item ("" any)
    std::vector<Int3> conveyors, cables;
    uint32_t refundTo = 0;                                // (#854) who removed machines in this change (their contents go to that player)
    std::vector<RotationPartMsg> rotation;                // (#1269) shafts, cogs, cranks and rotation machines
};
struct ItemMsg { int kind; uint64_t body; uint32_t plot; int part; int item; int count; int amount = 0; uint32_t to = 0; }; // amount: kind 8, how many (0 = all); kind 9, the slot // 0 insert, 1 extract, 2 seed, 3 take (one item from a chest), 4 take slot (item = slot, count: 0 input, 1 output, 2 fuel), 5 drain outputs, 6 insert into the fuel slot, 7 water (count = litres), 8 slot to the cursor, 9 insert into one chest slot (#862)

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
struct ItemOut { int item = 0, count = 0; uint32_t to = 0; }; // (#332) to: the player it goes to (0: this game's own inventory)
std::vector<ItemOut> results;     // sim -> main: (item, count) back to the inventory
std::vector<ItemOut> cursorResults; // sim -> main: (item, count) onto the cursor (#192)
void (*cursorSink)(const std::string&, int, float) = nullptr;
std::vector<ItemOut> refunds;     // sim thread: contents of removed machines (to results)

std::mutex viewMutex;                         // snapshot for the HUD
std::tuple<uint64_t, uint32_t, int> watchKey{ 0, 0, -1 };
MachineView watchView;
std::tuple<uint64_t, uint32_t, int, int, int> cableKey{ 0, 0, 0, 0, 0 }; // body, plot, cell (#570)
CableView cableView;
std::map<std::tuple<uint64_t, uint32_t, int>, std::vector<std::pair<int, int>>> chestView; // chests of plots near the player: (item, count)
std::map<std::tuple<uint64_t, uint32_t, int>, int> pumpView;     // water pumps (every plot): water units in their outputs
std::map<std::tuple<uint64_t, uint32_t, int>, int> pumpPending;  // pump water drained (queued) but not yet out of pumpView (viewMutex)
std::map<std::tuple<uint64_t, uint32_t, int>, int> tankView;     // boilers (every plot): litres of room in their tanks (#229 review)
std::map<std::tuple<uint64_t, uint32_t, int>, int> tankPending;  // litres sent by AddWater but not yet in a tank (viewMutex)
// A 64-bit mix of (body, plot, part) (#180: the old XOR of shifted fields collided easily).
inline uint64_t WorkingKey(uint64_t body, uint32_t plot, int part) {
    auto mix = [](uint64_t x) { x += 0x9e3779b97f4a7c15ull; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull; x = (x ^ (x >> 27)) * 0x94d049bb133111ebull; return x ^ (x >> 31); };
    return mix(mix(body) ^ (((uint64_t)plot << 32) | (uint32_t)part));
}
bool viewsRebuild = false; // (#1207) a moved plot: rebuild the near-plot views from scratch
std::unordered_set<uint64_t> workingView;                          // machines working on near plots (viewMutex)
std::map<std::pair<uint64_t, uint32_t>, std::vector<int>> workingByPlot; // (#332) the same, by plot
// (#1269) every rotation part of the plots near the player, as of the last tick (viewMutex): models, gauges, tooltips
struct RotPub { int part = -1; int64_t milli = 0; Rotation::Status status = Rotation::Status::Stopped; Rotation::NetResult net; };
std::map<std::pair<uint64_t, uint32_t>, std::vector<RotPub>> rotView;
std::vector<Factory::MachineKey> watchMany;                        // (#332) a server's open machines (viewMutex)
std::map<Factory::MachineKey, MachineView> watchManyViews;
std::vector<Vector3> nearPoints;                                    // (#332) main thread: a server's players
bool remoteMode = false;                                            // (#332) a client: views from the server
std::map<Factory::MachineKey, MachineView> remoteViews;
std::unordered_set<uint64_t> remoteWorking;
std::map<std::pair<uint64_t, uint32_t>, std::vector<int>> remoteWorkingByPlot;
FactoryStats stats;

std::thread simThread;
std::atomic<bool> running{ false };
std::atomic<long long> tickCount{ 0 };
std::atomic<int> simSpeed{ 1 };
std::condition_variable tickCv;
std::mutex tickCvMutex;
// The sim clock (#315): the main thread grants ticks from its simulation steps (Factory::Advance) and WaitTicks; the
// thread runs one per credit, so the factory pauses with the game and keeps its pace at any frame rate.
std::atomic<long long> tickCredit{ 0 };
std::condition_variable creditCv;
std::mutex creditMutex;
double creditAccum = 0.0; // main thread: fractional ticks owed
void GrantTicks(long long n, long long cap) {
    {
        std::lock_guard<std::mutex> lk(creditMutex);
        long long c = tickCredit.load() + n;
        if (cap > 0 && c > cap) c = std::max(tickCredit.load(), cap); // don't race to catch up after a stall (never drops WaitTicks credit)
        tickCredit = c;
    }
    creditCv.notify_all();
}

std::unordered_map<std::string, uint32_t> sentVersions; // main thread: layouts already sent

// (#1207) FindPlot through a hash of (body, plot) -> index, rebuilt when the plot list changes (a new plot, a reset, a
// moved plot: plotIndexDirty); every hit is checked against the plot itself. Sim thread (and callers under simMutex).
std::unordered_map<uint64_t, size_t> plotIndex;
bool plotIndexDirty = true;
uint64_t PlotHash(uint64_t body, uint32_t plot) { return body * 0x9E3779B97F4A7C15ull ^ (uint64_t)plot; }
PlotSim* FindPlot(uint64_t body, uint32_t plot) {
    if (plotIndexDirty) {
        plotIndex.clear();
        for (size_t i = 0; i < plots.size(); ++i) plotIndex[PlotHash(plots[i].body, plots[i].plot)] = i;
        plotIndexDirty = false;
    }
    auto it = plotIndex.find(PlotHash(body, plot));
    if (it != plotIndex.end() && it->second < plots.size() && plots[it->second].body == body && plots[it->second].plot == plot) return &plots[it->second];
    for (auto& p : plots) if (p.body == body && p.plot == plot) return &p; // (a hash collision: never in practice)
    return nullptr;
}

// Adds items to a stack array; returns the count that did not fit.
int AddItems(Stack* s, int n, int item, int count) {
    const int cap = tables.cap(item);
    for (int i = 0; i < n && count > 0; ++i)
        if (s[i].item == item && s[i].count < cap) { const int k = std::min(count, cap - s[i].count); s[i].count += k; count -= k; }
    for (int i = 0; i < n && count > 0; ++i)
        if (s[i].count == 0) { s[i].item = item; const int k = std::min(count, cap); s[i].count = k; count -= k; }
    return count;
}
bool RoomFor(const Stack* s, int n, int item) { for (int i = 0; i < n; ++i) if (s[i].count == 0 || (s[i].item == item && s[i].count < tables.cap(item))) return true; return false; }
int CountItems(const Stack* s, int n, int item) { int c = 0; for (int i = 0; i < n; ++i) if (s[i].item == item) c += s[i].count; return c; }
void TakeItems(Stack* s, int n, int item, int count) {
    for (int i = 0; i < n && count > 0; ++i)
        if (s[i].item == item) { const int k = std::min(count, s[i].count); s[i].count -= k; count -= k; if (!s[i].count) s[i].item = -1; }
}
int InSlots(const Machine& m) { return m.type == PartType::Chest ? kChestSlots : kSlots; }

bool IsFuel(const MachineDef& def, int item) { for (int f : def.fuels) if (f == item) return true; return false; }
bool RecipeWants(PartType type, int item) {
    const std::vector<uint8_t>& w = wantsOf[(int)type]; // (#1207)
    return item >= 0 && item < (int)w.size() && w[item];
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
            k = std::min(count, tables.cap(item) - m.fuel.count);
            if (k > 0) { m.fuel.item = item; m.fuel.count += k; }
        }
        count -= k;
        if (count <= 0 || toFuel || !dual) return count;
    }
    if (toFuel || !Accepts(m.type, item)) return count;
    const int left = AddItems(m.in.data(), InSlots(m), item, count);
    if (left < count && m.fire == kFireOut) m.fire = kFireAuto; // (#1103) something to cook: a put-out campfire lights by itself again
    return left;
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
            m.dirty = true; // (#1214) room in its fuel slot: its conveyor network may move more
            if (gCount) gCount->items[f].used.add(1); // (#1199) fuel burnt is used
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
// it cools. Returns the EU made. (#1205) A fuel generator keeps the EU of a burnt tick that nothing drew (m.fuelEU, at most
// one tick's output) and serves the next demand from it first: fuel burns in proportion to the EU drawn (a whole tick per
// tick under a 4 EU/t load wasted 88% of it).
long long GeneratorStep(Machine& m, const MachineDef& def, int dt, long long needEU) {
    const bool steam = def.waterPerTick > 0;
    long long fromBuffer = 0;
    if (!steam) { fromBuffer = std::min<long long>(std::max(0, m.fuelEU), std::max(0LL, needEU)); m.fuelEU -= (int)fromBuffer; needEU -= fromBuffer; }
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
        const long long made = (long long)burnt * def.euOut;
        const long long drawn = std::min(made, std::max(0LL, needEU));
        m.fuelEU = (int)std::min<long long>(def.euOut, m.fuelEU + (made - drawn)); // (#1205) what was not drawn waits for the next demand
        eu = fromBuffer + drawn;
    }
    const bool noWater = steam && m.water < def.waterPerTick;
    const bool noFuel = m.energy <= 0 && m.fuel.count <= 0;
    if (steam && (noWater || noFuel)) m.heat = std::max(0, m.heat - dt);             // the fire is out: it cools
    m.status = noWater ? MachineStatus::NoWater : (noFuel && eu == 0) ? MachineStatus::NoFuel
             : eu > 0 || (steam && burnt > 0 && !banked) ? MachineStatus::Working : MachineStatus::NoLoad;
    return eu;
}

void BuildNets(PlotSim& ps, const LayoutMsg& L) {
    ps.itemNets.clear();
    ps.powerNets.clear();
    ps.cableNet.clear();
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
                if (power) ps.cableNet[c] = (int)nets.size(); // this network's id (pushed below)
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
    if (!ps) { plots.push_back({}); ps = &plots.back(); ps->body = L.body; ps->plot = L.plot; plotIndexDirty = true; }
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
        m.dirty = true; m.asleep = false; // (#1214) a new layout: everything looks again
        m.in.resize(InSlots(m)); // sized by type (a new machine starts empty)
        m.variant = i < L.machineVariant.size() && !L.machineVariant[i].empty() ? tables.item(L.machineVariant[i]) : -1;
        // A loaded save's contents arrive before the layout: apply them once.
        auto pend = pendingStates.find({ L.body, L.plot, part });
        if (pend != pendingStates.end()) {
            const auto& s = pend->second; // recipe progress energy fuelEU, then "in" item count pairs, then "out" pairs
            size_t k = 0;
            auto num = [&]() { return k < s.size() ? atoi(s[k++].c_str()) : 0; };
            m.recipe = num(); m.progress = num(); m.energy = num(); m.fuelEU = num();
            const int nin = num();
            for (int j = 0; j < nin && k + 1 < s.size(); ++j) { const int item = tables.item(s[k]); const int c = atoi(s[k + 1].c_str()); k += 2; AddItems(m.in.data(), InSlots(m), item, c); }
            const int nout = num();
            for (int j = 0; j < nout && k + 1 < s.size(); ++j) { const int item = tables.item(s[k]); const int c = atoi(s[k + 1].c_str()); k += 2; AddItems(m.out, kSlots, item, c); }
            // Then tagged extras: "f item count" (the fuel slot), "w litres", "h heat".
            while (k < s.size()) {
                const std::string tag = s[k++];
                if (tag == "f" && k + 1 < s.size()) { m.fuel.item = tables.item(s[k]); m.fuel.count = atoi(s[k + 1].c_str()); k += 2; }
                else if (tag == "w") m.water = num();
                else if (tag == "h") m.heat = num();
                else if (tag == "fire") m.fire = (uint8_t)std::clamp(num(), 0, 2); // (#1103)
            }
            // Saves from before the fuel slot (#194) kept fuel in the inputs: move fuel no recipe of this machine wants
            // into the fuel slot (#229 review), so an old boiler or furnace still burns it.
            if (HasFuelSlot(m.type))
                for (int j = 0; j < InSlots(m); ++j) {
                    Stack& st = m.in[j];
                    if (st.count <= 0 || !IsFuel(Def(m.type), st.item) || RecipeWants(m.type, st.item)) continue;
                    if (m.fuel.count > 0 && m.fuel.item != st.item) continue;
                    const int k = std::min(st.count, tables.cap(st.item) - m.fuel.count);
                    m.fuel.item = st.item; m.fuel.count += k; st.count -= k;
                    if (!st.count) st.item = -1;
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
        for (int i = 0; i < InSlots(old); ++i) if (old.in[i].count > 0) refunds.push_back({ old.in[i].item, old.in[i].count, L.refundTo });
        for (int i = 0; i < kSlots; ++i) if (old.out[i].count > 0) refunds.push_back({ old.out[i].item, old.out[i].count, L.refundTo });
        if (old.fuel.count > 0) refunds.push_back({ old.fuel.item, old.fuel.count, L.refundTo });
    }
    ps->machines = std::move(next);
    ps->byPart = std::move(byPart);
    BuildNets(*ps, L);
    ps->rot.Set(L.rotation); // (#1269)
}

// (#1269) Rotation before power: rotation machines with a recipe load their network, then one solve per plot.
void RotationStep(PlotSim& ps, int dt) {
    if (ps.rot.empty()) return;
    for (const Machine& m : ps.machines) if (Def(m.type).power == Power::Rotation) ps.rot.SetRunning(m.part, m.recipe >= 0);
    ps.rot.Step(dt);
}

// One step of dt ticks for a plot.
void TickPlot(PlotSim& ps, int dt, long long& made, int& active) {
    gCount = &ps.stats; // (#1199) BurnFuel counts into this plot
    RotationStep(ps, dt); // (#1269)
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
            if (def.power != Power::Generator) continue;
            const long long eu = GeneratorStep(m, def, dt, std::max(0LL, demand + batteryRoom - supply));
            supply += eu;
            if (eu > 0) ps.stats.euMade[(int)m.type].add(eu); // (#1199)
        }
        ps.stats.euDemand.add(demand);
        // Share: generation first, then the batteries; a shortfall slows every consumer equally.
        const long long fromBatteries = std::max(0LL, std::min(batteryAvail, demand - supply));
        const long long avail = supply + fromBatteries;
        for (int mi : net.machines) {
            Machine& m = ps.machines[mi];
            const MachineDef& def = Def(m.type);
            if (def.power == Power::Eu && m.recipe >= 0 && def.volts >= net.voltage) {
                const long long need = (long long)tables.recipes[m.recipe].eu * dt;
                m.received = (int)(demand <= avail ? need : need * avail / std::max(1LL, demand));
                if (m.received > 0) ps.stats.euUsed[(int)m.type].add(m.received); // (#1199)
            }
        }
        net.supply = supply + fromBatteries; net.demand = demand; net.ticks = std::max(1, dt); // for the cable tooltip (#570)
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
        // (#1214) asleep (nothing could start, nothing changed since): the same answer as a search, without one
        if (m.dirty) m.asleep = false;
        if (m.asleep && sleepOn) { m.status = MachineStatus::Idle; continue; }
        // Inputs are only checked when a craft starts and taken when it finishes (#193): take an input out
        // meanwhile and the craft is cancelled, its progress lost. Fuel has its own slot (#194).
        auto inputsHere = [&](const Recipe& r) {
            for (const Stack& s : r.in) if (CountItems(m.in.data(), kSlots, s.item) < s.count) return false;
            return true;
        };
        if (m.recipe >= 0 && !inputsHere(tables.recipes[m.recipe])) { m.recipe = -1; m.progress = 0; }
        if (m.fire == kFireOut) { m.status = MachineStatus::Idle; continue; } // (#1103) put out: no burn, no cooking
        if (m.recipe < 0) {
            auto usable = [&](const Recipe& r) {
                return r.machine == m.type && (m.variant < 0 || (!r.out.empty() && r.out[0].item == m.variant)) && inputsHere(r);
            };
            // Left to right as the slots are shown (#195): the first stack that some recipe can use picks it; the file
            // order chose before, so a stack further right could go first.
            // (#1207) only this machine type's recipes (file order kept, so the same one is picked)
            const std::vector<int>& mine = recipesOf[(int)m.type];
            for (int si = 0; si < InSlots(m) && m.recipe < 0; ++si) {
                if (m.in[si].count <= 0 || !RecipeWants(m.type, m.in[si].item)) continue;
                for (int ri : mine) {
                    const Recipe& r = tables.recipes[ri];
                    bool uses = false;
                    for (const Stack& s : r.in) uses |= s.item == m.in[si].item;
                    if (uses && usable(r)) { m.recipe = ri; m.progress = 0; break; }
                }
            }
            for (size_t k = 0; k < mine.size() && m.recipe < 0; ++k)
                if (usable(tables.recipes[mine[k]])) { m.recipe = mine[k]; m.progress = 0; }
            if (m.recipe < 0 && m.fire == kFireLit) { // (#1103) lit with nothing to cook: it burns its fuel (light, warmth)
                const int burnt = BurnFuel(m, def, dt);
                if (burnt > 0) { m.status = MachineStatus::Working; ++active; }
                else { m.status = MachineStatus::NoFuel; m.fire = kFireAuto; } // burnt out: it lights again by itself
                continue;
            }
            if (m.recipe < 0) { m.status = MachineStatus::Idle; m.asleep = true; continue; } // asleep until items arrive (#1214: really now)
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
            } else if (def.power == Power::Rotation) { // (#1269) at the shaft's speed: permille of a tick per tick, carried over
                const int rate = ps.rot.ProgressPermille(m.part), acc = m.rotAcc + rate * dt;
                ticks = acc / 1000; m.rotAcc = acc % 1000;
                if (rate > 0 && ticks == 0) { m.status = MachineStatus::Working; ++active; continue; } // turning, under a tick yet
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
            for (const Stack& s : r.in) TakeItems(m.in.data(), kSlots, s.item, s.count); // consumed now, at the end (#193)
            std::copy(trial, trial + kSlots, m.out);
            for (const Stack& s : r.out) { made += s.count; ps.stats.items[s.item].made.add(s.count); } // (#1199)
            for (const Stack& s : r.in) ps.stats.items[s.item].used.add(s.count);
            m.recipe = -1; m.progress = 0;
            m.dirty = true; // (#1214) outputs to move, inputs taken
        }
    }
    // (#1199) Time per status by machine kind: every machine, every tick of this step.
    for (const Machine& m : ps.machines) ps.stats.statusTicks[(int)m.type][(int)m.status & 7] += dt;
    // (#1214) A machine whose contents changed wakes its conveyor network (the change is seen; the flag is spent).
    for (Machine& m : ps.machines) {
        if (!m.dirty) continue;
        if (m.itemNet >= 0 && m.itemNet < (int)ps.itemNets.size()) ps.itemNets[m.itemNet].asleep = false;
        m.dirty = false;
    }
    // Items: each conveyor network moves items from outputs (and chests) to machines that take them.
    for (Net& net : ps.itemNets) {
        net.credit = std::min(net.credit + kConveyorItemsPerSecond * dt, kConveyorItemsPerSecond * kTicksPerSecond * 2);
        const int n = (int)net.machines.size();
        // (#1214) asleep: its last pass moved nothing and nothing changed since, so a pass would move nothing again; it banks
        // its credit and steps its round-robin start as that pass did (the order of later moves stays the same).
        if (net.asleep && sleepOn) { if (net.credit >= kTicksPerSecond && n > 1) net.rr = (net.rr + 1) % n; continue; }
        while (net.credit >= kTicksPerSecond && n > 1) {
            bool moved = false;
            for (int s = 0; s < n && !moved; ++s) {
                Machine& src = ps.machines[net.machines[(net.rr + s) % n]];
                const bool chest = src.type == PartType::Chest;
                Stack* from = chest ? src.in.data() : src.out;
                const int fromN = chest ? kChestSlots : kSlots;
                for (int i = 0; i < fromN && !moved; ++i) {
                    if (from[i].count <= 0) continue;
                    const int item = from[i].item;
                    for (int d = 1; d < n && !moved; ++d) {
                        Machine& dst = ps.machines[net.machines[(net.rr + s + d) % n]];
                        if (&dst == &src || (chest && dst.type == PartType::Chest) || !Accepts(dst.type, item)) continue;
                        if (dst.type == PartType::Chest && !chest) { if (AddItems(dst.in.data(), kChestSlots, item, 1) == 0) moved = true; }
                        else if (dst.type != PartType::Chest && CountItems(dst.in.data(), kSlots, item) + (dst.fuel.item == item ? dst.fuel.count : 0) < kFeedMax && // a buffer, fuel too (#566)
                                 !(item == WaterBucketItem() && Def(dst.type).waterPerTick > 0 && !RoomFor(dst.out, kSlots, BucketItem()))) { // the empty bucket needs an output slot (#229 review)
                            int buckets = 0; // fuel to the fuel slot, water buckets into a boiler's tank (its empty bucket to its outputs)
                            if (InsertInto(dst, item, 1, false, buckets) == 0) { moved = true; if (buckets > 0) AddItems(dst.out, kSlots, BucketItem(), buckets); }
                        }
                        if (moved) { TakeItems(from, fromN, item, 1); dst.dirty = true; } // (#1214) new input: it looks again next tick
                    }
                }
            }
            net.rr = (net.rr + 1) % std::max(1, n);
            if (!moved) { net.asleep = true; break; } // (#1214) until one of its machines changes
            net.credit -= kTicksPerSecond;
        }
    }
}

// (#997, user decision 2026-11-21) Parallel plots: plots share no conveyors, cables or machines, so the plots due this tick
// tick on up to kMaxWorkers threads (the tick thread is one of them), each writing only its own plots; their counts are
// summed in plot order after the join, so the result is the same as one thread's (unit test: 1 vs 4 workers). Default:
// one worker per 4 hardware threads, 1 to 4 (the bench won: 4 workers -35 to -40% of the tick, no frame cost seen; it
// leaves most cores to chunk building); RAYLIBGAME_FACTORY_WORKERS=N overrides. A fork/join costs ~40-60 us (measured),
// so a tick with little to do stays on the tick thread.
constexpr int kMaxWorkers = 4;
constexpr int kHardwareThreadsPerWorker = 4;
constexpr int kParallelMinMachines = 800; // machines due this tick below which one thread is quicker (guess: to tune)
int parallelMinMachines = kParallelMinMachines; // (the unit test lowers it)
int workerCount = [] {
    if (const char* e = getenv("RAYLIBGAME_FACTORY_WORKERS")) return std::clamp(atoi(e), 1, kMaxWorkers);
    return std::clamp((int)std::thread::hardware_concurrency() / kHardwareThreadsPerWorker, 1, kMaxWorkers);
}();

struct PlotJob { PlotSim* ps; int dt; long long made = 0; int active = 0; };
class WorkerPool {
public:
    ~WorkerPool() { stopAll(); }
    // Runs jobs[b] for b in buckets of `parts` (begin index per part, the last part ends at jobs.size()); part 0 here.
    void run(std::vector<PlotJob>& jobs, const std::vector<size_t>& begins) {
        const int parts = (int)begins.size();
        ensure(parts - 1);
        {
            std::lock_guard<std::mutex> lk(m);
            work = &jobs; starts = &begins; left = parts - 1; ++gen;
        }
        go.notify_all();
        runPart(0);
        std::unique_lock<std::mutex> lk(m);
        done.wait(lk, [&] { return left == 0; });
        work = nullptr;
    }
    void stopAll() {
        { std::lock_guard<std::mutex> lk(m); stop = true; }
        go.notify_all();
        for (auto& t : threads) if (t.joinable()) t.join();
        threads.clear();
        stop = false;
    }
private:
    void ensure(int n) {
        while ((int)threads.size() < n) {
            const int index = (int)threads.size() + 1; // part index this thread runs
            threads.emplace_back([this, index] {
                Log_SetThreadName("factory-w");
                long long seen = 0;
                for (;;) {
                    {
                        std::unique_lock<std::mutex> lk(m);
                        go.wait(lk, [&] { return stop || gen != seen; });
                        if (stop) return;
                        seen = gen;
                        if (!work || index >= (int)starts->size()) continue; // (fewer parts this tick)
                    }
                    runPart(index);
                    std::lock_guard<std::mutex> lk(m);
                    if (--left == 0) done.notify_one();
                }
            });
        }
    }
    void runPart(int p) {
        const size_t b = (*starts)[p], e = p + 1 < (int)starts->size() ? (*starts)[p + 1] : work->size();
        for (size_t i = b; i < e; ++i) { PlotJob& j = (*work)[i]; TickPlot(*j.ps, j.dt, j.made, j.active); }
        gCount = nullptr;
    }
    std::vector<std::thread> threads;
    std::mutex m;
    std::condition_variable go, done;
    long long gen = 0;
    int left = 0;
    bool stop = false;
    std::vector<PlotJob>* work = nullptr;
    const std::vector<size_t>* starts = nullptr;
};
WorkerPool workerPool;

// Ticks the jobs (each a plot and its step) on `workers` threads, splitting them into runs of about equal machine counts
// in their order; then adds their counts up in that order.
void TickPlots(std::vector<PlotJob>& jobs, int workers, long long& made, int& active) {
    long long machines = 0;
    for (const PlotJob& j : jobs) machines += (long long)j.ps->machines.size();
    const int parts = (int)std::min<long long>({ (long long)workers, (long long)jobs.size(), std::max(1LL, machines / std::max(1, parallelMinMachines / 2)) });
    if (parts <= 1 || machines < parallelMinMachines) {
        for (PlotJob& j : jobs) TickPlot(*j.ps, j.dt, j.made, j.active);
        gCount = nullptr;
    } else {
        std::vector<size_t> begins{ 0 };
        long long acc = 0;
        for (size_t i = 0; i < jobs.size() && (int)begins.size() < parts; ++i) {
            acc += (long long)jobs[i].ps->machines.size();
            if (acc * parts >= machines * (long long)begins.size() && i + 1 < jobs.size()) begins.push_back(i + 1);
        }
        workerPool.run(jobs, begins);
    }
    for (const PlotJob& j : jobs) { made += j.made; active += j.active; }
}

// The watched cable's network (#570), under viewMutex.
void FillCableView() {
    const auto& [b, pl, x, y, z] = cableKey;
    cableView = CableView{};
    const PlotSim* p = FindPlot(b, pl);
    if (!p) return;
    cableView.valid = true;
    auto it = p->cableNet.find(Int3{ x, y, z });
    if (it == p->cableNet.end() || it->second < 0 || it->second >= (int)p->powerNets.size()) return;
    const Net& net = p->powerNets[it->second];
    cableView.voltage = net.voltage;
    for (int mi : net.machines) {
        const Power pw = Def(p->machines[mi].type).power;
        cableView.generators += pw == Power::Generator; cableView.batteries += pw == Power::Battery; cableView.consumers += pw == Power::Eu;
    }
    cableView.supply = (float)net.supply / (float)net.ticks;
    cableView.demand = (float)net.demand / (float)net.ticks;
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
    v.hotOutputs = HeatsItems(m.type);
    v.canLight = m.type == PartType::Campfire;
    v.fire = m.fire;
    if (m.fuel.count > 0) v.fuelSlot = { tables.name(m.fuel.item), m.fuel.count };
    if (Def(m.type).waterPerTick > 0) { v.water = m.water; v.waterMax = kWaterTank; v.heat = m.heat; v.heatMax = kHeatMax; }
    auto list = [](const Stack* s, int n) { // totals per item in slot order (#567: a chest of coal is "coal x3200", not 50 entries)
        std::vector<std::pair<int, int>> total;
        for (int i = 0; i < n; ++i) {
            if (s[i].count <= 0) continue;
            auto it = std::find_if(total.begin(), total.end(), [&](const auto& t) { return t.first == s[i].item; });
            if (it == total.end()) total.push_back({ s[i].item, s[i].count }); else it->second += s[i].count;
        }
        std::string out;
        for (const auto& t : total) { if (!out.empty()) out += ", "; out += tables.name(t.first) + " x" + std::to_string(t.second); }
        return out.empty() ? std::string("-") : out;
    };
    v.input = list(m.in.data(), InSlots(m));
    v.output = list(m.out, kSlots);
    for (int i = 0; i < InSlots(m); ++i) v.inSlots.push_back({ m.in[i].count > 0 ? tables.name(m.in[i].item) : std::string(), m.in[i].count });
    if (m.type != PartType::Chest) for (int i = 0; i < kSlots; ++i) v.outSlots.push_back({ m.out[i].count > 0 ? tables.name(m.out[i].item) : std::string(), m.out[i].count });
    if (Def(m.type).power == Power::Rotation) { // (#1269)
        v.rotation = true;
        v.rpm = (float)std::abs(ps.rot.SpeedMilli(part)) / 1000.0f;
        v.rotStatus = Rotation::StatusName(ps.rot.StatusOf(part));
        if (const Rotation::NetResult* n = ps.rot.NetOf(part)) { v.speedPercent = n->ratio / 10; v.stressPercent = Rotation::StressPercent(*n); }
    }
}

// (#1269) Each near plot's rotation parts as of this tick (viewMutex held); far plots, and plots without any, publish none.
void PublishRotation(bool reset) {
    if (reset) rotView.clear();
    for (const PlotSim& ps : plots) {
        const std::pair<uint64_t, uint32_t> key{ ps.body, ps.plot };
        if (!ps.nearPlayer || ps.rot.empty()) { rotView.erase(key); continue; }
        std::vector<RotPub>& out = rotView[key];
        out.clear();
        for (const RotationPartMsg& rp : ps.rot.parts) {
            RotPub p;
            p.part = rp.part; p.milli = ps.rot.SpeedMilli(rp.part); p.status = ps.rot.StatusOf(rp.part);
            if (const Rotation::NetResult* n = ps.rot.NetOf(rp.part)) p.net = *n;
            out.push_back(p);
        }
    }
}

void SimLoop() {
    Log_SetThreadName("factory");
    using clk = std::chrono::steady_clock;
    double emaMs = 0.0, maxMs = 0.0;
    while (running) {
        {   // One tick per credit from the sim clock (#315).
            std::unique_lock<std::mutex> lk(creditMutex);
            creditCv.wait_for(lk, std::chrono::milliseconds(100), [] { return tickCredit.load() > 0 || !running; });
            if (!running) break;
            if (tickCredit.load() <= 0) continue;
            --tickCredit;
        }
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
            if (reset) { plots.clear(); plotIndexDirty = true; emaMs = maxMs = 0.0; statSeconds = 0; tickRingAt = 0; }
            // Moved plots first, so the new plot's first layout finds its machines (same part indices).
            for (const MoveMsg& mv : moves)
                if (PlotSim* p = FindPlot(mv.fromBody, mv.fromPlot); p && !FindPlot(mv.toBody, mv.toPlot)) { p->body = mv.toBody; p->plot = mv.toPlot; plotIndexDirty = true; viewsRebuild = true; }
            for (const LayoutMsg& L : layouts) ApplyLayout(L);
            for (const auto& n : nears) if (PlotSim* p = FindPlot(n.first, n.second.first)) p->nearPlayer = n.second.second;
            std::vector<ItemOut> toCursor; // kind 8 (#192)
            std::vector<ItemOut> back;
            back.swap(refunds); // removed machines' contents (from this tick's layouts)
            std::vector<std::pair<std::tuple<uint64_t, uint32_t, int>, int>> drained; // pump drains applied this tick (pumpPending)
            std::vector<std::pair<std::tuple<uint64_t, uint32_t, int>, int>> filled;  // tank fills applied this tick (tankPending)
            for (const ItemMsg& im : items) {
                PlotSim* p = FindPlot(im.body, im.plot);
                if (im.kind == 11) { if (p) p->rot.Hold(im.part); continue; } // (#1269) a player turns a hand crank
                if (im.kind == 12) { if (p) p->rot.SetSourceRpm(im.part, im.count); continue; } // (#1270) a water wheel's or windmill's sampled speed
                auto it = p ? p->byPart.find(im.part) : decltype(p->byPart.end()){};
                if (!p || it == p->byPart.end()) {
                    if (im.kind == 0 || im.kind == 9) back.push_back({ im.item, im.count, im.to });
                    if (im.kind == 5) drained.push_back({ std::make_tuple(im.body, im.plot, im.part), im.count });
                    if (im.kind == 7) filled.push_back({ std::make_tuple(im.body, im.plot, im.part), im.count });
                    continue;
                }
                Machine& m = p->machines[it->second];
                m.dirty = true; // (#1214) any command on a machine wakes it and its network
                if (im.kind == 9 && m.type == PartType::Chest) { // (#862) into one storage slot: empty, or the same item up to its cap
                    int left = im.count;
                    if (im.amount >= 0 && im.amount < kChestSlots && im.item >= 0) {
                        Stack& s = m.in[im.amount];
                        if (s.count <= 0 || s.item == im.item) {
                            const int k = std::min(left, tables.cap(im.item) - std::max(0, s.count));
                            if (k > 0) { if (s.count <= 0) { s = Stack{}; s.item = im.item; } s.count += k; left -= k; }
                        }
                    }
                    if (im.count - left > 0) p->stats.items[im.item].handIn.add(im.count - left); // (#1199)
                    if (left > 0) back.push_back({ im.item, left, im.to });
                } else if (im.kind == 0 || im.kind == 2 || im.kind == 6 || im.kind == 9) { // (9 on a machine: routed, as Insert)
                    int buckets = 0;
                    const int left = InsertInto(m, im.item, im.count, im.kind == 6, buckets);
                    if (im.count - left > 0) p->stats.items[im.item].handIn.add(im.count - left); // (#1199)
                    if (left > 0 && im.kind != 2) back.push_back({ im.item, left, im.to });
                    if (buckets > 0 && im.kind != 2) back.push_back({ BucketItem(), buckets, im.to }); // the emptied buckets
                } else if (im.kind == 10) { // (#1103) light or put out a campfire
                    if (m.type == PartType::Campfire) m.fire = im.count ? kFireLit : kFireOut;
                } else if (im.kind == 7) {
                    if (Def(m.type).waterPerTick > 0) m.water = std::min(kWaterTank, m.water + im.count); // Waterworks sends only the room (TankRoom)
                    filled.push_back({ std::make_tuple(im.body, im.plot, im.part), im.count });
                } else if (im.kind == 3) {
                    if (m.type == PartType::Chest) { const int k = std::min(im.count, CountItems(m.in.data(), kChestSlots, im.item)); TakeItems(m.in.data(), kChestSlots, im.item, k); if (k > 0) { back.push_back({ im.item, k, im.to }); p->stats.items[im.item].handOut.add(k); } }
                } else if (im.kind == 5) {
                    TakeItems(m.out, kSlots, im.item, im.count); // a pump's water leaves through its pipes (nothing to the inventory)
                    drained.push_back({ std::make_tuple(im.body, im.plot, im.part), im.count });
                } else if (im.kind == 4) {
                    Stack* s = im.count == 2 ? &m.fuel : im.count ? m.out : m.in.data();
                    const int n = im.count == 2 ? 1 : im.count ? kSlots : InSlots(m);
                    const int hot = s == m.out && HeatsItems(m.type) ? kHotItem : 0; // fresh out of a furnace (#580)
                    if (im.item >= 0 && im.item < n && s[im.item].count > 0) { p->stats.items[s[im.item].item].handOut.add(s[im.item].count); back.push_back({ s[im.item].item | hot, s[im.item].count, im.to }); s[im.item] = Stack{}; }
                } else if (im.kind == 8) { // part of a slot onto the cursor (#192)
                    Stack* s = im.count == 2 ? &m.fuel : im.count ? m.out : m.in.data();
                    const int n = im.count == 2 ? 1 : im.count ? kSlots : InSlots(m);
                    const int hot = s == m.out && HeatsItems(m.type) ? kHotItem : 0;
                    if (im.item >= 0 && im.item < n && s[im.item].count > 0) {
                        const int k = im.amount > 0 ? std::min(im.amount, s[im.item].count) : s[im.item].count;
                        toCursor.push_back({ s[im.item].item | hot, k, im.to });
                        p->stats.items[s[im.item].item].handOut.add(k); // (#1199)
                        s[im.item].count -= k;
                        if (s[im.item].count <= 0) s[im.item] = Stack{};
                    }
                } else {
                    Stack* s = m.type == PartType::Chest ? m.in.data() : m.out;
                    const int n = m.type == PartType::Chest ? kChestSlots : kSlots;
                    const int hot = s == m.out && HeatsItems(m.type) ? kHotItem : 0; // #580
                    for (int i = 0; i < n; ++i) if (s[i].count > 0) { p->stats.items[s[i].item].handOut.add(s[i].count); back.push_back({ s[i].item | hot, s[i].count, im.to }); s[i] = Stack{}; }
                }
            }
            const auto t1 = clk::now(); // (#1207) commands done
            // Tick: near plots every tick; far plots once a second with a 20x step.
            const long long t = tickCount.load();
            long long made = 0;
            int active = 0, machines = 0, inets = 0, pnets = 0;
            static std::vector<PlotJob> jobs; // (#997) the plots due this tick, in plot order (parallel when it pays)
            jobs.clear();
            for (size_t i = 0; i < plots.size(); ++i) {
                PlotSim& ps = plots[i];
                machines += (int)ps.machines.size(); inets += (int)ps.itemNets.size(); pnets += (int)ps.powerNets.size();
                if (ps.nearPlayer) jobs.push_back({ &ps, 1 });
                else if ((t + (long long)i) % kTicksPerSecond == 0) jobs.push_back({ &ps, kTicksPerSecond });
            }
            TickPlots(jobs, workerCount, made, active);
            gCount = nullptr;
            if ((t + 1) % kTicksPerSecond == 0) { for (PlotSim& ps : plots) RollStats(ps.stats); ++statSeconds; } // (#1199) a second of factory time
            const auto t2 = clk::now(); // (#1207) plots done
            const double ms = std::chrono::duration<double, std::milli>(clk::now() - t0).count();
            emaMs = emaMs == 0.0 ? ms : emaMs * 0.95 + ms * 0.05;
            maxMs = std::max(maxMs, ms);
            {
                std::lock_guard<std::mutex> v(viewMutex);
                stats.machines = machines; stats.activeMachines = active; stats.itemNetworks = inets; stats.powerNetworks = pnets;
                stats.avgTickMs = emaMs; stats.maxTickMs = maxMs; stats.ticks = t + 1; stats.itemsMade += made;
                const auto& [b, pl, part] = watchKey;
                if (PlotSim* p = FindPlot(b, pl)) FillView(*p, part, watchView); else watchView = MachineView{};
                for (const auto& k : watchMany) { const auto& [wb, wp, wpart] = k; if (PlotSim* p = FindPlot(wb, wp)) FillView(*p, wpart, watchManyViews[k]); else watchManyViews[k] = MachineView{}; } // (#332)
                FillCableView();
                for (auto& r : back) results.push_back(r);
                for (auto& r : toCursor) cursorResults.push_back(r);
                pumpView.clear();
                for (const auto& [k, n] : drained) if (auto pit = pumpPending.find(k); pit != pumpPending.end() && (pit->second -= n) <= 0) pumpPending.erase(pit);
                if (reset) pumpPending.clear();
                tankView.clear();
                for (const auto& [k, n] : filled) if (auto tit = tankPending.find(k); tit != tankPending.end() && (tit->second -= n) <= 0) tankPending.erase(tit);
                if (reset) tankPending.clear();
                // Boilers' tank room and pumps' water (every plot; few machines): one pass.
                for (const PlotSim& ps : plots)
                    for (const Machine& m : ps.machines) {
                        if (Def(m.type).waterPerTick > 0) tankView[std::make_tuple(ps.body, ps.plot, m.part)] = kWaterTank - m.water;
                        else if (m.type == PartType::WaterPump) { int w = 0; for (int i = 0; i < kSlots; ++i) if (m.out[i].count > 0) w += m.out[i].count; pumpView[std::make_tuple(ps.body, ps.plot, m.part)] = w; }
                    }
                PublishRotation(reset || viewsRebuild); // (#1269)
                // (#1207) Machines working and chest contents on the plots near the player (models, smoke, chest screens): each
                // plot republishes only what changed since it last did (they were rebuilt for every machine every tick: ~60%
                // of a 5000-machine tick). A plot leaving "near" takes its entries away; a reset or a moved plot rebuilds all.
                if (reset || viewsRebuild) {
                    workingView.clear(); workingByPlot.clear(); chestView.clear();
                    for (PlotSim& ps : plots) ps.pub = PlotSim::Published{};
                    viewsRebuild = false;
                }
                static std::vector<int> working, chests;
                for (PlotSim& ps : plots) {
                    PlotSim::Published& pub = ps.pub;
                    auto unpublishWorking = [&] { for (int part : pub.working) workingView.erase(WorkingKey(ps.body, ps.plot, part)); workingByPlot.erase({ ps.body, ps.plot }); pub.working.clear(); };
                    auto unpublishChests = [&] { for (int part : pub.chests) chestView.erase(std::make_tuple(ps.body, ps.plot, part)); pub.chests.clear(); };
                    if (!ps.nearPlayer) {
                        if (pub.near) { unpublishWorking(); unpublishChests(); pub.near = false; }
                        continue;
                    }
                    working.clear(); chests.clear();
                    uint64_t h = 1469598103934665603ull; // FNV-1a over the chests' parts and slots
                    auto mix = [&h](uint64_t v) { h = (h ^ v) * 1099511628211ull; };
                    for (const Machine& m : ps.machines) {
                        if (m.status == MachineStatus::Working) working.push_back(m.part);
                        if (m.type != PartType::Chest) continue;
                        chests.push_back(m.part);
                        mix((uint64_t)m.part);
                        for (int i = 0; i < kChestSlots; ++i) if (m.in[i].count > 0) mix(((uint64_t)(uint32_t)m.in[i].item << 32) ^ (uint64_t)(uint32_t)(m.in[i].count * 64 + i));
                    }
                    if (!pub.near || working != pub.working) {
                        unpublishWorking();
                        for (int part : working) workingView.insert(WorkingKey(ps.body, ps.plot, part));
                        if (!working.empty()) workingByPlot[{ ps.body, ps.plot }] = working;
                        pub.working = working;
                    }
                    if (!pub.near || h != pub.chestHash || chests != pub.chests) {
                        unpublishChests();
                        for (const Machine& m : ps.machines) {
                            if (m.type != PartType::Chest) continue;
                            auto& v = chestView[std::make_tuple(ps.body, ps.plot, m.part)];
                            for (int i = 0; i < kChestSlots; ++i) if (m.in[i].count > 0) v.push_back({ m.in[i].item, m.in[i].count });
                        }
                        pub.chests = chests;
                        pub.chestHash = h;
                    }
                    pub.near = true;
                }
            }
            // (#1207) The whole tick (commands, plots, views) into the last minute's ring; its percentiles once a second.
            const auto t3 = clk::now();
            const float cmdMs = std::chrono::duration<float, std::milli>(t1 - t0).count();
            const float tickMs = std::chrono::duration<float, std::milli>(t2 - t1).count();
            const float viewMs = std::chrono::duration<float, std::milli>(t3 - t2).count();
            tickRing[tickRingAt++ % kTickRing] = { cmdMs + tickMs + viewMs, cmdMs, tickMs, viewMs };
            if ((t + 1) % kTicksPerSecond == 0) {
                const size_t n = std::min<size_t>(tickRingAt, kTickRing);
                std::vector<float> all(n);
                TickTimes sum{};
                for (size_t i = 0; i < n; ++i) { all[i] = tickRing[i].total; sum.cmd += tickRing[i].cmd; sum.tick += tickRing[i].tick; sum.view += tickRing[i].view; }
                std::sort(all.begin(), all.end());
                std::lock_guard<std::mutex> v(viewMutex);
                stats.recentTicks = (int)n;
                stats.meanTickMs = n ? (sum.cmd + sum.tick + sum.view) / n : 0.0;
                stats.p99TickMs = n ? all[std::min(n - 1, (size_t)(n * 0.99))] : 0.0;
                stats.worstRecentMs = n ? all[n - 1] : 0.0;
                stats.cmdMs = n ? sum.cmd / n : 0.0; stats.tickOnlyMs = n ? sum.tick / n : 0.0; stats.viewMs = n ? sum.view / n : 0.0;
            }
        }
        tickCount++;
        tickCv.notify_all();
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
    workerPool.stopAll(); // (#997)
}

void Factory::Reset() {
    std::lock_guard<std::mutex> q(queueMutex);
    resetRequested = true;
    layoutQueue.clear();
    itemQueue.clear();
    sentVersions.clear();
    RotationSources_Reset(); // (#1270)
}

void Factory::Sync(Scene& world) {
    if (!running) Start();
    static double lastNear = 0.0;
    const bool nearCheck = SimTime() - lastNear > 0.5 || SimTime() < lastNear;
    if (nearCheck) lastNear = SimTime();
    for (auto& obj : world.objects) {
        auto* body = dynamic_cast<Planetoid*>(obj.get());
        if (!body || !body->plots) continue;
        for (auto& pl : body->plots->plots) {
            const std::string key = PlotKey(body->id, pl->id);
            auto sv = sentVersions.find(key);
            if (sv == sentVersions.end() || sv->second != pl->version) {
                LayoutMsg L{ body->id, pl->id, pl->version, {}, {}, {}, {}, {} };
                L.refundTo = pl->refundTo; pl->refundTo = 0; // (#854)
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
                    if (RotationPartMsg rm; RotationPlot_Describe(part, (int)i, rm)) L.rotation.push_back(std::move(rm)); // (#1269)
                }
                std::lock_guard<std::mutex> q(queueMutex);
                layoutQueue.push_back(std::move(L));
                sentVersions[key] = pl->version;
            }
            if (nearCheck) {
                const Vector3 o = body->toWorld({ (float)pl->origin.x, (float)pl->origin.y, (float)pl->origin.z });
                bool near = Vector3Distance(o, world.camera.position) < 1000.0f;
                for (const Vector3& np : nearPoints) near = near || Vector3Distance(o, np) < 1000.0f; // (#332) a server: its players
                std::lock_guard<std::mutex> q(queueMutex);
                nearQueue.push_back({ body->id, { pl->id, near } });
            }
        }
    }
    RotationSources_Update(world); // (#1270) water wheels and windmills: their river or wind, twice a second
    // Items back from the simulation.
    std::vector<ItemOut> back;
    {
        std::lock_guard<std::mutex> v(viewMutex);
        back.swap(results);
    }
    std::vector<ItemOut> cur;
    {
        std::lock_guard<std::mutex> v(viewMutex);
        cur.swap(cursorResults);
    }
    // kHotItem marks what came out of a heating machine's output (#580): it arrives hot, in the inventory or on the cursor.
    for (auto& r : back) { const std::string id = tables.name(r.item & ~kHotItem); const float h = (r.item & kHotItem) ? 1.0f : 0.0f; if (r.to && gFactoryDeliver) gFactoryDeliver(r.to, id, r.count, h); else inventory.addAmount(id, r.count, h); } // (#332) a server: to that player
    for (auto& r : cur) { const std::string id = tables.name(r.item & ~kHotItem); const float h = (r.item & kHotItem) ? 1.0f : 0.0f; if (r.to && gFactoryDeliver) gFactoryDeliver(r.to, id, r.count, h); else if (cursorSink) cursorSink(id, r.count, h); else inventory.addAmount(id, r.count, h); }
}

void Factory::TakeSlotToCursor(uint64_t body, uint32_t plot, int part, int slot, int which, int count, uint32_t to) {
    std::lock_guard<std::mutex> q(queueMutex);
    ItemMsg m{ 8, body, plot, part, which == 2 ? 0 : slot, which };
    m.amount = count; m.to = to; // (#854) a server: that player's inventory (its client puts it on its cursor)
    itemQueue.push_back(m);
}
void Factory::InsertSlot(uint64_t body, uint32_t plot, int part, const std::string& item, int count, int slot, uint32_t to) {
    std::lock_guard<std::mutex> q(queueMutex);
    ItemMsg m{ 9, body, plot, part, tables.item(item), count };
    m.amount = slot; m.to = to;
    itemQueue.push_back(m);
}
void Factory::SetCursorSink(void (*sink)(const std::string&, int, float)) { cursorSink = sink; }

MachineView Factory::View(uint64_t body, uint32_t plot, int part) {
    std::lock_guard<std::mutex> v(viewMutex);
    const auto key = std::make_tuple(body, plot, part);
    if (remoteMode) { // (#332) a client: the server's view (asked for when the key changes)
        if (watchKey != key) { watchKey = key; if (gFactoryRemoteWatch) gFactoryRemoteWatch(body, plot, part); }
        auto it = remoteViews.find(key); return it == remoteViews.end() ? MachineView{} : it->second;
    }
    if (watchKey != key) { watchKey = key; return MachineView{}; } // filled from the next tick on
    return watchView;
}

CableView Factory::Cable(uint64_t body, uint32_t plot, Int3 cell) {
    std::lock_guard<std::mutex> v(viewMutex);
    const auto key = std::make_tuple(body, plot, cell.x, cell.y, cell.z);
    if (cableKey != key) { cableKey = key; return CableView{}; } // filled from the next tick on
    return cableView;
}

FactoryStats Factory::Stats() {
    std::lock_guard<std::mutex> v(viewMutex);
    return stats;
}

FactoryStatsView Factory::StatsView(StatWindow window, bool onePlot, uint64_t body, uint32_t plot) {
    FactoryStatsView out;
    std::lock_guard<std::mutex> lock(simMutex);
    // The window shown: a longer one with no finished bucket yet falls back to the shorter.
    const long long secs = statSeconds;
    if (window == StatWindow::TenHours && secs < 600) window = StatWindow::Hour;
    if (window == StatWindow::Hour && secs < 60) window = StatWindow::Minute;
    out.window = window;
    const int filled = (int)std::min<long long>(kStatBuckets, window == StatWindow::Minute ? secs : window == StatWindow::Hour ? secs / 60 : secs / 600);
    const double bucketMin = window == StatWindow::Minute ? 1.0 / 60.0 : window == StatWindow::Hour ? 1.0 : 10.0; // minutes per bucket
    out.coveredMinutes = filled * bucketMin;
    const long long newest = window == StatWindow::Minute ? secs - 1 : window == StatWindow::Hour ? secs / 60 - 1 : secs / 600 - 1;
    auto bucket = [&](const Series& s, size_t i) -> long long { return window == StatWindow::Minute ? s.sec[i] : window == StatWindow::Hour ? s.min[i] : s.ten[i]; };
    auto sum = [&](const Series& s) { long long n = 0; for (size_t i = 0; i < kStatBuckets; ++i) n += bucket(s, i); return n; };
    auto perMin = [&](const Series& s) { return out.coveredMinutes > 0 ? (double)sum(s) / out.coveredMinutes : 0.0; };
    auto perTick = [&](const Series& s) { return out.coveredMinutes > 0 ? (double)sum(s) / (out.coveredMinutes * 60.0 * kTicksPerSecond) : 0.0; };
    auto spark = [&](const Series& s, std::vector<float>& v) { // oldest first
        v.resize(filled);
        for (int i = 0; i < filled; ++i) v[i] = (float)bucket(s, (size_t)(((newest - (filled - 1) + i) % kStatBuckets + kStatBuckets) % kStatBuckets));
    };
    // Sum the plots in scope into one set of series (bucket by bucket).
    std::map<int, ItemSeries> items;
    std::map<int, Series> euMade, euUsed;
    Series demand;
    std::map<int, StatKind> kinds;
    auto addSeries = [](Series& into, const Series& s) {
        for (int i = 0; i < kStatBuckets; ++i) { into.sec[i] += s.sec[i]; into.min[i] += s.min[i]; into.ten[i] += s.ten[i]; }
        into.total += s.total + s.acc; // (live: this second's count too)
    };
    for (const PlotSim& ps : plots) {
        if (onePlot && (ps.body != body || ps.plot != plot)) continue;
        ++out.plots;
        for (const auto& [id, is] : ps.stats.items) {
            ItemSeries& t = items[id];
            addSeries(t.made, is.made); addSeries(t.used, is.used); addSeries(t.handIn, is.handIn); addSeries(t.handOut, is.handOut);
        }
        for (const auto& [k, s] : ps.stats.euMade) addSeries(euMade[k], s);
        for (const auto& [k, s] : ps.stats.euUsed) addSeries(euUsed[k], s);
        addSeries(demand, ps.stats.euDemand);
        for (const auto& [k, t] : ps.stats.statusTicks) { StatKind& sk = kinds[k]; sk.type = k; for (int i = 0; i < 8; ++i) sk.ticks[i] += t[i]; }
        for (const Machine& m : ps.machines) {
            StatKind& sk = kinds[(int)m.type]; sk.type = (int)m.type; ++sk.count; ++sk.now[(int)m.status & 7];
            if (m.type == PartType::Battery) { out.batteryEnergy += m.energy; out.batteryCapacity += kBatteryCapacity; }
        }
    }
    for (const auto& [id, is] : items) {
        StatLine l;
        l.item = tables.name(id);
        l.madePerMin = perMin(is.made); l.usedPerMin = perMin(is.used); l.handInPerMin = perMin(is.handIn); l.handOutPerMin = perMin(is.handOut);
        l.made = is.made.total; l.used = is.used.total; l.handIn = is.handIn.total; l.handOut = is.handOut.total;
        spark(is.made, l.madeSpark); spark(is.used, l.usedSpark);
        out.items.push_back(std::move(l));
    }
    std::sort(out.items.begin(), out.items.end(), [](const StatLine& a, const StatLine& b) { return a.item < b.item; });
    for (auto& [k, sk] : kinds) {
        if (auto it = euMade.find(k); it != euMade.end()) { sk.euMadePerTick = perTick(it->second); out.euMadePerTick += sk.euMadePerTick; }
        if (auto it = euUsed.find(k); it != euUsed.end()) { sk.euUsedPerTick = perTick(it->second); out.euUsedPerTick += sk.euUsedPerTick; }
        out.kinds.push_back(sk);
    }
    out.euDemandPerTick = perTick(demand);
    return out;
}

void Factory::Insert(uint64_t body, uint32_t plot, int part, const std::string& item, int count, uint32_t to) {
    std::lock_guard<std::mutex> q(queueMutex);
    ItemMsg m{ 0, body, plot, part, tables.item(item), count }; m.to = to; itemQueue.push_back(m);
}
void Factory::Extract(uint64_t body, uint32_t plot, int part, uint32_t to) {
    std::lock_guard<std::mutex> q(queueMutex);
    ItemMsg m{ 1, body, plot, part, -1, 0 }; m.to = to; itemQueue.push_back(m);
}
void Factory::TakeSlot(uint64_t body, uint32_t plot, int part, int slot, int which, uint32_t to) {
    std::lock_guard<std::mutex> q(queueMutex);
    ItemMsg m{ 4, body, plot, part, which == 2 ? 0 : slot, which }; m.to = to; itemQueue.push_back(m);
}
void Factory::InsertFuel(uint64_t body, uint32_t plot, int part, const std::string& item, int count, uint32_t to) {
    std::lock_guard<std::mutex> q(queueMutex);
    ItemMsg m{ 6, body, plot, part, tables.item(item), count }; m.to = to; itemQueue.push_back(m);
}
void Factory::SetFire(uint64_t body, uint32_t plot, int part, bool lit) {
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back(ItemMsg{ 10, body, plot, part, -1, lit ? 1 : 0 });
}
void Factory::Crank(uint64_t body, uint32_t plot, int part) { // (#1269)
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back(ItemMsg{ 11, body, plot, part, -1, 0 });
}
void Factory::SetSourceRpm(uint64_t body, uint32_t plot, int part, int rpm) { // (#1270)
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back(ItemMsg{ 12, body, plot, part, -1, rpm });
}
RotationView Factory::RotationAt(uint64_t body, uint32_t plot, int part) {
    RotationView v;
    std::lock_guard<std::mutex> l(viewMutex);
    auto it = rotView.find({ body, plot });
    if (it == rotView.end()) return v;
    for (const RotPub& p : it->second) {
        if (p.part != part) continue;
        v.valid = true;
        v.rpm = (float)p.milli / 1000.0f;
        v.status = Rotation::StatusName(p.status);
        v.jam = Rotation::JamName(p.net.jam);
        v.stressPercent = Rotation::StressPercent(p.net);
        v.speedPercent = p.net.ratio / 10;
        v.capacity = (int)p.net.capacity;
        v.demand = (int)(p.net.demandMilli / 1000);
        break;
    }
    return v;
}
void Factory::RotationSpeeds(uint64_t body, uint32_t plot, std::vector<std::pair<int, float>>& out) {
    out.clear();
    std::lock_guard<std::mutex> l(viewMutex);
    auto it = rotView.find({ body, plot });
    if (it == rotView.end()) return;
    for (const RotPub& p : it->second) out.push_back({ p.part, (float)p.milli / 1000.0f });
}
int Factory::TankRoom(uint64_t body, uint32_t plot, int part) {
    std::lock_guard<std::mutex> v(viewMutex);
    const auto key = std::make_tuple(body, plot, part);
    auto it = tankView.find(key);
    if (it == tankView.end()) return 0; // not simulated yet: no room until the sim has seen it
    auto pit = tankPending.find(key);
    return std::max(0, it->second - (pit == tankPending.end() ? 0 : pit->second)); // less what is already on its way in
}
void Factory::AddWater(uint64_t body, uint32_t plot, int part, int litres) {
    if (litres <= 0) return;
    { std::lock_guard<std::mutex> v(viewMutex); tankPending[std::make_tuple(body, plot, part)] += litres; } // TankRoom won't offer it twice
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back({ 7, body, plot, part, -1, litres });
}
bool Factory::TakesWater(PartType type) { LoadRecipes(); return Def(type).waterPerTick > 0; }
std::vector<Factory::RecipeInfo> Factory::RecipeList() {
    LoadRecipes();
    std::vector<RecipeInfo> out;
    for (const auto& r : tables.recipes) { // (the factory's own recipe type: "Recipe" alone is ambiguous with the hand recipes here)
        RecipeInfo ri;
        ri.machine = r.machine; ri.ticks = r.ticks; ri.eu = r.eu;
        for (const Stack& s : r.in) ri.in.push_back({ tables.name(s.item), s.count });
        for (const Stack& s : r.out) ri.out.push_back({ tables.name(s.item), s.count });
        out.push_back(std::move(ri));
    }
    return out;
}
std::vector<std::string> Factory::FuelsOf(PartType type) {
    LoadRecipes();
    std::vector<std::string> out;
    for (int f : Def(type).fuels) out.push_back(tables.name(f));
    return out;
}
void Factory::Seed(uint64_t body, uint32_t plot, int part, const std::string& item, int count) {
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back({ 2, body, plot, part, tables.item(item), count });
}

void Factory::Advance(float simDt) {
    if (remoteMode) return; // (#332) a client: the server simulates
    const int speed = std::max(1, simSpeed.load());
    creditAccum += (double)simDt * kTicksPerSecond * speed;
    const long long n = (long long)creditAccum;
    creditAccum -= (double)n;
    if (n > 0) GrantTicks(n, std::max(2, kTicksPerSecond * speed / 4)); // at most 250 ms behind
}

void Factory::WaitTicks(int n) {
    const long long target = tickCount.load() + n;
    GrantTicks(n, 0); // tests wait on the main thread, which grants no steps meanwhile

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
            if (m.fire != kFireAuto) fprintf(f, " fire %d", m.fire);
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
    // Storage takes anything (#283): looked up first, items no recipe or fuel list had interned (food, tools, gems,
    // dyes...) were refused before the storage rule was reached. Insert interns the name.
    if (Def(type).power == Power::Storage) return !item.empty();
    const int id = tables.find(item);
    return id >= 0 && ::Accepts(type, id);
}

void Factory::SetSpeed(int ticksPerTick) { simSpeed = std::max(1, ticksPerTick); }

void Factory::Take(uint64_t body, uint32_t plot, int part, const std::string& item, int count, uint32_t to) {
    std::lock_guard<std::mutex> q(queueMutex);
    ItemMsg m{ 3, body, plot, part, tables.item(item), count }; m.to = to; itemQueue.push_back(m);
}

std::vector<std::pair<std::string, int>> Factory::ChestContents(uint64_t body, uint32_t plot, int part) {
    std::vector<std::pair<std::string, int>> out;
    std::lock_guard<std::mutex> v(viewMutex);
    if (remoteMode) { // (#854) a client: the chest's view the server sends for the chests near us
        auto rv = remoteViews.find(std::make_tuple(body, plot, part));
        if (rv != remoteViews.end()) for (const auto& s : rv->second.inSlots) if (s.second > 0) {
            bool merged = false;
            for (auto& o : out) if (o.first == s.first) { o.second += s.second; merged = true; }
            if (!merged) out.push_back(s);
        }
        return out;
    }
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
    if (remoteMode) return remoteWorking.count(WorkingKey(body, plot, part)) > 0; // (#332)
    return workingView.count(WorkingKey(body, plot, part)) > 0;
}

// ---- Multiplayer (#332) ----
void Factory::WatchMany(const std::vector<MachineKey>& keys) {
    std::lock_guard<std::mutex> v(viewMutex);
    watchMany = keys;
    for (auto it = watchManyViews.begin(); it != watchManyViews.end();)
        it = std::find(keys.begin(), keys.end(), it->first) == keys.end() ? watchManyViews.erase(it) : std::next(it);
}
MachineView Factory::Watched(const MachineKey& key) {
    std::lock_guard<std::mutex> v(viewMutex);
    auto it = watchManyViews.find(key);
    return it == watchManyViews.end() ? MachineView{} : it->second;
}
void Factory::SetNearPoints(const std::vector<Vector3>& points) { nearPoints = points; }
std::vector<int> Factory::WorkingParts(uint64_t body, uint32_t plot) {
    std::lock_guard<std::mutex> v(viewMutex);
    auto it = workingByPlot.find({ body, plot });
    return it == workingByPlot.end() ? std::vector<int>{} : it->second;
}
void Factory::SetRemote(bool remote) {
    std::lock_guard<std::mutex> v(viewMutex);
    remoteMode = remote;
    remoteViews.clear(); remoteWorking.clear(); remoteWorkingByPlot.clear();
}
bool Factory::Remote() { return remoteMode; }
void Factory::SetRemoteView(const MachineKey& key, const MachineView& view) {
    std::lock_guard<std::mutex> v(viewMutex);
    remoteViews[key] = view;
}
void Factory::SetRemoteWorking(uint64_t body, uint32_t plot, const std::vector<int>& parts) {
    std::lock_guard<std::mutex> v(viewMutex);
    auto& was = remoteWorkingByPlot[{ body, plot }];
    for (int p : was) remoteWorking.erase(WorkingKey(body, plot, p));
    was = parts;
    for (int p : parts) remoteWorking.insert(WorkingKey(body, plot, p));
}
void (*Factory::CursorSink())(const std::string&, int, float) { return cursorSink; }

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1214): sleeping machines and conveyor networks change nothing but the time ----
#include "unit_helpers.h"

namespace {
// `groups` generator + 2 macerators + chest groups on one plot (cables and conveyors as perf/factory-scale lays them).
LayoutMsg TestChains(uint32_t plotId, int groups, std::vector<std::array<int, 4>>& parts, std::array<int, 4>* rot = nullptr) {
    LayoutMsg L{};
    L.body = 7; L.plot = plotId; L.version = 1;
    int part = 0;
    auto machine = [&](PartType t, int x, int z) {
        std::vector<Int3> cells;
        for (int dz = 0; dz < 2; ++dz) for (int dy = 0; dy < 2; ++dy) for (int dx = 0; dx < 2; ++dx) cells.push_back({ x + dx, dy, z + dz });
        L.machines.push_back({ part, t }); L.machineCells.push_back(cells); L.machineVariant.push_back("");
        return part++;
    };
    for (int g = 0; g < groups; ++g) {
        const int bx = (g % 5) * 12, bz = (g / 5) * 5;
        const int gen = machine(PartType::Generator, bx, bz), m1 = machine(PartType::Macerator, bx + 3, bz);
        const int m2 = machine(PartType::Macerator, bx + 6, bz), ch = machine(PartType::Chest, bx + 9, bz);
        for (int x = 0; x < 8; ++x) L.cables.push_back({ bx + x, 0, bz - 1 });
        for (int x = 3; x < 11; ++x) L.conveyors.push_back({ bx + x, 0, bz + 2 });
        parts.push_back({ gen, m1, m2, ch });
    }
    if (rot) { // (#1269) three millstones stacked under a hand crank on a vertical shaft: 384 SU on 256, slowed to 66%
        const int up = RotationPlot_OrientFor({ 0, 1, 0 });
        auto add = [&](PartType t, Int3 at, int orient) {
            RotationPartMsg rm;
            const Part pt = RotationPlot_PartAt(t, at, orient);
            RotationPlot_Describe(pt, part, rm);
            L.rotation.push_back(rm);
            if (IsMachine(t)) {
                std::vector<Int3> cells;
                for (const Rotation::Cell& c : rm.geom.cells) cells.push_back({ c.x, c.y, c.z });
                L.machines.push_back({ part, t }); L.machineCells.push_back(cells); L.machineVariant.push_back("");
            }
            return part++;
        };
        const int a = add(PartType::Millstone, { 70, 0, 0 }, 0), b = add(PartType::Millstone, { 70, 2, 0 }, 0), c = add(PartType::Millstone, { 70, 4, 0 }, 0);
        add(PartType::Shaft, { 70, 6, 0 }, up);
        *rot = { add(PartType::HandCrank, { 70, 7, 0 }, up), a, b, c };
    }
    return L;
}

uint64_t StateHash(const PlotSim& ps) {
    uint64_t h = 1469598103934665603ull;
    auto mix = [&h](long long v) { h = (h ^ (uint64_t)v) * 1099511628211ull; };
    for (const Machine& m : ps.machines) {
        mix(m.part); mix((int)m.type); mix(m.recipe); mix(m.progress); mix(m.energy); mix(m.fuelEU); mix((int)m.status);
        mix(m.water); mix(m.heat); mix(m.fuel.item); mix(m.fuel.count);
        for (const Stack& s : m.in) { mix(s.item); mix(s.count); }
        for (const Stack& s : m.out) { mix(s.item); mix(s.count); }
    }
    for (const Net& n : ps.itemNets) { mix(n.credit); mix(n.rr); }
    for (const Machine& m : ps.machines) mix(m.rotAcc);
    mix((long long)ps.rot.Hash()); // (#1269) the shaft networks' speeds, statuses and cranks
    return h;
}

struct ChainsRun { uint64_t hash = 0; long long made = 0; int asleepMachines = 0, asleepNets = 0, sandInLateChest = 0, gravelMilled = 0; };

// 60 groups on six plots for 2400 ticks (2 min): every third group starved (coal, no stone: idle macerators), every
// seventh without coal (unpowered); at tick 600 a starved group gets stone (it must wake and work). `sleep` on or off;
// the plots ticked by TickPlots on `workers` threads (#997: the parallel tick must give the same result).
ChainsRun RunChains(bool sleep, int workers = 1) {
    LoadRecipes();
    const bool was = sleepOn;
    const int minWas = parallelMinMachines;
    sleepOn = sleep;
    parallelMinMachines = 0; // (these plots are small: split them anyway)
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    constexpr int kPlots = 6;
    std::vector<std::vector<std::array<int, 4>>> parts(kPlots);
    std::vector<std::array<int, 4>> rot(kPlots); // (#1269) crank, 3 millstones per plot
    for (int p = 0; p < kPlots; ++p) ApplyLayout(TestChains(1 + p, 10, parts[p], &rot[p]));
    std::vector<PlotSim*> all;
    for (int p = 0; p < kPlots; ++p) all.push_back(FindPlot(7, 1 + p));
    PlotSim* a = all[0];
    const std::vector<std::array<int, 4>>& partsA = parts[0];
    auto feed = [](PlotSim& ps, const std::vector<std::array<int, 4>>& parts) {
        for (size_t g = 0; g < parts.size(); ++g) {
            int buckets = 0;
            Machine& gen = ps.machines[ps.byPart.at(parts[g][0])];
            if (g % 7 != 6) InsertInto(gen, tables.item("coal"), 16, false, buckets);
            if (g % 3 == 2) continue; // starved
            InsertInto(ps.machines[ps.byPart.at(parts[g][1])], tables.item("stone"), 24, false, buckets);
            InsertInto(ps.machines[ps.byPart.at(parts[g][2])], tables.item("stone"), 24, false, buckets);
        }
    };
    for (int p = 0; p < kPlots; ++p) feed(*all[p], parts[p]);
    for (int p = 0; p < kPlots; ++p) // (#1269) stone for the millstones (the third only on even plots: two or three loads)
        for (int k = 1; k <= 3; ++k) { int buckets = 0; if (k < 3 || p % 2 == 0) InsertInto(all[p]->machines[all[p]->byPart.at(rot[p][k])], tables.item("stone"), 24, false, buckets); }
    ChainsRun r;
    long long made = 0; int active = 0;
    std::vector<PlotJob> jobs;
    for (int t = 0; t < 2400; ++t) {
        if (t == 600) { // a player puts stone into a starved group's macerator (Insert: the command wakes it)
            Machine& m = a->machines[a->byPart.at(partsA[2][1])];
            m.dirty = true;
            int buckets = 0;
            InsertInto(m, tables.item("stone"), 8, false, buckets);
        }
        for (int p = 0; p < kPlots; ++p) if ((t / 150 + p) % 3 != 0) all[p]->rot.Hold(rot[p][0]); // (#1269) cranked on and off
        jobs.clear();
        for (PlotSim* ps : all) jobs.push_back({ ps, 1 });
        TickPlots(jobs, workers, made, active);
    }
    for (int p = 0; p < kPlots; ++p)
        for (int k = 1; k <= 3; ++k) r.gravelMilled += CountItems(all[p]->machines[all[p]->byPart.at(rot[p][k])].out, kSlots, tables.item("gravel"));
    for (const PlotSim* ps : all) r.hash = r.hash * 31 + StateHash(*ps);
    r.made = made;
    for (const PlotSim* ps : all) {
        for (const Machine& m : ps->machines) r.asleepMachines += m.asleep;
        for (const Net& n : ps->itemNets) r.asleepNets += n.asleep;
    }
    const Machine& lateChest = a->machines[a->byPart.at(partsA[2][3])];
    r.sandInLateChest = CountItems(lateChest.in.data(), kChestSlots, tables.item("sand"));
    sleepOn = was;
    parallelMinMachines = minWas;
    workerPool.stopAll();
    plots.clear(); plotIndexDirty = true;
    return r;
}
} // namespace

// (#997, user decision) the determinism gate for parallel plots: 1 and 4 workers give the same factory, bit for bit.
TEST(Factory, PlotsTickedOnFourWorkersMatchOneWorker) {
    const ChainsRun one = RunChains(true, 1), four = RunChains(true, 4);
    EXPECT_GT(one.made, 0);
    EXPECT_EQ(four.made, one.made);
    EXPECT_EQ(four.hash, one.hash) << "every machine and network the same whatever the worker count";
    EXPECT_GT(one.gravelMilled, 0) << "(#1269) the cranked millstones ground stone, slowed";
    EXPECT_EQ(four.gravelMilled, one.gravelMilled) << "rotation is per plot: the same on any worker";
}

TEST(Factory, SleepingMachinesAndNetworksChangeNothingButTheTime) {
    const ChainsRun awake = RunChains(false), sleeping = RunChains(true);
    EXPECT_GT(awake.made, 0) << "the chains work";
    EXPECT_EQ(sleeping.made, awake.made);
    EXPECT_EQ(sleeping.hash, awake.hash) << "every machine's slots, recipe, progress, fuel and status, and every network's credit and round-robin, are the same";
    EXPECT_GT(sleeping.asleepMachines, 0) << "the starved macerators sleep";
    EXPECT_GT(sleeping.asleepNets, 0) << "their conveyor networks sleep";
    EXPECT_GT(sleeping.sandInLateChest, 0) << "the starved group woke when stone came and its sand reached its chest";
}
#endif
