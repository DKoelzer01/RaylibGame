#include "factory.h"
#include "inventory.h"
#include "plots.h"
#include "machinetiers.h"
#include "ores.h"
#include "forms.h"
#include "disassembly.h"
#include "composed.h"
#include "itemdata.h"
#include "rotationplot.h"
#include "rotationsources.h"
#include "conduitnet.h"
#include "multiblock.h"
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
#include <filesystem>
#include <map>
#include <mutex>
#include <sstream>
#include <set>
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

// (#1102) centi: hundredths of one more item (0-99), the player inventory's resolution: chests keep fractions (0.4 dirt);
// a slot holding only a fraction is not empty. Machines craft and conveyors move whole items only.
struct Stack { int item = -1, count = 0, centi = 0; std::string data; }; // (#1446 2b) data: one item's own record (a chest holds it, one a slot; automation never moves it)
inline bool Plain(const Stack& s) { return s.data.empty(); }
inline bool Empty(const Stack& s) { return s.count <= 0 && s.centi <= 0; }
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
    // (#1734) Outputs by chance ("iron_ore 8%" in machine_recipes.txt): count = hundredths of an item per craft. A machine
    // carries the fractions (Machine::disCarry, saved), so every machine gives exactly its share, deterministically: the
    // 13th sieved gravel of a machine always gives its first iron ore. Single machines only (a multiblock's jobs drop them).
    std::vector<Stack> chance;
    int ticks = 0, eu = 0;   // generators: ticks of burn per fuel item, EU per tick produced
    int tier = -1;           // (#1427 L7a) its voltage tier (machine_recipes.txt's sixth column; else its machine's first)
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
        if (id < kMaxCapItems) fluids[id] = FluidName(n); // (#1254)
        itemIds[n] = id;
        itemNames.push_back(n);
        return id;
    }
    // (#1254) Fluid items (alchemy's essences): names with a registered prefix. Lock-free reads like cap().
    bool fluid(int id) const { return id >= 0 && id < kMaxCapItems && fluids[id]; }
    void addFluidPrefix(const std::string& prefix) {
        std::lock_guard<std::mutex> lock(mutex);
        if (prefix.empty() || std::find(fluidPrefixes.begin(), fluidPrefixes.end(), prefix) != fluidPrefixes.end()) return;
        fluidPrefixes.push_back(prefix);
        for (int id = 0; id < (int)itemNames.size() && id < kMaxCapItems; ++id) fluids[id] = FluidName(itemNames[id]);
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
    bool fluids[kMaxCapItems] = {};
    std::vector<std::string> fluidPrefixes;
    bool FluidName(const std::string& n) const { for (const std::string& p : fluidPrefixes) if (n.compare(0, p.size(), p) == 0) return true; return false; }
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
enum class Power : uint8_t { None, Fuel, Eu, Generator, Storage, Battery, Rotation, Steam, Engine, Dynamo, Motor, Manual }; // Rotation (#1269): at its shaft's speed; Steam..Motor (#1271 L3): a boiler for engines, a steam engine, a dynamo, an electric motor; Manual (#1271 L2): worked by a player's strikes (the anvil)
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
                : p == "storage" ? Power::Storage : p == "battery" ? Power::Battery : p == "rotation" ? Power::Rotation
                : p == "steam" ? Power::Steam : p == "engine" ? Power::Engine : p == "dynamo" ? Power::Dynamo : p == "motor" ? Power::Motor : p == "manual" ? Power::Manual : Power::None;
        d.euOut = atoi(cols[3].c_str());
        if (cols.size() > 4) { std::stringstream fs(cols[4]); std::string fuel; while (fs >> fuel) d.fuels.push_back(tables.item(fuel)); }
        if (cols.size() > 5) d.waterPerTick = std::max(0, atoi(cols[5].c_str()));
    }
}

// (#1207) built once by LoadRecipes
std::array<std::vector<int>, 256> recipesOf;      // by PartType: indices into tables.recipes
std::array<std::vector<uint8_t>, 256> wantsOf;   // by PartType, by item id: a recipe of it takes the item
std::vector<uint8_t> disassemblable;              // (#1427 L7c) by item id: the disassembly index takes it apart
// (#1214) sleeping machines and networks: on unless RAYLIBGAME_FACTORY_SLEEP=0 (A/B; the unit test runs both)
bool sleepOn = [] { const char* e = getenv("RAYLIBGAME_FACTORY_SLEEP"); return !(e && e[0] == '0'); }();

void LoadRecipes() {
    if (!tables.recipes.empty()) return;
    LoadMachineDefs();
    std::ifstream f("resources/data/machine_recipes.txt");
    std::string line;
    auto parseStacks = [](const std::string& s, std::vector<Stack>& out, std::vector<Stack>* chance = nullptr) {
        std::stringstream ss(s);
        std::string part;
        while (std::getline(ss, part, ',')) {
            std::stringstream ps(part);
            std::string name, n;
            if (!(ps >> name >> n) || name == "-") continue;
            if (chance && n.back() == '%') { if (atoi(n.c_str()) > 0) chance->push_back({ tables.item(name), atoi(n.c_str()) }); } // (#1734) "8%": hundredths per craft
            else out.push_back({ tables.item(name), atoi(n.c_str()) });
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
        parseStacks(cols[1], r.out, &r.chance);
        parseStacks(cols[2], r.in);
        r.ticks = std::max(1, atoi(cols[3].c_str()));
        r.eu = std::max(0, atoi(cols[4].c_str()));
        if (cols.size() >= 6 && !cols[5].empty()) r.tier = MachineTiers_ByName(cols[5]); // (#1427 L7a)
        tables.recipes.push_back(r);
    }
    { // (#1614) The form table's machine recipes (forms.txt): every metal's plates in the press and the compressor, ...,
      // where machine_recipes.txt has none for that machine and output
        std::set<std::pair<int, int>> have;
        for (const Recipe& r : tables.recipes) for (const Stack& o : r.out) have.insert({ (int)r.machine, o.item });
        for (const FormRecipe& f : Forms_All()) {
            if (f.machine.empty()) continue;
            Recipe r;
            r.machine = MachineByName(f.machine);
            if (r.machine == PartType::Count) continue;
            const int out = tables.item(f.output);
            if (out < 0 || have.count({ (int)r.machine, out })) continue;
            r.out.push_back({ out, f.outCount });
            bool ok = true;
            for (const auto& [name, n] : f.inputs) { const int id = tables.item(name); ok = ok && id >= 0; r.in.push_back({ id, n }); }
            if (!ok) continue;
            r.ticks = std::max(1, f.ticks);
            r.eu = std::max(0, f.eu);
            tables.recipes.push_back(r);
            have.insert({ (int)r.machine, out });
        }
    }
    if (gFactoryExtraRecipes) { // (#1254) recipes the game makes from its own data (alchemy's distillery and cauldron)
        std::vector<Factory::RecipeInfo> extra;
        gFactoryExtraRecipes(extra);
        for (const Factory::RecipeInfo& e : extra) {
            if (e.machine >= PartType::Count || e.in.empty() || e.out.empty()) continue;
            Recipe r;
            r.machine = e.machine;
            for (const auto& [name, n] : e.in) if (n > 0) r.in.push_back({ tables.item(name), n });
            for (const auto& [name, n] : e.out) if (n > 0) r.out.push_back({ tables.item(name), n });
            r.ticks = std::max(1, e.ticks);
            r.eu = std::max(0, e.eu);
            tables.recipes.push_back(r);
        }
    }
    // (#1427 L7a) A recipe with no tier is its machine's first: a tiered family's first tier, else its machines.txt voltage's.
    for (Recipe& r : tables.recipes)
        if (r.tier < 0) r.tier = MachineTiers_Tiered(r.machine) ? MachineTiers_First(r.machine) : MachineTiers_ForVolts(Def(r.machine).volts);
    // (#1207) Per machine type: its recipes in file order, and which items they take (a recipe start and Accepts looked
    // through every recipe of every machine, each tick).
    for (size_t ri = 0; ri < tables.recipes.size(); ++ri) {
        const Recipe& r = tables.recipes[ri];
        recipesOf[(int)r.machine].push_back((int)ri);
        std::vector<uint8_t>& w = wantsOf[(int)r.machine];
        for (const Stack& s : r.in) { if (s.item >= (int)w.size()) w.resize(s.item + 1, 0); w[s.item] = 1; }
    }
    // (#1427 L7c) The disassembly index, from the forward recipes: hand recipes, then the assembler's (one output each)
    std::vector<DisassemblyRecipe> forward;
    for (const ::Recipe& h : itemDB.recipes()) if (!h.formShape) forward.push_back({ h.output, h.outputCount, h.inputs, 0, 0 }); // (#1614: form shapes are raw material to disassembly)
    for (const Recipe& r : tables.recipes) {
        if (r.machine != PartType::LvAssembler || r.out.size() != 1) continue;
        DisassemblyRecipe d{ tables.name(r.out[0].item), r.out[0].count, {}, r.ticks, r.eu };
        for (const Stack& s : r.in) d.in.push_back({ tables.name(s.item), s.count });
        forward.push_back(std::move(d));
    }
    Composed_Is(""); // (its tables load here, on the loading thread: the sim reads them)
    Disassembly_Build(forward, gDisassemblyLadderTier);
    for (const std::string& n : Disassembly_Items()) { const int id = tables.item(n); if (id >= (int)disassemblable.size()) disassemblable.resize(id + 1, 0); disassemblable[id] = 1; }
    LOG_INFO("Factory", "%zu machine recipes, %zu items to take apart\n", tables.recipes.size(), Disassembly_Items().size());
}

// ---------------------------------------------------------------------------------------------
// Simulation state (owned by the simulation thread, guarded by simMutex)
// ---------------------------------------------------------------------------------------------

constexpr int kTierLV = 32;                 // volts
constexpr int kBatteryCapacity = 100000;    // EU
constexpr int kBatteryRate = 64;            // EU per tick in or out
constexpr int kConveyorItemsPerSecond = 8;
constexpr int kStrikeTicks = 40;          // (#1271 L2) one strike on an anvil: 2 s of its recipe's work (a bloom takes 3)

struct Machine {
    int part = -1;
    // (#1214) asleep: its recipe search found nothing and its contents have not changed since; dirty: they changed (an insert,
    // a take, fuel taken, a craft, a variant, fire or layout change) - it wakes, and its conveyor network wakes
    bool asleep = false, dirty = true;
    // (#1261) the summary sim's measurement: crafts per recipe since the plot's window began, and the caps they gave when
    // the plot entered Summary (crafts per kCapTicks of each recipe it ran)
    std::map<int, int> crafts;
    std::map<int, int> capPerWindow;
    long long nearEuCur = 0, nearEuPrev = 0, nearEuCap = 0; // (#1261) a generator's EU in them
    std::map<int, int> nearCur, nearPrev, nearCap; // (#1261) crafts per recipe while near, in two 30 s buckets; caps from them on leaving
    // (#1403) the recipe and progress at the start of the measuring window and of each near bucket: a window's work is
    // its completions less what the first was already done before it, plus what the last unfinished one has done since
    int winRecipe = -1, winProgress = 0, nearCurRecipe = -1, nearCurProgress = 0, nearPrevRecipe = -1, nearPrevProgress = 0;
    int euCarry = 0;                            // (#1266) EU received and not yet a whole tick of its recipe
    long long euWindow = 0, euCap = 0;          // (#1261) a generator: EU it delivered in the window, and per kCapTicks
    std::map<int, long long> carry;             // (#1261) part-made crafts per recipe between advances (fixed point)
    PartType type = PartType::Count;
    std::vector<Int3> cells;
    std::vector<Stack> in;                  // InSlots(): kSlots inputs, or a chest's kChestSlots of storage (#495: only chests pay for 50)
    Stack out[kSlots];
    Stack fuel;                             // fuel-burning machines and generators: their own fuel slot (#194)
    int water = 0, heat = 0;                // boilers (#229): litres in the tank; heat 0..kHeatMax
    uint8_t fire = kFireAuto;               // campfires (#1103): lit or put out by a player, else by itself
    uint8_t mode = 0;                       // (#1271 L6c) an assembler: 0 assembles, 1 disassembles (#1427; until then it idles)
    uint8_t anvilJob = 0;                   // (#1219) an anvil's job on a composed copy: 1 repair, 2 take a gem out (progress: its work)
    int pick = -1;                          // (#1271 L6c) the recipe output a player chose for a machine of several recipes (-1: any)
    int quality = 100;                      // (#1446 2c) the craft in progress' quality % for a composed output (T2: the anvil's bar sets it; back to 100 after)
    int recipe = -1, progress = 0;         // generators use progress as the banked-fire tick counter
    int energy = 0;                         // battery charge / generator fuel ticks left
    int fuelEU = 0;                         // generator: EU per tick of the current fuel
    int tier = 0;                           // (#1427 L7a) its voltage tier (machine_tiers.txt: 0 LV, 1 MV ...; -1 under LV), from its part
    int variant = -1;                       // item id: a miner only makes the ore of the deposit under it
    MachineStatus status = MachineStatus::Idle;
    int powerNet = -1, itemNet = -1;
    int received = 0;                       // EU this tick (consumers)
    std::vector<int> powerNodes;            // (#1620 C2b) its nodes in the plot's power conduit graph (where cables touch it), BuildNets
    long long reach = -1, canSend = -1;     // (#1620 C2b) this step: what the cables can bring a consumer / take from a generator or dynamo; -1 = no cap
    int rotAcc = 0;                        // (#1269) rotation machines: permille of a tick of progress carried over
    int boiler = -1;                        // (#1271 L3) a steam engine: the boiler it draws steam from (one touching it), else -1
    int bellows = -1;                       // (#1271 L2) a bloomery or blast furnace: the bellows touching it (blowing while its shaft turns), else -1
    int strike = 0;                         // (#1271 L2) a manual machine (the anvil): ticks of work a player's strikes have put in, not yet used
    // (#1427 L7c) an assembler set to disassemble: the input slot it is taking apart (-1: none), that job's time and EU a
    // tick (overclocked), the raw parts' fractions carried (hundredths, by item id) and the chance rolls' state (both saved)
    int disSlot = -1, disTicks = 0, disEu = 0;
    std::map<int, int> disCarry;
    uint32_t disRng = 0;
    // (#1429 PR 4) Multiblocks. A hatch (MbHatch) is a machine of its own: an item input hatch a buffer conveyors and hands fill
    // (its inputs), an item output hatch one they empty (its outputs), an energy hatch a consumer on its cables drawing what its
    // structures ask (hatchDemand, EU a tick) and passing them their share. Every block may be shared (the user, 2026-11-23).
    uint8_t hatch = 0;                      // kHatch*: its kind (0 not a hatch)
    long long hatchDemand = 0;              // an energy hatch: its users' EU a tick this step
    std::map<int, int> reserved;            // an item output hatch: item -> count running recipes will put in it (their room is kept)
    // A controller (MbController): formed (the plot says so and it validates), the machine whose recipes it runs, its jobs
    // (parallel recipes: inputs claimed at their start, outputs' room reserved), its hatches (indices into the plot's machines)
    struct MbJob { int recipe = -1, progress = 0, euCarry = 0; std::vector<std::pair<int, Stack>> res; }; // res: hatch part, its reserved stack
    bool mbFormed = false;
    PartType mbMachine = PartType::Count;
    int mbParallel = 1;
    std::vector<int> mbHatches;
    std::vector<MbJob> jobs;
    long long mbEu = 0;                     // EU its energy hatches passed it this step
};
constexpr uint8_t kHatchNone = 0, kHatchItemIn = 1, kHatchItemOut = 2, kHatchEnergy = 3, kHatchOther = 4; // other: fluid, rotation (not in the factory yet)
uint8_t HatchOf(const std::string& kind) { return kind == "item_in" ? kHatchItemIn : kind == "item_out" ? kHatchItemOut : kind == "energy" ? kHatchEnergy : kHatchOther; }
// (#1427 L7a) A recipe in a machine: overclocked by the tiers it is under the machine (machine_tiers.txt), its time down
// and its EU a tick up; and the voltage a machine takes (a tiered variant's, else machines.txt's).
static int StepsOf(const Machine& m, const Recipe& r) { return std::max(0, MachineTiers_Steps(m.tier, r.tier)); }
static int TicksOf(const Machine& m, const Recipe& r) { return MachineTiers_Ticks(r.ticks, StepsOf(m, r)); }
static int EuOf(const Machine& m, const Recipe& r) { return MachineTiers_Eu(r.eu, StepsOf(m, r)); }
static int VoltsOf(const Machine& m) {
    if (m.hatch == kHatchEnergy) return MachineTiers_Volts(m.tier); // (#1429 PR 4) its voltage tier
    const int v = Def(m.type).volts; return MachineTiers_Tiered(m.type) ? std::max(v, MachineTiers_Volts(m.tier)) : v;
}
// (#1429 PR 4) How a machine is powered: an energy hatch draws from its cables (an EU consumer); else its type's.
static Power PowerOf(const Machine& m) { return m.hatch == kHatchEnergy ? Power::Eu : Def(m.type).power; }
// EU a tick its job (a recipe, or a disassembly: #1427 L7c) asks for; -1 with no job
static int JobEu(const Machine& m);

struct Net {
    std::vector<int> machines; int credit = 0, rr = 0, voltage = kTierLV;
    bool asleep = false; // (#1214) its last pass moved nothing and none of its machines changed since
    long long supply = 0, demand = 0; int ticks = 1; // the last step's EU made and wanted, over its ticks (cable tooltips, #570)
    long long carried = 0;                           // (#1620 C2b) of that demand, what the cables could carry ("carries X of Y")
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
    ConduitGraph powerGraph;                // (#1620 C2b) the cables as 1x1 copper runs (capacity per cable path), empty without the copper material
    PlotStats stats;                        // (#1199)
    RotationPlot rot;                       // (#1269) its shaft networks (per plot: ticks on any worker)
    // (#1207) what this plot last put in the near-plot views (working machines, chests and a hash of their contents)
    struct Published { bool near = false; std::vector<int> working, chests, stalled; uint64_t chestHash = 0; std::vector<std::pair<int, int>> anvils; } pub; // anvils (#1613): (part, workpiece item)
    // (#1261) the summary sim's state (doc 2.1): Near ticks every tick; Measuring takes today's 20-tick step for 60 s and
    // counts each machine's crafts; Summary is advanced by rates (S2) until the plot comes near again.
    enum class Mode : uint8_t { Near, Measuring, Summary } mode = Mode::Near;
    long long modeSince = 0; // the tick the current mode (and the measuring window) began
    long long lastAdvance = 0; // (Summary) the tick its buffers were last advanced to
    bool loadedSummary = false; // (#1302) loaded from a save while in the summary sim: far, it goes on there with its saved caps
    long long bucketSince = 0; bool bucketFull = false; // (#1261) the near crafts' current 30 s bucket began; a whole one before it
    std::vector<std::pair<int, int>> completed;         // (#1728) (item, count) of recipes this plot really finished since the last hand-over (TickPlot only: never the summary sim); emptied under viewMutex
    bool hasStructures = false;                         // (#1429 PR 4) a formed multiblock: far, it keeps the 20-tick step (not in the summary sim yet: #1763)
    int mbRR = 0;                                       // (#1429 PR 4) the controller (machine index) that tries to start a recipe first next
};

// (#1261) The summary sim's plot states. RAYLIBGAME_FACTORY_SUMMARY=0: far plots keep today's 20-tick step (A/B).
constexpr long long kMeasureTicks = 60 * kTicksPerSecond; // 60 s of measuring (decided)
constexpr long long kCapTicks = 60 * kTicksPerSecond;     // a cap is crafts per this many ticks
bool summaryOn = [] { const char* e = getenv("RAYLIBGAME_FACTORY_SUMMARY"); return !(e && e[0] == '0'); }();
// RAYLIBGAME_FACTORY_SUMMARY_LOG=1: each plot's mode changes (with its caps) and each summary advance, to the log.
bool summaryLog = [] { const char* e = getenv("RAYLIBGAME_FACTORY_SUMMARY_LOG"); return e && e[0] == '1'; }();

void AdvanceSummary(PlotSim& ps, long long t0, long long t1); // (below)
void BeginMeasuring(PlotSim& ps, long long t) {
    if (summaryLog) LOG_INFO("Factory", "summary: plot %u measures from tick %lld (was mode %d)\n", ps.plot, t, (int)ps.mode);
    if (ps.mode == PlotSim::Mode::Summary) AdvanceSummary(ps, ps.lastAdvance, t); // (its buffers to now first)
    ps.mode = PlotSim::Mode::Measuring;
    ps.modeSince = t;
    for (Machine& m : ps.machines) { m.crafts.clear(); m.euWindow = 0; m.winRecipe = m.recipe; m.winProgress = m.progress; } // (#1403: where it stands)
}
// While near, crafts and EU go into 30 s buckets: the last whole bucket and the current one (30-60 s) are what the plot
// really did near the player.
constexpr long long kNearBucketTicks = 30 * kTicksPerSecond;
void RotateNearBuckets(PlotSim& ps, long long t) {
    if (t - ps.bucketSince < kNearBucketTicks) return;
    for (Machine& m : ps.machines) {
        m.nearPrev.swap(m.nearCur); m.nearCur.clear(); m.nearEuPrev = m.nearEuCur; m.nearEuCur = 0;
        m.nearPrevRecipe = m.nearCurRecipe; m.nearPrevProgress = m.nearCurProgress; m.nearCurRecipe = m.recipe; m.nearCurProgress = m.progress; // (#1403)
    }
    ps.bucketSince = t;
    ps.bucketFull = true;
}
// (#1403) A cap from a window: recipe r's work in it, in ticks of the recipe - n completions, less the progress the first
// already had when the window began (r0, p0), plus the progress of the one still going - per kCapTicks, rounded down.
// Counting whole completions made one craft more a window when crafts are long (a 15 s craft, 2-4 in a 30-60 s window): a
// craft nearly done before the window began counted whole, and the summary made more than the full sim.
int CapOf(const Machine& m, int r, int n, int r0, int p0, long long window) {
    if (r < 0 || r >= (int)tables.recipes.size() || window <= 0) return 0;
    const long long R = std::max(1, TicksOf(m, tables.recipes[r])); // (#1427 L7a) overclocked, as its progress counts
    long long work = (long long)n * R + (m.recipe == r ? m.progress : 0) - (r0 == r ? p0 : 0);
    if (work < 0) work = 0;
    return (int)(work * kCapTicks / (R * window));
}
// The player (or a server's players) came near or went away. Leaving after 30 s or more near, each machine keeps the caps
// of its last 30-60 s (what it really did near the player: the 20-tick measuring step loses part of a step at each craft, and a
// window right after a start-up is low); measuring then may only raise them.
void SetNear(PlotSim& ps, bool near, long long t) {
    if (summaryLog) LOG_INFO("Factory", "summary: plot %u near %d at tick %lld (mode %d)\n", ps.plot, near ? 1 : 0, t, (int)ps.mode);
    ps.nearPlayer = near;
    if (near) {
        if (ps.mode == PlotSim::Mode::Summary) AdvanceSummary(ps, ps.lastAdvance, t); // its buffers to now, then it ticks
        ps.mode = PlotSim::Mode::Near; ps.modeSince = t;
        ps.bucketSince = t; ps.bucketFull = false;
        for (Machine& m : ps.machines) { m.nearCur.clear(); m.nearPrev.clear(); m.nearEuCur = m.nearEuPrev = 0; m.nearCurRecipe = m.recipe; m.nearCurProgress = m.progress; m.nearPrevRecipe = -1; m.nearPrevProgress = 0; }
        return;
    }
    if (ps.mode != PlotSim::Mode::Near) return;
    // (#1302) Loaded from a save in the summary sim and far before 30 s near: it goes on in the summary with its saved caps.
    const bool loaded = ps.loadedSummary;
    ps.loadedSummary = false;
    if (loaded && !ps.bucketFull && summaryOn && !ps.hasStructures) {
        for (Machine& m : ps.machines) { m.capPerWindow = m.nearCap; m.euCap = m.nearEuCap; }
        ps.mode = PlotSim::Mode::Summary;
        ps.lastAdvance = ps.modeSince = t;
        if (summaryLog) LOG_INFO("Factory", "summary: plot %u goes on in Summary from the save at tick %lld\n", ps.plot, t);
        return;
    }
    const long long window = ps.bucketFull ? kNearBucketTicks + (t - ps.bucketSince) : 0;
    for (Machine& m : ps.machines) {
        m.nearCap.clear(); m.nearEuCap = 0;
        if (!ps.bucketFull) continue; // (under 30 s near: measuring alone; else the last 30-60 s near)
        std::map<int, int> sum = m.nearPrev;
        for (const auto& [r, n] : m.nearCur) sum[r] += n;
        for (const auto& [r, n] : sum) m.nearCap[r] = CapOf(m, r, n, m.nearPrevRecipe, m.nearPrevProgress, window); // (#1403)
        m.nearEuCap = (m.nearEuPrev + m.nearEuCur) * kCapTicks / window;
    }
    BeginMeasuring(ps, t);
}
// After 60 s of measuring, a far plot enters Summary with each machine's caps: its crafts per kCapTicks of each recipe it
// ran (in the measuring window, or the minute before it left near if that was more), rounded down (no margin, decided:
// the oracle test holds the summary at or under the full sim).
void UpdateMode(PlotSim& ps, long long t) {
    if (!summaryOn || ps.mode != PlotSim::Mode::Measuring || t - ps.modeSince < kMeasureTicks) return;
    if (ps.hasStructures) return; // (#1429 PR 4) not in the summary sim yet: far, it keeps the 20-tick step
    const long long window = t - ps.modeSince;
    for (Machine& m : ps.machines) {
        m.capPerWindow = m.nearCap;
        for (const auto& [recipe, n] : m.crafts) { int& c = m.capPerWindow[recipe]; c = std::max(c, CapOf(m, recipe, n, m.winRecipe, m.winProgress, window)); } // (#1403)
        m.euCap = std::max(m.nearEuCap, m.euWindow * kCapTicks / window);
        if (summaryLog && (!m.capPerWindow.empty() || m.euCap > 0)) {
            std::string caps;
            for (const auto& [r, c] : m.capPerWindow) caps += TextFormat(" r%d:%d(near %d, measured %d)", r, c, m.nearCap.count(r) ? m.nearCap.at(r) : 0, m.crafts.count(r) ? (int)((long long)m.crafts.at(r) * kCapTicks / window) : 0);
            LOG_INFO("Factory", "summary: plot %u part %d caps per %lld ticks:%s eu %lld (near %lld)\n", ps.plot, m.part, kCapTicks, caps.c_str(), m.euCap, m.nearEuCap);
        }
    }
    if (summaryLog) LOG_INFO("Factory", "summary: plot %u enters Summary at tick %lld after %lld ticks of measuring\n", ps.plot, t, window);
    ps.mode = PlotSim::Mode::Summary;
    ps.lastAdvance = t;
    ps.modeSince = t;
}
// The step plot i takes at tick t: 1 near, 20 once a second while measuring (staggered by plot), 0 in Summary (S1: a
// summary plot stands still; S2 advances it by rates).
int StepOf(const PlotSim& ps, size_t i, long long t) {
    if (ps.nearPlayer) return 1;
    if (summaryOn && ps.mode == PlotSim::Mode::Summary) return 0;
    return (t + (long long)i) % kTicksPerSecond == 0 ? kTicksPerSecond : 0;
}
thread_local PlotStats* gCount = nullptr;  // the plot ticking now on this thread (BurnFuel counts the fuel it takes into it; thread_local: plots may tick on several threads later, #997)
long long statSeconds = 0;    // seconds rolled since the reset (how much of each window is filled)
// (#1207) The last minute of factory ticks' times (sim thread), for the mean, p99 and worst and the phase split.
struct TickTimes { float total = 0, cmd = 0, tick = 0, view = 0; };
constexpr size_t kTickRing = 1200;
std::array<TickTimes, kTickRing> tickRing{};
size_t tickRingAt = 0;
std::atomic<bool> tickRingClear{ false }; // ResetTickStats: a timed window starts (the scene's own setup ticks leave the ring)

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
    std::vector<int> machineTier;                          // (#1427 L7a) per machine: its voltage tier (a tiered machine's variant)
    std::vector<Int3> conveyors, cables;
    uint32_t refundTo = 0;                                // (#854) who removed machines in this change (their contents go to that player)
    std::vector<RotationPartMsg> rotation;                // (#1269) shafts, cogs, cranks and rotation machines
    // (#1429 PR 4) per machine its hatch kind (kHatch*; empty: none); and the formed structures: controller part, the machine
    // whose recipes it runs, its parallel, its voltage tier, its hatches' parts
    std::vector<uint8_t> machineHatch;
    struct Structure { int controller = -1; PartType machine = PartType::Count; int parallel = 1, tier = -1; std::vector<int> hatches; };
    std::vector<Structure> structures;
};
struct ItemMsg { int kind; uint64_t body; uint32_t plot; int part; int item; int count; int amount = 0; uint32_t to = 0; int centi = 0; std::string data; }; // (#1446 2b) data: one item's record (kinds 0, 9: count 1) // (#1102) centi: a fraction (hundredths) more, for a chest // amount: kind 8, how many (0 = all); kind 9, the slot // 0 insert, 1 extract, 2 seed, 3 take (one item from a chest), 4 take slot (item = slot, count: 0 input, 1 output, 2 fuel), 5 drain outputs, 6 insert into the fuel slot, 7 water (count = litres), 8 slot to the cursor, 9 insert into one chest slot (#862)

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
// (#1271 L6b) machine recipes the research has not unlocked (by output item): queued from the main thread, applied at the
// command phase (lockedSim is read by the plots' ticks, never written while they run)
std::vector<int> lockedQueue; bool lockedQueued = false;
std::unordered_set<int> lockedSim;
bool RecipeLocked(const Recipe& r) { if (lockedSim.empty()) return false; for (const Stack& s : r.out) if (lockedSim.count(s.item)) return true; return false; }
struct ItemOut { int item = 0, count = 0; uint32_t to = 0; int centi = 0; std::string data; }; // (#1446 2b) data: one item's record (count 1) // (#1102) centi: a fraction more (from a chest) // (#332) to: the player it goes to (0: this game's own inventory)
std::vector<ItemOut> results;     // sim -> main: (item, count) back to the inventory
std::vector<ItemOut> cursorResults; // sim -> main: (item, count) onto the cursor (#192)
struct MadeOut { uint64_t body = 0; uint32_t plot = 0; int item = 0, count = 0; };
std::vector<MadeOut> madeResults; // (#1728) sim -> main: recipes machines finished, in plot order (a hand-over to gFactoryMade by Sync)
void (*cursorSink)(const std::string&, int, float, const std::string&) = nullptr; // (#1446) data: one item's record
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
std::map<std::tuple<uint64_t, uint32_t, int>, int> anvilView;     // (#1613) the workpiece on each anvil of a near plot (item; viewMutex)
std::map<std::pair<uint64_t, uint32_t>, std::vector<int>> workingByPlot; // (#332) the same, by plot
std::map<std::pair<uint64_t, uint32_t>, std::vector<int>> stalledByPlot; // (#1243) machines stopped for a reason (no power, output full, no fuel, no water) on near plots (viewMutex)
// (#1269) every rotation part of the plots near the player, as of the last tick (viewMutex): models, gauges, tooltips
struct RotPub { int part = -1; int64_t milli = 0; Rotation::Status status = Rotation::Status::Stopped; Rotation::NetResult net; };
std::map<std::pair<uint64_t, uint32_t>, std::vector<RotPub>> rotView;
std::vector<Factory::MachineKey> watchMany;                        // (#332) a server's open machines (viewMutex)
std::map<Factory::MachineKey, MachineView> watchManyViews;
std::vector<Vector3> nearPoints;                                    // (#332) main thread: a server's players
bool cameraNear = true;                                              // (#1302) main thread: the camera makes plots near (not on a dedicated server)
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
std::unordered_map<std::string, uint64_t> sentFormed;    // (#1429 PR 4) main thread: the formed list each layout was sent with (a hash)

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
        if (s[i].item == item && Plain(s[i]) && s[i].count < cap) { const int k = std::min(count, cap - s[i].count); s[i].count += k; count -= k; }
    for (int i = 0; i < n && count > 0; ++i)
        if (Empty(s[i])) { s[i] = Stack{}; s[i].item = item; const int k = std::min(count, cap); s[i].count = k; count -= k; }
    return count;
}
// (#1102) A fraction of an item (hundredths, 0-99) into these slots: onto a stack of it with room (100 hundredths make one
// more), else an empty slot. Returns the hundredths that did not fit.
int AddCenti(Stack* s, int n, int item, int centi) {
    if (centi <= 0) return 0;
    const int cap = tables.cap(item);
    for (int i = 0; i < n; ++i)
        if (s[i].item == item && Plain(s[i]) && !Empty(s[i]) && s[i].count < cap) {
            s[i].centi += centi;
            if (s[i].centi >= 100) { s[i].centi -= 100; ++s[i].count; }
            return 0;
        }
    for (int i = 0; i < n; ++i)
        if (Empty(s[i])) { s[i] = Stack{}; s[i].item = item; s[i].centi = centi; return 0; }
    return centi;
}
bool RoomFor(const Stack* s, int n, int item) { for (int i = 0; i < n; ++i) if (Empty(s[i]) || (s[i].item == item && Plain(s[i]) && s[i].count < tables.cap(item))) return true; return false; }
int CountItems(const Stack* s, int n, int item) { int c = 0; for (int i = 0; i < n; ++i) if (s[i].item == item && Plain(s[i])) c += s[i].count; return c; }
void TakeItems(Stack* s, int n, int item, int count) {
    for (int i = 0; i < n && count > 0; ++i)
        if (s[i].item == item && Plain(s[i])) { const int k = std::min(count, s[i].count); s[i].count -= k; count -= k; if (Empty(s[i])) s[i] = Stack{}; } // (#1102: a fraction keeps its item)
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

// (#1219) A unit an anvil repairs some composed copy with: any stone, or a part material's ingot (bronze_ingot, iron_ingot...).
bool IsRepairMaterial(int item) {
    if (item < 0) return false;
    const std::string n = tables.name(item);
    if (ItemMatches("any_stone", n)) return true;
    const size_t at = n.rfind("_ingot");
    return at != std::string::npos && at + 6 == n.size() && Composed_Material(n.substr(0, at)) != nullptr;
}
bool Accepts(PartType type, int item) {
    if (tables.fluid(item)) return false; // (#1254) fluid moves only through pipes and phials (MoveFluids), never as items
    if (Def(type).power == Power::Storage) return true;
    if (Def(type).waterPerTick > 0 && item == WaterBucketItem()) return true;
    if (type == PartType::Anvil && IsRepairMaterial(item)) return true; // (#1219) a unit to repair a composed copy with
    return IsFuel(Def(type), item) || RecipeWants(type, item);
}
// (#1427 L7c) An assembler set to disassemble takes what it can take apart (`disassemblable`), and nothing else.
bool Disassembles(int item) { return item >= 0 && item < (int)disassemblable.size() && disassemblable[item] && !tables.fluid(item); }
bool AcceptsIn(const Machine& m, int item) {
    if (m.type == PartType::MbHatch) return m.hatch == kHatchItemIn && item >= 0 && !tables.fluid(item); // (#1429 PR 4) anything a structure might use
    if (m.type == PartType::MbController) return false;
    return m.mode == 1 && m.type == PartType::LvAssembler ? Disassembles(item) : Accepts(m.type, item);
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
    if (toFuel || !AcceptsIn(m, item)) return count;
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
    // (#1620 C2b) The cables as a conduit graph: each cable cell is a centred 1x1 copper run (the design's conversion), so a
    // network's cables carry at most 32 EU/t a path. Runs join neighbouring cells; a machine's port is the face slot of a
    // cable cell touching it (only for the network the machine joined). Without the copper material: no graph, no capping.
    ps.powerGraph = ConduitGraph();
    for (Machine& m : ps.machines) m.powerNodes.clear();
    if (Conduit_Materials().find("copper")) {
        std::vector<Int3> cs = L.cables;
        std::sort(cs.begin(), cs.end(), [](const Int3& a, const Int3& b) { return a.x != b.x ? a.x < b.x : a.y != b.y ? a.y < b.y : a.z < b.z; });
        cs.erase(std::unique(cs.begin(), cs.end(), [](const Int3& a, const Int3& b) { return a.x == b.x && a.y == b.y && a.z == b.z; }), cs.end());
        std::unordered_set<Int3> cableCells(cs.begin(), cs.end());
        auto centre = [](const Int3& c) { return CInt3{ c.x * 8 + 4, c.y * 8 + 4, c.z * 8 + 4 }; };
        std::vector<ConduitRun> runs;
        auto addRun = [&](CInt3 a, CInt3 b) {
            ConduitRun r; r.id = (int)runs.size() + 1; r.kind = ConduitKind::Cable; r.material = "copper"; r.size = 1; r.points = { a, b };
            runs.push_back(std::move(r));
        };
        for (const Int3& c : cs) {
            if (cableCells.count({ c.x + 1, c.y, c.z })) addRun(centre(c), centre({ c.x + 1, c.y, c.z }));
            if (cableCells.count({ c.x, c.y + 1, c.z })) addRun(centre(c), centre({ c.x, c.y + 1, c.z }));
            if (cableCells.count({ c.x, c.y, c.z + 1 })) addRun(centre(c), centre({ c.x, c.y, c.z + 1 }));
        }
        std::vector<ConduitPortHit> hits;
        const Int3 dirs[6] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
        for (const Int3& c : cs) {
            const auto cn = ps.cableNet.find(c);
            if (cn == ps.cableNet.end()) continue;
            for (const Int3& d : dirs) {
                const auto mt = machineAt.find({ c.x + d.x, c.y + d.y, c.z + d.z });
                if (mt == machineAt.end() || ps.machines[mt->second].powerNet != cn->second) continue;
                const CInt3 cc = centre(c);
                const CInt3 f = { cc.x + d.x * 4, cc.y + d.y * 4, cc.z + d.z * 4 };
                addRun(cc, f);
                hits.push_back({ f, mt->second, 0 });
            }
        }
        ps.powerGraph = Conduit_BuildGraph(runs, Conduit_Materials(), hits);
        for (size_t i = 0; i < hits.size(); ++i) {
            const int n = ps.powerGraph.portNode[i];
            if (n < 0) continue;
            std::vector<int>& pn = ps.machines[hits[i].owner].powerNodes;
            if (std::find(pn.begin(), pn.end(), n) == pn.end()) pn.push_back(n);
        }
        // A machine touching the cables at several places passes power inside itself.
        for (const Machine& m : ps.machines)
            for (size_t k = 1; k < m.powerNodes.size(); ++k) {
                ConduitEdge e; e.a = m.powerNodes[0]; e.b = m.powerNodes[k]; e.runId = 0; e.capacity = 1LL << 40; e.length = 0;
                ps.powerGraph.edges.push_back(e);
            }
    }
    // A network runs at the highest voltage of its generators and batteries (lower machines refuse it).
    for (Net& net : ps.powerNets) {
        net.voltage = 0;
        for (int m : net.machines) {
            const MachineDef& d = Def(ps.machines[m].type);
            if (d.power == Power::Generator || d.power == Power::Battery || d.power == Power::Dynamo) net.voltage = std::max(net.voltage, d.volts);
        }
    }
    // (#1271 L3) Each steam engine draws from a boiler touching it (the lowest-numbered one).
    std::unordered_map<Int3, int> boilerAt;
    for (size_t i = 0; i < ps.machines.size(); ++i)
        if (Def(ps.machines[i].type).power == Power::Steam) for (const Int3& c : ps.machines[i].cells) boilerAt[c] = (int)i;
    const Int3 dirs[6] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
    for (Machine& m : ps.machines) {
        m.boiler = -1;
        if (Def(m.type).power != Power::Engine) continue;
        for (const Int3& c : m.cells)
            for (const Int3& d : dirs) {
                auto it = boilerAt.find({ c.x + d.x, c.y + d.y, c.z + d.z });
                if (it != boilerAt.end() && (m.boiler < 0 || it->second < m.boiler)) m.boiler = it->second;
            }
    }
    // (#1271 L2) Each bloomery and blast furnace takes air from bellows touching it (the lowest-numbered).
    std::unordered_map<Int3, int> bellowsAt;
    for (size_t i = 0; i < ps.machines.size(); ++i)
        if (ps.machines[i].type == PartType::Bellows) for (const Int3& c : ps.machines[i].cells) bellowsAt[c] = (int)i;
    for (Machine& m : ps.machines) {
        m.bellows = -1;
        if (m.type != PartType::Bloomery && m.type != PartType::BlastFurnace) continue;
        for (const Int3& c : m.cells)
            for (const Int3& d : dirs) {
                auto it = bellowsAt.find({ c.x + d.x, c.y + d.y, c.z + d.z });
                if (it != bellowsAt.end() && (m.bellows < 0 || it->second < m.bellows)) m.bellows = it->second;
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
        m.tier = i < L.machineTier.size() ? L.machineTier[i] : MachineTiers_PartTier(type, 0, Def(type).volts); // (#1427 L7a)
        m.hatch = i < L.machineHatch.size() ? L.machineHatch[i] : kHatchNone; // (#1429 PR 4)
        // A loaded save's contents arrive before the layout: apply them once.
        auto pend = pendingStates.find({ L.body, L.plot, part });
        if (pend != pendingStates.end()) {
            const auto& s = pend->second; // recipe progress energy fuelEU, then "in" item count pairs, then "out" pairs
            size_t k = 0;
            auto num = [&]() { return k < s.size() ? atoi(s[k++].c_str()) : 0; };
            m.recipe = num(); m.progress = num(); m.energy = num(); m.fuelEU = num();
            const int nin = num();
            for (int j = 0; j < nin && k + 1 < s.size(); ++j) {
                const std::string tok = s[k]; const size_t bar = tok.find('|'); // (#1446 2b) "name|record": one item with its own data
                const int item = tables.item(bar == std::string::npos ? tok : tok.substr(0, bar)); const double a = atof(s[k + 1].c_str()); k += 2;
                if (bar != std::string::npos) {
                    for (int q = 0; q < InSlots(m); ++q) if (Empty(m.in[q])) { m.in[q] = Stack{}; m.in[q].item = item; m.in[q].count = 1; m.in[q].data = tok.substr(bar + 1); break; }
                    continue;
                }
                const int c = (int)a, centi = (int)lround((a - c) * 100.0); AddItems(m.in.data(), InSlots(m), item, c); if (centi > 0) AddCenti(m.in.data(), InSlots(m), item, std::min(centi, 99)); // (#1102) "2.50": a fraction
            }
            const int nout = num();
            for (int j = 0; j < nout && k + 1 < s.size(); ++j) {
                const std::string tok = s[k]; const size_t bar = tok.find('|'); // (#1446 2c) "name|record" in an output too
                const int item = tables.item(bar == std::string::npos ? tok : tok.substr(0, bar)); const int c = atoi(s[k + 1].c_str()); k += 2;
                if (bar != std::string::npos) { for (int q = 0; q < kSlots; ++q) if (m.out[q].count <= 0) { m.out[q] = Stack{}; m.out[q].item = item; m.out[q].count = 1; m.out[q].data = tok.substr(bar + 1); break; } continue; }
                AddItems(m.out, kSlots, item, c);
            }
            // Then tagged extras: "f item count" (the fuel slot), "w litres", "h heat".
            while (k < s.size()) {
                const std::string tag = s[k++];
                if (tag == "f" && k + 1 < s.size()) { m.fuel.item = tables.item(s[k]); m.fuel.count = atoi(s[k + 1].c_str()); k += 2; }
                else if (tag == "w") m.water = num();
                else if (tag == "h") m.heat = num();
                else if (tag == "fire") m.fire = (uint8_t)std::clamp(num(), 0, 2); // (#1103)
                else if (tag == "mode") m.mode = (uint8_t)std::clamp(num(), 0, 1); // (#1271 L6c)
                else if (tag == "pick" && k < s.size()) m.pick = tables.item(s[k++]); // (#1271 L6c) a player's recipe choice
                else if (tag == "dis") { m.disSlot = std::clamp(num(), -1, kSlots - 1); m.progress = std::max(0, num()); m.recipe = -1; } // (#1427 L7c)
                else if (tag == "drng" && k < s.size()) m.disRng = (uint32_t)strtoul(s[k++].c_str(), nullptr, 10);
                else if (tag == "dc" && k + 1 < s.size()) { const int id = tables.item(s[k]); m.disCarry[id] = std::clamp(atoi(s[k + 1].c_str()), 0, 99); k += 2; }
                else if (tag == "sum") ps->loadedSummary = true; // (#1302) it was in the summary sim: its caps follow
                else if (tag == "cap" && k + 2 < s.size()) { // a cap, kept only if that recipe is still this machine's with that output
                    const int r = atoi(s[k].c_str()); const std::string out = s[k + 1]; const int c = atoi(s[k + 2].c_str()); k += 3;
                    if (r >= 0 && r < (int)tables.recipes.size() && tables.recipes[r].machine == type && !tables.recipes[r].out.empty() && tables.name(tables.recipes[r].out[0].item) == out) m.nearCap[r] = c;
                }
                else if (tag == "ecap" && k < s.size()) m.nearEuCap = atoll(s[k++].c_str());
                else if (tag == "mbjob" && k + 2 < s.size()) { // (#1429 PR 4) a controller's job: recipe (checked by its first output), progress
                    const int r = atoi(s[k].c_str()); const std::string out = s[k + 1]; const int p = atoi(s[k + 2].c_str()); k += 3;
                    if (r >= 0 && r < (int)tables.recipes.size() && !tables.recipes[r].out.empty() && tables.name(tables.recipes[r].out[0].item) == out) { Machine::MbJob j; j.recipe = r; j.progress = std::max(0, p); m.jobs.push_back(j); }
                }
            }
            // Saves from before the fuel slot (#194) kept fuel in the inputs: move fuel no recipe of this machine wants
            // into the fuel slot (#229 review), so an old boiler or furnace still burns it.
            if (HasFuelSlot(m.type))
                for (int j = 0; j < InSlots(m); ++j) {
                    Stack& st = m.in[j];
                    if (st.count <= 0 || !Plain(st) || !IsFuel(Def(m.type), st.item) || RecipeWants(m.type, st.item)) continue;
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
        // (#1254) fluid (essence) is not an item: it is lost with its vessel. (#1629) Every stack comes back exactly: a fraction
        // (0.5 dirt in a chest, #1102: it was skipped, and 2.5 came back as 2) and an item's own record (#1446).
        auto refund = [&](const Stack& s) {
            if (Empty(s) || tables.fluid(s.item)) return;
            ItemOut r; r.item = s.item; r.count = s.count; r.to = L.refundTo; r.centi = s.centi; r.data = s.data;
            refunds.push_back(r);
        };
        for (int i = 0; i < InSlots(old); ++i) refund(old.in[i]);
        for (int i = 0; i < kSlots; ++i) refund(old.out[i]);
        refund(old.fuel);
        for (const Machine::MbJob& j : old.jobs) // (#1429 PR 4) a controller's running recipes: the inputs they claimed
            if (j.recipe >= 0) for (const Stack& s : tables.recipes[j.recipe].in) refund(s);
    }
    ps->machines = std::move(next);
    // (#1429 PR 4) The formed structures: each controller's machine, parallel, tier and hatches (machine indices); the rest
    // unformed (they keep their jobs, which wait). The output hatches' reservations are their users' jobs' (a job whose
    // hatch is gone reserves again when it finishes).
    ps->hasStructures = !L.structures.empty();
    for (Machine& m : ps->machines) { m.mbFormed = false; m.mbHatches.clear(); m.reserved.clear(); }
    for (const LayoutMsg::Structure& s : L.structures) {
        const auto c = byPart.find(s.controller);
        if (c == byPart.end()) continue;
        Machine& m = ps->machines[c->second];
        if (m.type != PartType::MbController) continue;
        m.mbFormed = true; m.mbMachine = s.machine; m.mbParallel = std::max(1, s.parallel); m.tier = s.tier;
        for (int hp : s.hatches) { const auto h = byPart.find(hp); if (h != byPart.end() && ps->machines[h->second].type == PartType::MbHatch) m.mbHatches.push_back(h->second); }
    }
    for (Machine& m : ps->machines)
        for (Machine::MbJob& j : m.jobs) {
            std::vector<std::pair<int, Stack>> kept;
            for (auto& [hp, st] : j.res) {
                const auto h = byPart.find(hp);
                if (h == byPart.end() || ps->machines[h->second].hatch != kHatchItemOut) continue;
                Machine& hatch = ps->machines[h->second];
                hatch.reserved[st.item] += st.count;
                kept.push_back({ hp, st });
            }
            j.res = std::move(kept);
        }
    if (!ps->nearPlayer) BeginMeasuring(*ps, tickCount.load()); // (#1261) a layout change while far measures again
    ps->byPart = std::move(byPart);
    BuildNets(*ps, L);
    ps->rot.Set(L.rotation); // (#1269)
}

// (#1269) Rotation before power: rotation machines with a recipe load their network, then one solve per plot.
void RotationStep(PlotSim& ps, int dt) {
    // (#1271 L3) Boilers: one with steam engines burns for the engines whose networks have work (shared equally), the rest
    // warm up and bank their fire; each engine turns at the share of a full boiler's steam it gets.
    std::vector<int> wanting(ps.machines.size(), 0);
    for (const Machine& m : ps.machines)
        if (Def(m.type).power == Power::Engine && m.boiler >= 0 && ps.rot.NetHasLoad(m.part)) ++wanting[m.boiler];
    std::vector<long long> steamOf(ps.machines.size(), 0);
    for (size_t i = 0; i < ps.machines.size(); ++i) {
        Machine& b = ps.machines[i];
        const MachineDef& def = Def(b.type);
        if (def.power != Power::Steam) continue;
        steamOf[i] = GeneratorStep(b, def, dt, (long long)wanting[i] * Rotation::kEngineSteamPerTick * dt);
    }
    for (Machine& m : ps.machines) {
        if (Def(m.type).power != Power::Engine) continue;
        const bool wants = m.boiler >= 0 && ps.rot.NetHasLoad(m.part);
        const long long steam = wants ? steamOf[m.boiler] / std::max(1, wanting[m.boiler]) : 0;
        const int rpm = (int)(steam * Rotation::kSteamEngine.rpm / std::max(1, Rotation::kEngineSteamPerTick * dt));
        ps.rot.SetSourceRpm(m.part, rpm);
        m.status = rpm > 0 ? MachineStatus::Working : m.boiler < 0 ? MachineStatus::NoFuel : !wants ? MachineStatus::NoLoad
                 : ps.machines[m.boiler].status == MachineStatus::NoWater ? MachineStatus::NoWater : MachineStatus::NoFuel;
    }
    if (ps.rot.empty()) return;
    for (const Machine& m : ps.machines) {
        const Power p = Def(m.type).power;
        if (p == Power::Rotation) ps.rot.SetRunning(m.part, m.recipe >= 0);
        // (#1271 L3) a dynamo loads its shaft while its cable network wanted power last step
        else if (p == Power::Dynamo) ps.rot.SetRunning(m.part, m.powerNet >= 0 && m.powerNet < (int)ps.powerNets.size() && ps.powerNets[m.powerNet].demand > 0);
    }
    for (const Machine& m : ps.machines) // (#1271 L2) bellows blow (load their shaft) while their furnace has a recipe
        if (m.bellows >= 0 && m.bellows < (int)ps.machines.size() && m.recipe >= 0) ps.rot.SetRunning(ps.machines[m.bellows].part, true);
    ps.rot.Step(dt);
}

// (#1271 L2) A furnace's extra ticks of work from its bellows this step: up to one more a tick at the bellows' rated speed
// (doc: the bellows speed it up, not the yield), the fraction carried over in the bellows' own accumulator.
int BellowsTicks(PlotSim& ps, Machine& m, int dt) {
    if (m.bellows < 0 || m.bellows >= (int)ps.machines.size()) return 0;
    Machine& b = ps.machines[m.bellows];
    const int perm = std::min(1000, ps.rot.ProgressPermille(b.part)), acc = b.rotAcc + perm * dt;
    b.rotAcc = acc % 1000;
    b.status = perm > 0 ? MachineStatus::Working : MachineStatus::NoPower;
    return acc / 1000;
}

// (#1427 L7c) Disassembly: what taking a stack's item apart gives back (disassembly.h): a composed item by its record, else
// the index's entry. False if it is not taken apart.
static bool DisEntryOf(const Stack& s, DisassemblyEntry& e) {
    if (s.count <= 0 || s.item < 0 || tables.fluid(s.item)) return false;
    const std::string name = tables.name(s.item);
    if (!s.data.empty()) return Disassembly_Composed(name, s.data, e);
    const DisassemblyEntry* f = Disassembly_Find(name);
    if (!f) return false;
    e = *f;
    return true;
}
// Its forward recipe's overclock in this machine: from the item's tier, at least the assembler's first (a hand recipe or a
// tier-1 part is an LV job, not one 4 tiers under it)
static int DisSteps(const Machine& m, const DisassemblyEntry& e) { return std::max(0, MachineTiers_Steps(m.tier, std::max(e.tier, MachineTiers_First(m.type)))); }
static int JobEu(const Machine& m) {
    if (m.hatch == kHatchEnergy) return m.hatchDemand > 0 ? (int)std::min<long long>(m.hatchDemand, INT32_MAX) : -1; // (#1429 PR 4)
    if (m.type == PartType::MbController) return -1; // (its energy hatches draw for it)
    if (m.mode == 1) return m.disSlot >= 0 ? m.disEu : -1;
    return m.recipe >= 0 && m.progress < TicksOf(m, tables.recipes[m.recipe]) ? EuOf(m, tables.recipes[m.recipe]) : -1; // (#1847) a done craft waiting on full outputs draws nothing
}
// An assembler set to disassemble, one step: the first input it can take apart starts a job (the forward recipe's time and EU,
// overclocked); a finished job takes one item and gives back its parts at the tier yield, raw parts as hundredths carried in
// the machine, components rolled per unit. Every output must fit, or it waits (nothing rolled or taken).
// The first input it can take apart starts a job: its slot, time and EU. False: none.
static bool DisStartJob(Machine& m, DisassemblyEntry& e) {
    for (int si = 0; si < InSlots(m); ++si)
        if (DisEntryOf(m.in[si], e) && MachineTiers_YieldPercent(m.tier, e.tier) > 0) {
            const int st = DisSteps(m, e);
            m.disSlot = si; m.progress = 0; m.disTicks = MachineTiers_Ticks(e.ticks, st); m.disEu = MachineTiers_Eu(e.eu, st);
            return true;
        }
    return false;
}
static void DisassembleStep(PlotSim& ps, Machine& m, int dt, long long& made, int& active) {
    DisassemblyEntry e;
    if (m.disSlot >= 0 && (m.disSlot >= InSlots(m) || !DisEntryOf(m.in[m.disSlot], e))) { m.disSlot = -1; m.progress = 0; } // taken out meanwhile
    if (m.disSlot < 0 && !DisStartJob(m, e)) { m.status = MachineStatus::Idle; m.asleep = true; return; } // nothing it takes apart: until items change
    if (m.disTicks <= 0) { const int st = DisSteps(m, e); m.disTicks = MachineTiers_Ticks(e.ticks, st); m.disEu = MachineTiers_Eu(e.eu, st); } // (after a load)
    if (m.progress < m.disTicks) {
        if (m.powerNet < 0) { m.status = MachineStatus::NoPower; return; }
        if (VoltsOf(m) < ps.powerNets[m.powerNet].voltage) { m.status = MachineStatus::Overvoltage; return; }
        int ticks = dt;
        if (m.disEu > 0) { const int eu = m.received + m.euCarry; ticks = eu / m.disEu; m.euCarry = eu % m.disEu; }
        m.progress = std::min(m.disTicks, m.progress + ticks);
        m.status = ticks > 0 ? MachineStatus::Working : MachineStatus::NoPower;
        if (ticks > 0) ++active;
    }
    if (m.progress < m.disTicks) return;
    if (m.disRng == 0) m.disRng = (uint32_t)WorkingKey(ps.body, ps.plot, m.part) | 1u; // seeded by the machine: the same run every time
    uint32_t rng = m.disRng;
    std::vector<DisassemblyGive> give;
    Disassembly_Give(e, MachineTiers_YieldPercent(m.tier, e.tier), rng, give);
    Stack trial[kSlots];
    std::copy(m.out, m.out + kSlots, trial);
    std::map<int, int> carry = m.disCarry;
    std::vector<std::pair<int, int>> whole;
    for (const DisassemblyGive& g : give) {
        const int id = tables.item(g.item);
        const int c = carry[id] + g.centi;
        carry[id] = c % 100;
        if (c >= 100) {
            if (AddItems(trial, kSlots, id, c / 100) > 0) { m.status = MachineStatus::OutputFull; return; }
            whole.push_back({ id, c / 100 });
        }
    }
    for (auto it = carry.begin(); it != carry.end();) it = it->second == 0 ? carry.erase(it) : std::next(it);
    Stack& src = m.in[m.disSlot];
    ps.stats.items[src.item].used.add(1);
    if (--src.count <= 0) src = Stack{};
    std::copy(trial, trial + kSlots, m.out);
    for (const auto& [id, n] : whole) { made += n; ps.stats.items[id].made.add(n); }
    m.disCarry = std::move(carry);
    m.disRng = rng;
    m.disSlot = -1; m.progress = 0; m.disTicks = 0; m.euCarry = 0;
    // The next job now, so the next power split counts its demand (else one tick of "no power" between jobs)
    DisStartJob(m, e);
    m.dirty = true;
}

// One step of dt ticks for a plot.
// (#1219, user 2026-11-22) The anvil's own jobs on a composed copy in its inputs, hammered as a recipe is (a player's strikes):
//   repair: with a unit of its head's material beside it (its ingot; any stone for a stone head): its wear down by half its
//           durability (GUESS, Progression tunes the costs in T2), the unit used;
//   a gem out: without that material, the last gem set comes out beside it, or breaks (25%, GUESS).
// The copy (and the gem) go to empty output slots. True if a copy is on it (its job ran, or waits): no recipe this tick.
constexpr int kAnvilJobTicks = 80;           // two strikes (GUESS)
constexpr int kGemBreakPercent = 25;         // GUESS
bool RepairMaterial(const std::string& head, int item) {
    const std::string n = tables.name(item);
    if (head == "stone") return ItemMatches("any_stone", n);
    return n == head + "_ingot" || n == head;
}
bool AnvilJob(PlotSim& ps, Machine& m, int& active) {
    int di = -1;
    for (int si = 0; si < InSlots(m); ++si) if (m.in[si].count > 0 && !Plain(m.in[si])) { di = si; break; }
    if (di < 0) { m.anvilJob = 0; return false; }
    const std::string cls = tables.name(m.in[di].item);
    const ComposedStats st = Composed_Stats(cls, m.in[di].data);
    int mi = -1;
    if (st.valid && st.wear > 0)
        for (int si = 0; si < InSlots(m); ++si) if (si != di && m.in[si].count > 0 && Plain(m.in[si]) && RepairMaterial(st.head, m.in[si].item)) { mi = si; break; }
    const int job = mi >= 0 ? 1 : (!st.gems.empty() ? 2 : 0);
    if (job == 0) { m.anvilJob = 0; m.progress = 0; m.status = MachineStatus::Idle; return true; } // it waits to be taken back
    if (m.anvilJob != job) { m.anvilJob = (uint8_t)job; m.progress = 0; }
    const int ticks = std::min(m.strike, kAnvilJobTicks - m.progress);
    m.strike -= ticks; m.progress += ticks;
    if (m.progress < kAnvilJobTicks) { m.status = ticks > 0 ? MachineStatus::Working : MachineStatus::Idle; if (ticks > 0) ++active; return true; }
    int q1 = -1, q2 = -1;
    for (int q = 0; q < kSlots; ++q) if (m.out[q].count <= 0) { if (q1 < 0) q1 = q; else if (q2 < 0) q2 = q; }
    if (q1 < 0 || (job == 2 && q2 < 0)) { m.status = MachineStatus::OutputFull; return true; }
    std::string d = m.in[di].data;
    if (job == 1) {
        d = ItemData::SetInt(d, "w", std::max(0, st.wear - st.durability / 2));
        if (--m.in[mi].count <= 0 && m.in[mi].centi <= 0) m.in[mi] = Stack{};
    } else {
        std::string j;
        for (size_t g = 0; g + 1 < st.gems.size(); ++g) j += (j.empty() ? "" : ",") + st.gems[g];
        d = ItemData::Set(d, "j", j);
        // A roll from the plot, the part and the tick (deterministic: tests and replays agree)
        const uint32_t h = (uint32_t)(ps.plot * 2654435761u) ^ (uint32_t)(m.part * 40503u) ^ (uint32_t)(tickCount.load() * 2246822519u);
        if ((h >> 7) % 100 >= (uint32_t)kGemBreakPercent) { m.out[q2] = Stack{}; m.out[q2].item = tables.item(st.gems.back()); m.out[q2].count = 1; }
    }
    m.out[q1] = Stack{}; m.out[q1].item = m.in[di].item; m.out[q1].count = 1; m.out[q1].data = d;
    m.in[di] = Stack{};
    m.anvilJob = 0; m.progress = 0;
    m.status = MachineStatus::Working; ++active;
    m.dirty = true;
    return true;
}

// ---- (#1429 PR 4) Multiblocks: a formed structure is one machine, keyed by its controller ----
// It runs its machine's recipes (multiblocks.txt `machine`), up to its `parallel` at once, at its tier (overclocked as a single
// machine). A recipe claims all its inputs from the pooled item input hatches when it starts (never half-claimed) and
// reserves its outputs' room in the output hatches then, so a running recipe always finishes (the user, 2026-11-23). An
// electric one draws through its energy hatches: a shortfall slows its jobs evenly, and a shared hatch splits what it gets
// by its users' demand. Shared hatches: the controllers take turns to start, one recipe each a round, the one after the
// last to start going first (machine index order, the turn kept in the plot: deterministic).
static bool MbElectric(const Machine& c) { return Def(c.mbMachine).power == Power::Eu; }
static long long MbNeed(const Machine& c) { // EU a tick its unfinished jobs ask
    long long n = 0;
    for (const Machine::MbJob& j : c.jobs) { const Recipe& r = tables.recipes[j.recipe]; if (j.progress < TicksOf(c, r)) n += EuOf(c, r); }
    return n;
}
static void MbHatchesOf(const PlotSim& ps, const Machine& c, uint8_t kind, std::vector<int>& out) {
    out.clear();
    for (int h : c.mbHatches) if (ps.machines[h].hatch == kind) out.push_back(h);
}
// Its share of each energy hatch's demand: its need split evenly over its energy hatches (the remainder to the first).
static long long MbShareAt(const Machine& c, const std::vector<int>& energy, size_t k) {
    const long long need = MbNeed(c), n = (long long)energy.size();
    return need / n + ((long long)k < need % n ? 1 : 0);
}
// Before the power step: what each energy hatch asks its cables for.
void MbDemand(PlotSim& ps) {
    std::vector<int> eh;
    for (Machine& m : ps.machines) if (m.hatch == kHatchEnergy) { m.hatchDemand = 0; m.status = MachineStatus::Idle; }
    for (Machine& c : ps.machines) {
        c.mbEu = 0;
        if (c.type != PartType::MbController || !c.mbFormed || !MbElectric(c)) continue;
        MbHatchesOf(ps, c, kHatchEnergy, eh);
        for (size_t k = 0; k < eh.size(); ++k) ps.machines[eh[k]].hatchDemand += MbShareAt(c, eh, k);
    }
}
// After it: each user's share of what its energy hatches received (by its part of each hatch's demand).
void MbTakePower(PlotSim& ps) {
    std::vector<int> eh;
    for (Machine& c : ps.machines) {
        if (c.type != PartType::MbController || !c.mbFormed || !MbElectric(c)) continue;
        MbHatchesOf(ps, c, kHatchEnergy, eh);
        for (size_t k = 0; k < eh.size(); ++k) {
            const Machine& h = ps.machines[eh[k]];
            if (h.hatchDemand > 0 && h.received > 0) c.mbEu += (long long)h.received * MbShareAt(c, eh, k) / h.hatchDemand;
        }
    }
}
// Reserve room for recipe r's outputs in c's output hatches (an output may spread over several), against what is in them and
// already reserved. False (nothing reserved) when it does not all fit.
static bool MbReserve(PlotSim& ps, const Machine& c, const Recipe& r, Machine::MbJob& job) {
    std::vector<int> outs;
    MbHatchesOf(ps, c, kHatchItemOut, outs);
    std::map<int, std::map<int, int>> res; // hatch machine -> its reservations with this job's
    std::vector<std::pair<int, Stack>> plan;
    for (const Stack& s : r.out) {
        int left = s.count;
        for (int h : outs) {
            if (left <= 0) break;
            auto rh = res.find(h);
            if (rh == res.end()) rh = res.emplace(h, ps.machines[h].reserved).first;
            Stack trial[kSlots];
            std::copy(ps.machines[h].out, ps.machines[h].out + kSlots, trial);
            bool over = false;
            for (const auto& [item, n] : rh->second) over |= AddItems(trial, kSlots, item, n) > 0;
            if (over) continue;
            const int after = AddItems(trial, kSlots, s.item, left);
            if (after == left) continue;
            Stack put; put.item = s.item; put.count = left - after;
            rh->second[s.item] += put.count;
            plan.push_back({ ps.machines[h].part, put });
            left = after;
        }
        if (left > 0) return false;
    }
    for (auto& [h, m] : res) ps.machines[h].reserved = m;
    job.res = std::move(plan);
    return true;
}
// Start one recipe in c if its pooled inputs allow it and its outputs' room can be reserved. outFull: one had its inputs but
// no room.
static bool MbStart(PlotSim& ps, Machine& c, bool& outFull) {
    std::vector<int> ins;
    MbHatchesOf(ps, c, kHatchItemIn, ins);
    auto pooled = [&](int item) { int n = 0; for (int h : ins) n += CountItems(ps.machines[h].in.data(), kSlots, item); return n; };
    for (int ri : recipesOf[(int)c.mbMachine]) {
        const Recipe& r = tables.recipes[ri];
        if (MachineTiers_Steps(c.tier, r.tier) < 0 || RecipeLocked(r) || r.in.empty()) continue;
        bool have = true;
        for (const Stack& s : r.in) have &= pooled(s.item) >= s.count;
        if (!have) continue;
        Machine::MbJob job; job.recipe = ri;
        if (!MbReserve(ps, c, r, job)) { outFull = true; continue; }
        for (const Stack& s : r.in) { // claimed now, from the hatches in order
            int need = s.count;
            for (int h : ins) {
                if (need <= 0) break;
                Machine& hm = ps.machines[h];
                const int k = std::min(need, CountItems(hm.in.data(), kSlots, s.item));
                if (k <= 0) continue;
                TakeItems(hm.in.data(), kSlots, s.item, k);
                need -= k; hm.dirty = true;
            }
        }
        c.jobs.push_back(std::move(job));
        return true;
    }
    return false;
}
// The structures' work this step: jobs progress and finish into their reserved room, then the controllers take turns to start.
void MbStep(PlotSim& ps, int dt, long long& made, int& active) {
    std::vector<int> ctl;
    for (size_t i = 0; i < ps.machines.size(); ++i) if (ps.machines[i].type == PartType::MbController) ctl.push_back((int)i);
    if (ctl.empty()) return;
    std::vector<bool> worked(ctl.size(), false), full(ctl.size(), false);
    for (size_t q = 0; q < ctl.size(); ++q) {
        Machine& c = ps.machines[ctl[q]];
        if (!c.mbFormed) continue;
        int unfinished = 0;
        for (const Machine::MbJob& j : c.jobs) unfinished += j.progress < TicksOf(c, tables.recipes[j.recipe]);
        int k = 0;
        for (Machine::MbJob& j : c.jobs) {
            const Recipe& r = tables.recipes[j.recipe];
            const int rTicks = TicksOf(c, r), rEu = EuOf(c, r);
            if (j.progress >= rTicks) continue;
            int ticks = dt;
            if (MbElectric(c) && rEu > 0) { // its even share of the EU its hatches passed it (the remainder to the first)
                const long long share = c.mbEu / unfinished + (k == 0 ? c.mbEu % unfinished : 0);
                ++k;
                const long long eu = share + j.euCarry;
                ticks = (int)std::min<long long>(eu / rEu, rTicks);
                j.euCarry = (int)(eu % rEu);
                if (share > 0) worked[q] = true; // (EU short of a whole tick still works: slowed, not stopped)
            }
            j.progress = std::min(rTicks, j.progress + ticks);
            if (ticks > 0) worked[q] = true;
        }
        for (size_t ji = 0; ji < c.jobs.size();) {
            Machine::MbJob& j = c.jobs[ji];
            const Recipe& r = tables.recipes[j.recipe];
            if (j.progress < TicksOf(c, r)) { ++ji; continue; }
            // A job loaded from a save, or whose hatch was removed, reserves its room now (or waits for it).
            bool held = !j.res.empty();
            for (const auto& [hp, st] : j.res) { const auto h = ps.byPart.find(hp); held &= h != ps.byPart.end() && ps.machines[h->second].hatch == kHatchItemOut; }
            if (!held) {
                for (const auto& [hp, st] : j.res) { const auto h = ps.byPart.find(hp); if (h != ps.byPart.end()) { auto& rv = ps.machines[h->second].reserved; if ((rv[st.item] -= st.count) <= 0) rv.erase(st.item); } }
                j.res.clear();
                if (!MbReserve(ps, c, r, j)) { full[q] = true; ++ji; continue; }
            }
            for (const auto& [hp, st] : j.res) {
                Machine& h = ps.machines[ps.byPart.at(hp)];
                if ((h.reserved[st.item] -= st.count) <= 0) h.reserved.erase(st.item);
                AddItems(h.out, kSlots, st.item, st.count);
                h.dirty = true;
            }
            for (const Stack& s : r.out) { made += s.count; ps.stats.items[s.item].made.add(s.count); if (ps.completed.size() < 4096) ps.completed.emplace_back(s.item, s.count); }
            for (const Stack& s : r.in) ps.stats.items[s.item].used.add(s.count);
            ++c.crafts[j.recipe];
            c.jobs.erase(c.jobs.begin() + (std::ptrdiff_t)ji);
        }
    }
    // Starts: rounds of one recipe per controller, from the one after the last to start.
    for (bool any = true; any;) {
        any = false;
        const int n = (int)ctl.size(), first = ps.mbRR % n;
        for (int s = 0; s < n; ++s) {
            const int q = (first + s) % n;
            Machine& c = ps.machines[ctl[q]];
            if (!c.mbFormed || (int)c.jobs.size() >= c.mbParallel) continue;
            bool outFull = false;
            if (MbStart(ps, c, outFull)) { any = true; ps.mbRR = (q + 1) % n; }
            else if (outFull) full[q] = true;
        }
    }
    for (size_t q = 0; q < ctl.size(); ++q) {
        Machine& c = ps.machines[ctl[q]];
        c.status = !c.mbFormed ? MachineStatus::Idle : worked[q] ? MachineStatus::Working : full[q] ? MachineStatus::OutputFull
                 : MbElectric(c) && !c.jobs.empty() ? MachineStatus::NoPower : MachineStatus::Idle;
        if (worked[q]) ++active;
    }
}

void TickPlot(PlotSim& ps, int dt, long long& made, int& active) {
    gCount = &ps.stats; // (#1199) BurnFuel counts into this plot
    MbDemand(ps); // (#1429 PR 4) what the energy hatches ask
    RotationStep(ps, dt); // (#1269)
    // Power: generators burn fuel into EU at their voltage; batteries buffer; consumers on a network
    // above their voltage refuse it.
    for (Net& net : ps.powerNets) {
        long long supply = 0, demand = 0, batteryAvail = 0, batteryRoom = 0;
        for (int mi : net.machines) {
            Machine& m = ps.machines[mi];
            const MachineDef& def = Def(m.type);
            m.received = 0; m.reach = -1; m.canSend = -1;
            if (def.power == Power::Battery) {
                batteryAvail += std::min(m.energy, kBatteryRate * dt);
                batteryRoom += std::min(kBatteryCapacity - m.energy, kBatteryRate * dt);
            } else if (PowerOf(m) == Power::Eu) {
                if (VoltsOf(m) < net.voltage) { m.status = MachineStatus::Overvoltage; continue; }
                if (JobEu(m) >= 0) demand += (long long)JobEu(m) * dt;
            } else if (def.power == Power::Motor) { // (#1271 L3) a motor draws while its shaft network has work
                if (VoltsOf(m) < net.voltage) { m.status = MachineStatus::Overvoltage; continue; }
                if (ps.rot.NetHasLoad(m.part)) demand += (long long)Rotation::kMotorEuPerTick * dt;
            }
        }
        // (#1620 C2b) What the cables can carry: sources (generators, dynamos, batteries) to sinks (the consumers that want
        // power), a max-flow over the plot's cable graph. A sink's reach and a source's canSend cap them below; batteries'
        // charging is not capped (only what they give out counts as a source).
        long long capped = demand;
        if (!ps.powerGraph.nodes.empty()) {
            std::vector<ConduitFlowNode> sources, sinks;
            std::vector<int> srcMachine, sinkMachine;
            for (int mi : net.machines) {
                Machine& m = ps.machines[mi];
                if (m.powerNodes.empty()) continue;
                const MachineDef& def = Def(m.type);
                const int node = m.powerNodes[0];
                if (def.power == Power::Generator) { // (out of fuel, nothing banked: it can't send, so it mustn't take the cable's path)
                    const bool empty = m.energy <= 0 && m.fuel.count <= 0 && m.fuelEU <= 0;
                    sources.push_back({ node, empty ? 0LL : (long long)def.euOut * dt }); srcMachine.push_back(mi);
                }
                else if (def.power == Power::Dynamo) {
                    const long long can = (long long)Rotation::kDynamoEuPerTick * dt * std::min(1000, ps.rot.ProgressPermille(m.part)) / 1000;
                    sources.push_back({ node, can }); srcMachine.push_back(mi);
                } else if (def.power == Power::Battery) sources.push_back({ node, (long long)std::min(m.energy, kBatteryRate * dt) }), srcMachine.push_back(-1);
                else if (PowerOf(m) == Power::Eu && VoltsOf(m) >= net.voltage && JobEu(m) >= 0) { sinks.push_back({ node, (long long)JobEu(m) * dt }); sinkMachine.push_back(mi); }
                else if (def.power == Power::Motor && VoltsOf(m) >= net.voltage && ps.rot.NetHasLoad(m.part)) { sinks.push_back({ node, (long long)Rotation::kMotorEuPerTick * dt }); sinkMachine.push_back(mi); }
            }
            std::vector<int64_t> srcFlow;
            std::vector<int64_t> reach = Conduit_MaxFlow(ps.powerGraph, sources, sinks, dt, &srcFlow);
            long long want = 0, got = 0;
            for (size_t i = 0; i < sinks.size(); ++i) { want += sinks[i].amount; got += reach[i]; }
            if (got < want) {
                // A max-flow fills sinks in whatever order its paths come: some machines would run full, others get nothing.
                // The rule is that a shortfall slows everyone evenly: the largest equal share of every sink's need (in 1/1024)
                // the cables can deliver in full, found by halving (10 solves, only on a shortfall).
                std::vector<ConduitFlowNode> scaled = sinks;
                std::vector<int64_t> bestFlow(sources.size(), 0);
                std::vector<int64_t> bestReach(sinks.size(), 0);
                int lo = 0, hi = 1024;
                while (lo < hi) {
                    const int mid = (lo + hi + 1) / 2;
                    long long ask = 0, gave = 0;
                    for (size_t i = 0; i < sinks.size(); ++i) { scaled[i].amount = sinks[i].amount * mid / 1024; ask += scaled[i].amount; }
                    std::vector<int64_t> f;
                    const std::vector<int64_t> r = Conduit_MaxFlow(ps.powerGraph, sources, scaled, dt, &f);
                    for (int64_t v : r) gave += v;
                    if (gave >= ask) { lo = mid; bestFlow = f; bestReach = r; } else hi = mid - 1;
                }
                reach = bestReach;
                srcFlow = bestFlow;
            }
            capped = 0;
            for (size_t i = 0; i < sinks.size(); ++i) { ps.machines[sinkMachine[i]].reach = reach[i]; capped += reach[i]; }
            for (size_t i = 0; i < sources.size(); ++i) if (srcMachine[i] >= 0) ps.machines[srcMachine[i]].canSend = srcFlow[i];
            // (sinks without a cable node stay uncapped: add what they want back)
            for (int mi : net.machines) {
                const Machine& m = ps.machines[mi];
                if (!m.powerNodes.empty()) continue;
                const MachineDef& def = Def(m.type);
                if (PowerOf(m) == Power::Eu && VoltsOf(m) >= net.voltage && JobEu(m) >= 0) capped += (long long)JobEu(m) * dt;
                else if (def.power == Power::Motor && VoltsOf(m) >= net.voltage && ps.rot.NetHasLoad(m.part)) capped += (long long)Rotation::kMotorEuPerTick * dt;
            }
        }
        // Generators burn for the consumers and the batteries' room, in a fixed order (#229: not whatever the load).
        for (int mi : net.machines) {
            Machine& m = ps.machines[mi];
            const MachineDef& def = Def(m.type);
            if (def.power == Power::Dynamo) { // (#1271 L3) EU from its shaft's speed (full at its rated speed), as much as is wanted
                const long long can = (long long)Rotation::kDynamoEuPerTick * dt * std::min(1000, ps.rot.ProgressPermille(m.part)) / 1000;
                long long eu = std::min(can, std::max(0LL, capped + batteryRoom - supply));
                if (m.canSend >= 0) eu = std::min(eu, m.canSend);
                supply += eu;
                if (eu > 0) ps.stats.euMade[(int)m.type].add(eu);
                m.status = eu > 0 ? MachineStatus::Working : can > 0 ? MachineStatus::NoLoad : MachineStatus::NoPower;
                continue;
            }
            if (def.power != Power::Generator) continue;
            long long need = std::max(0LL, capped + batteryRoom - supply);
            if (m.canSend >= 0) need = std::min(need, m.canSend); // (#1620 C2b) it burns only for what its cables carry away
            const long long eu = GeneratorStep(m, def, dt, need);
            supply += eu;
            if (eu > 0) ps.stats.euMade[(int)m.type].add(eu); // (#1199)
            m.euWindow += eu; if (ps.mode == PlotSim::Mode::Near) m.nearEuCur += eu; // (#1261)
        }
        ps.stats.euDemand.add(demand);
        // Share: generation first, then the batteries; a shortfall slows every consumer equally.
        const long long fromBatteries = std::max(0LL, std::min(batteryAvail, capped - supply));
        const long long avail = supply + fromBatteries;
        for (int mi : net.machines) {
            Machine& m = ps.machines[mi];
            const MachineDef& def = Def(m.type);
            if (PowerOf(m) == Power::Eu && JobEu(m) >= 0 && VoltsOf(m) >= net.voltage) {
                const long long need = (long long)JobEu(m) * dt;
                const long long cap = m.reach >= 0 ? std::min(need, m.reach) : need; // (#1620 C2b) what its cables bring it
                m.received = (int)(capped <= avail ? cap : cap * avail / std::max(1LL, capped));
                if (m.received > 0) ps.stats.euUsed[(int)m.type].add(m.received); // (#1199)
            } else if (def.power == Power::Motor && VoltsOf(m) >= net.voltage && ps.rot.NetHasLoad(m.part)) { // (#1271 L3) its share
                const long long need = (long long)Rotation::kMotorEuPerTick * dt;
                const long long cap = m.reach >= 0 ? std::min(need, m.reach) : need;
                m.received = (int)(capped <= avail ? cap : cap * avail / std::max(1LL, capped));
                if (m.received > 0) ps.stats.euUsed[(int)m.type].add(m.received);
            }
        }
        net.supply = supply + fromBatteries; net.demand = demand; net.carried = capped; net.ticks = std::max(1, dt); // for the cable tooltip (#570)
        long long surplus = supply - std::min(supply, capped);
        long long drain = fromBatteries;
        for (int mi : net.machines) {
            Machine& m = ps.machines[mi];
            if (Def(m.type).power != Power::Battery) continue;
            if (drain > 0) { const long long k = std::min<long long>(drain, std::min(m.energy, kBatteryRate * dt)); m.energy -= (int)k; drain -= k; }
            else if (surplus > 0) { const long long k = std::min<long long>(surplus, std::min(kBatteryCapacity - m.energy, kBatteryRate * dt)); m.energy += (int)k; surplus -= k; }
            m.status = m.energy > 0 ? MachineStatus::Working : MachineStatus::Idle;
        }
    }
    MbTakePower(ps); // (#1429 PR 4) the structures' share of their energy hatches'
    MbStep(ps, dt, made, active);
    // Machines: start recipes, progress (EU received, fuel burnt, or by themselves), finish into
    // the outputs.
    for (Machine& m : ps.machines) {
        if (m.type == PartType::MbController || m.type == PartType::MbHatch) { // (#1429 PR 4) stepped above; a hatch is a buffer
            if (m.hatch == kHatchEnergy && m.status != MachineStatus::Overvoltage) m.status = m.received > 0 ? MachineStatus::Working : MachineStatus::Idle;
            continue;
        }
        const MachineDef& def = Def(m.type);
        if (def.power == Power::Generator || def.power == Power::Battery || def.power == Power::Storage) {
            // A generator on no network still lights: a boiler warms up and banks its fire (#229).
            if (def.power == Power::Generator && m.powerNet < 0) GeneratorStep(m, def, dt, 0);
            if (def.power == Power::Generator && m.status == MachineStatus::Working) ++active;
            continue;
        }
        if (def.power == Power::Steam || def.power == Power::Engine || def.power == Power::Dynamo) { // (#1271 L3) stepped above
            if (m.status == MachineStatus::Working) ++active;
            continue;
        }
        if (def.power == Power::Motor) { // (#1271 L3) its shaft turns at the share of full power it got (next tick's solve)
            const int rpm = (int)((long long)m.received * Rotation::kMotor.rpm / std::max(1, Rotation::kMotorEuPerTick * dt));
            ps.rot.SetSourceRpm(m.part, rpm);
            m.status = rpm > 0 ? MachineStatus::Working : m.powerNet < 0 ? MachineStatus::NoPower : ps.rot.NetHasLoad(m.part) ? MachineStatus::NoPower : MachineStatus::NoLoad;
            if (rpm > 0) ++active;
            if (m.powerNet < 0) m.received = 0;
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
        if (m.mode == 1) { m.recipe = -1; DisassembleStep(ps, m, dt, made, active); continue; } // (#1427 L7c) set to disassemble
        if (m.type == PartType::Anvil && m.recipe < 0 && AnvilJob(ps, m, active)) continue; // (#1219) a composed copy on it: repair / a gem out
        if (m.recipe < 0) {
            auto usable = [&](const Recipe& r) {
                return r.machine == m.type && (m.variant < 0 || (!r.out.empty() && r.out[0].item == m.variant)) && (m.pick < 0 || (!r.out.empty() && r.out[0].item == m.pick)) && MachineTiers_Steps(m.tier, r.tier) >= 0 && inputsHere(r) && !RecipeLocked(r); // (#1271 L6b) unlocked by research
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
        const int rTicks = TicksOf(m, r), rEu = EuOf(m, r); // (#1427 L7a) overclocked under its machine's tier
        if (m.progress < rTicks) {
            int ticks = 0;
            if (def.power == Power::Eu) {
                if (m.powerNet < 0) { m.status = MachineStatus::NoPower; continue; }
                if (VoltsOf(m) < ps.powerNets[m.powerNet].voltage) { m.status = MachineStatus::Overvoltage; continue; }
                // (#1266) the EU under one tick's worth waits for the next step (an overloaded net slows its machines, it does not
                // stall them: a 5-of-8 EU share made 0 ticks every tick near the player)
                if (rEu > 0) { const int eu = m.received + m.euCarry; ticks = eu / rEu; m.euCarry = eu % rEu; } else ticks = dt;
            } else if (def.power == Power::Manual) { // (#1271 L2) as far as the player's strikes go
                ticks = std::min(m.strike, rTicks - m.progress);
                m.strike -= ticks;
                if (ticks == 0) { m.status = MachineStatus::Idle; continue; }
            } else if (def.power == Power::Fuel) {
                ticks = BurnFuel(m, def, std::min(dt + BellowsTicks(ps, m, dt), rTicks - m.progress)); // (#1271 L2) bellows: up to twice as fast
            } else if (def.power == Power::Rotation) { // (#1269) at the shaft's speed: permille of a tick per tick, carried over
                const int rate = ps.rot.ProgressPermille(m.part), acc = m.rotAcc + rate * dt;
                ticks = acc / 1000; m.rotAcc = acc % 1000;
                if (rate > 0 && ticks == 0) { m.status = MachineStatus::Working; ++active; continue; } // turning, under a tick yet
            } else {
                ticks = dt + BellowsTicks(ps, m, dt); // runs by itself (coke oven, bloomery; #1271 L2: with bellows up to twice as fast)
            }
            m.progress = std::min(rTicks, m.progress + ticks);
            m.status = ticks > 0 ? MachineStatus::Working : def.power == Power::Fuel ? MachineStatus::NoFuel : MachineStatus::NoPower;
            if (ticks > 0) ++active;
        }
        if (m.progress >= rTicks) {
            // Finish only if every output fits.
            Stack trial[kSlots];
            std::copy(m.out, m.out + kSlots, trial);
            bool fits = true;
            for (const Stack& s : r.out) {
                if (Composed_Is(tables.name(s.item))) { // (#1446 2c) a composed output: one copy an empty slot, its record from the parts used
                    std::vector<std::string> parts; for (const Stack& p : r.in) parts.push_back(tables.name(p.item));
                    const std::string data = Composed_FromParts(tables.name(s.item), parts, m.quality);
                    if (!data.empty()) {
                        for (int u = 0; u < std::max(1, s.count) && fits; ++u) {
                            int q = 0; while (q < kSlots && trial[q].count > 0) ++q;
                            if (q == kSlots) { fits = false; break; }
                            trial[q] = Stack{}; trial[q].item = s.item; trial[q].count = 1; trial[q].data = data;
                        }
                        if (!fits) break;
                        continue;
                    }
                }
                if (AddItems(trial, kSlots, s.item, s.count) > 0) { fits = false; break; }
            }
            std::map<int, int> carry; // (#1734) outputs by chance: the machine's carried hundredths, kept only if all fits
            std::vector<std::pair<int, int>> lucky;
            if (fits && !r.chance.empty()) {
                carry = m.disCarry;
                for (const Stack& s : r.chance) {
                    const int c = carry[s.item] + s.count;
                    carry[s.item] = c % 100;
                    if (c >= 100) { lucky.push_back({ s.item, c / 100 }); if (AddItems(trial, kSlots, s.item, c / 100) > 0) { fits = false; break; } }
                }
            }
            if (!fits) { m.status = MachineStatus::OutputFull; continue; }
            if (!r.chance.empty()) { for (auto it = carry.begin(); it != carry.end();) it = it->second == 0 ? carry.erase(it) : std::next(it); m.disCarry = std::move(carry); }
            for (const auto& [id, n] : lucky) { made += n; ps.stats.items[id].made.add(n); if (ps.completed.size() < 4096) ps.completed.emplace_back(id, n); }
            m.quality = 100; // (#1446 2c) used
            for (const Stack& s : r.in) TakeItems(m.in.data(), kSlots, s.item, s.count); // consumed now, at the end (#193)
            std::copy(trial, trial + kSlots, m.out);
            for (const Stack& s : r.out) { made += s.count; ps.stats.items[s.item].made.add(s.count); if (ps.completed.size() < 4096) ps.completed.emplace_back(s.item, s.count); } // (#1199; #1728 completed: counts for the owner's quests)
            for (const Stack& s : r.in) ps.stats.items[s.item].used.add(s.count);
            ++m.crafts[m.recipe]; // (#1261) the measuring window
            if (ps.mode == PlotSim::Mode::Near) ++m.nearCur[m.recipe]; // (and the near minute)
            m.euCarry = 0; // (#1266)
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
                    if (from[i].count <= 0 || !Plain(from[i])) continue; // (#1446) conveyors never move an item with data (TakeItems skips it: it would copy)
                    const int item = from[i].item;
                    for (int d = 1; d < n && !moved; ++d) {
                        Machine& dst = ps.machines[net.machines[(net.rr + s + d) % n]];
                        if (&dst == &src || (chest && dst.type == PartType::Chest) || !AcceptsIn(dst, item)) continue;
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

// ---- (#1261) The summary sim (the factory doc, 2.2-2.5): a far plot advanced by rates, never above the full sim ----
// Every machine runs the recipes it ran while measuring, at most at its measured cap; a monotone fixpoint (rates only
// fall) lowers them for power, inputs, output room and the conveyor budget; time advances in segments that end at the
// next buffer event (a stock empties, a room fills). Buffers: chests, and each machine's own input, fuel and output
// slots. Integer fixed point, fixed orders, every rounding against the summary.
constexpr long long kFix = 1000000;                         // rates: millionths per tick
constexpr int kSolveRounds = 64;                            // (guess)
constexpr int kMaxSegments = 32;                            // (guess) events per advance; time beyond them makes nothing
constexpr long long kAdvanceEvery = 10 * kTicksPerSecond;   // (guess) an advance every 10 s per plot (staggered)

struct SumNode { int mi, recipe; long long rate; int fuelItem = -1; }; // crafts per tick x kFix; a fuel machine's fuel
struct SumGen { int mi; long long eu; long long scale = kFix; int fuelItem = -1; long long drawn = 0; }; // EU/tick x kFix it may give

long long Scaled(long long v, long long num, long long den) { return den <= 0 ? 0 : (long long)((__int128)v * num / den); }

void AdvanceSummary(PlotSim& ps, long long t0, long long t1) {
    if (t1 <= t0) return;
    const int nm = (int)ps.machines.size(), nets = (int)ps.itemNets.size() + nm; // a machine on no conveyor is a net of its own
    auto netOf = [&](int mi) { const int n = ps.machines[mi].itemNet; return n >= 0 && n < (int)ps.itemNets.size() ? n : (int)ps.itemNets.size() + mi; };
    auto fuelOf = [&](int fuel) { auto it = fuelTicks.find(fuel); return it != fuelTicks.end() ? it->second : 200; };
    long long t = t0, logCrafts = 0;
    int segment = 0;
    for (; segment < kMaxSegments && t < t1; ++segment) {
        // Stocks and rooms now: per net and item (chests), per machine and item (own slots).
        std::map<std::pair<int, int>, long long> chestStock, chestRoom;
        std::vector<std::map<int, long long>> own(nm), ownRoom(nm);
        for (int mi = 0; mi < nm; ++mi) {
            const Machine& m = ps.machines[mi];
            if (m.type == PartType::Chest) {
                for (int i = 0; i < kChestSlots; ++i) if (m.in[i].count > 0 && Plain(m.in[i])) chestStock[{ netOf(mi), m.in[i].item }] += m.in[i].count; // (#1446) never an item with data
                continue;
            }
            for (int i = 0; i < InSlots(m); ++i) if (m.in[i].count > 0 && Plain(m.in[i])) own[mi][m.in[i].item] += m.in[i].count;
            if (m.fuel.count > 0) own[mi][m.fuel.item] += m.fuel.count;
        }
        auto roomFor = [&](const Stack* s, int n, int item) { long long r = 0; for (int i = 0; i < n; ++i) r += s[i].count <= 0 ? tables.cap(item) : s[i].item == item && Plain(s[i]) ? std::max(0, tables.cap(item) - s[i].count) : 0; return r; };
        auto chestRoomFor = [&](int net, int item) {
            auto key = std::make_pair(net, item);
            auto it = chestRoom.find(key);
            if (it != chestRoom.end()) return it->second;
            long long r = 0;
            for (int mi = 0; mi < nm; ++mi) if (ps.machines[mi].type == PartType::Chest && netOf(mi) == net) r += roomFor(ps.machines[mi].in.data(), kChestSlots, item);
            return chestRoom[key] = r;
        };
        // The nodes and the generators, at their caps.
        std::vector<SumNode> nodes;
        std::vector<SumGen> gens;
        for (int mi = 0; mi < nm; ++mi) {
            const Machine& m = ps.machines[mi];
            const MachineDef& def = Def(m.type);
            if (def.power == Power::Generator) {
                if (m.euCap <= 0) continue;
                SumGen g{ mi, m.euCap * kFix / kCapTicks };
                g.fuelItem = m.fuel.count > 0 ? m.fuel.item : -1;
                for (int f : def.fuels) if (g.fuelItem < 0 && chestStock.count({ netOf(mi), f })) g.fuelItem = f;
                if (g.fuelItem < 0) continue; // nothing to burn (water-only boilers come with the rotation step)
                gens.push_back(g);
                continue;
            }
            if (m.fire == kFireOut) continue;
            for (const auto& [recipe, cap] : m.capPerWindow) {
                if (cap <= 0) continue;
                const Recipe& r = tables.recipes[recipe];
                SumNode n{ mi, recipe, std::min<long long>((long long)cap * kFix / kCapTicks, kFix / std::max(1, TicksOf(m, r))) };
                if (def.power == Power::Fuel) {
                    n.fuelItem = m.fuel.count > 0 ? m.fuel.item : -1;
                    for (int f : def.fuels) if (n.fuelItem < 0 && chestStock.count({ netOf(mi), f })) n.fuelItem = f;
                    if (n.fuelItem < 0) continue;
                }
                if (def.power == Power::Eu && (m.powerNet < 0 || VoltsOf(m) < ps.powerNets[m.powerNet].voltage)) continue;
                nodes.push_back(n);
            }
        }
        // Per tick x kFix: what a node takes and makes of each item; what a generator burns.
        auto nodeNeeds = [&](const SumNode& n, std::vector<std::pair<int, long long>>& out) {
            const Recipe& r = tables.recipes[n.recipe];
            out.clear();
            for (const Stack& s : r.in) out.push_back({ s.item, n.rate * s.count });
            if (n.fuelItem >= 0) out.push_back({ n.fuelItem, n.rate * r.ticks / fuelOf(n.fuelItem) + 1 }); // (+1: rounding against the summary)
        };
        auto genBurn = [&](const SumGen& g) { const int euOut = std::max(1, Def(ps.machines[g.mi].type).euOut); return g.drawn / ((long long)euOut * fuelOf(g.fuelItem)) + (g.drawn > 0 ? 1 : 0); };
        std::vector<std::pair<int, long long>> needs;
        for (int round = 0; round < kSolveRounds; ++round) {
            bool changed = false;
            auto scaleNode = [&](SumNode& n, long long num, long long den) { const long long v = Scaled(n.rate, num, den); if (v < n.rate) { n.rate = v; changed = true; } };
            // Power: a net's consumers share its generators' supply (the TickPlot rule: everyone slows equally).
            for (size_t pn = 0; pn < ps.powerNets.size(); ++pn) {
                long long demand = 0, supply = 0;
                for (const SumNode& n : nodes) if (ps.machines[n.mi].powerNet == (int)pn && Def(ps.machines[n.mi].type).power == Power::Eu) { const Recipe& r = tables.recipes[n.recipe]; demand += n.rate * EuOf(ps.machines[n.mi], r) * TicksOf(ps.machines[n.mi], r); } // (EU a tick while it crafts x the ticks of a craft)
                for (const SumGen& g : gens) if (ps.machines[g.mi].powerNet == (int)pn) supply += Scaled(g.eu, g.scale, kFix);
                if (demand > supply) for (SumNode& n : nodes) if (ps.machines[n.mi].powerNet == (int)pn && Def(ps.machines[n.mi].type).power == Power::Eu) scaleNode(n, supply, demand);
                const long long used = std::min(demand, supply);
                for (SumGen& g : gens) if (ps.machines[g.mi].powerNet == (int)pn) g.drawn = Scaled(Scaled(g.eu, g.scale, kFix), used, std::max(1LL, supply));
            }
            // Items: per net and item, what producers make (P) against what consumers without their own stock want (D').
            std::map<std::pair<int, int>, long long> made, want;
            for (const SumNode& n : nodes) {
                const Recipe& r = tables.recipes[n.recipe];
                for (const Stack& s : r.out) made[{ netOf(n.mi), s.item }] += n.rate * s.count;
                for (const Stack& s : r.chance) made[{ netOf(n.mi), s.item }] += n.rate * s.count / 100; // (#1734) hundredths a craft
                nodeNeeds(n, needs);
                for (const auto& [item, need] : needs) if (!own[n.mi].count(item)) want[{ netOf(n.mi), item }] += need;
            }
            for (const SumGen& g : gens) if (!own[g.mi].count(g.fuelItem)) want[{ netOf(g.mi), g.fuelItem }] += genBurn(g);
            for (const auto& [key, w] : want) {
                const long long p = made.count(key) ? made[key] : 0;
                if (w <= p || chestStock.count(key)) continue; // enough made, or a chest covers the rest (until it empties: an event)
                for (SumNode& n : nodes) {
                    if (netOf(n.mi) != key.first || own[n.mi].count(key.second)) continue;
                    nodeNeeds(n, needs);
                    for (const auto& [item, need] : needs) if (item == key.second && need > 0) { scaleNode(n, p, w); break; }
                }
                for (SumGen& g : gens) if (netOf(g.mi) == key.first && g.fuelItem == key.second && !own[g.mi].count(g.fuelItem)) { const long long v = Scaled(g.scale, p, w); if (v < g.scale) { g.scale = v; changed = true; } }
            }
            // Outputs: what is made and not taken needs a chest's room, or the maker's own room; else the maker slows to the take.
            for (const auto& [key, p] : made) {
                const long long w = want.count(key) ? want[key] : 0;
                if (p <= w || chestRoomFor(key.first, key.second) > 0) continue;
                long long fromRoomless = 0;
                for (const SumNode& n : nodes) {
                    if (netOf(n.mi) != key.first) continue;
                    for (const Stack& s : tables.recipes[n.recipe].out) if (s.item == key.second && roomFor(ps.machines[n.mi].out, kSlots, s.item) <= 0) fromRoomless += n.rate * s.count;
                }
                if (fromRoomless <= 0) continue;
                const long long take = std::max(0LL, w - (p - fromRoomless));
                for (SumNode& n : nodes) {
                    if (netOf(n.mi) != key.first) continue;
                    for (const Stack& s : tables.recipes[n.recipe].out) if (s.item == key.second && roomFor(ps.machines[n.mi].out, kSlots, s.item) <= 0) { scaleNode(n, take, fromRoomless); break; }
                }
            }
            // The conveyor budget: 8 items/s a network.
            const long long budget = (long long)kConveyorItemsPerSecond * kFix / kTicksPerSecond;
            for (int net = 0; net < (int)ps.itemNets.size(); ++net) {
                long long flow = 0;
                for (const auto& [key, p] : made) if (key.first == net) flow += std::max(p, want.count(key) ? want[key] : 0LL);
                for (const auto& [key, w] : want) if (key.first == net && !made.count(key)) flow += w;
                if (flow <= budget) continue;
                for (SumNode& n : nodes) if (netOf(n.mi) == net) scaleNode(n, budget, flow);
            }
            if (!changed) break;
        }
        // The next event: a stock that empties or a room that fills at these rates (ticks, rounded down; at least 1).
        long long dt = t1 - t;
        std::map<std::pair<int, int>, long long> made, want;
        std::vector<std::map<int, long long>> ownUse(nm), ownMake(nm);
        for (const SumNode& n : nodes) {
            const Recipe& r = tables.recipes[n.recipe];
            for (const Stack& s : r.out) { made[{ netOf(n.mi), s.item }] += n.rate * s.count; ownMake[n.mi][s.item] += n.rate * s.count; }
            nodeNeeds(n, needs);
            for (const auto& [item, need] : needs) { if (own[n.mi].count(item)) ownUse[n.mi][item] += need; else want[{ netOf(n.mi), item }] += need; }
        }
        for (const SumGen& g : gens) { const long long b = genBurn(g); if (own[g.mi].count(g.fuelItem)) ownUse[g.mi][g.fuelItem] += b; else want[{ netOf(g.mi), g.fuelItem }] += b; }
        for (int mi = 0; mi < nm; ++mi) for (const auto& [item, use] : ownUse[mi]) if (use > 0) dt = std::min(dt, own[mi][item] * kFix / use);
        for (const auto& [key, w] : want) { const long long p = made.count(key) ? made[key] : 0; if (w > p && chestStock.count(key)) dt = std::min(dt, chestStock[key] * kFix / (w - p)); }
        for (const auto& [key, p] : made) {
            const long long w = want.count(key) ? want[key] : 0;
            if (p <= w) continue;
            const long long room = chestRoomFor(key.first, key.second);
            if (room > 0) { dt = std::min(dt, room * kFix / (p - w)); continue; }
            for (int mi = 0; mi < nm; ++mi) if (netOf(mi) == key.first && ownMake[mi].count(key.second) && ownMake[mi][key.second] > 0)
                dt = std::min(dt, roomFor(ps.machines[mi].out, kSlots, key.second) * kFix / ownMake[mi][key.second]);
        }
        dt = std::max(1LL, dt);
        // Apply the segment: crafts (whole ones, the rest carried), their inputs (own slots first, then what the net made,
        // then chests) and outputs (to what the net used, then chests, then the maker's own slots).
        std::vector<long long> crafts(nodes.size());
        for (size_t k = 0; k < nodes.size(); ++k) {
            long long& carry = ps.machines[nodes[k].mi].carry[nodes[k].recipe];
            const long long total = nodes[k].rate * dt + carry;
            crafts[k] = total / kFix;
            carry = total % kFix;
        }
        // Whole items: a shortfall from rounding (an input that is not there) takes crafts away, never adds any.
        for (int pass = 0; pass < 8; ++pass) {
            std::map<std::pair<int, int>, long long> makeN, needN;
            std::vector<std::map<int, long long>> ownNeed(nm);
            for (size_t k = 0; k < nodes.size(); ++k) {
                const Recipe& r = tables.recipes[nodes[k].recipe];
                const int mi = nodes[k].mi;
                for (const Stack& s : r.out) makeN[{ netOf(mi), s.item }] += crafts[k] * s.count;
                for (const Stack& s : r.chance) makeN[{ netOf(mi), s.item }] += crafts[k] * s.count / 100; // (#1734; the carry may add one)
                for (const Stack& s : r.in) ownNeed[mi][s.item] += crafts[k] * s.count;
                if (nodes[k].fuelItem >= 0) ownNeed[mi][nodes[k].fuelItem] += (crafts[k] * r.ticks + fuelOf(nodes[k].fuelItem) - 1) / fuelOf(nodes[k].fuelItem);
            }
            for (int mi = 0; mi < nm; ++mi) for (const auto& [item, need] : ownNeed[mi]) { const long long fromOwn = std::min(need, own[mi].count(item) ? own[mi][item] : 0LL); if (need > fromOwn) needN[{ netOf(mi), item }] += need - fromOwn; }
            bool short_ = false;
            for (const auto& [key, need] : needN) {
                const long long have = (makeN.count(key) ? makeN[key] : 0) + (chestStock.count(key) ? chestStock[key] : 0);
                if (need <= have) continue;
                short_ = true;
                for (size_t k = 0; k < nodes.size(); ++k) {
                    if (netOf(nodes[k].mi) != key.first || crafts[k] <= 0) continue;
                    bool uses = nodes[k].fuelItem == key.second;
                    for (const Stack& s : tables.recipes[nodes[k].recipe].in) uses |= s.item == key.second;
                    if (uses) crafts[k] = crafts[k] * have / need;
                }
            }
            if (!short_) break;
        }
        // Generators: the fuel for the EU their consumers drew in whole crafts (rounded up).
        std::vector<long long> genFuel(gens.size());
        for (size_t g = 0; g < gens.size(); ++g) {
            long long eu = 0, euAll = 0, supply = 0;
            const int pn = ps.machines[gens[g].mi].powerNet;
            for (size_t k = 0; k < nodes.size(); ++k) if (Def(ps.machines[nodes[k].mi].type).power == Power::Eu && ps.machines[nodes[k].mi].powerNet == pn) { const Recipe& r = tables.recipes[nodes[k].recipe]; euAll += crafts[k] * EuOf(ps.machines[nodes[k].mi], r) * TicksOf(ps.machines[nodes[k].mi], r); }
            for (const SumGen& o : gens) if (ps.machines[o.mi].powerNet == pn) supply += Scaled(o.eu, o.scale, kFix);
            eu = Scaled(euAll, Scaled(gens[g].eu, gens[g].scale, kFix), std::max(1LL, supply));
            const long long per = (long long)std::max(1, Def(ps.machines[gens[g].mi].type).euOut) * fuelOf(gens[g].fuelItem);
            genFuel[g] = (eu + per - 1) / per;
        }
        // Take the inputs.
        std::map<std::pair<int, int>, long long> fromNet, makeN;
        auto takeFrom = [&](int mi, int item, long long n) { // own slots (inputs, then the fuel slot); returns what is left
            Machine& m = ps.machines[mi];
            const long long a = std::min<long long>(n, CountItems(m.in.data(), InSlots(m), item));
            if (a > 0) TakeItems(m.in.data(), InSlots(m), item, (int)a);
            n -= a;
            if (n > 0 && m.fuel.item == item && m.fuel.count > 0) { const long long b = std::min<long long>(n, m.fuel.count); m.fuel.count -= (int)b; if (!m.fuel.count) m.fuel.item = -1; n -= b; }
            return n;
        };
        for (size_t k = 0; k < nodes.size(); ++k) {
            if (crafts[k] <= 0) continue;
            const Recipe& r = tables.recipes[nodes[k].recipe];
            const int mi = nodes[k].mi;
            for (const Stack& s : r.in) { const long long left = takeFrom(mi, s.item, crafts[k] * s.count); if (left) fromNet[{ netOf(mi), s.item }] += left; ps.stats.items[s.item].used.add(crafts[k] * s.count); }
            if (nodes[k].fuelItem >= 0) { const long long nf = (crafts[k] * r.ticks + fuelOf(nodes[k].fuelItem) - 1) / fuelOf(nodes[k].fuelItem); const long long left = takeFrom(mi, nodes[k].fuelItem, nf); if (left) fromNet[{ netOf(mi), nodes[k].fuelItem }] += left; ps.stats.items[nodes[k].fuelItem].used.add(nf); }
            for (const Stack& s : r.out) { makeN[{ netOf(mi), s.item }] += crafts[k] * s.count; ps.stats.items[s.item].made.add(crafts[k] * s.count); }
            for (const Stack& s : r.chance) { // (#1734) outputs by chance: the machine's carry, as when it runs near
                int& c = ps.machines[mi].disCarry[s.item];
                const long long all = c + crafts[k] * s.count;
                c = (int)(all % 100);
                if (all >= 100) { makeN[{ netOf(mi), s.item }] += all / 100; ps.stats.items[s.item].made.add(all / 100); }
            }
        }
        for (size_t g = 0; g < gens.size(); ++g) if (genFuel[g] > 0) { const long long left = takeFrom(gens[g].mi, gens[g].fuelItem, genFuel[g]); if (left) fromNet[{ netOf(gens[g].mi), gens[g].fuelItem }] += left; ps.stats.items[gens[g].fuelItem].used.add(genFuel[g]); }
        // What the net made covers what it needed; the rest of the need comes out of its chests; the rest of the making goes
        // into its chests, then its makers' own output slots.
        for (const auto& [key, need] : fromNet) {
            const long long fromMade = std::min(need, makeN.count(key) ? makeN[key] : 0LL);
            if (fromMade) makeN[key] -= fromMade;
            long long left = need - fromMade;
            for (int mi = 0; mi < nm && left > 0; ++mi) {
                Machine& c = ps.machines[mi];
                if (c.type != PartType::Chest || netOf(mi) != key.first) continue;
                const long long a = std::min<long long>(left, CountItems(c.in.data(), kChestSlots, key.second));
                if (a > 0) { TakeItems(c.in.data(), kChestSlots, key.second, (int)a); left -= a; c.dirty = true; }
            }
        }
        for (auto& [key, extra] : makeN) {
            long long left = extra;
            for (int mi = 0; mi < nm && left > 0; ++mi) {
                Machine& c = ps.machines[mi];
                if (c.type != PartType::Chest || netOf(mi) != key.first) continue;
                left = AddItems(c.in.data(), kChestSlots, key.second, (int)left);
                c.dirty = true;
            }
            for (size_t k = 0; k < nodes.size() && left > 0; ++k) {
                Machine& m = ps.machines[nodes[k].mi];
                if (netOf(nodes[k].mi) != key.first) continue;
                for (const Stack& s : tables.recipes[nodes[k].recipe].out) if (s.item == key.second && left > 0) left = AddItems(m.out, kSlots, key.second, (int)left);
            }
        }
        for (Machine& m : ps.machines) { m.dirty = true; m.asleep = false; } // (back near, everything looks again)
        for (size_t k = 0; k < nodes.size(); ++k) ps.machines[nodes[k].mi].status = crafts[k] > 0 ? MachineStatus::Working : MachineStatus::Idle;
        for (const Machine& m : ps.machines) ps.stats.statusTicks[(int)m.type][(int)m.status & 7] += dt; // (#1199 counters: the bottleneck finder sees far plots too)
        if (summaryLog) {
            long long c = 0, rates = 0;
            for (size_t k = 0; k < nodes.size(); ++k) { c += crafts[k]; rates += nodes[k].rate; }
            logCrafts += c;
            LOG_INFO("Factory", "summary: plot %u segment %d: %lld ticks, %zu nodes at %lld crafts/Mtick in all, %lld crafts\n", ps.plot, segment, dt, nodes.size(), rates, c);
        }
        t += dt;
    }
    if (summaryLog) LOG_INFO("Factory", "summary: plot %u advanced %lld -> %lld (%d segments, %lld crafts)\n", ps.plot, t0, t, segment, logCrafts);
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
        cableView.generators += pw == Power::Generator || pw == Power::Dynamo; cableView.batteries += pw == Power::Battery; cableView.consumers += pw == Power::Eu || pw == Power::Motor;
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
    if (m.recipe >= 0) { v.progress = m.progress; v.duration = TicksOf(m, tables.recipes[m.recipe]); }
    else if (m.mode == 1 && m.disSlot >= 0 && m.disTicks > 0) { v.progress = m.progress; v.duration = m.disTicks; } // (#1427 L7c)
    else if (!m.jobs.empty()) { // (#1429 PR 4) a structure: its first job's progress
        const Machine::MbJob& j = m.jobs[0];
        v.progress = j.progress; v.duration = TicksOf(m, tables.recipes[j.recipe]);
    }
    v.energy = m.energy;
    v.energyMax = m.type == PartType::Battery ? kBatteryCapacity : 0;
    v.hasFuelSlot = HasFuelSlot(m.type);
    v.hotOutputs = HeatsItems(m.type);
    v.canLight = m.type == PartType::Campfire;
    v.fire = m.fire;
    v.mode = m.mode; // (#1271 L6c)
    v.pick = m.pick >= 0 ? tables.name(m.pick) : std::string(); // (#1271 L6c)
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
    for (int i = 0; i < InSlots(m); ++i) { v.inSlots.push_back({ !Empty(m.in[i]) ? tables.name(m.in[i].item) : std::string(), m.in[i].count }); v.inCenti.push_back(m.in[i].centi); v.inData.push_back(m.in[i].data); } // (#1102; #1446 data)
    if (m.type != PartType::Chest) for (int i = 0; i < kSlots; ++i) { v.outSlots.push_back({ m.out[i].count > 0 ? tables.name(m.out[i].item) : std::string(), m.out[i].count }); v.outData.push_back(m.out[i].count > 0 ? m.out[i].data : std::string()); } // (#1446 2c)
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
                if (lockedQueued) { lockedSim = std::unordered_set<int>(lockedQueue.begin(), lockedQueue.end()); lockedQueued = false; } // (#1271 L6b)
            }
            if (reset) { plots.clear(); plotIndexDirty = true; emaMs = maxMs = 0.0; statSeconds = 0; tickRingAt = 0; }
            // Moved plots first, so the new plot's first layout finds its machines (same part indices).
            for (const MoveMsg& mv : moves)
                if (PlotSim* p = FindPlot(mv.fromBody, mv.fromPlot); p && !FindPlot(mv.toBody, mv.toPlot)) { p->body = mv.toBody; p->plot = mv.toPlot; plotIndexDirty = true; viewsRebuild = true; }
            for (const LayoutMsg& L : layouts) ApplyLayout(L);
            for (const auto& n : nears) if (PlotSim* p = FindPlot(n.first, n.second.first); p && p->nearPlayer != n.second.second) SetNear(*p, n.second.second, tickCount.load()); // (#1261) leaving near starts measuring
            std::vector<ItemOut> toCursor; // kind 8 (#192)
            std::vector<ItemOut> back;
            back.swap(refunds); // removed machines' contents (from this tick's layouts)
            std::vector<std::pair<std::tuple<uint64_t, uint32_t, int>, int>> drained; // pump drains applied this tick (pumpPending)
            std::vector<std::pair<std::tuple<uint64_t, uint32_t, int>, int>> filled;  // tank fills applied this tick (tankPending)
            for (const ItemMsg& im : items) {
                PlotSim* p = FindPlot(im.body, im.plot);
                if (im.kind == 11) { if (p) p->rot.Hold(im.part); continue; } // (#1269) a player turns a hand crank
                if (im.kind == 14) { if (p) if (auto sit = p->byPart.find(im.part); sit != p->byPart.end()) { Machine& sm = p->machines[sit->second]; if (im.count < 0) { sm.pick = im.item; if (sm.recipe >= 0 && sm.pick >= 0 && (tables.recipes[sm.recipe].out.empty() || tables.recipes[sm.recipe].out[0].item != sm.pick)) { sm.recipe = -1; sm.progress = 0; } } else { sm.mode = (uint8_t)std::clamp(im.count, 0, 1); sm.recipe = -1; sm.progress = 0; sm.disSlot = -1; sm.disTicks = 0; } sm.dirty = true; } continue; } // (#1271 L6c) an assembler's mode (count 0/1), or a machine's recipe pick (count -1: item, -1 any)
                if (im.kind == 13) { // (#1271 L2) a strike on an anvil: its work, up to two strikes ahead
                    if (p) if (auto sit = p->byPart.find(im.part); sit != p->byPart.end()) { Machine& sm = p->machines[sit->second]; sm.strike = std::min(sm.strike + kStrikeTicks, 2 * kStrikeTicks); sm.dirty = true; }
                    continue;
                }
                if (im.kind == 12) { if (p) p->rot.SetSourceRpm(im.part, im.count); continue; } // (#1270) a water wheel's or windmill's sampled speed
                auto it = p ? p->byPart.find(im.part) : decltype(p->byPart.end()){};
                if (!p || it == p->byPart.end()) {
                    if (im.kind == 0 || im.kind == 9) back.push_back({ im.item, im.count, im.to, im.centi }); // (#1849) and its fraction
                    if (im.kind == 5) drained.push_back({ std::make_tuple(im.body, im.plot, im.part), im.count });
                    if (im.kind == 7) filled.push_back({ std::make_tuple(im.body, im.plot, im.part), im.count });
                    continue;
                }
                Machine& m = p->machines[it->second];
                m.dirty = true; // (#1214) any command on a machine wakes it and its network
                if (!p->nearPlayer && im.kind != 5 && im.kind != 7) BeginMeasuring(*p, tickCount.load()); // (#1261) a player's insert or take while far (not pipe flows: pumps, tanks)
                if (!im.data.empty() && (im.kind == 0 || im.kind == 2 || im.kind == 6 || im.kind == 9)) { // (#1446 2b) one item with its record
                    bool placed = false;
                    Stack* s = nullptr; int n = 0;
                    if (m.type == PartType::Chest) { s = m.in.data(); n = kChestSlots; }
                    // (#1271 L7c hook) an assembler set to disassemble takes a composed item into an input slot, record intact;
                    // what it gives back is L7c's (Progression). Every other machine and automation path refuses items with data.
                    else if (m.mode == 1 && im.kind != 6) { s = m.in.data(); n = InSlots(m); }
                    // (#1219, user 2026-11-22) the anvil takes a composed copy to repair it or take a gem out (AnvilJob)
                    else if (m.type == PartType::Anvil && im.kind != 6 && Composed_Is(tables.name(im.item))) { s = m.in.data(); n = InSlots(m); }
                    if (s) {
                        const int want = im.kind == 9 && m.type == PartType::Chest ? im.amount : -1;
                        for (int i = 0; i < n && !placed; ++i) {
                            if (want >= 0 && i != want) continue;
                            if (!Empty(s[i])) continue;
                            s[i] = Stack{}; s[i].item = im.item; s[i].count = 1; s[i].data = im.data; placed = true;
                        }
                    }
                    if (placed) p->stats.items[im.item].handIn.add(1);
                    else if (im.kind != 2) { ItemOut o{ im.item, 1, im.to }; o.data = im.data; back.push_back(o); }
                } else if (im.kind == 9 && m.type == PartType::Chest) { // (#862) into one storage slot: empty, or the same item up to its cap
                    int left = im.count, leftCenti = im.centi;
                    if (im.amount >= 0 && im.amount < kChestSlots && im.item >= 0) {
                        Stack& s = m.in[im.amount];
                        if (Empty(s) || (s.item == im.item && Plain(s))) {
                            if (Empty(s)) { s = Stack{}; s.item = im.item; }
                            const int k = std::min(left, tables.cap(im.item) - std::max(0, s.count));
                            if (k > 0) { s.count += k; left -= k; }
                            if (leftCenti > 0 && s.count < tables.cap(im.item)) { // (#1102) the fraction too, while there is room
                                s.centi += leftCenti; leftCenti = 0;
                                if (s.centi >= 100) { s.centi -= 100; ++s.count; }
                            }
                            if (Empty(s)) s = Stack{};
                        }
                    }
                    if (im.count - left > 0) p->stats.items[im.item].handIn.add(im.count - left); // (#1199)
                    if (left > 0 || leftCenti > 0) back.push_back({ im.item, left, im.to, leftCenti });
                } else if (im.kind == 0 || im.kind == 2 || im.kind == 6 || im.kind == 9) { // (9 on a machine: routed, as Insert)
                    int buckets = 0;
                    const int left = InsertInto(m, im.item, im.count, im.kind == 6, buckets);
                    if (im.count - left > 0) p->stats.items[im.item].handIn.add(im.count - left); // (#1199)
                    if (left > 0 && im.kind != 2) back.push_back({ im.item, left, im.to });
                    if (buckets > 0 && im.kind != 2) back.push_back({ BucketItem(), buckets, im.to }); // the emptied buckets
                    if (im.centi > 0) { // (#1102) a fraction: a chest keeps it; a machine crafts whole items, so it goes back
                        const int c = m.type == PartType::Chest ? AddCenti(m.in.data(), kChestSlots, im.item, im.centi) : im.centi;
                        if (c > 0 && im.kind != 2) back.push_back({ im.item, 0, im.to, c });
                    }
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
                    if (im.item >= 0 && im.item < n && !Empty(s[im.item]) && !tables.fluid(s[im.item].item)) { p->stats.items[s[im.item].item].handOut.add(s[im.item].count); ItemOut o{ s[im.item].item | hot, s[im.item].count, im.to, s[im.item].centi }; o.data = s[im.item].data; back.push_back(o); s[im.item] = Stack{}; } // (#1254) never fluid; (#1446) the record goes with it
                } else if (im.kind == 8) { // part of a slot onto the cursor (#192)
                    Stack* s = im.count == 2 ? &m.fuel : im.count ? m.out : m.in.data();
                    const int n = im.count == 2 ? 1 : im.count ? kSlots : InSlots(m);
                    const int hot = s == m.out && HeatsItems(m.type) ? kHotItem : 0;
                    if (im.item >= 0 && im.item < n && !Empty(s[im.item]) && !tables.fluid(s[im.item].item)) { // (#1254) never fluid
                        const int k = im.amount > 0 ? std::min(im.amount, s[im.item].count) : s[im.item].count;
                        if (k > 0) { ItemOut o{ s[im.item].item | hot, k, im.to }; o.data = s[im.item].data; toCursor.push_back(o); } // (#1446) with its record
                        p->stats.items[s[im.item].item].handOut.add(k); // (#1199)
                        s[im.item].count -= k;
                        // (#1102) the whole slot taken: its fraction goes to the inventory (the cursor holds whole items from here)
                        if (s[im.item].count <= 0 && s[im.item].centi > 0) { back.push_back({ s[im.item].item, 0, im.to, s[im.item].centi }); s[im.item].centi = 0; }
                        if (Empty(s[im.item])) s[im.item] = Stack{};
                    }
                } else {
                    Stack* s = m.type == PartType::Chest ? m.in.data() : m.out;
                    const int n = m.type == PartType::Chest ? kChestSlots : kSlots;
                    const int hot = s == m.out && HeatsItems(m.type) ? kHotItem : 0; // #580
                    for (int i = 0; i < n; ++i) if (!Empty(s[i]) && !tables.fluid(s[i].item)) { p->stats.items[s[i].item].handOut.add(s[i].count); ItemOut o{ s[i].item | hot, s[i].count, im.to, s[i].centi }; o.data = s[i].data; back.push_back(o); s[i] = Stack{}; } // (#1446) with its record
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
                if (ps.mode == PlotSim::Mode::Near) RotateNearBuckets(ps, t); // (#1261) the near minute
                UpdateMode(ps, t); // (#1261) measuring done: Summary
                if (summaryOn && ps.mode == PlotSim::Mode::Summary && !ps.nearPlayer && (t + (long long)i) % kAdvanceEvery == 0) { AdvanceSummary(ps, ps.lastAdvance, t); ps.lastAdvance = t; } // (#1261) every 10 s, staggered
                if (const int step = StepOf(ps, i, t)) jobs.push_back({ &ps, step });
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
                for (PlotSim& ps : plots) { // (#1728) finished recipes, in plot order (deterministic); a full queue drops the newest
                    for (const auto& [item, count] : ps.completed) if (madeResults.size() < 65536) madeResults.push_back({ ps.body, ps.plot, item, count });
                    ps.completed.clear();
                }
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
                    workingView.clear(); workingByPlot.clear(); chestView.clear(); stalledByPlot.clear(); anvilView.clear();
                    for (PlotSim& ps : plots) ps.pub = PlotSim::Published{};
                    viewsRebuild = false;
                }
                static std::vector<int> working, chests, stalled;
                static std::vector<std::pair<int, int>> anvils;
                for (PlotSim& ps : plots) {
                    PlotSim::Published& pub = ps.pub;
                    auto unpublishWorking = [&] { for (int part : pub.working) workingView.erase(WorkingKey(ps.body, ps.plot, part)); workingByPlot.erase({ ps.body, ps.plot }); pub.working.clear();
                                                  stalledByPlot.erase({ ps.body, ps.plot }); pub.stalled.clear(); };
                    auto unpublishChests = [&] { for (int part : pub.chests) chestView.erase(std::make_tuple(ps.body, ps.plot, part)); pub.chests.clear(); };
                    auto unpublishAnvils = [&] { for (const auto& a : pub.anvils) anvilView.erase(std::make_tuple(ps.body, ps.plot, a.first)); pub.anvils.clear(); };
                    if (!ps.nearPlayer) {
                        if (pub.near) { unpublishWorking(); unpublishChests(); unpublishAnvils(); pub.near = false; }
                        continue;
                    }
                    working.clear(); chests.clear(); stalled.clear(); anvils.clear();
                    uint64_t h = 1469598103934665603ull; // FNV-1a over the chests' parts and slots
                    auto mix = [&h](uint64_t v) { h = (h ^ v) * 1099511628211ull; };
                    for (const Machine& m : ps.machines) {
                        if (m.status == MachineStatus::Working) working.push_back(m.part);
                        else if (m.status != MachineStatus::Idle && m.status != MachineStatus::NoLoad) stalled.push_back(m.part); // (#1243)
                        if (m.type == PartType::Anvil) { // (#1613) its workpiece, drawn on it: the finished piece first, else the last put on
                            int item = -1;
                            for (int i = 0; i < kSlots && item < 0; ++i) if (!Empty(m.out[i])) item = m.out[i].item;
                            for (int i = InSlots(m) - 1; i >= 0 && item < 0; --i) if (!Empty(m.in[i])) item = m.in[i].item;
                            if (item >= 0) anvils.push_back({ m.part, item });
                        }
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
                    if (!pub.near || stalled != pub.stalled) { // (#1243) machine sounds' stall cues
                        if (stalled.empty()) stalledByPlot.erase({ ps.body, ps.plot }); else stalledByPlot[{ ps.body, ps.plot }] = stalled;
                        pub.stalled = stalled;
                    }
                    if (!pub.near || h != pub.chestHash || chests != pub.chests) {
                        unpublishChests();
                        for (const Machine& m : ps.machines) {
                            if (m.type != PartType::Chest) continue;
                            auto& v = chestView[std::make_tuple(ps.body, ps.plot, m.part)];
                            for (int i = 0; i < kChestSlots; ++i) if (m.in[i].count > 0 && Plain(m.in[i])) v.push_back({ m.in[i].item, m.in[i].count }); // (#1446) crafting from chests never takes an item with data
                        }
                        pub.chests = chests;
                        pub.chestHash = h;
                    }
                    if (!pub.near || anvils != pub.anvils) {
                        unpublishAnvils();
                        for (const auto& a : anvils) anvilView[std::make_tuple(ps.body, ps.plot, a.first)] = a.second;
                        pub.anvils = anvils;
                    }
                    pub.near = true;
                }
            }
            // (#1207) The whole tick (commands, plots, views) into the last minute's ring; its percentiles once a second.
            const auto t3 = clk::now();
            const float cmdMs = std::chrono::duration<float, std::milli>(t1 - t0).count();
            const float tickMs = std::chrono::duration<float, std::milli>(t2 - t1).count();
            const float viewMs = std::chrono::duration<float, std::milli>(t3 - t2).count();
            if (tickRingClear.exchange(false)) tickRingAt = 0;
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
        { // (#1524) under the waiters' mutex: a WaitTicks between its check and its sleep would miss the notify (30 s stall)
            std::lock_guard<std::mutex> lk(tickCvMutex);
            tickCount++;
        }
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
    sentFormed.clear();
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
            // (#1429 PR 4) forming or unforming a structure changes no part (no version): its list resends the layout too
            uint64_t formedHash = 1469598103934665603ULL;
            for (int f : pl->formed) formedHash = (formedHash ^ (uint64_t)(uint32_t)f) * 1099511628211ULL;
            if (sv == sentVersions.end() || sv->second != pl->version || sentFormed[key] != formedHash) {
                sentFormed[key] = formedHash;
                LayoutMsg L{ body->id, pl->id, pl->version, {}, {}, {}, {}, {} };
                L.refundTo = pl->refundTo; pl->refundTo = 0; // (#854)
                std::vector<std::pair<Int3, CellShape>> fp;
                for (int c : pl->formed) { // (#1429 PR 4) each formed structure that still validates: what it runs and its hatches
                    MbMachine mb;
                    if (!Multiblock_Members(*pl, c, mb) || mb.def->machine.empty()) continue;
                    LayoutMsg::Structure s;
                    s.controller = c; s.machine = MachineByName(mb.def->machine); s.parallel = mb.def->parallel;
                    if (s.machine == PartType::Count) continue;
                    // Its tier: the lowest voltage of its energy hatches and casings that have one (#1429 design, decided);
                    // none (a steam-age structure): its machine's own tier, as a single block of it.
                    s.tier = MachineTiers_PartTier(s.machine, 0, Def(s.machine).volts);
                    bool capped = false;
                    for (const std::string& v : mb.volts) { const int t = MachineTiers_ByName(v); if (t >= 0) { s.tier = capped ? std::min(s.tier, t) : t; capped = true; } }
                    for (const auto& h : mb.hatches) s.hatches.push_back(h.first);
                    L.structures.push_back(std::move(s));
                }
                for (size_t i = 0; i < pl->parts.size(); ++i) {
                    const Part& part = pl->parts[i];
                    if (part.type == PartType::Count) continue;
                    if (part.type == PartType::MbController || part.type == PartType::MbHatch) { // (#1429 PR 4) machines of the factory
                        Plot::Footprint(part, fp);
                        std::vector<Int3> cells;
                        for (auto& c : fp) cells.push_back(c.first);
                        const std::string hid = part.type == PartType::MbHatch ? Multiblock_PartItem(part.type, part.material) : "";
                        const auto hv = Multiblock_All().volts.find(hid);
                        L.machines.push_back({ (int)i, part.type });
                        L.machineCells.push_back(std::move(cells));
                        L.machineVariant.push_back("");
                        L.machineTier.push_back(hv != Multiblock_All().volts.end() ? MachineTiers_ByName(hv->second) : 0);
                        L.machineHatch.resize(L.machines.size() - 1, kHatchNone);
                        L.machineHatch.push_back(part.type == PartType::MbHatch ? HatchOf(Multiblock_HatchKind(part.material)) : kHatchNone);
                        continue;
                    }
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
                        if ((part.type >= PartType::LvLab && part.type <= PartType::HvLab) || part.type == PartType::SteamLab) variant = "research_none"; // (#1193) idle until the game sets its team's node
                        if (part.type == PartType::ManaSiphon) variant = "jar_none"; // (#1252) idle until the game finds the ley node it stands on
                        L.machineVariant.push_back(variant);
                        L.machineTier.push_back(MachineTiers_PartTier(part.type, part.material, Def(part.type).volts)); // (#1427 L7a) its variant
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
                bool near = cameraNear && Vector3Distance(o, world.camera.position) < 1000.0f; // (#1302: not on a dedicated server)
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
    std::vector<MadeOut> finished; // (#1728) recipes the machines completed since the last Sync: the owner's quests count them (nobody need take them)
    {
        std::lock_guard<std::mutex> v(viewMutex);
        finished.swap(madeResults);
    }
    if (gFactoryMade) for (const MadeOut& f : finished) gFactoryMade(f.body, f.plot, tables.name(f.item & ~kHotItem), f.count);
    // kHotItem marks what came out of a heating machine's output (#580): it arrives hot, in the inventory or on the cursor.
    for (auto& r : back) { const std::string id = tables.name(r.item & ~kHotItem); const float h = (r.item & kHotItem) ? 1.0f : 0.0f; const double amount = r.count + r.centi / 100.0; /* (#1102) */ if (r.to && gFactoryDeliver) gFactoryDeliver(r.to, id, amount, h, r.data); else if (!r.data.empty()) inventory.addInstance(id, r.data, h); else inventory.addAmount(id, amount, h); } // (#332) a server: to that player
    for (auto& r : cur) { const std::string id = tables.name(r.item & ~kHotItem); const float h = (r.item & kHotItem) ? 1.0f : 0.0f; if (r.to && gFactoryDeliver) gFactoryDeliver(r.to, id, r.count, h, r.data); else if (cursorSink) cursorSink(id, r.count, h, r.data); else if (!r.data.empty()) inventory.addInstance(id, r.data, h); else inventory.addAmount(id, r.count, h); }
}

void Factory::TakeSlotToCursor(uint64_t body, uint32_t plot, int part, int slot, int which, int count, uint32_t to) {
    std::lock_guard<std::mutex> q(queueMutex);
    ItemMsg m{ 8, body, plot, part, which == 2 ? 0 : slot, which };
    m.amount = count; m.to = to; // (#854) a server: that player's inventory (its client puts it on its cursor)
    itemQueue.push_back(m);
}
void Factory::InsertSlot(uint64_t body, uint32_t plot, int part, const std::string& item, int count, int slot, uint32_t to, int centi, const std::string& data) {
    std::lock_guard<std::mutex> q(queueMutex);
    ItemMsg m{ 9, body, plot, part, tables.item(item), count };
    m.amount = slot; m.to = to; m.centi = data.empty() ? std::clamp(centi, 0, 99) : 0; m.data = data; if (!data.empty()) m.count = 1; // (#1102; #1446 data: one item)
    itemQueue.push_back(m);
}
void Factory::SetCursorSink(void (*sink)(const std::string&, int, float, const std::string&)) { cursorSink = sink; }

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

void Factory::ResetTickStats() { tickRingClear = true; }

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
    std::set<int> estimated; // (#1261) items counted on plots in the summary sim
    for (const PlotSim& ps : plots) {
        if (onePlot && (ps.body != body || ps.plot != plot)) continue;
        ++out.plots;
        const bool summary = ps.mode == PlotSim::Mode::Summary;
        out.estimatedPlots += summary;
        for (const auto& [id, is] : ps.stats.items) {
            if (summary) estimated.insert(id);
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
        l.estimated = estimated.count(id) > 0;
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

void Factory::Insert(uint64_t body, uint32_t plot, int part, const std::string& item, int count, uint32_t to, int centi, const std::string& data) {
    std::lock_guard<std::mutex> q(queueMutex);
    ItemMsg m{ 0, body, plot, part, tables.item(item), data.empty() ? count : 1 }; m.to = to; m.centi = data.empty() ? std::clamp(centi, 0, 99) : 0; m.data = data; itemQueue.push_back(m); // (#1102; #1446 data: one item)
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
// (#1193) Between ticks (the sim holds simMutex for a whole tick): a machine's one output, and its outputs taken by name.
bool Factory::SetVariant(uint64_t body, uint32_t plot, int part, const std::string& item) {
    std::lock_guard<std::mutex> lock(simMutex);
    PlotSim* p = FindPlot(body, plot);
    auto it = p ? p->byPart.find(part) : decltype(p->byPart.end()){};
    if (!p || it == p->byPart.end()) return false;
    Machine& m = p->machines[it->second];
    const int v = item.empty() ? -1 : tables.item(item);
    if (m.variant == v) return true;
    m.variant = v;
    m.dirty = true; // (#1214) a sleeping lab wakes: its new node's recipe may start now
    // (#1302) a far lab's caps (or a measuring window) were for its old node's recipe: its plot measures again, as after
    // a layout change (the summary can't run a recipe it never measured)
    if (!p->nearPlayer) BeginMeasuring(*p, tickCount.load());
    // A recipe of another output stops (its inputs are only taken when it finishes, #193: nothing is lost).
    if (m.recipe >= 0 && (v < 0 || tables.recipes[m.recipe].out.empty() || tables.recipes[m.recipe].out[0].item != v)) { m.recipe = -1; m.progress = 0; }
    return true;
}
std::vector<std::pair<std::string, int>> Factory::DrainOutputs(uint64_t body, uint32_t plot, int part, const std::string& prefix) {
    std::vector<std::pair<std::string, int>> out;
    std::lock_guard<std::mutex> lock(simMutex);
    PlotSim* p = FindPlot(body, plot);
    auto it = p ? p->byPart.find(part) : decltype(p->byPart.end()){};
    if (!p || it == p->byPart.end()) return out;
    for (Stack& s : p->machines[it->second].out) {
        if (s.count <= 0) continue;
        const std::string name = tables.name(s.item);
        if (name.compare(0, prefix.size(), prefix) != 0) continue;
        out.push_back({ name, s.count });
        s = Stack{};
        p->machines[it->second].dirty = true; // (#1214) room in its outputs: a lab that waited on them goes on
    }
    return out;
}
// (#1254) Fluid items in machine slots (alchemy's essences).
void Factory::AddFluidPrefix(const std::string& prefix) { tables.addFluidPrefix(prefix); }
bool Factory::IsFluidItem(const std::string& item) { return tables.fluid(tables.item(item)); }
std::vector<Factory::FluidSlots> Factory::ReadFluids(const std::vector<MachineKey>& keys) {
    std::vector<FluidSlots> out(keys.size());
    std::lock_guard<std::mutex> lock(simMutex);
    for (size_t k = 0; k < keys.size(); ++k) {
        const auto& [body, plot, part] = keys[k];
        PlotSim* p = FindPlot(body, plot);
        auto it = p ? p->byPart.find(part) : decltype(p->byPart.end()){};
        if (!p || it == p->byPart.end()) continue;
        const Machine& m = p->machines[it->second];
        FluidSlots& f = out[k];
        f.valid = true;
        auto read = [&](const Stack* s, int n, std::vector<std::pair<std::string, int>>& list, int& freeSlots) {
            for (int i = 0; i < n; ++i) {
                if (Empty(s[i])) { ++freeSlots; continue; }
                if (!tables.fluid(s[i].item)) continue;
                const std::string name = tables.name(s[i].item);
                bool merged = false;
                for (auto& e : list) if (e.first == name) { e.second += s[i].count; merged = true; }
                if (!merged) list.push_back({ name, s[i].count });
            }
        };
        read(m.in.data(), InSlots(m), f.in, f.freeIn);
        read(m.out, kSlots, f.out, f.freeOut);
    }
    return out;
}
std::vector<int> Factory::MoveFluids(const std::vector<FluidMove>& moves) {
    std::vector<int> done(moves.size(), 0);
    std::lock_guard<std::mutex> lock(simMutex);
    auto machineOf = [&](const MachineKey& key) -> Machine* {
        const auto& [body, plot, part] = key;
        PlotSim* p = FindPlot(body, plot);
        auto it = p ? p->byPart.find(part) : decltype(p->byPart.end()){};
        return !p || it == p->byPart.end() ? nullptr : &p->machines[it->second];
    };
    for (size_t k = 0; k < moves.size(); ++k) {
        const FluidMove& mv = moves[k];
        if (mv.count <= 0 || mv.item.empty()) continue;
        const int item = tables.item(mv.item);
        if (!tables.fluid(item)) continue; // only fluid moves this way
        Machine* from = mv.hasFrom ? machineOf(mv.from) : nullptr;
        Machine* to = mv.hasTo ? machineOf(mv.to) : nullptr;
        if ((mv.hasFrom && !from) || (mv.hasTo && !to)) continue;
        int n = mv.count;
        Stack* fs = nullptr; int fn = 0;
        if (from) { fs = mv.fromWhich ? from->out : from->in.data(); fn = mv.fromWhich ? kSlots : InSlots(*from); n = std::min(n, CountItems(fs, fn, item)); }
        Stack* ts = nullptr; int tn = 0;
        if (to) {
            ts = mv.toWhich ? to->out : to->in.data(); tn = mv.toWhich ? kSlots : InSlots(*to);
            if (mv.cap > 0) n = std::min(n, std::max(0, mv.cap - CountItems(ts, tn, item)));
            Stack trial[kChestSlots];
            std::copy(ts, ts + tn, trial);
            n -= AddItems(trial, tn, item, n); // what fits
        }
        if (n <= 0) continue;
        if (from) { TakeItems(fs, fn, item, n); from->dirty = true; } // (#1214) room in its outputs: it goes on
        if (to) { AddItems(ts, tn, item, n); to->dirty = true; }       // new input: it looks again next tick
        done[k] = n;
    }
    return done;
}
void Factory::SetFire(uint64_t body, uint32_t plot, int part, bool lit) {
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back(ItemMsg{ 10, body, plot, part, -1, lit ? 1 : 0 });
}
void Factory::SetLockedOutputs(const std::vector<std::string>& outputs) { // (#1271 L6b)
    LoadRecipes();
    std::vector<int> ids;
    for (const std::string& o : outputs) ids.push_back(tables.item(o));
    std::lock_guard<std::mutex> q(queueMutex);
    lockedQueue = std::move(ids); lockedQueued = true;
}
void Factory::SetMode(uint64_t body, uint32_t plot, int part, int mode) { // (#1271 L6c)
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back(ItemMsg{ 14, body, plot, part, -1, mode });
}
void Factory::SetPick(uint64_t body, uint32_t plot, int part, const std::string& item) { // (#1271 L6c)
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back(ItemMsg{ 14, body, plot, part, item.empty() ? -1 : tables.item(item), -1 });
}
void Factory::Strike(uint64_t body, uint32_t plot, int part) { // (#1271 L2)
    std::lock_guard<std::mutex> q(queueMutex);
    itemQueue.push_back(ItemMsg{ 13, body, plot, part, -1, 0 });
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
        ri.machine = r.machine; ri.ticks = r.ticks; ri.eu = r.eu; ri.tier = r.tier; // (#1427 L7d: tier)
        for (const Stack& s : r.in) ri.in.push_back({ tables.name(s.item), s.count });
        for (const Stack& s : r.out) ri.out.push_back({ tables.name(s.item), s.count });
        for (const Stack& s : r.chance) ri.chance.push_back({ tables.name(s.item), s.count });
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
    // (#1302) Far plots in the summary sim first come up to now, so the save holds their chests as they are.
    const long long now = tickCount.load();
    for (PlotSim& ps : plots) if (ps.body == bodyId && ps.mode == PlotSim::Mode::Summary && !ps.nearPlayer) { AdvanceSummary(ps, ps.lastAdvance, now); ps.lastAdvance = now; }
    int n = 0;
    for (const PlotSim& ps : plots) if (ps.body == bodyId) n += (int)ps.machines.size();
    // (#1881) Loaded machines wait in pendingStates until their plot reaches the sim (the next Sync): a save before that
    // wrote "machines 0" and lost every chest and machine on the body. Those go back out as they were read.
    for (const auto& [key, rest] : pendingStates) if (std::get<0>(key) == bodyId) ++n;
    fprintf(f, "machines %d\n", n);
    for (const auto& [key, rest] : pendingStates) {
        if (std::get<0>(key) != bodyId) continue;
        fprintf(f, "m %u %d", std::get<1>(key), std::get<2>(key));
        for (const std::string& w : rest) fprintf(f, " %s", w.c_str());
        fprintf(f, "\n");
    }
    for (const PlotSim& ps : plots) {
        if (ps.body != bodyId) continue;
        const bool summary = ps.mode == PlotSim::Mode::Summary;
        for (const Machine& m : ps.machines) {
            fprintf(f, "m %u %d %d %d %d %d", ps.plot, m.part, m.recipe, m.progress, m.energy, m.fuelEU);
            int nin = 0; for (int i = 0; i < InSlots(m); ++i) if (!Empty(m.in[i])) ++nin;
            fprintf(f, " %d", nin);
            for (int i = 0; i < InSlots(m); ++i) if (!Empty(m.in[i])) { if (!m.in[i].data.empty()) fprintf(f, " %s|%s 1", tables.name(m.in[i].item).c_str(), m.in[i].data.c_str()); /* (#1446) its record */ else if (m.in[i].centi > 0) fprintf(f, " %s %d.%02d", tables.name(m.in[i].item).c_str(), m.in[i].count, m.in[i].centi); else fprintf(f, " %s %d", tables.name(m.in[i].item).c_str(), m.in[i].count); } // (#1102) a fraction: "2.50"
            int nout = 0; for (int i = 0; i < kSlots; ++i) if (m.out[i].count > 0) ++nout;
            fprintf(f, " %d", nout);
            for (int i = 0; i < kSlots; ++i) if (m.out[i].count > 0) { if (!m.out[i].data.empty()) fprintf(f, " %s|%s 1", tables.name(m.out[i].item).c_str(), m.out[i].data.c_str()); /* (#1446 2c) */ else fprintf(f, " %s %d", tables.name(m.out[i].item).c_str(), m.out[i].count); }
            if (m.fuel.count > 0) fprintf(f, " f %s %d", tables.name(m.fuel.item).c_str(), m.fuel.count);
            if (m.water > 0) fprintf(f, " w %d", m.water);
            if (m.heat > 0) fprintf(f, " h %d", m.heat);
            if (m.fire != kFireAuto) fprintf(f, " fire %d", m.fire);
            if (m.mode != 0) fprintf(f, " mode %d", m.mode); // (#1271 L6c)
            if (m.pick >= 0) fprintf(f, " pick %s", tables.name(m.pick).c_str()); // (#1271 L6c)
            // (#1427 L7c) a disassembly: its job (slot, progress), the chance rolls' state, the fractions carried
            if (m.disSlot >= 0) fprintf(f, " dis %d %d", m.disSlot, m.progress);
            if (m.disRng != 0) fprintf(f, " drng %u", m.disRng);
            for (const auto& [id, c] : m.disCarry) if (c > 0) fprintf(f, " dc %s %d", tables.name(id).c_str(), c);
            for (const Machine::MbJob& j : m.jobs) // (#1429 PR 4) a structure's running recipes (their inputs are claimed)
                if (j.recipe >= 0 && !tables.recipes[j.recipe].out.empty()) fprintf(f, " mbjob %d %s %d", j.recipe, tables.name(tables.recipes[j.recipe].out[0].item).c_str(), j.progress);
            // (#1302) a plot in the summary sim: its caps (recipe index with its first output, to check it on load), so loaded
            // far it goes on in the summary at once instead of measuring for a minute
            if (summary) {
                fprintf(f, " sum");
                for (const auto& [r, c] : m.capPerWindow)
                    if (c > 0 && r >= 0 && r < (int)tables.recipes.size() && !tables.recipes[r].out.empty()) fprintf(f, " cap %d %s %d", r, tables.name(tables.recipes[r].out[0].item).c_str(), c);
                if (m.euCap > 0) fprintf(f, " ecap %lld", m.euCap);
            }
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
        while (ss >> w) { if (const char* now = RetiredItem(w); *now) w = now; rest.push_back(w); } // #1090: item ids only match it
        pendingStates[{ bodyId, plot, part }] = rest;
    }
}

void Factory::Load(FILE* f, uint64_t bodyId) { LoadSaved(bodyId, ReadSaved(f)); }

void Factory::MovePlot(uint64_t fromBody, uint32_t fromPlot, uint64_t toBody, uint32_t toPlot) {
    std::lock_guard<std::mutex> q(queueMutex);
    moveQueue.push_back({ fromBody, fromPlot, toBody, toPlot });
}

bool Factory::Accepts(PartType type, const std::string& item, int mode) {
    LoadRecipes();
    if (!item.empty() && tables.fluid(tables.item(item))) return false; // (#1254) fluid (essence) never goes in as an item, not even into a chest
    // Storage takes anything (#283): looked up first, items no recipe or fuel list had interned (food, tools, gems,
    // dyes...) were refused before the storage rule was reached. Insert interns the name.
    if (Def(type).power == Power::Storage) return !item.empty();
    const int id = tables.find(item);
    if (mode == 1 && type == PartType::LvAssembler) return Disassembles(id); // (#1427 L7c)
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

std::string Factory::Workpiece(uint64_t body, uint32_t plot, int part) {
    std::lock_guard<std::mutex> v(viewMutex);
    if (remoteMode) return std::string(); // (#1613) not replicated yet: a joined client sees a bare anvil (the MP pass, #1691)
    auto it = anvilView.find(std::make_tuple(body, plot, part));
    return it == anvilView.end() ? std::string() : tables.name(it->second);
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
void Factory::SetCameraNear(bool counts) { cameraNear = counts; } // (#1302)
std::vector<int> Factory::WorkingParts(uint64_t body, uint32_t plot) {
    std::lock_guard<std::mutex> v(viewMutex);
    auto it = workingByPlot.find({ body, plot });
    return it == workingByPlot.end() ? std::vector<int>{} : it->second;
}
std::vector<int> Factory::StalledParts(uint64_t body, uint32_t plot) {
    std::lock_guard<std::mutex> v(viewMutex);
    if (remoteMode) return {}; // (not replicated yet: a joined client hears only working and idle)
    auto it = stalledByPlot.find({ body, plot });
    return it == stalledByPlot.end() ? std::vector<int>{} : it->second;
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
void (*Factory::CursorSink())(const std::string&, int, float, const std::string&) { return cursorSink; }

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

// (#997, user decision) the determinism gate for parallel plots: 1, 2 and 4 workers give the same factory, bit for bit
// (#1015: the issue's 1/2/N).
TEST(Factory, PlotsTickedOnFourWorkersMatchOneWorker) {
    const ChainsRun one = RunChains(true, 1), two = RunChains(true, 2), four = RunChains(true, 4);
    EXPECT_GT(one.made, 0);
    EXPECT_EQ(four.made, one.made);
    EXPECT_EQ(two.made, one.made);
    EXPECT_EQ(four.hash, one.hash) << "every machine and network the same whatever the worker count";
    EXPECT_EQ(two.hash, one.hash) << "every machine and network the same whatever the worker count";
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

// (#1271 L3) Steam on rotation and the bridge to LV, on one plot:
//  - a boiler, a steam engine against it, a shaft to a mechanical press (iron ingots) and on to a dynamo cabled to an electric
//    furnace (iron ore): the boiler warms up, the engine turns the press and the dynamo, the dynamo powers the furnace;
//  - an LV generator cabled to an electric motor turning a shaft into a second press: EU back into rotation.
TEST(Factory, ASteamEngineDrivesAPressAndADynamoAndAMotorTurnsAShaft) {
    LoadRecipes();
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    LayoutMsg L{};
    L.body = 7; L.plot = 90; L.version = 1;
    int part = 0;
    const int xAxis = RotationPlot_OrientFor({ 1, 0, 0 });
    auto add = [&](PartType t, Int3 at, int orient = 0) {
        const Part pt = RotationPlot_PartAt(t, at, orient);
        std::vector<std::pair<Int3, CellShape>> fp;
        Plot::Footprint(pt, fp);
        if (IsMachine(t)) {
            std::vector<Int3> cells;
            for (const auto& c : fp) cells.push_back(c.first);
            L.machines.push_back({ part, t }); L.machineCells.push_back(cells); L.machineVariant.push_back("");
        }
        RotationPartMsg rm;
        if (RotationPlot_Describe(pt, part, rm)) L.rotation.push_back(rm);
        return part++;
    };
    const int boiler = add(PartType::Boiler, { 0, 0, 0 }), engine = add(PartType::SteamEngine, { 2, 0, 0 });
    add(PartType::Shaft, { 4, 0, 0 }, xAxis);
    const int press = add(PartType::Press, { 5, 0, 0 }), dynamo = add(PartType::Dynamo, { 7, 0, 0 });
    L.cables = { { 9, 0, 0 }, { 10, 0, 0 } };
    const int furnace = add(PartType::Furnace, { 11, 0, 0 });
    const int gen = add(PartType::Generator, { 0, 0, 10 });
    L.cables.push_back({ 2, 0, 10 });
    const int motor = add(PartType::Motor, { 3, 0, 10 });
    add(PartType::Shaft, { 5, 0, 10 }, xAxis);
    const int press2 = add(PartType::Press, { 6, 0, 10 });
    // (#1614) a lathe on its own motor (one motor drives one machine at full speed)
    const int gen3 = add(PartType::Generator, { 0, 0, 20 });
    L.cables.push_back({ 2, 0, 20 });
    add(PartType::Motor, { 3, 0, 20 });
    add(PartType::Shaft, { 5, 0, 20 }, xAxis);
    const int lathe = add(PartType::Lathe, { 6, 0, 20 });
    ApplyLayout(L);
    PlotSim& ps = *FindPlot(7, 90);
    auto M = [&](int p) -> Machine& { return ps.machines[ps.byPart.at(p)]; };
    int buckets = 0;
    InsertInto(M(boiler), tables.item("coal"), 64, true, buckets);
    M(boiler).water = kWaterTank;
    InsertInto(M(press), tables.item("iron_ingot"), 8, false, buckets);
    InsertInto(M(furnace), tables.item("iron_ore"), 8, false, buckets);
    InsertInto(M(gen), tables.item("coal"), 16, true, buckets);
    InsertInto(M(press2), tables.item("iron_ingot"), 8, false, buckets);
    InsertInto(M(lathe), tables.item("iron_ingot"), 8, false, buckets);
    InsertInto(M(gen3), tables.item("coal"), 16, true, buckets);
    EXPECT_EQ(M(engine).boiler, ps.byPart.at(boiler)) << "the engine finds the boiler it touches";
    long long made = 0; int active = 0;
    bool engineTurned = false, dynamoWorked = false, motorTurned = false;
    for (int t = 0; t < 3000; ++t) {
        TickPlot(ps, 1, made, active);
        engineTurned |= ps.rot.SpeedMilli(engine) != 0;
        dynamoWorked |= M(dynamo).status == MachineStatus::Working;
        motorTurned |= ps.rot.SpeedMilli(motor) != 0;
    }
    EXPECT_TRUE(engineTurned) << "the boiler's steam turned the engine";
    EXPECT_TRUE(dynamoWorked) << "the dynamo made EU for the furnace";
    EXPECT_TRUE(motorTurned) << "the generator's EU turned the motor";
    EXPECT_GT(CountItems(M(press).out, kSlots, tables.item("iron_plate")), 0) << "the steam-driven press pressed plates";
    EXPECT_GE(CountItems(M(lathe).out, kSlots, tables.item("iron_rod")), 3) << "(#1614) the lathe on the motor's shaft turned an ingot into 3 rods";
    EXPECT_GT(CountItems(M(furnace).out, kSlots, tables.item("iron_ingot")), 0) << "the dynamo-powered furnace smelted";
    EXPECT_GT(CountItems(M(press2).out, kSlots, tables.item("iron_plate")), 0) << "the motor-driven press pressed plates";
}

// (#1734) The rotation-age machines on shaft power: the trip hammer hammers blooms into ingots, the mixer stirs copper and
// tin (dust or crushed ore) 3 to 1 into bronze dust. Each on its own motor (64 RPM, over the mixer's 32 RPM minimum).
TEST(Factory, ATripHammerWorksBloomsAndAMixerMakesBronzeDust) {
    LoadRecipes();
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    LayoutMsg L{};
    L.body = 7; L.plot = 91; L.version = 1;
    int part = 0;
    const int xAxis = RotationPlot_OrientFor({ 1, 0, 0 });
    auto add = [&](PartType t, Int3 at, int orient = 0) {
        const Part pt = RotationPlot_PartAt(t, at, orient);
        std::vector<std::pair<Int3, CellShape>> fp;
        Plot::Footprint(pt, fp);
        if (IsMachine(t)) {
            std::vector<Int3> cells;
            for (const auto& c : fp) cells.push_back(c.first);
            L.machines.push_back({ part, t }); L.machineCells.push_back(cells); L.machineVariant.push_back("");
        }
        RotationPartMsg rm;
        if (RotationPlot_Describe(pt, part, rm)) L.rotation.push_back(rm);
        return part++;
    };
    int machine[3], gen[3];
    const PartType kinds[3] = { PartType::TripHammer, PartType::Mixer, PartType::Mixer };
    for (int i = 0; i < 3; ++i) { // one line each: generator, cable, motor, shaft, machine
        const int z = 10 * i;
        gen[i] = add(PartType::Generator, { 0, 0, z });
        L.cables.push_back({ 2, 0, z });
        add(PartType::Motor, { 3, 0, z });
        add(PartType::Shaft, { 5, 0, z }, xAxis);
        machine[i] = add(kinds[i], { 6, 0, z });
    }
    ApplyLayout(L);
    PlotSim& ps = *FindPlot(7, 91);
    auto M = [&](int p) -> Machine& { return ps.machines[ps.byPart.at(p)]; };
    int buckets = 0;
    for (int i = 0; i < 3; ++i) InsertInto(M(gen[i]), tables.item("coal"), 16, true, buckets);
    InsertInto(M(machine[0]), tables.item("iron_bloom"), 4, false, buckets);
    InsertInto(M(machine[1]), tables.item("copper_dust"), 3, false, buckets);
    InsertInto(M(machine[1]), tables.item("tin_dust"), 1, false, buckets);
    InsertInto(M(machine[2]), tables.item("copper_crushed"), 3, false, buckets);
    InsertInto(M(machine[2]), tables.item("tin_crushed"), 1, false, buckets);
    long long made = 0; int active = 0;
    for (int t = 0; t < 1200; ++t) TickPlot(ps, 1, made, active);
    EXPECT_EQ(CountItems(M(machine[0]).out, kSlots, tables.item("iron_ingot")), 4) << "the trip hammer worked all 4 blooms (80 ticks each)";
    EXPECT_EQ(CountItems(M(machine[1]).out, kSlots, tables.item("bronze_dust")), 4) << "copper and tin dust 3:1 -> 4 bronze dust";
    EXPECT_EQ(CountItems(M(machine[2]).out, kSlots, tables.item("bronze_dust")), 4) << "crushed copper and tin 3:1 -> 4 bronze dust";
}

// (#1734 PR 2) The rotary sieve's outputs by chance ("iron_ore 8%") come from the machine's carried hundredths: 12 gravel
// carry 96 hundredths of an iron ore and give none, the 13th gives the first; every 100 gravel give exactly 8 iron, 6 copper
// and 4 tin ore with their 100 sand, and nothing is left carried.
TEST(Factory, ARotarySieveGivesItsChanceOutputsExactly) {
    LoadRecipes();
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    LayoutMsg L{};
    L.body = 7; L.plot = 92; L.version = 1;
    int part = 0;
    auto add = [&](PartType t, Int3 at, int orient = 0) {
        const Part pt = RotationPlot_PartAt(t, at, orient);
        std::vector<std::pair<Int3, CellShape>> fp;
        Plot::Footprint(pt, fp);
        if (IsMachine(t)) {
            std::vector<Int3> cells;
            for (const auto& c : fp) cells.push_back(c.first);
            L.machines.push_back({ part, t }); L.machineCells.push_back(cells); L.machineVariant.push_back("");
        }
        RotationPartMsg rm;
        if (RotationPlot_Describe(pt, part, rm)) L.rotation.push_back(rm);
        return part++;
    };
    const int gen = add(PartType::Generator, { 0, 0, 0 });
    L.cables.push_back({ 2, 0, 0 });
    add(PartType::Motor, { 3, 0, 0 });
    add(PartType::Shaft, { 5, 0, 0 }, RotationPlot_OrientFor({ 1, 0, 0 }));
    const int sieve = add(PartType::RotarySieve, { 6, 0, 0 });
    ApplyLayout(L);
    PlotSim& ps = *FindPlot(7, 92);
    Machine& m = ps.machines[ps.byPart.at(sieve)];
    int buckets = 0;
    InsertInto(ps.machines[ps.byPart.at(gen)], tables.item("coal"), 64, true, buckets);
    const int gravel = tables.item("gravel"), sand = tables.item("sand"), iron = tables.item("iron_ore");
    long long made = 0; int active = 0;
    auto runUntil = [&](int sandWanted) { for (int t = 0; t < 20000 && CountItems(m.out, kSlots, sand) < sandWanted; ++t) TickPlot(ps, 1, made, active); };
    InsertInto(m, gravel, 12, false, buckets);
    runUntil(12);
    ASSERT_EQ(CountItems(m.out, kSlots, sand), 12) << "12 gravel sifted";
    EXPECT_EQ(CountItems(m.out, kSlots, iron), 0) << "12 x 8% = 96 hundredths: no iron ore yet";
    EXPECT_EQ(m.disCarry[iron], 96);
    InsertInto(m, gravel, 88, false, buckets);
    runUntil(100);
    ASSERT_EQ(CountItems(m.out, kSlots, sand), 100);
    EXPECT_EQ(CountItems(m.out, kSlots, iron), 8) << "8% of 100 gravel";
    EXPECT_EQ(CountItems(m.out, kSlots, tables.item("copper_ore")), 6);
    EXPECT_EQ(CountItems(m.out, kSlots, tables.item("tin_ore")), 4);
    EXPECT_TRUE(m.disCarry.empty()) << "100 gravel leave nothing carried";
}

// (#1261, S1) A plot that leaves "near" measures for 60 s on the 20-tick step, then enters Summary with caps from what each
// machine crafted; a player's insert while far measures again; coming near ticks every tick; SUMMARY=0 never leaves the
// 20-tick step.
TEST(Factory, AFarPlotMeasuresThenEntersSummaryWithItsCaps) {
    LoadRecipes();
    const bool was = summaryOn;
    for (const bool on : { true, false }) {
        summaryOn = on;
        plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
        std::vector<std::array<int, 4>> parts;
        ApplyLayout(TestChains(1, 10, parts));
        PlotSim& ps = *FindPlot(7, 1);
        for (size_t g = 0; g < parts.size(); ++g) { // fed groups, every third starved (as RunChains)
            int buckets = 0;
            InsertInto(ps.machines[ps.byPart.at(parts[g][0])], tables.item("coal"), 16, false, buckets);
            if (g % 3 == 2) continue;
            InsertInto(ps.machines[ps.byPart.at(parts[g][1])], tables.item("stone"), 32, false, buckets);
        }
        long long made = 0; int active = 0, steps = 0;
        auto tick = [&](long long t) { UpdateMode(ps, t); if (const int s = StepOf(ps, 0, t)) { TickPlot(ps, s, made, active); ++steps; } };
        SetNear(ps, false, 0);
        EXPECT_EQ(ps.mode, PlotSim::Mode::Measuring);
        for (long long t = 0; t < kMeasureTicks; ++t) tick(t);
        EXPECT_EQ(steps, (int)(kMeasureTicks / kTicksPerSecond)) << "measuring: one 20-tick step a second";
        tick(kMeasureTicks);
        if (!on) { EXPECT_EQ(ps.mode, PlotSim::Mode::Measuring) << "SUMMARY=0: today's far step forever"; continue; }
        ASSERT_EQ(ps.mode, PlotSim::Mode::Summary);
        const Machine& fed = ps.machines[ps.byPart.at(parts[0][1])];
        const Machine& starved = ps.machines[ps.byPart.at(parts[2][1])];
        ASSERT_EQ(fed.capPerWindow.size(), 1u);
        const int cap = fed.capPerWindow.begin()->second;
        EXPECT_GE(cap, 10) << "a fed macerator: 100 ticks a craft, ~12 in a minute";
        EXPECT_LE(cap, 12);
        EXPECT_TRUE(starved.capPerWindow.empty()) << "a starved one made nothing, so it has no cap";
        // (this loop does not advance it: only the sim loop and the oracle do) - nothing ticks it either
        const uint64_t before = StateHash(ps);
        const int stepsBefore = steps;
        for (long long t = kMeasureTicks + 1; t < kMeasureTicks + 200; ++t) tick(t);
        EXPECT_EQ(steps, stepsBefore);
        EXPECT_EQ(StateHash(ps), before);
        BeginMeasuring(ps, kMeasureTicks + 200); // what a player's insert while far does
        EXPECT_EQ(ps.mode, PlotSim::Mode::Measuring);
        EXPECT_TRUE(ps.machines[ps.byPart.at(parts[0][1])].crafts.empty()) << "a new window";
        SetNear(ps, true, kMeasureTicks + 201);
        EXPECT_EQ(ps.mode, PlotSim::Mode::Near);
        EXPECT_EQ(StepOf(ps, 0, kMeasureTicks + 202), 1);
    }
    summaryOn = was;
    plots.clear(); plotIndexDirty = true;
}

// ---- (#1261, S2) The oracle (the factory doc, section 3): the summary never makes more than the full sim ----
namespace {
struct OraclePart { PartType type; int x, z; };
LayoutMsg OracleLayout(uint32_t plotId, const std::vector<OraclePart>& machines, const std::vector<Int3>& conveyors, const std::vector<Int3>& cables) {
    LayoutMsg L{};
    L.body = 9; L.plot = plotId; L.version = 1;
    for (size_t i = 0; i < machines.size(); ++i) {
        std::vector<Int3> cells;
        for (int dz = 0; dz < 2; ++dz) for (int dy = 0; dy < 2; ++dy) for (int dx = 0; dx < 2; ++dx) cells.push_back({ machines[i].x + dx, dy, machines[i].z + dz });
        L.machines.push_back({ (int)i, machines[i].type }); L.machineCells.push_back(cells); L.machineVariant.push_back("");
    }
    L.conveyors = conveyors; L.cables = cables;
    return L;
}
std::vector<Int3> Row(int x0, int x1, int z) { std::vector<Int3> r; for (int x = x0; x <= x1; ++x) r.push_back({ x, 0, z }); return r; }

// Four plots: (1) a steady chain, chest of stone -> macerator -> chest; (2) a power-limited line, one generator for six
// macerators wanting 48 EU/t of its 32, their ore in their own slots (which run dry: an event); (3) a chest with only 20
// stone (it runs dry: an event); (4) a hand-fed furnace with no conveyor. Totals made per plot and item over `ticks`.
// far: the plots run near for 2 minutes, then leave "near", measure for 60 s, and advance by the summary every 10 s.
// maceratorTicks > 0 (#1403): every macerator recipe takes that many ticks (long crafts against the measuring window)
std::vector<std::map<std::string, long long>> OracleRun(bool far, long long ticks, int maceratorTicks = 0) {
    LoadRecipes();
    std::vector<std::pair<Recipe*, int>> restore; // (the tables outlive the test: LoadRecipes doesn't reload them)
    if (maceratorTicks > 0) for (Recipe& r : tables.recipes) if (r.machine == PartType::Macerator) { restore.push_back({ &r, r.ticks }); r.ticks = maceratorTicks; }
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    const std::vector<Int3> cable12 = Row(0, 12, -1);
    ApplyLayout(OracleLayout(1, { { PartType::Generator, 0, 0 }, { PartType::Chest, 3, 0 }, { PartType::Macerator, 6, 0 }, { PartType::Chest, 9, 0 } }, Row(3, 10, 2), Row(0, 7, -1)));
    ApplyLayout(OracleLayout(2, { { PartType::Generator, 0, 0 }, { PartType::Macerator, 3, 0 }, { PartType::Macerator, 6, 0 }, { PartType::Macerator, 9, 0 },
                                  { PartType::Macerator, 12, 0 }, { PartType::Macerator, 15, 0 }, { PartType::Macerator, 18, 0 } }, {}, Row(0, 19, -1)));
    ApplyLayout(OracleLayout(3, { { PartType::Generator, 0, 0 }, { PartType::Chest, 3, 0 }, { PartType::Macerator, 6, 0 }, { PartType::Chest, 9, 0 } }, Row(3, 10, 2), Row(0, 7, -1)));
    ApplyLayout(OracleLayout(4, { { PartType::Generator, 0, 0 }, { PartType::Furnace, 3, 0 } }, {}, Row(0, 4, -1)));
    std::vector<PlotSim*> all;
    for (uint32_t p = 1; p <= 4; ++p) all.push_back(FindPlot(9, p));
    int b = 0;
    auto put = [&](PlotSim& ps, int part, const char* item, int n, bool fuel = false) { InsertInto(ps.machines[ps.byPart.at(part)], tables.item(item), n, fuel, b); };
    put(*all[0], 0, "coal", 64, true); put(*all[0], 1, "stone", 640);
    put(*all[1], 0, "coal", 64, true); for (int k = 1; k <= 6; ++k) put(*all[1], k, "iron_ore", 64);
    put(*all[2], 0, "coal", 64, true); put(*all[2], 1, "stone", 20);
    put(*all[3], 0, "coal", 64, true); put(*all[3], 1, "sand", 64);
    long long madeAll = 0; int active = 0;
    constexpr long long kWarm = (2 * 60 + 7) * kTicksPerSecond; // B runs near for 2 min 7 s (not on a 30 s bucket boundary), then the player leaves
    for (long long t = 0; t < ticks; ++t)
        for (size_t i = 0; i < all.size(); ++i) {
            PlotSim& ps = *all[i];
            if (far && t == kWarm) SetNear(ps, false, t);
            if (ps.mode == PlotSim::Mode::Near) RotateNearBuckets(ps, t);
            UpdateMode(ps, t);
            if (ps.mode == PlotSim::Mode::Summary && (t + (long long)i) % kAdvanceEvery == 0) { AdvanceSummary(ps, ps.lastAdvance, t); ps.lastAdvance = t; }
            if (const int s = StepOf(ps, i, t)) TickPlot(ps, s, madeAll, active);
        }
    std::vector<std::map<std::string, long long>> made(all.size());
    for (size_t i = 0; i < all.size(); ++i) {
        PlotSim& ps = *all[i];
        if (ps.mode == PlotSim::Mode::Summary) AdvanceSummary(ps, ps.lastAdvance, ticks);
        for (const auto& [item, s] : ps.stats.items) if (s.made.total + s.made.acc > 0) made[i][tables.name(item)] = s.made.total + s.made.acc;
    }
    plots.clear(); plotIndexDirty = true;
    for (const auto& [r, t] : restore) r->ticks = t;
    return made;
}
} // namespace

// (#1261) A plot that leaves "near" at any time after 30 s keeps the caps of its last 30-60 s near (the first version kept
// them only when leaving exactly on a 30 s boundary: the in-game run used the measuring window alone and made 82.5%).
TEST(Factory, LeavingNearAtAnyTimeKeepsTheNearCaps) {
    LoadRecipes();
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    ApplyLayout(OracleLayout(1, { { PartType::Generator, 0, 0 }, { PartType::Chest, 3, 0 }, { PartType::Macerator, 6, 0 }, { PartType::Chest, 9, 0 } }, Row(3, 10, 2), Row(0, 7, -1)));
    PlotSim& ps = *FindPlot(9, 1);
    int b = 0;
    InsertInto(ps.machines[ps.byPart.at(0)], tables.item("coal"), 64, true, b);
    InsertInto(ps.machines[ps.byPart.at(1)], tables.item("stone"), 640, false, b);
    long long made = 0; int active = 0;
    for (long long t = 0; t < (60 + 13) * kTicksPerSecond; ++t) { RotateNearBuckets(ps, t); TickPlot(ps, 1, made, active); } // 1 min 13 s near
    SetNear(ps, false, (60 + 13) * kTicksPerSecond);
    const Machine& mac = ps.machines[ps.byPart.at(2)];
    ASSERT_FALSE(mac.nearCap.empty()) << "the macerator's near caps are kept";
    EXPECT_GE(mac.nearCap.begin()->second, 11) << "12 crafts a minute near (100 ticks a craft), at least 11";
    EXPECT_GT(ps.machines[ps.byPart.at(0)].nearEuCap, 0) << "and the generator's EU";
    plots.clear(); plotIndexDirty = true;
}

TEST(Factory, TheSummaryNeverMakesMoreThanTheFullSimAndNotMuchLess) {
    const bool was = summaryOn;
    summaryOn = true;
    constexpr long long kTicks = 10 * 60 * kTicksPerSecond; // 10 minutes (the doc's oracle)
    const auto full = OracleRun(false, kTicks), summary = OracleRun(true, kTicks);
    summaryOn = was;
    const char* names[] = { "steady chain", "power-limited line", "dry input chest", "hand-fed furnace" };
    for (size_t p = 0; p < full.size(); ++p) {
        EXPECT_FALSE(full[p].empty()) << names[p] << ": the full sim makes something";
        for (const auto& [item, a] : full[p]) {
            const long long bsum = summary[p].count(item) ? summary[p].at(item) : 0;
            printf("oracle: %-18s %-14s full %5lld  summary %5lld  (%.1f%%)\n", names[p], item.c_str(), a, bsum, 100.0 * bsum / std::max(1LL, a));
            EXPECT_LE(bsum, a) << names[p] << ", " << item << ": the summary may never make more than the full sim";
            EXPECT_GE(bsum * 100, a * 85) << names[p] << ", " << item << ": full " << a << ", summary " << bsum << " (at least 85%)";
        }
        for (const auto& [item, bsum] : summary[p]) EXPECT_TRUE(full[p].count(item)) << names[p] << ": the summary made " << item << " the full sim never did";
    }
}

// (#1403) The same oracle with 300-tick macerator crafts (15 s: two to four in a 30-60 s window), as Progression's ladder
// had: the summary made one craft more per machine than the full sim on the power-limited line (450 vs 432). It may never.
TEST(Factory, TheSummaryNeverMakesMoreThanTheFullSimWithLongCrafts) {
    const bool was = summaryOn;
    summaryOn = true;
    constexpr long long kTicks = 10 * 60 * kTicksPerSecond;
    const auto full = OracleRun(false, kTicks, 300), summary = OracleRun(true, kTicks, 300);
    summaryOn = was;
    const char* names[] = { "steady chain", "power-limited line", "dry input chest", "hand-fed furnace" };
    for (size_t p = 0; p < full.size(); ++p)
        for (const auto& [item, a] : full[p]) {
            const long long bsum = summary[p].count(item) ? summary[p].at(item) : 0;
            printf("oracle (300-tick crafts): %-18s %-14s full %5lld  summary %5lld\n", names[p], item.c_str(), a, bsum);
            EXPECT_LE(bsum, a) << names[p] << ", " << item << ": the summary may never make more than the full sim";
        }
}

// (#1266) Six 8 EU/t macerators on one 32 EU/t generator, near the player (1-tick steps): each gets 5 EU of 8 a tick. It
// stalled (5 / 8 = 0 ticks, every tick); now the share adds up, and the line runs at about 5/8 of its full speed.
TEST(Factory, AnOverloadedPowerNetSlowsItsMachinesNearThePlayerAndDoesNotStallThem) {
    LoadRecipes();
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    ApplyLayout(OracleLayout(2, { { PartType::Generator, 0, 0 }, { PartType::Macerator, 3, 0 }, { PartType::Macerator, 6, 0 }, { PartType::Macerator, 9, 0 },
                                  { PartType::Macerator, 12, 0 }, { PartType::Macerator, 15, 0 }, { PartType::Macerator, 18, 0 } }, {}, Row(0, 19, -1)));
    PlotSim& ps = *FindPlot(9, 2);
    int b = 0;
    InsertInto(ps.machines[ps.byPart.at(0)], tables.item("coal"), 64, true, b);
    for (int k = 1; k <= 6; ++k) InsertInto(ps.machines[ps.byPart.at(k)], tables.item("iron_ore"), 64, false, b);
    long long made = 0; int active = 0;
    constexpr long long kTicks = 2 * 60 * kTicksPerSecond;
    for (long long t = 0; t < kTicks; ++t) TickPlot(ps, 1, made, active);
    // (#996, ladder L4) the macerator's iron recipe from the table (2 ore -> 3 crushed in 150 ticks; it was 1 -> 2), unlimited power
    int craftTicks = 150, perCraft = 2;
    for (const auto& r : tables.recipes) if (r.machine == PartType::Macerator && !r.in.empty() && r.in[0].item == tables.item("iron_ore") && !r.out.empty()) { craftTicks = r.ticks; perCraft = r.out[0].count; break; }
    const long long full = kTicks / craftTicks * 6 * perCraft;
    // (each gets floor(8 x 32 / 48) = 5 EU of 8: 5/8 of full speed; a craft's last part-tick of EU is not carried over)
    EXPECT_GE(made * 100, full * 5 / 8 * 85) << "about 5/8 of the full rate (at least 85% of it): made " << made << " of " << full;
    EXPECT_LE(made, full * 32 / 48) << "and no more than the share allows";
    plots.clear(); plotIndexDirty = true;
}

// (#1620 C2b) Cables carry what their copper allows (32 EU/t a path). Two LV generators (64 EU/t) and eight macerators wanting
// 64 EU/t in all: a straight line of cable carries 32 of it.
namespace {
struct CableRun { long long used = 0, made = 0; };
CableRun RunCableLine(uint32_t plotId, const std::vector<Int3>& cables, int ticks) {
    LoadRecipes();
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    std::vector<OraclePart> parts = { { PartType::Generator, 1, 0 }, { PartType::Generator, 4, 0 } };
    for (int k = 0; k < 8; ++k) parts.push_back({ PartType::Macerator, 7 + 3 * k, 0 });
    ApplyLayout(OracleLayout(plotId, parts, {}, cables));
    PlotSim& ps = *FindPlot(9, plotId);
    int b = 0;
    for (int k = 0; k < 2; ++k) InsertInto(ps.machines[ps.byPart.at(k)], tables.item("coal"), 64, true, b);
    for (int k = 2; k < 10; ++k) InsertInto(ps.machines[ps.byPart.at(k)], tables.item("iron_ore"), 64, false, b);
    long long made = 0; int active = 0;
    for (int t = 0; t < ticks; ++t) TickPlot(ps, 1, made, active);
    CableRun r;
    for (const auto& kv : ps.stats.euUsed) r.used += kv.second.total + kv.second.acc;
    for (const auto& kv : ps.stats.euMade) r.made += kv.second.total + kv.second.acc;
    plots.clear(); plotIndexDirty = true;
    return r;
}
} // namespace

TEST(Factory, ACableCarriesItsCapacity) {
    constexpr int kT = 10 * kTicksPerSecond;
    const CableRun r = RunCableLine(5, Row(0, 31, -1), kT);
    EXPECT_LE(r.used, 32LL * kT + 64) << "one cable carries 32 EU/t, not the 64 the generators could make: used " << r.used;
    EXPECT_GE(r.used, 32LL * kT * 9 / 10) << "and it does carry that: used " << r.used;
    EXPECT_LE(r.made, 32LL * kT + 64) << "the generators burn only for what is carried: made " << r.made;
}

TEST(Factory, ParallelCablesCarryMore) {
    constexpr int kT = 10 * kTicksPerSecond;
    std::vector<Int3> cables = Row(-1, 31, -1);
    for (const Int3& c : Row(-1, 31, 2)) cables.push_back(c);
    for (int z = 0; z <= 1; ++z) { cables.push_back({ -1, 0, z }); cables.push_back({ 31, 0, z }); }
    const CableRun r = RunCableLine(6, cables, kT);
    EXPECT_GT(r.used, 32LL * kT * 3 / 2) << "two lines carry about 64 EU/t: used " << r.used;
    EXPECT_LE(r.used, 64LL * kT + 64) << "but never more than the generators make: used " << r.used;
}

// (#1847) A machine whose craft is done but can't get out (outputs full) draws no power: its generator stops burning.
TEST(Factory, AMachineWithFullOutputsDrawsNoPower) {
    LoadRecipes();
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    ApplyLayout(OracleLayout(2, { { PartType::Generator, 0, 0 }, { PartType::Macerator, 3, 0 } }, {}, Row(0, 4, -1)));
    PlotSim& ps = *FindPlot(9, 2);
    int b = 0;
    InsertInto(ps.machines[ps.byPart.at(0)], tables.item("coal"), 64, true, b);
    Machine& mac = ps.machines[ps.byPart.at(1)];
    InsertInto(mac, tables.item("iron_ore"), 64, false, b);
    for (Stack& o : mac.out) { o = Stack{}; o.item = tables.item("dirt"); o.count = 64; } // no room for what it makes
    long long made = 0; int active = 0;
    for (int t = 0; t < 400; ++t) TickPlot(ps, 1, made, active); // the first craft finishes its work and waits
    ASSERT_GE(mac.recipe, 0);
    ASSERT_EQ(mac.progress, tables.recipes[mac.recipe].ticks) << "done, waiting on the outputs";
    const Machine& gen = ps.machines[ps.byPart.at(0)];
    const int fuel = gen.fuel.count, energy = gen.energy;
    for (int t = 0; t < 2000; ++t) TickPlot(ps, 1, made, active);
    EXPECT_EQ(made, 0) << "nothing fits";
    EXPECT_EQ(gen.fuel.count, fuel) << "no coal burnt for a machine that can't craft";
    EXPECT_LE(gen.energy, energy);
    EXPECT_EQ(mac.received, 0);
    plots.clear(); plotIndexDirty = true;
}

// (#1302) A plot saved in the summary sim: the save first brings its chests up to now, and loaded far it goes on in the
// summary at once with the caps it had (no minute of measuring).
TEST(Factory, APlotSavedInTheSummaryGoesOnInItWhenLoadedFar) {
    LoadRecipes();
    const bool was = summaryOn;
    summaryOn = true;
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    const LayoutMsg layout = OracleLayout(1, { { PartType::Generator, 0, 0 }, { PartType::Chest, 3, 0 }, { PartType::Macerator, 6, 0 }, { PartType::Chest, 9, 0 } }, Row(3, 10, 2), Row(0, 7, -1));
    ApplyLayout(layout);
    PlotSim* ps = FindPlot(9, 1);
    int b = 0;
    InsertInto(ps->machines[ps->byPart.at(0)], tables.item("coal"), 64, true, b);
    InsertInto(ps->machines[ps->byPart.at(1)], tables.item("stone"), 640, false, b);
    long long made = 0; int active = 0;
    const long long leave = 70 * kTicksPerSecond, summaryAt = leave + kMeasureTicks;
    for (long long t = 0; t < summaryAt + 400; ++t) {
        if (t == leave) SetNear(*ps, false, t);
        if (ps->mode == PlotSim::Mode::Near) RotateNearBuckets(*ps, t);
        UpdateMode(*ps, t);
        if (const int s = StepOf(*ps, 0, t)) TickPlot(*ps, s, made, active);
    }
    ASSERT_EQ(ps->mode, PlotSim::Mode::Summary);
    const std::map<int, int> caps = ps->machines[ps->byPart.at(2)].capPerWindow;
    const long long euCap = ps->machines[ps->byPart.at(0)].euCap;
    ASSERT_FALSE(caps.empty());
    auto chestSand = [&](PlotSim& p) { long long n = 0; for (int part : { 1, 3 }) { const Machine& c = p.machines[p.byPart.at(part)]; n += CountItems(c.in.data(), kChestSlots, tables.item("sand")); } return n; }; // (both chests: a conveyor puts sand in either)
    const long long sandBefore = chestSand(*ps);
    tickCount = summaryAt + 400 + 1200; // a minute later, the game saves
    const std::string path = (std::filesystem::temp_directory_path() / "raylibgame_factory_save_test.txt").string(); // (tmpfile() needs the drive root on Windows)
    FILE* f = fopen(path.c_str(), "w+");
    ASSERT_TRUE(f);
    Factory::Save(f, 9);
    const long long sandSaved = chestSand(*ps);
    EXPECT_GT(sandSaved, sandBefore) << "the save brought the far chain up to now first";
    rewind(f);
    const std::vector<std::string> lines = Factory::ReadSaved(f);
    fclose(f);
    std::remove(path.c_str());
    ASSERT_EQ(lines.size(), ps->machines.size());
    // Load into a fresh factory: the layout arrives, the plot is far from the start.
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    Factory::LoadSaved(9, lines);
    ApplyLayout(layout);
    ps = FindPlot(9, 1);
    ASSERT_TRUE(ps);
    EXPECT_EQ(chestSand(*ps), sandSaved) << "the chest as saved";
    SetNear(*ps, false, 5);
    EXPECT_EQ(ps->mode, PlotSim::Mode::Summary) << "it goes on in the summary at once";
    EXPECT_EQ(ps->machines[ps->byPart.at(2)].capPerWindow, caps) << "with the caps it had";
    EXPECT_EQ(ps->machines[ps->byPart.at(0)].euCap, euCap);
    // An old save (no caps) measures first, as before.
    std::vector<std::string> old = lines;
    for (std::string& l : old) { const size_t at = l.find(" sum"); if (at != std::string::npos) l.erase(at); }
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    Factory::LoadSaved(9, old);
    ApplyLayout(layout);
    ps = FindPlot(9, 1);
    SetNear(*ps, false, 5);
    EXPECT_EQ(ps->mode, PlotSim::Mode::Measuring) << "a save without caps measures first";
    summaryOn = was;
    tickCount = 0;
    plots.clear(); plotIndexDirty = true;
}

// (#1102) A chest keeps fractions as the player inventory does: 2.5 + 0.4 dirt is 2.9 in its slot; a slot holding only a
// fraction is not overwritten; a machine takes whole units and gives the fraction back; taking the chest's contents gives
// the fractions back. Through the factory thread's messages, as the game sends them.
TEST(Factory, AChestKeepsFractionsAndGivesThemBack) {
    LoadRecipes();
    Factory::Start();
    {
        LayoutMsg L{};
        L.body = 11; L.plot = 1; L.version = 1;
        const PartType types[2] = { PartType::Chest, PartType::Macerator };
        for (int i = 0; i < 2; ++i) {
            std::vector<Int3> cells;
            for (int dz = 0; dz < 2; ++dz) for (int dy = 0; dy < 2; ++dy) for (int dx = 0; dx < 2; ++dx) cells.push_back({ i * 3 + dx, dy, dz });
            L.machines.push_back({ i, types[i] }); L.machineCells.push_back(cells); L.machineVariant.push_back("");
        }
        std::lock_guard<std::mutex> q(queueMutex);
        layoutQueue.push_back(L);
    }
    Factory::WaitTicks(2);
    { std::lock_guard<std::mutex> v(viewMutex); results.clear(); }
    Factory::InsertSlot(11, 1, 0, "dirt", 2, 0, 0, 50);  // 2.5 dirt into slot 0
    Factory::Insert(11, 1, 0, "dirt", 0, 0, 40);          // 0.4 more (routed: onto the dirt)
    Factory::InsertSlot(11, 1, 0, "sand", 0, 1, 0, 30);   // only 0.3 sand in slot 1
    Factory::InsertSlot(11, 1, 0, "stone", 5, 1);         // 5 stone aimed at slot 1: it holds sand, so they come back
    Factory::Insert(11, 1, 1, "stone", 1, 0, 50);         // 1.5 stone into the macerator: it keeps 1, the 0.5 comes back
    Factory::View(11, 1, 0);
    Factory::WaitTicks(3);
    const MachineView v = Factory::View(11, 1, 0);
    ASSERT_GE(v.inSlots.size(), 2u);
    ASSERT_GE(v.inCenti.size(), 2u);
    EXPECT_EQ(v.inSlots[0].first, "dirt");
    EXPECT_EQ(v.inSlots[0].second * 100 + v.inCenti[0], 290) << "2.5 + 0.4 dirt";
    EXPECT_EQ(v.inSlots[1].first, "sand") << "a fraction-only slot is not overwritten";
    EXPECT_EQ(v.inSlots[1].second * 100 + v.inCenti[1], 30);
    auto back = [&](const char* item) { long long centi = 0; std::lock_guard<std::mutex> l(viewMutex); for (const ItemOut& r : results) if (tables.name(r.item & ~kHotItem) == item) centi += (long long)r.count * 100 + r.centi; return centi; };
    EXPECT_EQ(back("stone"), 500 + 50) << "the 5 stone the sand's slot refused, and the macerator's 0.5";
    { std::lock_guard<std::mutex> l(viewMutex); results.clear(); }
    Factory::Extract(11, 1, 0); // take everything out of the chest
    Factory::WaitTicks(2);
    EXPECT_EQ(back("dirt"), 290);
    EXPECT_EQ(back("sand"), 30);

    // (#1849) An insert whose plot is gone comes back whole, its fraction too.
    { std::lock_guard<std::mutex> l(viewMutex); results.clear(); }
    Factory::Insert(99, 7, 0, "dirt", 2, 0, 50);
    Factory::WaitTicks(2);
    EXPECT_EQ(back("dirt"), 250) << "2.5 dirt back";

    // (#1848) Removing the chest gives back its contents with their fractions (0.4 dirt is a slot with count 0).
    Factory::InsertSlot(11, 1, 0, "dirt", 0, 0, 0, 40);
    Factory::InsertSlot(11, 1, 0, "sand", 2, 1, 0, 50);
    Factory::WaitTicks(2);
    { std::lock_guard<std::mutex> l(viewMutex); results.clear(); }
    {
        LayoutMsg L{};
        L.body = 11; L.plot = 1; L.version = 2; // the macerator only: the chest was broken
        std::vector<Int3> cells;
        for (int dz = 0; dz < 2; ++dz) for (int dy = 0; dy < 2; ++dy) for (int dx = 0; dx < 2; ++dx) cells.push_back({ 3 + dx, dy, dz });
        L.machines.push_back({ 1, PartType::Macerator }); L.machineCells.push_back(cells); L.machineVariant.push_back("");
        std::lock_guard<std::mutex> q(queueMutex);
        layoutQueue.push_back(L);
    }
    Factory::WaitTicks(3);
    EXPECT_EQ(back("dirt"), 40) << "the chest's 0.4 dirt";
    EXPECT_EQ(back("sand"), 250) << "and its 2.5 sand";
    Factory::Stop();
    plots.clear(); plotIndexDirty = true;
    { std::lock_guard<std::mutex> l(viewMutex); results.clear(); }
}

// #1629 (playtest: "Chests with partial slots delete the item when a chest is broken"): a chest removed from the layout (broken)
// gives back every stack exactly - a fraction with its whole part (2.5 dirt), a fraction alone (0.5 sand), a whole stack (5
// stone) and an item with its own record (a composed pickaxe, #1446).
TEST(Factory, ABrokenChestGivesBackEveryStackExactly) {
    LoadRecipes();
    Factory::Start();
    LayoutMsg L{};
    L.body = 11; L.plot = 3; L.version = 1;
    std::vector<Int3> cells;
    for (int dz = 0; dz < 2; ++dz) for (int dy = 0; dy < 2; ++dy) for (int dx = 0; dx < 2; ++dx) cells.push_back({ dx, dy, dz });
    L.machines.push_back({ 0, PartType::Chest }); L.machineCells.push_back(cells); L.machineVariant.push_back("");
    { std::lock_guard<std::mutex> q(queueMutex); layoutQueue.push_back(L); }
    Factory::WaitTicks(2);
    const std::string record = "head=bronze;w=7";
    Factory::InsertSlot(11, 3, 0, "dirt", 2, 0, 0, 50);   // 2.5 dirt
    Factory::InsertSlot(11, 3, 0, "sand", 0, 1, 0, 50);   // 0.5 sand alone
    Factory::InsertSlot(11, 3, 0, "stone", 5, 2);         // 5 stone
    Factory::InsertSlot(11, 3, 0, "bronze_pickaxe", 1, 3, 0, 0, record);
    Factory::WaitTicks(3);
    { std::lock_guard<std::mutex> v(viewMutex); results.clear(); }
    LayoutMsg gone{}; // the plot without its chest: broken
    gone.body = 11; gone.plot = 3; gone.version = 2;
    { std::lock_guard<std::mutex> q(queueMutex); layoutQueue.push_back(gone); }
    Factory::WaitTicks(3);
    long long dirt = 0, sand = 0, stone = 0; int picks = 0; std::string pickRecord;
    {
        std::lock_guard<std::mutex> l(viewMutex);
        for (const ItemOut& r : results) {
            const std::string id = tables.name(r.item & ~kHotItem);
            const long long centi = (long long)r.count * 100 + r.centi;
            if (id == "dirt") dirt += centi; else if (id == "sand") sand += centi; else if (id == "stone") stone += centi;
            else if (id == "bronze_pickaxe") { picks += r.count; pickRecord = r.data; }
        }
    }
    EXPECT_EQ(dirt, 250) << "2.5 dirt (2 came back: the fraction was dropped)";
    EXPECT_EQ(sand, 50) << "0.5 sand (a fraction-only stack came back as nothing)";
    EXPECT_EQ(stone, 500);
    EXPECT_EQ(picks, 1);
    EXPECT_EQ(pickRecord, record) << "the pickaxe keeps its own record";
    Factory::Stop();
    plots.clear(); plotIndexDirty = true;
    { std::lock_guard<std::mutex> l(viewMutex); results.clear(); }
}

// (#1446 step 2b) An item with its own record (a composed tool): a chest keeps it in a slot of its own, never merged, and gives
// it back with the record; a machine refuses it (automation never moves one) unless it is an assembler set to disassemble (the
// L7c hook); the save keeps the record.
TEST(Factory, AChestKeepsAnItemsRecordAndOnlyADisassemblerTakesOne) {
    LoadRecipes();
    Factory::Start();
    LayoutMsg L{};
    L.body = 11; L.plot = 2; L.version = 1;
    const PartType types[2] = { PartType::Chest, PartType::LvAssembler };
    for (int i = 0; i < 2; ++i) {
        std::vector<Int3> cells;
        for (int dz = 0; dz < 2; ++dz) for (int dy = 0; dy < 2; ++dy) for (int dx = 0; dx < 2; ++dx) cells.push_back({ i * 3 + dx, dy, dz });
        L.machines.push_back({ i, types[i] }); L.machineCells.push_back(cells); L.machineVariant.push_back("");
    }
    { std::lock_guard<std::mutex> q(queueMutex); layoutQueue.push_back(L); }
    Factory::WaitTicks(2);
    { std::lock_guard<std::mutex> v(viewMutex); results.clear(); }
    const std::string d1 = "head=bronze;w=4", d2 = "head=bronze;w=9";
    Factory::InsertSlot(11, 2, 0, "bronze_pickaxe", 1, 3, 0, 0, d1); // into slot 3
    Factory::Insert(11, 2, 0, "bronze_pickaxe", 5, 0, 0, d2);        // routed: an empty slot of its own (count forced to 1)
    Factory::Insert(11, 2, 0, "bronze_pickaxe", 2);                  // plain ones: a slot of their own too, never onto a record
    Factory::Insert(11, 2, 1, "bronze_pickaxe", 1, 0, 0, d1);        // the assembler makes: it comes back
    Factory::WaitTicks(3);
    auto backData = [&]() { std::vector<std::string> d; std::lock_guard<std::mutex> l(viewMutex); for (const ItemOut& r : results) if (!r.data.empty()) d.push_back(r.data); return d; };
    EXPECT_EQ(backData(), std::vector<std::string>{ d1 }) << "an assembler making refuses an item with a record, and returns it whole";
    { std::lock_guard<std::mutex> l(viewMutex); results.clear(); }
    Factory::SetMode(11, 2, 1, 1);
    Factory::WaitTicks(2);
    Factory::Insert(11, 2, 1, "bronze_pickaxe", 1, 0, 0, d1);
    Factory::WaitTicks(3);
    EXPECT_TRUE(backData().empty()) << "set to disassemble it takes one (L7c)";
    Factory::View(11, 2, 1); // (the sim fills the one watched view: ask, wait, read)
    Factory::WaitTicks(3);
    const MachineView a = Factory::View(11, 2, 1);
    int inAsm = 0; for (const std::string& d : a.inData) inAsm += d == d1 ? 1 : 0;
    EXPECT_EQ(inAsm, 1);
    Factory::View(11, 2, 0);
    Factory::WaitTicks(3);
    const MachineView v = Factory::View(11, 2, 0);
    EXPECT_GE(v.inData.size(), 4u);
    if (v.inData.size() >= 4 && v.inSlots.size() >= 4) {
        EXPECT_EQ(v.inSlots[3].first, "bronze_pickaxe");
        EXPECT_EQ(v.inSlots[3].second, 1);
        EXPECT_EQ(v.inData[3], d1);
    }
    int withD2 = 0, plain = 0;
    for (size_t i = 0; i < v.inData.size(); ++i) {
        if (v.inData[i] == d2) { ++withD2; EXPECT_EQ(v.inSlots[i].second, 1) << "one item to a record"; }
        if (v.inData[i].empty() && v.inSlots[i].first == "bronze_pickaxe") plain += v.inSlots[i].second;
    }
    EXPECT_EQ(withD2, 1);
    EXPECT_EQ(plain, 2);
    // The save keeps the records; a fresh factory loads them back.
    Factory::Stop();
    const std::string path = (std::filesystem::temp_directory_path() / "raylibgame_factory_data_test.txt").string();
    FILE* f = fopen(path.c_str(), "w+");
    ASSERT_TRUE(f);
    Factory::Save(f, 11);
    rewind(f);
    const std::vector<std::string> lines = Factory::ReadSaved(f);
    fclose(f);
    std::remove(path.c_str());
    bool saved = false; for (const std::string& l : lines) saved |= l.find(" bronze_pickaxe|" + d1 + " 1") != std::string::npos;
    EXPECT_TRUE(saved) << "a record is saved as name|record";
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    Factory::LoadSaved(11, lines);
    ApplyLayout(L);
    PlotSim* ps = FindPlot(11, 2);
    ASSERT_TRUE(ps);
    const Machine& chest = ps->machines[ps->byPart.at(0)];
    int records = 0; for (int i = 0; i < kChestSlots; ++i) if (!chest.in[i].data.empty()) { ++records; EXPECT_EQ(chest.in[i].count, 1); }
    EXPECT_EQ(records, 2);
    EXPECT_EQ(CountItems(chest.in.data(), kChestSlots, tables.item("bronze_pickaxe")), 2) << "plain counts never include a record";
    const Machine& asm1 = ps->machines[ps->byPart.at(1)];
    EXPECT_EQ(asm1.mode, 1);
    bool asmKept = false; for (int i = 0; i < InSlots(asm1); ++i) asmKept |= asm1.in[i].data == d1;
    EXPECT_TRUE(asmKept);
    plots.clear(); plotIndexDirty = true; pendingStates.clear();
    { std::lock_guard<std::mutex> l(viewMutex); results.clear(); }
}

// (#996, ladder L4) The ore split, from machine_recipes.txt: the LV macerator alone gives half again as many ingots as the ore
// (crushed ore smelted 1:1); the MV chain (washer, centrifuge, dust smelted 1:1) gives twice as many, plus byproducts.
TEST(Factory, TheOreSplitGivesHalfAgainAtLvAndTwiceAtMv) {
    LoadRecipes();
    // Ratio out/in of `machine`'s recipe that takes `in` (the first such), counting its first output `out` (or its first).
    auto ratio = [](PartType machine, const std::string& in, const std::string& out) {
        const int inId = tables.item(in), outId = out.empty() ? -1 : tables.item(out);
        for (const auto& r : tables.recipes) {
            if (r.machine != machine || r.in.empty() || r.out.empty() || r.in[0].item != inId) continue;
            for (const auto& o : r.out) if (outId < 0 || o.item == outId) return (double)o.count / r.in[0].count;
        }
        return 0.0;
    };
    for (const char* metal : { "iron", "copper", "tin" }) {
        const std::string m = metal;
        const double crushed = ratio(PartType::Macerator, m + "_ore", m + "_crushed");
        const double basic = crushed * ratio(PartType::Furnace, m + "_crushed", m + "_ingot");
        const double full = crushed * ratio(PartType::OreWasher, m + "_crushed", m + "_purified") * ratio(PartType::Centrifuge, m + "_purified", m + "_dust")
                          * ratio(PartType::Furnace, m + "_dust", m + "_ingot");
        EXPECT_DOUBLE_EQ(basic, 1.5) << m << ": the LV macerator alone";
        EXPECT_DOUBLE_EQ(full, 2.0) << m << ": the MV chain";
    }
    EXPECT_GT(ratio(PartType::Centrifuge, "copper_purified", "gold_dust"), 0.0) << "gold still comes from the MV centrifuge";
}

// (#1271 L6) The steam lab: a rotation machine. Steered to tier 3 (the game's research labs set its variant), it turns steam
// packs into research units while a hand crank on its shaft turns, at the crank's speed; with nobody cranking it waits.
TEST(Factory, ASteamLabResearchesWhileItsShaftTurns) {
    LoadRecipes();
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    LayoutMsg L{};
    L.body = 7; L.plot = 91; L.version = 1;
    int part = 0;
    const int up = RotationPlot_OrientFor({ 0, 1, 0 });
    auto add = [&](PartType t, Int3 at, int orient = 0) {
        const Part pt = RotationPlot_PartAt(t, at, orient);
        std::vector<std::pair<Int3, CellShape>> fp;
        Plot::Footprint(pt, fp);
        if (IsMachine(t)) {
            std::vector<Int3> cells;
            for (const auto& c : fp) cells.push_back(c.first);
            L.machines.push_back({ part, t }); L.machineCells.push_back(cells); L.machineVariant.push_back(t == PartType::SteamLab ? "research_t3" : "");
        }
        RotationPartMsg rm;
        if (RotationPlot_Describe(pt, part, rm)) L.rotation.push_back(rm);
        return part++;
    };
    const int lab = add(PartType::SteamLab, { 0, 0, 0 });
    add(PartType::Shaft, { 0, 2, 0 }, up);
    const int crank = add(PartType::HandCrank, { 0, 3, 0 }, up);
    ApplyLayout(L);
    PlotSim& ps = *FindPlot(7, 91);
    Machine& m = ps.machines[ps.byPart.at(lab)];
    int buckets = 0;
    InsertInto(m, tables.item("steam_pack"), 4, false, buckets);
    long long made = 0; int active = 0;
    for (int t = 0; t < 300; ++t) TickPlot(ps, 1, made, active);
    EXPECT_EQ(CountItems(m.out, kSlots, tables.item("research_t3")), 0) << "nobody turns the crank: no research";
    for (int t = 0; t < 450; ++t) { ps.rot.Hold(crank); TickPlot(ps, 1, made, active); }
    EXPECT_EQ(CountItems(m.out, kSlots, tables.item("research_t3")), 2) << "450 ticks at the crank's 32 RPM (its rated speed): two 200-tick units";
    EXPECT_EQ(CountItems(m.in.data(), InSlots(m), tables.item("steam_pack")), 2);
    plots.clear(); plotIndexDirty = true;
}
// (#211, #1271 L2) Metallurgy: a bloomery makes iron blooms from ore and charcoal (400 ticks a bloom); bellows against it, their
// shaft cranked, make it up to twice as fast; an anvil hammers a bloom into an ingot only as far as a player's strikes go.
TEST(Factory, ABloomeryWithBellowsAndAnAnvilMakeIron) {
    LoadRecipes();
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    auto run = [](bool cranked) {
        LayoutMsg L{};
        L.body = 7; L.plot = cranked ? 93 : 92; L.version = 1;
        int part = 0;
        const int up = RotationPlot_OrientFor({ 0, 1, 0 });
        auto add = [&](PartType t, Int3 at, int orient = 0) {
            const Part pt = RotationPlot_PartAt(t, at, orient);
            std::vector<std::pair<Int3, CellShape>> fp;
            Plot::Footprint(pt, fp);
            if (IsMachine(t)) {
                std::vector<Int3> cells;
                for (const auto& c : fp) cells.push_back(c.first);
                L.machines.push_back({ part, t }); L.machineCells.push_back(cells); L.machineVariant.push_back("");
            }
            RotationPartMsg rm;
            if (RotationPlot_Describe(pt, part, rm)) L.rotation.push_back(rm);
            return part++;
        };
        const int bloomery = add(PartType::Bloomery, { 0, 0, 0 });
        add(PartType::Bellows, { 2, 0, 0 }, up);               // against the bloomery's side, its shaft up
        add(PartType::Shaft, { 2, 1, 0 }, up);
        const int crank = add(PartType::HandCrank, { 2, 2, 0 }, up);
        ApplyLayout(L);
        PlotSim& ps = *FindPlot(7, L.plot);
        Machine& m = ps.machines[ps.byPart.at(bloomery)];
        EXPECT_GE(m.bellows, 0) << "the bloomery finds the bellows touching it";
        int b = 0;
        InsertInto(m, tables.item("iron_ore"), 8, false, b);
        InsertInto(m, tables.item("charcoal"), 8, false, b);
        long long made = 0; int active = 0;
        for (int t = 0; t < 1210; ++t) { if (cranked) ps.rot.Hold(crank); TickPlot(ps, 1, made, active); }
        return CountItems(m.out, kSlots, tables.item("iron_bloom"));
    };
    EXPECT_EQ(run(false), 3) << "1210 ticks, 400 a bloom, no air";
    EXPECT_GE(run(true), 5) << "with the bellows blowing (the crank's 32 RPM, over their 16 rated): nearly twice as fast";

    // The anvil: a bloom needs 120 ticks of work, 40 a strike: two strikes leave it unfinished, the third finishes it.
    LayoutMsg L{};
    L.body = 7; L.plot = 94; L.version = 1;
    L.machines.push_back({ 0, PartType::Anvil }); L.machineCells.push_back({ { 0, 0, 0 } }); L.machineVariant.push_back("");
    ApplyLayout(L);
    PlotSim& ps = *FindPlot(7, 94);
    Machine& anvil = ps.machines[0];
    int b = 0;
    InsertInto(anvil, tables.item("iron_bloom"), 1, false, b);
    long long made = 0; int active = 0;
    for (int t = 0; t < 200; ++t) TickPlot(ps, 1, made, active);
    EXPECT_EQ(CountItems(anvil.out, kSlots, tables.item("iron_ingot")), 0) << "nobody hammers: nothing";
    for (int s = 0; s < 2; ++s) { anvil.strike = std::min(anvil.strike + kStrikeTicks, 2 * kStrikeTicks); for (int t = 0; t < 60; ++t) TickPlot(ps, 1, made, active); }
    EXPECT_EQ(CountItems(anvil.out, kSlots, tables.item("iron_ingot")), 0) << "two strikes: not yet";
    anvil.strike = std::min(anvil.strike + kStrikeTicks, 2 * kStrikeTicks);
    for (int t = 0; t < 60; ++t) TickPlot(ps, 1, made, active);
    EXPECT_EQ(CountItems(anvil.out, kSlots, tables.item("iron_ingot")), 1) << "the third strike: an ingot";
    plots.clear(); plotIndexDirty = true;
}
// (#1219, user 2026-11-22) The anvil's jobs on a composed copy, hammered: with a unit of its head's material, a repair (wear
// down by half its durability, the unit used); without one, its last gem comes out (or breaks); nothing to do: it waits.
TEST(Factory, AnAnvilRepairsACopyAndTakesAGemOut) {
    LoadRecipes();
    ASSERT_TRUE(Composed_Load());
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    LayoutMsg L{};
    L.body = 7; L.plot = 95; L.version = 1;
    L.machines.push_back({ 0, PartType::Anvil }); L.machineCells.push_back({ { 0, 0, 0 } }); L.machineVariant.push_back("");
    ApplyLayout(L);
    PlotSim& ps = *FindPlot(7, 95);
    Machine& anvil = ps.machines[0];
    long long made = 0; int active = 0;
    auto hammer = [&](int strikes) { for (int s = 0; s < strikes; ++s) { anvil.strike = std::min(anvil.strike + kStrikeTicks, 2 * kStrikeTicks); anvil.dirty = true; /* (as a strike does: it wakes the anvil) */ for (int t = 0; t < 60; ++t) TickPlot(ps, 1, made, active); } };
    auto outCopy = [&]() -> Stack* { for (int q = 0; q < kSlots; ++q) if (anvil.out[q].count > 0 && !anvil.out[q].data.empty()) return &anvil.out[q]; return nullptr; };
    EXPECT_TRUE(Accepts(PartType::Anvil, tables.item("bronze_ingot"))) << "a repair unit goes in";
    EXPECT_TRUE(Accepts(PartType::Anvil, tables.item("stone")));
    EXPECT_FALSE(Accepts(PartType::Anvil, tables.item("dirt")));
    // A worn bronze sword with a ruby, and a bronze ingot: repaired first (the ingot is there), the ruby kept.
    const std::string worn = ItemData::Set(ItemData::SetInt(Composed_MakeData("bronze"), "w", 300), "j", "ruby");
    const ComposedStats w0 = Composed_Stats("weapon_sword", worn);
    anvil.in[0] = Stack{}; anvil.in[0].item = tables.item("weapon_sword"); anvil.in[0].count = 1; anvil.in[0].data = worn;
    int b = 0;
    InsertInto(anvil, tables.item("bronze_ingot"), 1, false, b);
    for (int t = 0; t < 100; ++t) TickPlot(ps, 1, made, active);
    EXPECT_EQ(outCopy(), nullptr) << "nobody hammers: nothing";
    hammer(1);
    EXPECT_EQ(outCopy(), nullptr) << "one strike: not yet";
    hammer(1);
    Stack* r = outCopy();
    ASSERT_NE(r, nullptr) << "two strikes: repaired";
    EXPECT_EQ(ItemData::GetInt(r->data, "w", 0), std::max(0, w0.wear - w0.durability / 2)) << "wear down by half its durability";
    EXPECT_EQ(ItemData::Get(r->data, "j"), "ruby") << "a repair keeps the gem";
    EXPECT_EQ(CountItems(anvil.in.data(), InSlots(anvil), tables.item("bronze_ingot")), 0) << "the ingot used";
    // Back on the anvil with no ingot: its ruby comes out (or breaks), the sword without it.
    const std::string repaired = r->data;
    *r = Stack{};
    anvil.in[0] = Stack{}; anvil.in[0].item = tables.item("weapon_sword"); anvil.in[0].count = 1; anvil.in[0].data = repaired;
    hammer(2);
    r = outCopy();
    ASSERT_NE(r, nullptr);
    EXPECT_EQ(ItemData::Get(r->data, "j"), "") << "the gem is out";
    const int rubies = CountItems(anvil.out, kSlots, tables.item("ruby"));
    EXPECT_LE(rubies, 1) << "the ruby back, or broken";
    // Nothing to do (no wear to mend with no unit, no gem): it waits, no recipe runs.
    *r = Stack{};
    anvil.in[0] = Stack{}; anvil.in[0].item = tables.item("weapon_sword"); anvil.in[0].count = 1; anvil.in[0].data = Composed_MakeData("bronze");
    hammer(2);
    EXPECT_EQ(outCopy(), nullptr) << "nothing to do: it stays on the anvil";
    EXPECT_EQ(anvil.in[0].data, Composed_MakeData("bronze"));
    // Over many removals the break rate is about 25% (GUESS).
    int broke = 0;
    for (int k = 0; k < 200; ++k) {
        for (int q = 0; q < kSlots; ++q) anvil.out[q] = Stack{};
        anvil.in[0] = Stack{}; anvil.in[0].item = tables.item("weapon_sword"); anvil.in[0].count = 1; anvil.in[0].data = ItemData::Set(Composed_MakeData("bronze"), "j", "ruby");
        tickCount += 7;
        hammer(2);
        broke += CountItems(anvil.out, kSlots, tables.item("ruby")) == 0 ? 1 : 0;
    }
    EXPECT_GT(broke, 25);
    EXPECT_LT(broke, 80) << "about a quarter of 200";
    tickCount = 0;
    plots.clear(); plotIndexDirty = true;
}
// (#1271 L6c) The LV assembler: on a generator's power it assembles basic circuits from their parts; switched to
// disassembly (#1427 to come) it waits, and switched back it goes on.
// (#1446 step 2c) A machine making a composed class: each copy lands in an output slot of its own with its record built from
// the parts used and the machine's quality (T2's anvil sets it), never merged; plain counts don't see it; the save keeps it.
// (#1498) A conveyor never moves an item with data out of a chest (it used to copy it).
TEST(Factory, AMachineMakesComposedCopiesAndConveyorsLeaveRecordsAlone) {
    LoadRecipes();
    ASSERT_TRUE(Composed_Load());
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    const size_t nRecipes = tables.recipes.size();
    decltype(tables.recipes)::value_type r; r.machine = PartType::LvAssembler; r.ticks = 20; r.eu = 1;
    r.in = { Stack{ tables.item("bronze_pick_head"), 1 }, Stack{ tables.item("wood_handle"), 1 } };
    r.out = { Stack{ tables.item("tool_pickaxe"), 1 } };
    r.tier = MachineTiers_Tiered(r.machine) ? MachineTiers_First(r.machine) : MachineTiers_ForVolts(Def(r.machine).volts); // (as LoadRecipes sets it)
    tables.recipes.push_back(r);
    recipesOf[(int)r.machine].push_back((int)nRecipes); // (the caches LoadRecipes builds)
    std::vector<uint8_t>& wants = wantsOf[(int)r.machine];
    std::vector<int> wanted;
    for (const Stack& s : r.in) { if (s.item >= (int)wants.size()) wants.resize(s.item + 1, 0); if (!wants[s.item]) wanted.push_back(s.item); wants[s.item] = 1; }
    const LayoutMsg L = OracleLayout(7, { { PartType::Generator, 0, 0 }, { PartType::LvAssembler, 3, 0 }, { PartType::Chest, 6, 0 } }, Row(3, 7, 2), Row(0, 4, -1));
    ApplyLayout(L);
    PlotSim* ps = FindPlot(9, 7);
    ASSERT_TRUE(ps);
    Machine& gen = ps->machines[ps->byPart.at(0)];
    Machine& asm1 = ps->machines[ps->byPart.at(1)];
    Machine& chest = ps->machines[ps->byPart.at(2)];
    int b = 0;
    InsertInto(gen, tables.item("coal"), 16, true, b);
    // (#1498) a pick head with a record in the chest: the conveyor would offer it to the assembler (which takes pick heads)
    chest.in[0] = Stack{}; chest.in[0].item = tables.item("bronze_pick_head"); chest.in[0].count = 1; chest.in[0].data = "h=bronze;w=3";
    InsertInto(asm1, tables.item("bronze_pick_head"), 2, false, b);
    InsertInto(asm1, tables.item("wood_handle"), 2, false, b);
    asm1.pick = tables.item("tool_pickaxe"); asm1.quality = 106;
    long long made = 0; int active = 0;
    for (int t = 0; t < 200; ++t) TickPlot(*ps, 1, made, active);
    std::vector<std::string> recs;
    for (int i = 0; i < kSlots; ++i) if (asm1.out[i].count > 0 && asm1.out[i].item == tables.item("tool_pickaxe")) { EXPECT_EQ(asm1.out[i].count, 1) << "one copy a slot"; recs.push_back(asm1.out[i].data); }
    ASSERT_EQ(recs.size(), 2u) << "two crafts, two copies";
    EXPECT_EQ(ItemData::Get(recs[0], "h"), "bronze");
    EXPECT_EQ(ItemData::Get(recs[0], "g"), "wood");
    EXPECT_EQ(ItemData::GetInt(recs[0], "q", 100), 106) << "the first craft's quality";
    EXPECT_EQ(ItemData::GetInt(recs[1], "q", 100), 100) << "quality is used once";
    EXPECT_EQ(CountItems(asm1.out, kSlots, tables.item("tool_pickaxe")), 0) << "plain counts never see a copy";
    EXPECT_EQ(chest.in[0].data, "h=bronze;w=3") << "the chest keeps its pick head with a record";
    EXPECT_EQ(chest.in[0].count, 1);
    EXPECT_EQ(CountItems(asm1.in.data(), InSlots(asm1), tables.item("bronze_pick_head")), 0) << "no plain copy reached the assembler (#1498)";
    // The view and the save
    MachineView v; FillView(*ps, 1, v);
    int withData = 0; for (const std::string& d : v.outData) withData += d.empty() ? 0 : 1;
    EXPECT_EQ(withData, 2);
    const std::string path = (std::filesystem::temp_directory_path() / "raylibgame_factory_outdata_test.txt").string();
    FILE* f = fopen(path.c_str(), "w+");
    ASSERT_TRUE(f);
    Factory::Save(f, 9);
    rewind(f);
    const std::vector<std::string> lines = Factory::ReadSaved(f);
    fclose(f);
    std::remove(path.c_str());
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    Factory::LoadSaved(9, lines);
    ApplyLayout(L);
    ps = FindPlot(9, 7);
    ASSERT_TRUE(ps);
    const Machine& back = ps->machines[ps->byPart.at(1)];
    std::vector<std::string> loaded; for (int i = 0; i < kSlots; ++i) if (!back.out[i].data.empty()) loaded.push_back(back.out[i].data);
    EXPECT_EQ(loaded, recs) << "the save keeps each copy's record";
    tables.recipes.resize(nRecipes);
    recipesOf[(int)r.machine].pop_back();
    for (int item : wanted) wants[item] = 0;
    plots.clear(); plotIndexDirty = true; pendingStates.clear();
}

// (#1612) Every ore kind the world places (ores.txt: today's six and the vein families' new ones) is good for something: a
// machine recipe (or a hand recipe) takes it.
TEST(Factory, EveryOreKindHasARecipeThatTakesIt) {
    LoadRecipes();
    ASSERT_GE(Ores_Count(), 6);
    for (int k = 0; k < Ores_Count(); ++k) {
        const std::string item = Ores_Get(k)->item;
        bool used = false;
        for (const auto& r : tables.recipes) for (const Stack& s : r.in) used = used || (s.item >= 0 && tables.name(s.item) == item);
        for (const ::Recipe& r : itemDB.recipes()) for (const auto& in : r.inputs) used = used || ItemMatches(in.first, item);
        EXPECT_TRUE(used) << item << ": no recipe takes it";
    }
}

// (#1473) The material column's alloys (#1614): each has a machine recipe, made in a machine of its tier or earlier (the
// kiln is tier 1, the blast furnace tier 3, the arc furnace HV), and the ore chains give every metal they need.
// (#1614 PR 1, structure: fails CI, unlike the cost report) Every component shape of every metal (a materials.txt form that
// forms.txt makes) can be made by hand, and where a machine makes it too, the machine takes less metal than the hand.
TEST(Factory, EveryComponentShapeHasAHandRecipeAndMachinesSaveMetal) {
    UnitTest::EnsureItemDB();
    LoadRecipes();
    int shapes = 0;
    for (const ItemDef& d : itemDB.items()) {
        const size_t u = d.id.rfind('_');
        if (u == std::string::npos) continue;
        const std::string shape = d.id.substr(u + 1);
        if (std::find(Forms_Shapes().begin(), Forms_Shapes().end(), shape) == Forms_Shapes().end()) continue;
        ++shapes;
        // Its hand recipe, and the ingots (or rods, plates...) one item takes by hand
        double handPer = -1.0;
        for (const ::Recipe& r : itemDB.recipes())
            if (r.output == d.id) { int n = 0; for (const auto& in : r.inputs) n += in.second; handPer = (double)n / std::max(1, r.outputCount); break; }
        EXPECT_GE(handPer, 0.0) << d.id << ": no hand recipe (every step can be done by hand)";
        if (handPer < 0.0) continue;
        for (const auto& r : tables.recipes)
            for (const Stack& o : r.out) {
                if (o.item != tables.item(d.id) || o.count < 1) continue;
                int n = 0; for (const Stack& s : r.in) n += s.count;
                EXPECT_LT((double)n / o.count, handPer) << d.id << " in machine " << (int)r.machine << ": a machine takes less metal than the hand";
            }
    }
    EXPECT_GT(shapes, 20) << "the metals' plates, rods, wire, gears, screws, springs";
}

TEST(Factory, EachAlloyOfTheMaterialColumnIsMadeByItsTier) {
    LoadRecipes();
    const struct { const char* alloy; PartType latest; } kAlloys[] = {
        { "brass_ingot", PartType::Kiln }, { "solder_ingot", PartType::Kiln }, { "electrum_ingot", PartType::Kiln },
        { "invar_ingot", PartType::BlastFurnace }, { "stainless_steel_ingot", PartType::ArcFurnace }, { "tungstensteel_ingot", PartType::ArcFurnace } };
    for (const auto& a : kAlloys) {
        const int id = tables.item(a.alloy);
        ASSERT_GE(id, 0) << a.alloy << " is an item";
        bool made = false;
        for (const auto& r : tables.recipes)
            for (const Stack& o : r.out) if (o.item == id && r.machine == a.latest) made = true;
        EXPECT_TRUE(made) << a.alloy << ": made in its tier's machine";
    }
    // Their metals come from the ore chains: lead, zinc, silver, nickel, chromium and tungsten each have an ingot recipe
    for (const char* metal : { "lead_ingot", "zinc_ingot", "silver_ingot", "nickel_ingot", "chromium_ingot", "tungsten_ingot", "gold_ingot" }) {
        const int id = tables.item(metal);
        bool made = false;
        for (const auto& r : tables.recipes) for (const Stack& o : r.out) made = made || o.item == id;
        EXPECT_TRUE(made) << metal << ": an ore chain makes it";
    }
}

TEST(Factory, AnLvAssemblerMakesCircuitsAndWaitsWhenSetToDisassemble) {
    LoadRecipes();
    Disassembly_Build({}, nullptr); // (#1427 L7c) nothing to take apart here (other tests leave their own index)
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    ApplyLayout(OracleLayout(5, { { PartType::Generator, 0, 0 }, { PartType::LvAssembler, 3, 0 } }, {}, Row(0, 4, -1)));
    PlotSim& ps = *FindPlot(9, 5);
    Machine& gen = ps.machines[ps.byPart.at(0)];
    Machine& asm1 = ps.machines[ps.byPart.at(1)];
    int b = 0;
    InsertInto(gen, tables.item("coal"), 16, true, b);
    InsertInto(asm1, tables.item("vacuum_tube"), 8, false, b);
    InsertInto(asm1, tables.item("copper_wire"), 16, false, b);
    InsertInto(asm1, tables.item("steel_plate"), 4, false, b);
    int ticks = 160;
    for (const auto& r : tables.recipes) if (r.machine == PartType::LvAssembler && !r.out.empty() && r.out[0].item == tables.item("basic_circuit")) { ticks = r.ticks; break; }
    long long made = 0; int active = 0;
    asm1.mode = 1;
    for (int t = 0; t < 3 * ticks; ++t) TickPlot(ps, 1, made, active);
    EXPECT_EQ(CountItems(asm1.out, kSlots, tables.item("basic_circuit")), 0) << "set to disassemble: it waits (nothing it holds is taken apart)";
    EXPECT_EQ(asm1.status, MachineStatus::Idle);
    asm1.mode = 0; asm1.dirty = true; // (as SetMode does: it wakes)
    for (int t = 0; t < 2 * ticks + 40; ++t) TickPlot(ps, 1, made, active);
    EXPECT_EQ(CountItems(asm1.out, kSlots, tables.item("basic_circuit")), 2) << "two circuits in two crafts' time";
    EXPECT_EQ(CountItems(asm1.in.data(), InSlots(asm1), tables.item("vacuum_tube")), 4);
    // Its recipe picked (the screen's "Makes"): LV packs from the same parts, never a circuit
    asm1.pick = tables.item("lv_pack"); asm1.recipe = -1; asm1.progress = 0; asm1.dirty = true; // (as SetPick does: the circuit it began is dropped)
    int packTicks = 200;
    for (const auto& r : tables.recipes) if (r.machine == PartType::LvAssembler && !r.out.empty() && r.out[0].item == asm1.pick) { packTicks = r.ticks; break; }
    for (int t = 0; t < packTicks + 40; ++t) TickPlot(ps, 1, made, active);
    EXPECT_EQ(CountItems(asm1.out, kSlots, tables.item("basic_circuit")), 2) << "picked LV packs: no more circuits";
    EXPECT_EQ(CountItems(asm1.out, kSlots, tables.item("lv_pack")), 2) << "one craft of 2 LV packs";
    plots.clear(); plotIndexDirty = true;
}
// (#1427 L7c) Set to disassemble, an LV assembler takes LV circuits apart at 50%: the raw steel plate as hundredths carried
// (3 circuits: 1 plate and half of one waiting), the made parts (wires, tubes) rolled per unit; each in the forward recipe's
// time; the job and the carry survive a save; an MV assembler takes them apart at 75%, overclocked.
TEST(Factory, AnAssemblerSetToDisassembleGivesBackPartsAtItsTiersYield) {
    LoadRecipes();
    // The game's ladder tiers (the hook): circuits unlock at tier 4 (LV); their parts are LV jobs too, by their ingredients
    std::vector<DisassemblyRecipe> forward;
    for (const auto& r : tables.recipes) {
        if (r.machine != PartType::LvAssembler || r.out.size() != 1) continue;
        DisassemblyRecipe d{ tables.name(r.out[0].item), r.out[0].count, {}, r.ticks, r.eu };
        for (const Stack& s : r.in) d.in.push_back({ tables.name(s.item), s.count });
        forward.push_back(d);
    }
    forward.push_back({ "vacuum_tube", 1, { { "glass", 1 }, { "copper_wire", 1 } }, 0, 0 });
    forward.push_back({ "copper_wire", 2, { { "copper_ingot", 1 } }, 0, 0 });
    Disassembly_Build(forward, [](const std::string& item) { return item == "basic_circuit" ? 4 : -1; });
    for (const std::string& n : Disassembly_Items()) { const int id = tables.item(n); if (id >= (int)disassemblable.size()) disassemblable.resize(id + 1, 0); disassemblable[id] = 1; }
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    LayoutMsg L = OracleLayout(7, { { PartType::Generator, 0, 0 }, { PartType::LvAssembler, 3, 0 }, { PartType::LvAssembler, 6, 0 } }, {}, Row(0, 7, -1));
    L.machineTier = { MachineTiers_ByName("LV"), MachineTiers_ByName("LV"), MachineTiers_ByName("MV") };
    ApplyLayout(L);
    PlotSim* ps = FindPlot(9, 7);
    auto M = [&](int part) -> Machine& { return ps->machines[ps->byPart.at(part)]; };
    int b = 0;
    InsertInto(M(0), tables.item("coal"), 64, true, b);
    const int circuit = tables.item("basic_circuit"), steel = tables.item("steel_plate"), wire = tables.item("copper_wire"), tube = tables.item("vacuum_tube");
    M(1).mode = 1; M(2).mode = 1;
    EXPECT_EQ(InsertInto(M(1), tables.item("steel_plate"), 1, false, b), 1) << "set to disassemble it takes only what it can take apart";
    EXPECT_EQ(InsertInto(M(1), circuit, 3, false, b), 0);
    int ticks = 160;
    for (const auto& r : tables.recipes) if (r.machine == PartType::LvAssembler && !r.out.empty() && r.out[0].item == circuit) { ticks = r.ticks; break; }
    long long made = 0; int active = 0;
    int noPower = 0;
    for (int t = 0; t < ticks + ticks / 2; ++t) { TickPlot(*ps, 1, made, active); if (t > 0 && M(1).status == MachineStatus::NoPower) ++noPower; }
    EXPECT_EQ(noPower, 0) << "powered all along: no \"no power\" tick between one job and the next";
    EXPECT_EQ(CountItems(M(1).in.data(), kSlots, circuit), 2) << "one circuit in one forward recipe's time";
    EXPECT_EQ(M(1).disSlot, 0) << "the second one under way";
    EXPECT_EQ(M(1).disCarry[steel], 50) << "half a steel plate carried";
    // Save and load mid-job: the job, the carry and the rolls' state come back
    const uint32_t rngBefore = M(1).disRng;
    const int progressBefore = M(1).progress;
    const std::string path = (std::filesystem::temp_directory_path() / "raylibgame_disassembly_save_test.txt").string();
    FILE* f = fopen(path.c_str(), "w+");
    ASSERT_TRUE(f);
    Factory::Save(f, 9);
    rewind(f);
    const std::vector<std::string> lines = Factory::ReadSaved(f);
    fclose(f);
    std::remove(path.c_str());
    plots.clear(); plotIndexDirty = true; pendingStates.clear();
    Factory::LoadSaved(9, lines);
    ApplyLayout(L);
    ps = FindPlot(9, 7);
    EXPECT_EQ(M(1).mode, 1);
    EXPECT_EQ(M(1).disSlot, 0);
    EXPECT_EQ(M(1).progress, progressBefore);
    EXPECT_EQ(M(1).disRng, rngBefore);
    EXPECT_EQ(M(1).disCarry[steel], 50);
    for (int t = 0; t < 4 * ticks; ++t) TickPlot(*ps, 1, made, active);
    EXPECT_EQ(CountItems(M(1).in.data(), kSlots, circuit), 0) << "all three taken apart";
    EXPECT_EQ(CountItems(M(1).out, kSlots, steel), 1) << "3 x 1 plate x 50% = 1.5: one out";
    EXPECT_EQ(M(1).disCarry[steel], 50) << "and half of one carried";
    const int wires = CountItems(M(1).out, kSlots, wire), tubes = CountItems(M(1).out, kSlots, tube);
    EXPECT_GT(wires, 0); EXPECT_LT(wires, 12) << "12 wires went in: rolled at 50% each, never all";
    EXPECT_LE(tubes, 6);
    // MV: 75% (one tier below), overclocked one step: half the time at 3x the EU a tick (48 on one 32 EU/t generator: 2/3 speed)
    EXPECT_EQ(InsertInto(M(2), circuit, 4, false, b), 0);
    M(2).dirty = true;
    for (int t = 0; t < 4 * ticks; ++t) TickPlot(*ps, 1, made, active);
    EXPECT_EQ(CountItems(M(2).in.data(), kSlots, circuit), 0) << "MV: four taken apart in 4 x 160 x 1/2 x 3/2 = 480 ticks";
    EXPECT_EQ(CountItems(M(2).out, kSlots, steel), 3) << "4 x 75% = 3 plates";
    plots.clear(); plotIndexDirty = true;
}
// (#1427 L7a) Tiers: an MV macerator runs the LV iron recipe overclocked one step (half the time, three times the EU a
// tick); a recipe a tier above a machine never starts in it.
TEST(Factory, AHigherTierMachineOverclocksALowerRecipeAndALowerOneRefusesAHigherRecipe) {
    LoadRecipes();
    plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear();
    LayoutMsg L = OracleLayout(6, { { PartType::Generator, 0, 0 }, { PartType::Macerator, 3, 0 }, { PartType::Macerator, 6, 0 } }, {}, Row(0, 7, -1));
    L.machineTier = { MachineTiers_ByName("LV"), MachineTiers_ByName("MV"), MachineTiers_ByName("LV") };
    ApplyLayout(L);
    PlotSim& ps = *FindPlot(9, 6);
    Machine& gen = ps.machines[ps.byPart.at(0)];
    Machine& mv = ps.machines[ps.byPart.at(1)];
    Machine& lv = ps.machines[ps.byPart.at(2)];
    EXPECT_EQ(mv.tier, 1);
    EXPECT_EQ(lv.tier, 0);
    decltype(&tables.recipes[0]) iron = nullptr;
    for (auto& r : tables.recipes) if (r.machine == PartType::Macerator && !r.in.empty() && r.in[0].item == tables.item("iron_ore")) { iron = &r; break; }
    ASSERT_NE(iron, nullptr);
    EXPECT_EQ(iron->tier, 0) << "no tier column: the macerator's first, LV";
    EXPECT_EQ(TicksOf(mv, *iron), iron->ticks / 2);
    EXPECT_EQ(EuOf(mv, *iron), iron->eu * 3);
    EXPECT_EQ(TicksOf(lv, *iron), iron->ticks);
    int b = 0;
    InsertInto(gen, tables.item("coal"), 32, true, b);
    InsertInto(mv, tables.item("iron_ore"), 8, false, b);
    long long made = 0; int active = 0;
    for (int t = 0; t < iron->ticks + 20; ++t) TickPlot(ps, 1, made, active);
    EXPECT_EQ(CountItems(mv.out, kSlots, tables.item("iron_crushed")), 2 * iron->out[0].count) << "two crafts in one LV craft's time";
    // The same recipe made MV: the LV macerator never starts it; the MV one runs it at its own tier (no overclock).
    const int was = iron->tier;
    iron->tier = 1;
    InsertInto(lv, tables.item("iron_ore"), 4, false, b);
    for (int t = 0; t < iron->ticks + 20; ++t) TickPlot(ps, 1, made, active);
    EXPECT_EQ(CountItems(lv.out, kSlots, tables.item("iron_crushed")), 0) << "an MV recipe in an LV machine: never";
    EXPECT_EQ(TicksOf(mv, *iron), iron->ticks) << "at its own tier: no overclock";
    iron->tier = was;
    plots.clear(); plotIndexDirty = true;
}

// ---- (#1429 PR 4) Multiblocks in the factory ----
namespace {
// A plot of one-cell machines in a row along x (part = index, cell (2i, 0, 0)), each with its hatch kind, and the structures.
LayoutMsg MbLayout(uint32_t plotId, const std::vector<std::pair<PartType, uint8_t>>& parts, const std::vector<LayoutMsg::Structure>& structures, const std::vector<Int3>& cables = {}) {
    LayoutMsg L{};
    L.body = 9; L.plot = plotId; L.version = 1;
    for (size_t i = 0; i < parts.size(); ++i) {
        L.machines.push_back({ (int)i, parts[i].first }); L.machineCells.push_back({ { 2 * (int)i, 0, 0 } }); L.machineVariant.push_back("");
        L.machineTier.push_back(parts[i].second == kHatchEnergy ? 0 : MachineTiers_PartTier(parts[i].first, 0, Def(parts[i].first).volts)); // (an LV energy hatch)
        L.machineHatch.push_back(parts[i].second);
    }
    L.cables = cables;
    L.structures = structures;
    return L;
}
LayoutMsg::Structure MbOf(int controller, PartType machine, int parallel, std::vector<int> hatches, int tier = -1) {
    LayoutMsg::Structure s; s.controller = controller; s.machine = machine; s.parallel = parallel; s.hatches = std::move(hatches);
    s.tier = tier >= 0 ? tier : MachineTiers_PartTier(machine, 0, Def(machine).volts);
    return s;
}
const Recipe* RecipeFor(PartType machine, const char* input) {
    for (const Recipe& r : tables.recipes) if (r.machine == machine && !r.in.empty() && r.in[0].item == tables.item(input) && !r.out.empty()) return &r;
    return nullptr;
}
void MbReset() { plots.clear(); plotIndexDirty = true; pendingStates.clear(); refunds.clear(); }
} // namespace

// A coke oven battery (parallel 2) takes coal from two input hatches as one pool: inputs claimed at the start, two cokes at
// once into the output hatch. Not formed, it does nothing.
TEST(Factory, AFormedStructureRunsParallelRecipesFromItsPooledHatches) {
    LoadRecipes();
    MbReset();
    const auto* coke = RecipeFor(PartType::CokeOven, "coal");
    ASSERT_TRUE(coke);
    const std::vector<std::pair<PartType, uint8_t>> parts = { { PartType::MbController, kHatchNone }, { PartType::MbHatch, kHatchItemIn }, { PartType::MbHatch, kHatchItemIn }, { PartType::MbHatch, kHatchItemOut } };
    ApplyLayout(MbLayout(1, parts, {}));
    PlotSim& ps = *FindPlot(9, 1);
    auto M = [&](int part) -> Machine& { return ps.machines[ps.byPart.at(part)]; };
    int b = 0;
    const int need = coke->in[0].count; // per craft, split over the two hatches
    InsertInto(M(1), tables.item("coal"), need, false, b);
    InsertInto(M(2), tables.item("coal"), need, false, b);
    long long made = 0; int active = 0;
    for (int t = 0; t < 40; ++t) TickPlot(ps, 1, made, active);
    EXPECT_TRUE(M(0).jobs.empty()) << "unformed: nothing starts";
    EXPECT_EQ(CountItems(M(1).in.data(), kSlots, tables.item("coal")), need);
    ApplyLayout(MbLayout(1, parts, { MbOf(0, PartType::CokeOven, 2, { 1, 2, 3 }) }));
    TickPlot(ps, 1, made, active);
    EXPECT_EQ((int)M(0).jobs.size(), 2) << "two recipes at once from the pooled coal";
    EXPECT_EQ(CountItems(M(1).in.data(), kSlots, tables.item("coal")) + CountItems(M(2).in.data(), kSlots, tables.item("coal")), 0) << "inputs claimed at the start";
    for (int t = 0; t < coke->ticks + 5; ++t) TickPlot(ps, 1, made, active);
    EXPECT_EQ(CountItems(M(3).out, kSlots, coke->out[0].item), 2 * coke->out[0].count) << "both finished into the output hatch in one recipe's time";
    EXPECT_TRUE(M(0).jobs.empty());
    MbReset();
}

// Two structures share one input hatch: with coal for one recipe at a time, they start in turn.
TEST(Factory, StructuresSharingAnInputHatchStartInTurn) {
    LoadRecipes();
    MbReset();
    const auto* coke = RecipeFor(PartType::CokeOven, "coal");
    ASSERT_TRUE(coke);
    ApplyLayout(MbLayout(2, { { PartType::MbController, kHatchNone }, { PartType::MbController, kHatchNone }, { PartType::MbHatch, kHatchItemIn },
                              { PartType::MbHatch, kHatchItemOut }, { PartType::MbHatch, kHatchItemOut } },
                         { MbOf(0, PartType::CokeOven, 1, { 2, 3 }), MbOf(1, PartType::CokeOven, 1, { 2, 4 }) }));
    PlotSim& ps = *FindPlot(9, 2);
    auto M = [&](int part) -> Machine& { return ps.machines[ps.byPart.at(part)]; };
    long long made = 0; int active = 0;
    std::vector<int> starters;
    for (int round = 0; round < 4; ++round) {
        int b = 0;
        InsertInto(M(2), tables.item("coal"), coke->in[0].count, false, b);
        TickPlot(ps, 1, made, active);
        starters.push_back(!M(0).jobs.empty() ? 0 : !M(1).jobs.empty() ? 1 : -1);
        for (int t = 0; t < coke->ticks + 5; ++t) TickPlot(ps, 1, made, active);
    }
    EXPECT_EQ(starters, (std::vector<int>{ 0, 1, 0, 1 })) << "round-robin per recipe start";
    EXPECT_EQ(CountItems(M(3).out, kSlots, coke->out[0].item), CountItems(M(4).out, kSlots, coke->out[0].item));
    MbReset();
}

// A recipe reserves its outputs' room when it starts: with room for one craft's coke, only one of two parallel recipes starts,
// and nothing put in the output hatch later takes the room the running one holds.
TEST(Factory, AStructureReservesItsOutputsRoomAtTheStart) {
    LoadRecipes();
    MbReset();
    const auto* coke = RecipeFor(PartType::CokeOven, "coal");
    ASSERT_TRUE(coke);
    ApplyLayout(MbLayout(3, { { PartType::MbController, kHatchNone }, { PartType::MbHatch, kHatchItemIn }, { PartType::MbHatch, kHatchItemOut } },
                         { MbOf(0, PartType::CokeOven, 2, { 1, 2 }) }));
    PlotSim& ps = *FindPlot(9, 3);
    auto M = [&](int part) -> Machine& { return ps.machines[ps.byPart.at(part)]; };
    const int cokeId = coke->out[0].item, cap = tables.cap(cokeId);
    // Three slots full of stone, the fourth holding coke with room for exactly one craft.
    for (int s = 0; s < kSlots - 1; ++s) { M(2).out[s].item = tables.item("stone"); M(2).out[s].count = tables.cap(tables.item("stone")); }
    M(2).out[kSlots - 1].item = cokeId; M(2).out[kSlots - 1].count = cap - coke->out[0].count;
    int b = 0;
    InsertInto(M(1), tables.item("coal"), 2 * coke->in[0].count, false, b);
    long long made = 0; int active = 0;
    TickPlot(ps, 1, made, active);
    EXPECT_EQ((int)M(0).jobs.size(), 1) << "room for one craft: one starts";
    EXPECT_EQ(CountItems(M(1).in.data(), kSlots, tables.item("coal")), coke->in[0].count) << "the other's coal stays in the hatch";
    for (int t = 0; t < coke->ticks + 5; ++t) TickPlot(ps, 1, made, active);
    EXPECT_EQ(M(2).out[kSlots - 1].count, cap) << "the first finished into its reserved room";
    EXPECT_EQ(M(0).status, MachineStatus::OutputFull) << "the second waits for room";
    MbReset();
}

// Two macerator structures behind one LV energy hatch, wanting twice what its generator makes: each gets half (evenly).
TEST(Factory, AnEnergyHatchSharedByTwoStructuresSlowsBothEvenly) {
    LoadRecipes();
    MbReset();
    const auto* iron = RecipeFor(PartType::Macerator, "iron_ore");
    ASSERT_TRUE(iron);
    const int parallel = std::max(1, 32 / std::max(1, iron->eu)); // each alone could use the generator's whole 32 EU/t
    ApplyLayout(MbLayout(4, { { PartType::Generator, kHatchNone }, { PartType::MbHatch, kHatchEnergy }, { PartType::MbController, kHatchNone }, { PartType::MbController, kHatchNone },
                              { PartType::MbHatch, kHatchItemIn }, { PartType::MbHatch, kHatchItemIn }, { PartType::MbHatch, kHatchItemOut }, { PartType::MbHatch, kHatchItemOut } },
                         { MbOf(2, PartType::Macerator, parallel, { 1, 4, 6 }, 0), MbOf(3, PartType::Macerator, parallel, { 1, 5, 7 }, 0) },
                         Row(0, 2, -1)));
    PlotSim& ps = *FindPlot(9, 4);
    auto M = [&](int part) -> Machine& { return ps.machines[ps.byPart.at(part)]; };
    int b = 0;
    InsertInto(M(0), tables.item("coal"), 64, true, b);
    InsertInto(M(4), tables.item("iron_ore"), 256, false, b);
    InsertInto(M(5), tables.item("iron_ore"), 256, false, b);
    long long made = 0; int active = 0;
    const int ticks = 4 * iron->ticks;
    for (int t = 0; t < ticks; ++t) TickPlot(ps, 1, made, active);
    auto crafts = [](const Machine& m) { int n = 0; for (const auto& kv : m.crafts) n += kv.second; return n; };
    const int a = crafts(M(2)), c = crafts(M(3));
    EXPECT_GT(a, 0);
    EXPECT_NEAR(a, c, parallel) << "evenly: within one round of crafts";
    // (32 EU/t buys this many crafts in the window; the last round of each job may finish just after it ends)
    EXPECT_GE(a + c, (int)(32LL * ticks / ((long long)iron->eu * iron->ticks)) - 2 * parallel) << "crafts " << a << " + " << c;
    EXPECT_EQ(M(2).status, MachineStatus::Working);
    const Series& used = ps.stats.euUsed[(int)PartType::MbHatch];
    const long long eu = used.total + used.acc;
    EXPECT_LE(eu, 32LL * ticks) << "the hatch passed no more than the generator's 32 EU/t";
    EXPECT_GE(eu, 32LL * ticks * 9 / 10) << "and nearly all of it (the demand is twice that)";
    MbReset();
}

// A structure's running recipes save with it and finish after a load; a removed controller gives back what its recipes claimed.
TEST(Factory, AStructuresJobsSaveAndARemovedControllerGivesBackItsClaimedInputs) {
    LoadRecipes();
    MbReset();
    const auto* coke = RecipeFor(PartType::CokeOven, "coal");
    ASSERT_TRUE(coke);
    const std::vector<std::pair<PartType, uint8_t>> parts = { { PartType::MbController, kHatchNone }, { PartType::MbHatch, kHatchItemIn }, { PartType::MbHatch, kHatchItemOut } };
    const LayoutMsg L = MbLayout(5, parts, { MbOf(0, PartType::CokeOven, 1, { 1, 2 }) });
    ApplyLayout(L);
    PlotSim* ps = FindPlot(9, 5);
    int b = 0;
    InsertInto(ps->machines[ps->byPart.at(1)], tables.item("coal"), coke->in[0].count, false, b);
    long long made = 0; int active = 0;
    for (int t = 0; t < coke->ticks / 2; ++t) TickPlot(*ps, 1, made, active);
    ASSERT_EQ((int)ps->machines[ps->byPart.at(0)].jobs.size(), 1);
    const int progress = ps->machines[ps->byPart.at(0)].jobs[0].progress;
    const std::string path = (std::filesystem::temp_directory_path() / "raylibgame_factory_mb_test.txt").string();
    FILE* f = fopen(path.c_str(), "w+");
    ASSERT_TRUE(f);
    Factory::Save(f, 9);
    rewind(f);
    const std::vector<std::string> lines = Factory::ReadSaved(f);
    fclose(f);
    std::remove(path.c_str());
    MbReset();
    Factory::LoadSaved(9, lines);
    ApplyLayout(L);
    ps = FindPlot(9, 5);
    ASSERT_TRUE(ps);
    Machine& c = ps->machines[ps->byPart.at(0)];
    ASSERT_EQ((int)c.jobs.size(), 1) << "the job came back";
    EXPECT_EQ(c.jobs[0].progress, progress);
    for (int t = 0; t < coke->ticks; ++t) TickPlot(*ps, 1, made, active);
    EXPECT_EQ(CountItems(ps->machines[ps->byPart.at(2)].out, kSlots, coke->out[0].item), coke->out[0].count) << "and finished (its room reserved again)";
    // Claim again, then break the controller: the claimed coal comes back.
    InsertInto(ps->machines[ps->byPart.at(1)], tables.item("coal"), coke->in[0].count, false, b);
    TickPlot(*ps, 1, made, active);
    ASSERT_EQ((int)ps->machines[ps->byPart.at(0)].jobs.size(), 1);
    refunds.clear();
    LayoutMsg gone = MbLayout(5, { { PartType::MbHatch, kHatchItemIn }, { PartType::MbHatch, kHatchItemOut } }, {});
    gone.machines = { { 1, PartType::MbHatch }, { 2, PartType::MbHatch } };
    ApplyLayout(gone);
    int back = 0;
    for (const ItemOut& r : refunds) if (r.item == tables.item("coal")) back += r.count;
    EXPECT_EQ(back, coke->in[0].count);
    MbReset();
}
#endif
