#include "inventory.h"
#include "composed.h"
#include "forms.h"
#include "utils/logger.h"
#include "world.h"
#include "voxel.h"
#include "plots.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <set>
#include <sstream>

ItemDB itemDB;
// (inventory: the local PlayerSession's, src/space/playersession.cpp, #313)

const char* RetiredItem(const std::string& id) { return id == "concrete" ? "stone" : ""; }

static std::string Trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

static std::vector<std::string> Split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string part;
    while (std::getline(ss, part, sep)) out.push_back(Trim(part));
    return out;
}

bool ItemDB::loadExtra(const char* itemsPath) {
    std::ifstream file(itemsPath);
    if (!file) { LOG_ERROR("Items", "Cannot open %s\n", itemsPath); return false; }
    std::stringstream text;
    text << file.rdbuf();
    return loadExtraText(text.str());
}

bool ItemDB::loadExtraText(const std::string& text) {
    std::istringstream items(text);
    std::string line;
    std::vector<std::string> recipeLines; // (#1252) "out n = in n, ..." lines: parsed once every item of the file is known
    while (std::getline(items, line)) {
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;
        if (line.find('|') == std::string::npos && line.find('=') != std::string::npos) { recipeLines.push_back(line); continue; }
        const auto f = Split(line, '|');
        if (f.size() < 5 || byId.count(f[0])) continue;
        ItemDef d;
        d.id = f[0];
        d.name = f[1];
        d.kind = f[2] == "tool" || f[2] == "1" ? ItemKind::Tool : f[2] == "block" ? ItemKind::Block : ItemKind::Resource;
        const auto rgb = Split(f[3], ',');
        if (rgb.size() == 3) d.color = { (unsigned char)atoi(rgb[0].c_str()), (unsigned char)atoi(rgb[1].c_str()), (unsigned char)atoi(rgb[2].c_str()), 255 };
        d.description = f[4];
        byId[d.id] = itemList.size();
        itemList.push_back(d);
    }
    for (const std::string& rl : recipeLines) {
        const size_t eq = rl.find('=');
        Recipe r;
        std::istringstream out(rl.substr(0, eq));
        out >> r.output >> r.outputCount;
        for (const std::string& in : Split(rl.substr(eq + 1), ',')) {
            std::istringstream is(in);
            std::string id; int n = 0;
            if (is >> id >> n && n > 0) r.inputs.push_back({ id, n });
        }
        bool valid = find(r.output) != nullptr && !r.inputs.empty();
        for (auto& in : r.inputs) valid = valid && find(in.first) != nullptr;
        if (valid) recipeList.push_back(r);
        else LOG_WARN("Items", "recipes: unknown item in '%s'\n", rl.c_str());
    }
    // Placeable parts named here are blocks, as load() makes them.
    for (int t = 0; t < (int)PartType::Count; ++t)
        if (const char* item = PartCostItem((PartType)t); item && !IsShape((PartType)t) && (PartType)t != PartType::Lever)
            if (auto it = byId.find(item); it != byId.end() && itemList[it->second].kind == ItemKind::Resource) itemList[it->second].kind = ItemKind::Block;
    return true;
}

bool ItemDB::load(const char* itemsPath, const char* recipesPath) {
    itemList.clear();
    byId.clear();
    recipeList.clear();

    std::ifstream items(itemsPath);
    if (!items) { LOG_ERROR("Items", "Cannot open %s\n", itemsPath); return false; }
    std::string line;
    while (std::getline(items, line)) {
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;
        const auto f = Split(line, '|');
        if (f.size() < 5) { LOG_WARN("Items", "items: bad line '%s'\n", line.c_str()); continue; }
        ItemDef d;
        d.id = f[0];
        d.name = f[1];
        // Kind: "tool", "block" or "resource"; a legacy stack-limit number of 1 means a tool.
        d.kind = f[2] == "tool" || f[2] == "1" ? ItemKind::Tool : f[2] == "block" ? ItemKind::Block : ItemKind::Resource;
        const auto rgb = Split(f[3], ',');
        if (rgb.size() == 3) d.color = { (unsigned char)atoi(rgb[0].c_str()), (unsigned char)atoi(rgb[1].c_str()), (unsigned char)atoi(rgb[2].c_str()), 255 };
        d.description = f[4];
        byId[d.id] = itemList.size();
        itemList.push_back(d);
    }

    // Material families (materials.txt next to items.txt): one line adds a material's ore, dust,
    // ingot, plate... items; ids are material_form, and items.txt entries win.
    {
        std::string dir = itemsPath;
        const size_t slash = dir.find_last_of("/\\");
        dir = slash == std::string::npos ? std::string() : dir.substr(0, slash + 1);
        std::ifstream mats(dir + "materials.txt");
        while (std::getline(mats, line)) {
            line = Trim(line);
            if (line.empty() || line[0] == '#') continue;
            const auto f = Split(line, '|');
            if (f.size() < 4) continue;
            const auto rgb = Split(f[2], ',');
            Color c = GRAY;
            if (rgb.size() == 3) c = { (unsigned char)atoi(rgb[0].c_str()), (unsigned char)atoi(rgb[1].c_str()), (unsigned char)atoi(rgb[2].c_str()), 255 };
            std::istringstream forms(f[3]);
            std::string form;
            while (forms >> form) {
                ItemDef d;
                d.id = f[0] + "_" + form;
                if (byId.count(d.id)) continue;
                // "Iron ore", "Crushed iron ore", "Purified iron ore", "Iron dust", "Iron ingot"...
                std::string lower = f[1];
                if (!lower.empty()) lower[0] = (char)tolower(lower[0]);
                d.name = form == "crushed" ? "Crushed " + lower + " ore" : form == "purified" ? "Purified " + lower + " ore" : f[1] + " " + form;
                d.color = form == "ore" ? Color{ (unsigned char)(c.r * 0.7f), (unsigned char)(c.g * 0.6f), (unsigned char)(c.b * 0.5f), 255 } : c;
                d.description = form == "ore" ? "Ore from deep veins." : form == "dust" ? "Crushed; smelts quickly." : form == "crushed" ? "Ore crushed in a macerator: wash it (ore washer) for more." : form == "purified" ? "Washed ore: a centrifuge splits off its byproduct." : "";
                byId[d.id] = itemList.size();
                itemList.push_back(d);
            }
        }
    }

    std::ifstream recipes(recipesPath);
    if (!recipes) { LOG_ERROR("Items", "Cannot open %s\n", recipesPath); return false; }
    while (std::getline(recipes, line)) {
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;
        const size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        Recipe r;
        std::istringstream out(line.substr(0, eq));
        out >> r.output >> r.outputCount;
        for (const std::string& in : Split(line.substr(eq + 1), ',')) {
            std::istringstream is(in);
            std::string id;
            int n = 0;
            if (is >> id >> n && n > 0) r.inputs.push_back({ id, n });
        }
        bool valid = find(r.output) != nullptr && !r.inputs.empty();
        for (auto& in : r.inputs) valid = valid && find(in.first) != nullptr;
        if (valid) recipeList.push_back(r);
        else LOG_WARN("Items", "recipes: unknown item in '%s'\n", line.c_str());
    }
    // (#1614) The form table's hand recipes (forms.txt): every metal's plates, rods, wire, gears, screws and springs by hand,
    // where recipes.txt has no recipe for that output
    {
        std::set<std::string> made;
        for (const Recipe& r : recipeList) made.insert(r.output);
        for (const FormRecipe& f : Forms_All()) {
            if (!f.machine.empty() || made.count(f.output) || !find(f.output)) continue;
            Recipe r;
            r.output = f.output; r.outputCount = f.outCount; r.inputs = f.inputs; r.formShape = true;
            bool valid = true;
            for (auto& in : r.inputs) valid = valid && find(in.first) != nullptr;
            if (valid) { recipeList.push_back(r); made.insert(r.output); }
        }
    }
    // Placeable parts (machines, conveyors, cables, rocket parts) are blocks: stacks of 256. Not the parts only the world
    // places (#812): a meteorite and dungeon traps name the material they break into (nickel-iron ore, stone), which must
    // stay a resource (stone as a block stacked to 256 and lost nothing on death).
    for (int t = 0; t < (int)PartType::Count; ++t)
        if (IsShape((PartType)t) || (PartType)t == PartType::Meteorite || (PartType)t == PartType::PressurePlate || (PartType)t == PartType::GasVent || (PartType)t == PartType::Lever) continue; else
        if (const char* item = PartCostItem((PartType)t)) {
            auto it = byId.find(item);
            if (it != byId.end() && itemList[it->second].kind == ItemKind::Resource) itemList[it->second].kind = ItemKind::Block;
        }
    LOG_INFO("Items", "%zu items, %zu recipes\n", itemList.size(), recipeList.size());
    return true;
}

const ItemDef* ItemDB::find(const std::string& id) const {
    auto it = byId.find(id);
    return it == byId.end() ? nullptr : &itemList[it->second];
}

const char* OreKindAt(Vector3 local, unsigned int seed) {
    // Deposits about 48 m across hold one kind of ore: iron 30%, copper 25%, tin 17%, coal 13%,
    // bauxite (aluminium) 10%, ilmenite (titanium) 5%.
    const int x = (int)floorf(local.x / 48.0f), y = (int)floorf(local.y / 48.0f), z = (int)floorf(local.z / 48.0f);
    uint32_t h = (uint32_t)x * 73856093u ^ (uint32_t)y * 19349663u ^ (uint32_t)z * 83492791u ^ seed * 2654435761u;
    h ^= h >> 16; h *= 0x7feb352du; h ^= h >> 15; h *= 0x846ca68bu; h ^= h >> 16;
    const unsigned r = h % 100;
    return r < 30 ? "iron_ore" : r < 55 ? "copper_ore" : r < 72 ? "tin_ore" : r < 85 ? "coal" : r < 95 ? "bauxite_ore" : "ilmenite_ore";
}

const char* ItemForMaterial(VoxelMaterial m, bool soil) {
    switch (m) {
        case VoxelMaterial::Regolith: return soil ? "dirt" : "regolith";
        case VoxelMaterial::Ice: return "ice";
        case VoxelMaterial::Ore: return "iron_ore";
        case VoxelMaterial::Sand: return "sand";
        case VoxelMaterial::Gravel: return "gravel";
        case VoxelMaterial::Grass: return soil ? "dirt" : "regolith"; // turf: soil
        case VoxelMaterial::Snow: return "ice";
        case VoxelMaterial::Clay: return "clay"; // #206
        case VoxelMaterial::Pebbles: return "gravel"; // river stones (#595)
        case VoxelMaterial::Silt: case VoxelMaterial::Mud: case VoxelMaterial::Farmland: return soil ? "dirt" : "regolith"; // floodplain silt, shore mud (#595), tilled dirt (#1085)
        case VoxelMaterial::Granite: return "granite"; // #587: a stratum gives its own rock
        case VoxelMaterial::Basalt: return "basalt";
        case VoxelMaterial::Sandstone: return "sandstone";
        case VoxelMaterial::Limestone: return "limestone";
        case VoxelMaterial::Marble: return "marble";
        case VoxelMaterial::Ash: return "gravel";  // (#1591) volcanic ash and scoria (GUESS: its own item, and sulfur, later)
        case VoxelMaterial::Lava: return "basalt"; // (#1591) the cooled crust over a crater's lava
        default: return "stone";
    }
}

ItemKind KindOf(const std::string& id) {
    const ItemDef* d = itemDB.find(id);
    return d ? d->kind : ItemKind::Resource;
}

// Recipe wildcards: "any_planks" takes planks of any wood, "any_log" any log, "any_stone" stone or any rock (#587).
bool ItemMatches(const std::string& want, const std::string& have) {
    if (want == have) return true;
    auto ends = [&](const char* suf) { const size_t n = strlen(suf); return have.size() > n && have.compare(have.size() - n, n, suf) == 0; };
    if (want == "any_planks") return ends("_planks");
    if (want == "any_log") return ends("_log");
    if (want == "any_stone") return have == "stone" || have == "granite" || have == "basalt" || have == "sandstone" || have == "limestone" || have == "marble"; // #587: any rock
    return false;
}

int StackLimit(const std::string& id) {
    switch (KindOf(id)) {
        case ItemKind::Tool: return 1;
        case ItemKind::Block: return 256;
        default: return 8192;
    }
}

static int64_t ToCenti(double amount) { return (int64_t)llround(amount * 100.0); }

int64_t Inventory::room(const std::string& id) const {
    const int64_t cap = (int64_t)StackLimit(id) * 100;
    int64_t r = 0;
    for (int i = 0; i < kSlots; ++i) { // the grid and the hotbar: nothing is put into the off hand but by hand (#165)
        const InvSlot& s = slots[i];
        if (s.empty()) r += cap;
        else if (s.id == id && !s.hasData()) r += std::max<int64_t>(0, cap - s.centi); // (#1446) never on top of an item with data
    }
    return r;
}

void Inventory::addFractional(const std::string& id, float amount) {
    if (amount <= 0.0f) return;
    // Hundredths go in; less than that waits until enough has been dug.
    float& p = partial[id];
    p += amount;
    const double hundredths = floor(p * 100.0);
    if (hundredths >= 1.0) { p -= (float)(hundredths / 100.0); addAmount(id, hundredths / 100.0); }
}

ItemGlow Item_GlowOf(const std::string& id, float heat) {
    // Constant glows (radioactive materials, #580): resources/data/item_glow.txt "id | r,g,b | strength", read once.
    static std::unordered_map<std::string, ItemGlow> table;
    static bool loaded = false;
    if (!loaded) {
        loaded = true;
        std::ifstream f("resources/data/item_glow.txt");
        std::string line;
        while (std::getline(f, line)) {
            line = Trim(line);
            if (line.empty() || line[0] == '#') continue;
            const auto parts = Split(line, '|');
            if (parts.size() < 3) continue;
            const auto rgb = Split(parts[1], ',');
            if (rgb.size() != 3) continue;
            table[parts[0]] = { Color{ (unsigned char)atoi(rgb[0].c_str()), (unsigned char)atoi(rgb[1].c_str()), (unsigned char)atoi(rgb[2].c_str()), 255 },
                                std::clamp((float)atof(parts[2].c_str()), 0.0f, 1.0f) };
        }
    }
    ItemGlow g;
    if (auto it = table.find(id); it != table.end()) g = it->second;
    if (heat > 0.0f) { // incandescence: dull red when nearly cool, orange, then yellow-white fresh out of the furnace
        const float h = std::clamp(heat, 0.0f, 1.0f);
        auto lerp = [](Color a, Color b, float t) { return Color{ (unsigned char)(a.r + (b.r - a.r) * t), (unsigned char)(a.g + (b.g - a.g) * t), (unsigned char)(a.b + (b.b - a.b) * t), 255 }; };
        const Color c = h < 0.5f ? lerp({ 150, 20, 5, 255 }, { 255, 110, 25, 255 }, h * 2.0f) : lerp({ 255, 110, 25, 255 }, { 255, 225, 150, 255 }, (h - 0.5f) * 2.0f);
        const float s = std::sqrt(h); // stays visible well into the cooling
        if (s > g.strength) g = { c, s };
    }
    return g;
}

// The heat of `a` hundredths at heat ha joined by `b` at hb (#580).
static float MixHeat(int64_t a, float ha, int64_t b, float hb) {
    return a + b > 0 ? (float)(((double)a * ha + (double)b * hb) / (double)(a + b)) : 0.0f;
}

void Inventory::cool(float dt) {
    if (dt <= 0.0f) return;
    for (InvSlot& s : slots) if (s.heat > 0.0f) s.heat = std::max(0.0f, s.heat - dt / kCoolSeconds);
}

double Inventory::addAmount(const std::string& id, double amount, float heat, bool credit) {
    if (amount <= 0.0 || id.empty()) return 0.0;
    if (IsWildcardItem(id)) return amount; // a recipe tag (any_planks), never an item you hold (#269)
    heat = std::clamp(heat, 0.0f, 1.0f);
    // Tools and blocks come in whole units.
    int64_t left = KindOf(id) == ItemKind::Resource ? ToCenti(amount) : (int64_t)floor(amount + 1e-6) * 100;
    const int64_t cap = (int64_t)StackLimit(id) * 100;
    for (int i = 0; i < kSlots; ++i) { // top up stacks of it first (not the off hand's, #165)
        InvSlot& s = slots[i];
        if (left <= 0) break;
        if (s.empty() || s.id != id || s.hasData()) continue; // (#1446) plain stacks only
        const int64_t put = std::min(left, cap - s.centi);
        if (put > 0) { s.heat = MixHeat(s.centi, s.heat, put, heat); s.centi += put; left -= put; }
    }
    // Then empty cells: tools first take an empty hotbar cell (they go to hand, as the old hotbar bound new tools), the
    // rest fill the main grid and only then the hotbar (#209).
    auto fill = [&](int from, int to) {
        for (int i = from; i < to && left > 0; ++i) {
            InvSlot& s = slots[i];
            if (!s.empty()) continue;
            const int64_t put = std::min(left, cap);
            s.id = id; s.centi = put; s.heat = heat; s.data.clear(); left -= put;
        }
    };
    const bool tool = KindOf(id) == ItemKind::Tool;
    if (tool) fill(kGrid, kSlots);
    fill(0, kGrid);
    if (!tool) fill(kGrid, kSlots);
    if (left > 0) overflowed = true;
    const double added = (KindOf(id) == ItemKind::Resource ? ToCenti(amount) : (int64_t)floor(amount + 1e-6) * 100) / 100.0 - left / 100.0;
    if (credit) obtained[id] += added;
    if (onAdded && muted <= 0 && added > 0.0) onAdded(id, added); // the pickup feed (#208)
    return left / 100.0;
}

void Inventory::noteObtained(const std::string& id, double amount, const std::string& data) {
    if (id.empty() || IsWildcardItem(id)) return;
    if (data.empty()) { if (amount > 0.0) obtained[id] += amount; return; }
    obtained[id] += 1.0; // as addInstance: one, and the old tool id a composed tool also counts as
    const std::string legacy = Composed_LegacyId(id, data);
    if (!legacy.empty()) obtained[legacy] += 1.0;
}

bool Inventory::removeAmount(const std::string& id, double amount) {
    const int64_t want = ToCenti(amount);
    if (want <= 0) return true;
    int64_t have = 0;
    for (const InvSlot& s : slots) if (!s.empty() && !s.hasData() && ItemMatches(id, s.id)) have += s.centi;
    if (have < want) return false;
    int64_t left = want;
    auto takeFrom = [&](int i) {
        InvSlot& s = slots[i];
        if (left <= 0 || s.empty() || s.hasData() || !ItemMatches(id, s.id)) return; // (#1446) never an item with data
        const int64_t take = std::min(left, s.centi);
        s.centi -= take; left -= take;
        if (s.centi <= 0) s = InvSlot{};
    };
    // The preferred cell first (the hotbar cell in hand: its tool wears, its blocks run out), then the main grid from
    // its last stacks (so the first stay put), then the other hotbar cells, the off hand last (#165).
    if (preferred >= 0 && preferred < kSlots) takeFrom(preferred);
    for (int i = kGrid - 1; i >= 0 && left > 0; --i) takeFrom(i);
    for (int i = kSlots - 1; i >= kGrid && left > 0; --i) takeFrom(i);
    takeFrom(kOffHand);
    return true;
}

double Inventory::amount(const std::string& id) const {
    int64_t c = 0;
    for (const InvSlot& s : slots) if (!s.empty() && !s.hasData() && ItemMatches(id, s.id)) c += s.centi;
    return c / 100.0;
}

int Inventory::count(const std::string& id) const {
    int64_t c = 0;
    for (const InvSlot& s : slots) if (!s.empty() && !s.hasData() && ItemMatches(id, s.id)) c += s.centi;
    return (int)(c / 100);
}

int Inventory::countInstances(const std::string& id) const {
    int n = 0;
    for (const InvSlot& s : slots) if (!s.empty() && s.hasData() && ItemMatches(id, s.id)) ++n;
    return n;
}

int Inventory::addInstance(const std::string& id, const std::string& data, float heat, bool credit) {
    if (id.empty() || IsWildcardItem(id)) return -1;
    if (data.empty()) { const bool fit = addAmount(id, 1, heat, credit) <= 0.0; return fit ? 0 : -1; }
    auto findFree = [&](int from, int to) { for (int i = from; i < to; ++i) if (slots[i].empty()) return i; return -1; };
    const bool tool = KindOf(id) == ItemKind::Tool;
    int cell = tool ? findFree(kGrid, kSlots) : -1;
    if (cell < 0) cell = findFree(0, kGrid);
    if (cell < 0 && !tool) cell = findFree(kGrid, kSlots);
    if (cell < 0) { overflowed = true; return -1; }
    InvSlot& s = slots[cell];
    s.id = id; s.centi = 100; s.heat = std::clamp(heat, 0.0f, 1.0f); s.data = data;
    if (journal) instanceLog.push_back({ id, +1, data, {} });
    if (credit) obtained[id] += 1.0;
    // (#1449 T1) a composed tool also counts as the old tool its head makes it (tasks say "get stone_pickaxe 1"); the
    // pickup line and a team's book (onAdded) get that name: "Stone pickaxe"
    const std::string legacy = Composed_LegacyId(id, data);
    if (credit && !legacy.empty()) obtained[legacy] += 1.0;
    if (onAdded && muted <= 0) onAdded(legacy.empty() ? id : legacy, 1.0);
    return cell;
}

InvSlot Inventory::takeInstance(int cell) {
    if (cell < 0 || cell >= kCells || slots[cell].empty()) return InvSlot{};
    InvSlot out = slots[cell];
    slots[cell] = InvSlot{};
    if (journal && out.hasData()) instanceLog.push_back({ out.id, -1, out.data, {} });
    return out;
}

bool Inventory::setInstanceData(int cell, const std::string& newData) {
    if (cell < 0 || cell >= kCells || slots[cell].empty() || !slots[cell].hasData() || newData.empty()) return false;
    InvSlot& s = slots[cell];
    if (s.data == newData) return true;
    if (journal) instanceLog.push_back({ s.id, 0, s.data, newData });
    s.data = newData;
    return true;
}

int Inventory::findInstance(const std::string& id, const std::string& data) const {
    for (int i = 0; i < kCells; ++i) if (!slots[i].empty() && slots[i].id == id && slots[i].data == data) return i;
    return -1;
}

std::unordered_map<std::string, int> Inventory::all() const {
    std::unordered_map<std::string, int64_t> c;
    for (const InvSlot& s : slots) if (!s.empty() && !s.hasData()) c[s.id] += s.centi; // (#1446) plain items
    std::unordered_map<std::string, int> out;
    for (auto& kv : c) if (kv.second >= 100) out[kv.first] = (int)(kv.second / 100);
    return out;
}

int Inventory::takeHeadParts(const std::string& id, int n) {
    if (n <= 0) return 100;
    std::vector<int> forged; // cells holding a forged head: best quality first, then the lowest cell
    for (int i = 0; i < kCells; ++i)
        if (!slots[i].empty() && slots[i].hasData() && slots[i].id == id) forged.push_back(i);
    std::stable_sort(forged.begin(), forged.end(), [&](int a, int b) {
        return Composed_PartQuality(slots[a].data) > Composed_PartQuality(slots[b].data);
    });
    const int fromForged = std::min(n, (int)forged.size());
    if (fromForged + count(id) < n) return -1;
    int sum = 100 * (n - fromForged); // plain heads count as 100 %
    for (int k = 0; k < fromForged; ++k) sum += Composed_PartQuality(slots[forged[k]].data);
    for (int k = 0; k < fromForged; ++k) takeInstance(forged[k]);
    if (n > fromForged) removeAmount(id, n - fromForged);
    return sum / n;
}

bool Inventory::canCraft(const Recipe& r) const {
    const bool composed = Composed_Is(r.output);
    for (auto& in : r.inputs) {
        const bool head = composed && Composed_IsHeadPartOf(r.output, in.first); // (#1613) a forged head counts too
        if ((head ? held(in.first) : count(in.first)) < in.second) return false;
    }
    // The output must fit (after the inputs are taken out), or crafting would destroy ingredients.
    Inventory after = *this;
    after.journal = false;
    for (auto& in : r.inputs) {
        if (composed && Composed_IsHeadPartOf(r.output, in.first)) after.takeHeadParts(in.first, in.second);
        else after.removeAmount(in.first, in.second);
    }
    if (composed) { // (#1449) one empty cell a copy
        int empty = 0; for (int k = 0; k < kGrid; ++k) empty += after.slots[k].empty() ? 1 : 0;
        return empty >= std::max(1, r.outputCount);
    }
    return after.room(r.output) >= (int64_t)r.outputCount * 100;
}

bool Inventory::craft(const Recipe& r) {
    if (!canCraft(r)) return false;
    const bool composed = Composed_Is(r.output);
    int quality = 100;
    for (auto& in : r.inputs) {
        if (composed && Composed_IsHeadPartOf(r.output, in.first)) quality = takeHeadParts(in.first, in.second); // (#1613)
        else removeAmount(in.first, in.second);
    }
    if (composed) { // (#1449) the assembly hook: the parts consumed become the copy's record
        std::vector<std::string> parts;
        for (auto& in : r.inputs) parts.push_back(in.first);
        const std::string data = Composed_FromParts(r.output, parts, std::max(quality, kComposedQualityMin));
        if (!data.empty()) { for (int k = 0; k < std::max(1, r.outputCount); ++k) addInstance(r.output, data); return true; }
    }
    addAmount(r.output, r.outputCount);
    return true;
}

void Inventory::swapSlots(int a, int b) {
    if (a == b || a < 0 || b < 0 || a >= kCells || b >= kCells) return;
    InvSlot& from = slots[a];
    InvSlot& to = slots[b];
    if (!from.empty() && !to.empty() && from.id == to.id && !from.hasData() && !to.hasData()) { // (#1446) items with data swap
        const int64_t cap = (int64_t)StackLimit(from.id) * 100;
        const int64_t put = std::min(from.centi, cap - to.centi);
        to.heat = MixHeat(to.centi, to.heat, put, from.heat); // #580
        to.centi += put; from.centi -= put;
        if (from.centi <= 0) from = InvSlot{};
        return;
    }
    std::swap(from, to);
}

void Inventory::sort() {
    std::unordered_map<std::string, int64_t> total;
    std::unordered_map<std::string, double> heatSum; // #580: each item keeps its average heat
    std::vector<InvSlot> own;                          // (#1446) items with data: each keeps its cell's contents, first in id order
    for (int k = 0; k < kGrid; ++k) {
        if (slots[k].empty()) continue;
        if (slots[k].hasData()) { own.push_back(slots[k]); continue; }
        total[slots[k].id] += slots[k].centi; heatSum[slots[k].id] += (double)slots[k].centi * slots[k].heat; // the hotbar stays as it is
    }
    std::stable_sort(own.begin(), own.end(), [](const InvSlot& x, const InvSlot& y) { return x.id < y.id; });
    std::vector<std::pair<std::string, int64_t>> items(total.begin(), total.end());
    std::sort(items.begin(), items.end(), [](const auto& x, const auto& y) {
        const ItemKind kx = KindOf(x.first), ky = KindOf(y.first);
        if (kx != ky) return (int)kx > (int)ky; // tools, then blocks, then resources
        const ItemDef* dx = itemDB.find(x.first);
        const ItemDef* dy = itemDB.find(y.first);
        return (dx ? dx->name : x.first) < (dy ? dy->name : y.first);
    });
    for (int k = 0; k < kGrid; ++k) slots[k] = InvSlot{};
    size_t i = 0;
    for (const InvSlot& o : own) if (i < (size_t)kGrid) slots[i++] = o;
    for (auto& it : items) {
        const int64_t cap = (int64_t)StackLimit(it.first) * 100;
        for (int64_t left = it.second; left > 0 && i < (size_t)kGrid; ++i) {
            const int64_t put = std::min(left, cap);
            slots[i].id = it.first; slots[i].centi = put; slots[i].heat = it.second > 0 ? (float)(heatSum[it.first] / (double)it.second) : 0.0f; left -= put;
        }
    }
}

bool HasSoil(const Planetoid& p) {
    const auto* f = dynamic_cast<const PlanetField*>(p.field.get());
    return f && !f->moon;
}

const char* OreKindOf(const Planetoid& body, Vector3 local) {
    if (const auto* pf = dynamic_cast<const PlanetField*>(body.field.get())) {
        const float d = Vector3Length(local);
        if (d > 1.0f) {
            const float depth = pf->radius + pf->elevation(Vector3Scale(local, 1.0f / d)) - d;
            const int k = pf->oreKind(local, depth);
            if (k >= 0) return PlanetField::OreKindName(k);
        }
        return "iron_ore"; // (ore material where oreKind finds none: only at a different detail level)
    }
    if (body.ringRock >= 0) {
        // Ring rocks (M24): metal-rich - nickel-iron most, then iron, platinum, a little iridium and copper
        // (deposits ~48 m across, as on asteroids).
        const int x = (int)floorf(local.x / 48.0f), y = (int)floorf(local.y / 48.0f), z = (int)floorf(local.z / 48.0f);
        uint32_t h = (uint32_t)x * 73856093u ^ (uint32_t)y * 19349663u ^ (uint32_t)z * 83492791u ^ (body.seed + 17u) * 2654435761u;
        h ^= h >> 16; h *= 0x7feb352du; h ^= h >> 15; h *= 0x846ca68bu; h ^= h >> 16;
        const unsigned r = h % 100;
        return r < 45 ? "nickel_iron_ore" : r < 65 ? "iron_ore" : r < 80 ? "platinum_ore" : r < 86 ? "iridium_ore" : "copper_ore";
    }
    return OreKindAt(local, body.seed);
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#869): stacking, merging, heat ----
#include "unit_helpers.h"

namespace {
struct InventoryTest : ::testing::Test {
    void SetUp() override { UnitTest::EnsureItemDB(); }
    static std::string FirstOfKind(ItemKind k) {
        for (const ItemDef& d : itemDB.items()) if (d.kind == k && !IsWildcardItem(d.id)) return d.id;
        return {};
    }
};
}

TEST_F(InventoryTest, KindsStackToTheirLimits) {
    EXPECT_EQ(KindOf("stone_pickaxe"), ItemKind::Tool);
    EXPECT_EQ(StackLimit("stone_pickaxe"), 1);
    EXPECT_EQ(KindOf("stone"), ItemKind::Resource) << "#812: stone is a resource, not a block";
    EXPECT_EQ(StackLimit("stone"), 8192);
    const std::string block = FirstOfKind(ItemKind::Block);
    ASSERT_FALSE(block.empty());
    EXPECT_EQ(StackLimit(block), 256);
}

TEST_F(InventoryTest, ResourcesKeepHundredths) {
    Inventory inv;
    EXPECT_EQ(inv.addAmount("iron_ore", 1.234), 0.0);
    EXPECT_DOUBLE_EQ(inv.amount("iron_ore"), 1.23);
    EXPECT_EQ(inv.count("iron_ore"), 1);
    inv.addFractional("iron_ore", 0.004f); // under a hundredth: waits
    EXPECT_DOUBLE_EQ(inv.amount("iron_ore"), 1.23);
    inv.addFractional("iron_ore", 0.007f); // now 0.011: one hundredth goes in
    EXPECT_DOUBLE_EQ(inv.amount("iron_ore"), 1.24);
}

TEST_F(InventoryTest, StacksFillUpBeforeNewCells) {
    Inventory inv;
    inv.add("stone", 8000);
    inv.add("stone", 500); // 192 tops up the first stack, 308 starts the second
    EXPECT_EQ(inv.grid()[0].centi, 8192 * 100);
    EXPECT_EQ(inv.grid()[1].centi, 308 * 100);
    EXPECT_EQ(inv.count("stone"), 8500);
}

TEST_F(InventoryTest, ToolsGoToTheHotbarAndNeverStack) {
    Inventory inv;
    inv.add("stone_pickaxe", 2);
    EXPECT_EQ(inv.grid()[Inventory::kGrid].id, "stone_pickaxe");
    EXPECT_EQ(inv.grid()[Inventory::kGrid + 1].id, "stone_pickaxe");
    EXPECT_EQ(inv.grid()[Inventory::kGrid].centi, 100);
    inv.add("stone", 5); // resources fill the main grid first
    EXPECT_EQ(inv.grid()[0].id, "stone");
}

TEST_F(InventoryTest, OverflowReturnsWhatDidNotFit) {
    Inventory inv;
    const double left = inv.addAmount("stone", 8192.0 * Inventory::kSlots + 7.0);
    EXPECT_DOUBLE_EQ(left, 7.0);
    EXPECT_TRUE(inv.overflowed);
    EXPECT_TRUE(inv.offHand().empty()) << "nothing goes into the off hand but by hand (#165)";
    EXPECT_DOUBLE_EQ(inv.obtained["stone"], 8192.0 * Inventory::kSlots);
}

TEST_F(InventoryTest, WildcardsAreTagsNotItems) {
    Inventory inv;
    EXPECT_DOUBLE_EQ(inv.addAmount("any_planks", 3), 3.0);
    EXPECT_TRUE(inv.all().empty());
    inv.add("oak_planks", 4);
    EXPECT_EQ(inv.count("any_planks"), 4);
    EXPECT_TRUE(ItemMatches("any_stone", "granite"));
    EXPECT_FALSE(ItemMatches("any_stone", "dirt"));
    EXPECT_TRUE(inv.remove("any_planks", 3));
    EXPECT_EQ(inv.count("oak_planks"), 1);
}

TEST_F(InventoryTest, RemoveTakesThePreferredCellFirstThenTheLastStacks) {
    Inventory inv;
    inv.grid()[0] = { "stone", 1000 };
    inv.grid()[5] = { "stone", 1000 };
    inv.grid()[Inventory::kGrid + 2] = { "stone", 1000 };
    inv.preferred = Inventory::kGrid + 2;
    ASSERT_TRUE(inv.removeAmount("stone", 15));
    EXPECT_TRUE(inv.grid()[Inventory::kGrid + 2].empty());
    EXPECT_EQ(inv.grid()[5].centi, 500);
    EXPECT_EQ(inv.grid()[0].centi, 1000);
    EXPECT_FALSE(inv.removeAmount("stone", 16)) << "not enough: nothing is taken";
    EXPECT_EQ(inv.count("stone"), 15);
}

TEST_F(InventoryTest, SwapMergesSameItemsUpToTheLimit) {
    Inventory inv;
    inv.grid()[0] = { "stone", 8000 * 100 };
    inv.grid()[1] = { "stone", 500 * 100 };
    inv.swapSlots(1, 0);
    EXPECT_EQ(inv.grid()[0].centi, 8192 * 100);
    EXPECT_EQ(inv.grid()[1].centi, 308 * 100);
    inv.grid()[2] = { "dirt", 100 };
    inv.swapSlots(2, 1); // different items swap
    EXPECT_EQ(inv.grid()[1].id, "dirt");
    EXPECT_EQ(inv.grid()[2].id, "stone");
}

TEST_F(InventoryTest, HeatAveragesByAmountAndCools) {
    Inventory inv;
    inv.addAmount("iron_ingot", 1, 1.0f);
    inv.addAmount("iron_ingot", 3, 0.0f);
    EXPECT_NEAR(inv.grid()[0].heat, 0.25f, 1e-6f);
    inv.grid()[1] = { "iron_ingot", 400, 1.0f };
    inv.swapSlots(1, 0);
    EXPECT_NEAR(inv.grid()[0].heat, (4 * 0.25f + 4 * 1.0f) / 8.0f, 1e-6f);
    inv.cool(Inventory::kCoolSeconds * 0.25f);
    EXPECT_NEAR(inv.grid()[0].heat, 0.625f - 0.25f, 1e-5f);
    inv.cool(Inventory::kCoolSeconds);
    EXPECT_EQ(inv.grid()[0].heat, 0.0f);
}

TEST_F(InventoryTest, SortMergesTheGridAndKeepsAverageHeat) {
    Inventory inv;
    inv.grid()[3] = { "stone", 300, 0.0f };
    inv.grid()[7] = { "stone", 100, 1.0f };
    inv.grid()[9] = { "stone_pickaxe", 100 };
    inv.grid()[Inventory::kGrid] = { "dirt", 100 }; // the hotbar stays as it is
    inv.sort();
    EXPECT_EQ(inv.grid()[0].id, "stone_pickaxe") << "tools first";
    EXPECT_EQ(inv.grid()[1].id, "stone");
    EXPECT_EQ(inv.grid()[1].centi, 400);
    EXPECT_NEAR(inv.grid()[1].heat, 0.25f, 1e-6f);
    EXPECT_TRUE(inv.grid()[2].empty());
    EXPECT_EQ(inv.grid()[Inventory::kGrid].id, "dirt");
}

// (#1219) Weapons from recipes.txt: hafting a head makes one composed copy whose record names the head's material; every
// weapon recipe names a head the class takes (so none crafts a plain, partless weapon).
TEST_F(InventoryTest, HaftingAHeadMakesAComposedWeapon) {
    ASSERT_TRUE(Composed_Load());
    const Recipe* sword = nullptr;
    for (const Recipe& r : itemDB.recipes()) {
        if (r.output.rfind("weapon_", 0) != 0) continue;
        std::vector<std::string> parts; for (const auto& in : r.inputs) parts.push_back(in.first);
        const std::string d = Composed_FromParts(r.output, parts);
        EXPECT_TRUE(Composed_Stats(r.output, d).valid) << r.output << " from " << (parts.empty() ? std::string("?") : parts[0]);
        for (const auto& in : r.inputs) EXPECT_TRUE(itemDB.find(in.first) || IsWildcardItem(in.first)) << in.first << " is an item";
        if (r.output == "weapon_sword" && !r.inputs.empty() && r.inputs[0].first == "bronze_blade") sword = &r;
    }
    ASSERT_NE(sword, nullptr) << "a bronze sword recipe";
    Inventory inv;
    for (const auto& in : sword->inputs) inv.add(ItemMatches(in.first, "oak_planks") ? "oak_planks" : in.first, in.second);
    ASSERT_TRUE(inv.craft(*sword));
    EXPECT_EQ(inv.countInstances("weapon_sword"), 1);
    EXPECT_EQ(inv.count("weapon_sword"), 0) << "no plain sword";
    EXPECT_EQ(inv.count("bronze_blade"), 0) << "the blades went into it";
    const int cell = inv.findInstance("weapon_sword", Composed_MakeData("bronze"));
    EXPECT_GE(cell, 0) << "a bronze head, the default handle and binding";

    // (#1613) Forged blades carry their quality into the sword (2 blades: the average): the best go first, then plain ones.
    ASSERT_EQ(sword->inputs[0].second, 2) << "a sword takes 2 blades";
    auto give = [&](Inventory& v) { for (const auto& in : sword->inputs) if (in.first != "bronze_blade") v.add(ItemMatches(in.first, "oak_planks") ? "oak_planks" : in.first, in.second); };
    Inventory forged;
    give(forged);
    forged.add("bronze_blade", 3);
    forged.addInstance("bronze_blade", "q=112");
    forged.addInstance("bronze_blade", "q=127");
    forged.addInstance("bronze_blade", "q=90");
    EXPECT_EQ(forged.held("bronze_blade"), 6);
    ASSERT_TRUE(forged.craft(*sword));
    EXPECT_GE(forged.findInstance("weapon_sword", Composed_MakeData("bronze", {}, {}, 119)), 0) << "127 and 112 went first: 119";
    EXPECT_EQ(forged.countInstances("bronze_blade"), 1);
    EXPECT_EQ(forged.count("bronze_blade"), 3) << "the plain blades stay";
    give(forged);
    ASSERT_TRUE(forged.craft(*sword));
    EXPECT_GE(forged.findInstance("weapon_sword", Composed_MakeData("bronze", {}, {}, 95)), 0) << "90 and a plain 100: 95";
    give(forged);
    ASSERT_TRUE(forged.craft(*sword));
    EXPECT_GE(forged.findInstance("weapon_sword", Composed_MakeData("bronze")), 0) << "plain blades: 100 %";
    give(forged);
    EXPECT_FALSE(forged.canCraft(*sword)) << "no blades left";

    // only forged blades (no plain ones): they still count
    Inventory only;
    give(only);
    only.addInstance("bronze_blade", "q=80");
    EXPECT_FALSE(only.canCraft(*sword)) << "one blade is not two";
    only.addInstance("bronze_blade", "q=84");
    ASSERT_TRUE(only.craft(*sword));
    EXPECT_GE(only.findInstance("weapon_sword", Composed_MakeData("bronze", {}, {}, 82)), 0) << "Crude blades make a Crude sword";
    EXPECT_EQ(only.held("bronze_blade"), 0);
}

// (#1449 T1) Tools from recipes.txt: no recipe makes an old per-material tool any more; every tool_ recipe hafts a head the
// class takes, and the copy counts as the old tool its head makes it for tasks ("get stone_pickaxe 1", "get bronze_pickaxe 1").
TEST_F(InventoryTest, HaftingAHeadMakesAComposedToolThatCountsAsTheOldOne) {
    ASSERT_TRUE(Composed_Load());
    std::map<std::string, int> heads; // class -> recipes
    const Recipe* stonePick = nullptr;
    for (const Recipe& r : itemDB.recipes()) {
        std::string cls, data;
        EXPECT_FALSE(Composed_FromLegacy(r.output, cls, data)) << r.output << ": an old tool id is never made";
        if (r.output.rfind("tool_", 0) != 0) continue;
        std::vector<std::string> parts; for (const auto& in : r.inputs) parts.push_back(in.first);
        const std::string d = Composed_FromParts(r.output, parts);
        EXPECT_TRUE(Composed_Stats(r.output, d).valid) << r.output << " from " << (parts.empty() ? std::string("?") : parts[0]);
        for (const auto& in : r.inputs) EXPECT_TRUE(itemDB.find(in.first) || IsWildcardItem(in.first)) << in.first << " is an item";
        ++heads[r.output];
        if (r.output == "tool_pickaxe" && !r.inputs.empty() && r.inputs[0].first == "stone_pick_head") stonePick = &r;
    }
    EXPECT_EQ(heads["tool_pickaxe"], 5) << "stone, iron, steel, titanium, bronze";
    EXPECT_EQ(heads["tool_hoe"], 4);
    ASSERT_NE(stonePick, nullptr);
    Inventory inv;
    std::string added;
    inv.onAdded = [&](const std::string& id, double) { added = id; };
    for (const auto& in : stonePick->inputs) inv.add(ItemMatches(in.first, "oak_planks") ? "oak_planks" : in.first, in.second);
    ASSERT_TRUE(inv.craft(*stonePick));
    EXPECT_EQ(inv.countInstances("tool_pickaxe"), 1);
    EXPECT_GE(inv.findInstance("tool_pickaxe", Composed_MakeData("stone")), 0);
    EXPECT_DOUBLE_EQ(inv.obtainedTotal("stone_pickaxe"), 1.0) << "the tier-0 task counts it";
    EXPECT_DOUBLE_EQ(inv.obtainedTotal("tool_pickaxe"), 1.0);
    EXPECT_EQ(added, "stone_pickaxe") << "the pickup line and a team's book get the old name";
    EXPECT_EQ(Composed_LegacyId("tool_axe", Composed_MakeData("bronze")), "bronze_axe");
    EXPECT_EQ(Composed_LegacyId("tool_hoe", Composed_MakeData("bronze")), "") << "there was no bronze hoe";
    EXPECT_EQ(Composed_LegacyId("tool_pickaxe", Composed_MakeData("flint")), "");
    EXPECT_EQ(Composed_LegacyId("weapon_sword", Composed_MakeData("steel")), "");
}

TEST_F(InventoryTest, NoteObtainedCreditsLikeAddAmountAndACreditedStackDoesNotDoubleOnReturn) { // (#1728)
    Inventory a, b;
    a.addAmount("stone", 3.0);
    b.noteObtained("stone", 3.0);
    EXPECT_DOUBLE_EQ(b.obtainedTotal("stone"), a.obtainedTotal("stone")) << "noteObtained credits as addAmount does";
    EXPECT_EQ(b.count("stone"), 0) << "and puts nothing in a cell";
    bool fed = false;
    b.onAdded = [&](const std::string&, double) { fed = true; };
    b.noteObtained("stone", 2.0);
    EXPECT_FALSE(fed) << "not a pickup";
    EXPECT_DOUBLE_EQ(b.obtainedTotal("stone"), 5.0);
    b.noteObtained("stone", 0.0); b.noteObtained("", 4.0);
    EXPECT_DOUBLE_EQ(b.obtainedTotal("stone"), 5.0) << "nothing for a zero amount or no item";
    b.addAmount("stone", 5.0, 0.0f, false); // the credited cursor stack coming back in
    EXPECT_DOUBLE_EQ(b.obtainedTotal("stone"), 5.0) << "credit=false: not counted twice";
    EXPECT_EQ(b.count("stone"), 5);
    b.addAmount("stone", 1.0);
    EXPECT_DOUBLE_EQ(b.obtainedTotal("stone"), 6.0) << "the default still credits";
}

TEST_F(InventoryTest, CraftNeedsRoomForTheOutput) {
    const Recipe* planks = nullptr;
    for (const Recipe& r : itemDB.recipes()) if (r.output == "oak_planks" && r.inputs.size() == 1 && ItemMatches(r.inputs[0].first, "oak_log")) planks = &r;
    ASSERT_NE(planks, nullptr) << "a recipe oak_log -> oak_planks";
    Inventory inv;
    inv.add("oak_log", planks->inputs[0].second);
    ASSERT_TRUE(inv.craft(*planks));
    EXPECT_EQ(inv.count("oak_planks"), planks->outputCount);
    EXPECT_EQ(inv.count("oak_log"), 0);
    // Every cell full of tools but one holding the logs: the logs' cell frees up, so the planks fit.
    Inventory full;
    for (int i = 0; i < Inventory::kSlots; ++i) full.grid()[i] = { "stone_pickaxe", 100 };
    full.grid()[0] = { "oak_log", (int64_t)planks->inputs[0].second * 100 };
    EXPECT_TRUE(full.canCraft(*planks));
    // With twice the logs the cell stays taken: no room, and nothing is used up.
    full.grid()[0].centi *= 2;
    EXPECT_FALSE(full.canCraft(*planks));
    EXPECT_FALSE(full.craft(*planks));
    EXPECT_EQ(full.count("oak_log"), planks->inputs[0].second * 2);
}

TEST_F(InventoryTest, SavedCellsMoveUpFromNineWideSaves) {
    EXPECT_EQ(Inventory::SavedCell(3, 9), 3);
    EXPECT_EQ(Inventory::SavedCell(Inventory::kOldGrid, 9), Inventory::kGrid);
    EXPECT_EQ(Inventory::SavedCell(Inventory::kOldGrid, Inventory::kCols), Inventory::kOldGrid);
}

TEST_F(InventoryTest, LunacreteIsGoneAndOldStacksBecomeStone) {
    EXPECT_EQ(itemDB.find("concrete"), nullptr) << "#1090: no lunacrete on the planet";
    EXPECT_STREQ(RetiredItem("concrete"), "stone");
    EXPECT_STREQ(RetiredItem("stone"), "") << "a current item keeps its id";
    for (const Recipe& r : itemDB.recipes())
        for (const auto& in : r.inputs) EXPECT_NE(in.first, "concrete") << r.output;
}

// (#1250) magic_items.txt adds its items after the table; items.txt's own stay as they were.
TEST(ItemDBExtra, MagicItemsAddWithoutReplacing) {
    ItemDB db;
    ASSERT_TRUE(db.load("resources/data/items.txt", "resources/data/recipes.txt"));
    const size_t before = db.items().size();
    EXPECT_EQ(db.find("first_sight_page"), nullptr) << "not in the base game";
    ASSERT_TRUE(db.loadExtra("resources/data/magic_items.txt"));
    ASSERT_NE(db.find("first_sight_page"), nullptr);
    ASSERT_NE(db.find("resonant_crystal"), nullptr);
    const size_t withMagic = db.items().size();
    EXPECT_GT(withMagic, before + 2);
    // (#1252) the magic machines are placeable blocks, jars are resources, and their recipes load with them
    ASSERT_NE(db.find("mana_siphon"), nullptr);
    EXPECT_EQ(db.find("mana_siphon")->kind, ItemKind::Block);
    EXPECT_EQ(db.find("research_table")->kind, ItemKind::Block);
    EXPECT_EQ(db.find("jar_storm")->kind, ItemKind::Resource);
    int magicRecipes = 0;
    for (const Recipe& r : db.recipes()) magicRecipes += r.output == "mana_siphon" || r.output == "mana_dynamo" || r.output == "research_table" || r.output == "mana_jar";
    EXPECT_EQ(magicRecipes, 4) << "every magic recipe's items are known";
    const std::string stoneName = db.find("stone")->name;
    EXPECT_TRUE(db.loadExtra("resources/data/items.txt")) << "loading known ids again";
    EXPECT_EQ(db.items().size(), withMagic) << "known ids are kept, not added twice";
    EXPECT_EQ(db.find("stone")->name, stoneName);
    EXPECT_FALSE(db.loadExtra("resources/data/no_such_file.txt"));
}

// (#1446) Items with their own data: one per cell, never merged, never taken to pay for anything.
TEST_F(InventoryTest, ItemsWithDataStaySeparateAndAreNeverSpent) {
    Inventory inv;
    const int a = inv.addInstance("stone_pickaxe", "h=flint;w=3");
    const int b = inv.addInstance("stone_pickaxe", "h=granite;w=0");
    ASSERT_GE(a, 0); ASSERT_GE(b, 0);
    EXPECT_NE(a, b) << "each in its own cell";
    EXPECT_EQ(inv.grid()[a].data, "h=flint;w=3");
    EXPECT_EQ(inv.grid()[a].centi, 100) << "exactly one unit";
    EXPECT_EQ(inv.countInstances("stone_pickaxe"), 2);
    // A plain one stacks on neither (pickaxes stack 1 anyway: use a resource for the merge rules).
    EXPECT_EQ(inv.count("stone_pickaxe"), 0) << "plain counts ignore items with data";
    EXPECT_FALSE(inv.remove("stone_pickaxe", 1)) << "a recipe can't take one";
    EXPECT_EQ(inv.countInstances("stone_pickaxe"), 2);

    // Resources: a stack with data keeps its unit, plain stone fills other cells.
    const int s = inv.addInstance("stone", "mark=1");
    inv.add("stone", 5);
    EXPECT_EQ(inv.grid()[s].centi, 100) << "nothing was added on top of it";
    EXPECT_EQ(inv.count("stone"), 5);
    ASSERT_TRUE(inv.remove("stone", 5));
    EXPECT_EQ(inv.grid()[s].data, "mark=1") << "paying 5 stone left it alone";
    // Two cells of the same id with data swap instead of merging.
    const int s2 = inv.addInstance("stone", "mark=2");
    inv.swapSlots(s, s2);
    EXPECT_EQ(inv.grid()[s].data, "mark=2");
    EXPECT_EQ(inv.grid()[s2].data, "mark=1");
    // Sorting keeps each one whole.
    inv.sort();
    EXPECT_EQ(inv.countInstances("stone"), 2);
    // Taking one hands back its record.
    int cell = -1;
    for (int i = 0; i < Inventory::kCells; ++i) if (inv.grid()[i].data == "mark=2") cell = i;
    ASSERT_GE(cell, 0);
    const InvSlot taken = inv.takeInstance(cell);
    EXPECT_EQ(taken.data, "mark=2");
    EXPECT_TRUE(inv.grid()[cell].empty());
    EXPECT_EQ(inv.countInstances("stone"), 1);
}

// (#1446 step 2) A server's journal of items with data: what the client is told, in order.
TEST_F(InventoryTest, InstanceJournalLogsEveryChange) {
    Inventory inv;
    inv.journal = true;
    const int a = inv.addInstance("stone_pickaxe", "h=flint;w=0");
    ASSERT_GE(a, 0);
    EXPECT_TRUE(inv.setInstanceData(a, "h=flint;w=1"));
    EXPECT_TRUE(inv.setInstanceData(a, "h=flint;w=1")) << "unchanged: no entry";
    inv.add("stone", 3); // plain: never in the journal
    const int b = inv.addInstance("stone", "mark=1");
    inv.takeInstance(a);
    inv.clear();
    ASSERT_EQ(inv.instanceLog.size(), 5u);
    EXPECT_EQ(inv.instanceLog[0].sign, +1); EXPECT_EQ(inv.instanceLog[0].data, "h=flint;w=0");
    EXPECT_EQ(inv.instanceLog[1].sign, 0); EXPECT_EQ(inv.instanceLog[1].data, "h=flint;w=0"); EXPECT_EQ(inv.instanceLog[1].newData, "h=flint;w=1");
    EXPECT_EQ(inv.instanceLog[2].sign, +1); EXPECT_EQ(inv.instanceLog[2].id, "stone");
    EXPECT_EQ(inv.instanceLog[3].sign, -1); EXPECT_EQ(inv.instanceLog[3].data, "h=flint;w=1");
    EXPECT_EQ(inv.instanceLog[4].sign, -1); EXPECT_EQ(inv.instanceLog[4].data, "mark=1") << "clear() logs what it removes";
    (void)b;
    Inventory quiet;
    quiet.addInstance("stone", "x=1");
    EXPECT_TRUE(quiet.instanceLog.empty()) << "no journal unless asked (single player, clients)";
}

TEST_F(InventoryTest, AddInstanceFailsWhenNoCellIsFree) {
    Inventory inv;
    for (int i = 0; i < Inventory::kSlots; ++i) ASSERT_GE(inv.addInstance("stone", "n=" + std::to_string(i)), 0);
    inv.overflowed = false;
    EXPECT_EQ(inv.addInstance("stone", "n=extra"), -1);
    EXPECT_TRUE(inv.overflowed);
    EXPECT_TRUE(inv.grid()[Inventory::kOffHand].empty()) << "never into the off hand";
}
#endif
