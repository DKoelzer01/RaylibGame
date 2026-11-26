#include "inventory.h"
#include "utils/logger.h"
#include "world.h"
#include "voxel.h"
#include "plots.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>

ItemDB itemDB;
// (inventory: the local PlayerSession's, src/space/playersession.cpp, #313)

namespace {
// A placeable part's item stacks as a block (256). Not shapes, the lever, or the parts only the world places (#812): a
// meteorite and dungeon traps name the material they break into (nickel-iron ore, stone), which must stay a resource.
bool PartItemIsBlock(PartType t) {
    return !IsShape(t) && t != PartType::Meteorite && t != PartType::PressurePlate && t != PartType::GasVent && t != PartType::Lever;
}
} // namespace

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
    // Placeable parts named here are blocks, as load() makes them (#1858 M16: the same parts skipped, or stone and nickel-iron
    // ore from the world-only parts became blocks stacking to 256).
    for (int t = 0; t < (int)PartType::Count; ++t)
        if (const char* item = PartCostItem((PartType)t); item && PartItemIsBlock((PartType)t))
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
    // Placeable parts (machines, conveyors, cables, rocket parts) are blocks: stacks of 256. Not the parts only the world
    // places (#812): a meteorite and dungeon traps name the material they break into (nickel-iron ore, stone), which must
    // stay a resource (stone as a block stacked to 256 and lost nothing on death).
    for (int t = 0; t < (int)PartType::Count; ++t)
        if (!PartItemIsBlock((PartType)t)) continue; else
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
        else if (s.id == id) r += std::max<int64_t>(0, cap - s.centi);
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

double Inventory::addAmount(const std::string& id, double amount, float heat) {
    if (amount <= 0.0 || id.empty()) return 0.0;
    if (IsWildcardItem(id)) return amount; // a recipe tag (any_planks), never an item you hold (#269)
    heat = std::clamp(heat, 0.0f, 1.0f);
    // Tools and blocks come in whole units.
    int64_t left = KindOf(id) == ItemKind::Resource ? ToCenti(amount) : (int64_t)floor(amount + 1e-6) * 100;
    const int64_t cap = (int64_t)StackLimit(id) * 100;
    for (int i = 0; i < kSlots; ++i) { // top up stacks of it first (not the off hand's, #165)
        InvSlot& s = slots[i];
        if (left <= 0) break;
        if (s.empty() || s.id != id) continue;
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
            s.id = id; s.centi = put; s.heat = heat; left -= put;
        }
    };
    const bool tool = KindOf(id) == ItemKind::Tool;
    if (tool) fill(kGrid, kSlots);
    fill(0, kGrid);
    if (!tool) fill(kGrid, kSlots);
    if (left > 0) overflowed = true;
    const double added = (KindOf(id) == ItemKind::Resource ? ToCenti(amount) : (int64_t)floor(amount + 1e-6) * 100) / 100.0 - left / 100.0;
    obtained[id] += added;
    if (onAdded && muted <= 0 && added > 0.0) onAdded(id, added); // the pickup feed (#208)
    return left / 100.0;
}

bool Inventory::removeAmount(const std::string& id, double amount) {
    const int64_t want = ToCenti(amount);
    if (want <= 0) return true;
    int64_t have = 0;
    for (const InvSlot& s : slots) if (!s.empty() && ItemMatches(id, s.id)) have += s.centi;
    if (have < want) return false;
    int64_t left = want;
    auto takeFrom = [&](int i) {
        InvSlot& s = slots[i];
        if (left <= 0 || s.empty() || !ItemMatches(id, s.id)) return;
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
    for (const InvSlot& s : slots) if (!s.empty() && ItemMatches(id, s.id)) c += s.centi;
    return c / 100.0;
}

int Inventory::count(const std::string& id) const {
    int64_t c = 0;
    for (const InvSlot& s : slots) if (!s.empty() && ItemMatches(id, s.id)) c += s.centi;
    return (int)(c / 100);
}

std::unordered_map<std::string, int> Inventory::all() const {
    std::unordered_map<std::string, int64_t> c;
    for (const InvSlot& s : slots) if (!s.empty()) c[s.id] += s.centi;
    std::unordered_map<std::string, int> out;
    for (auto& kv : c) if (kv.second >= 100) out[kv.first] = (int)(kv.second / 100);
    return out;
}

bool Inventory::canCraft(const Recipe& r) const {
    for (auto& in : r.inputs) if (count(in.first) < in.second) return false;
    // The output must fit (after the inputs are taken out), or crafting would destroy ingredients.
    Inventory after = *this;
    for (auto& in : r.inputs) after.removeAmount(in.first, in.second);
    return after.room(r.output) >= (int64_t)r.outputCount * 100;
}

bool Inventory::craft(const Recipe& r) {
    if (!canCraft(r)) return false;
    for (auto& in : r.inputs) removeAmount(in.first, in.second);
    addAmount(r.output, r.outputCount);
    return true;
}

void Inventory::swapSlots(int a, int b) {
    if (a == b || a < 0 || b < 0 || a >= kCells || b >= kCells) return;
    InvSlot& from = slots[a];
    InvSlot& to = slots[b];
    if (!from.empty() && !to.empty() && from.id == to.id) {
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
    for (int k = 0; k < kGrid; ++k) if (!slots[k].empty()) { total[slots[k].id] += slots[k].centi; heatSum[slots[k].id] += (double)slots[k].centi * slots[k].heat; } // the hotbar stays as it is
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
// (#1858 M16) Loading extra text keeps the world-only parts' items resources, as load() does: stone (a pressure plate's
// material) and nickel-iron ore (a meteorite's) must not become blocks stacking to 256.
TEST(ItemDBExtra, ExtraTextKeepsWorldPartItemsResources) {
    ItemDB db;
    ASSERT_TRUE(db.load("resources/data/items.txt", "resources/data/recipes.txt"));
    std::vector<std::pair<std::string, ItemKind>> before;
    for (int t = 0; t < (int)PartType::Count; ++t)
        if (const char* item = PartCostItem((PartType)t)) if (const ItemDef* it = db.find(item)) before.push_back({ item, it->kind });
    ASSERT_TRUE(db.loadExtraText(""));
    for (const auto& b : before) EXPECT_EQ(db.find(b.first)->kind, b.second) << b.first << " changed kind on an extra load";
    if (const char* meteor = PartCostItem(PartType::Meteorite)) EXPECT_EQ(db.find(meteor)->kind, ItemKind::Resource) << meteor;
}

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
#endif
