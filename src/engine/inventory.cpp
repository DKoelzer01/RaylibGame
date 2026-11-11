#include "inventory.h"
#include "utils/logger.h"
#include "world.h"
#include "voxel.h"
#include "plots.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

ItemDB itemDB;
Inventory inventory;

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

bool ItemDB::load(const char* itemsPath, const char* recipesPath) {
    itemList.clear();
    byId.clear();
    recipeList.clear();

    std::ifstream items(itemsPath);
    if (!items) { logger.logf("[ERROR] Cannot open %s\n", itemsPath); return false; }
    std::string line;
    while (std::getline(items, line)) {
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;
        const auto f = Split(line, '|');
        if (f.size() < 5) { logger.logf("[WARN] items: bad line '%s'\n", line.c_str()); continue; }
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
    if (!recipes) { logger.logf("[ERROR] Cannot open %s\n", recipesPath); return false; }
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
        else logger.logf("[WARN] recipes: unknown item in '%s'\n", line.c_str());
    }
    // Placeable parts (machines, conveyors, cables, rocket parts) are blocks: stacks of 256.
    for (int t = 0; t < (int)PartType::Count; ++t)
        if (IsShape((PartType)t)) continue; else
        if (const char* item = PartCostItem((PartType)t)) {
            auto it = byId.find(item);
            if (it != byId.end() && itemList[it->second].kind == ItemKind::Resource) itemList[it->second].kind = ItemKind::Block;
        }
    logger.logf("[Items] %zu items, %zu recipes\n", itemList.size(), recipeList.size());
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
        default: return "stone";
    }
}

ItemKind KindOf(const std::string& id) {
    const ItemDef* d = itemDB.find(id);
    return d ? d->kind : ItemKind::Resource;
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
    for (const InvSlot& s : slots) {
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

double Inventory::addAmount(const std::string& id, double amount) {
    if (amount <= 0.0 || id.empty()) return 0.0;
    // Tools and blocks come in whole units.
    int64_t left = KindOf(id) == ItemKind::Resource ? ToCenti(amount) : (int64_t)floor(amount + 1e-6) * 100;
    const int64_t cap = (int64_t)StackLimit(id) * 100;
    for (InvSlot& s : slots) {       // top up stacks of it first
        if (left <= 0) break;
        if (s.empty() || s.id != id) continue;
        const int64_t put = std::min(left, cap - s.centi);
        if (put > 0) { s.centi += put; left -= put; }
    }
    for (InvSlot& s : slots) {       // then fill empty cells in order
        if (left <= 0) break;
        if (!s.empty()) continue;
        const int64_t put = std::min(left, cap);
        s.id = id; s.centi = put; left -= put;
    }
    if (left > 0) overflowed = true;
    obtained[id] += (KindOf(id) == ItemKind::Resource ? ToCenti(amount) : (int64_t)floor(amount + 1e-6) * 100) / 100.0 - left / 100.0;
    return left / 100.0;
}

bool Inventory::removeAmount(const std::string& id, double amount) {
    const int64_t want = ToCenti(amount);
    if (want <= 0) return true;
    int64_t have = 0;
    for (const InvSlot& s : slots) if (!s.empty() && s.id == id) have += s.centi;
    if (have < want) return false;
    int64_t left = want;
    for (int i = kSlots - 1; i >= 0 && left > 0; --i) { // from the last stacks, so the first stay put
        InvSlot& s = slots[i];
        if (s.empty() || s.id != id) continue;
        const int64_t take = std::min(left, s.centi);
        s.centi -= take; left -= take;
        if (s.centi <= 0) s = InvSlot{};
    }
    return true;
}

double Inventory::amount(const std::string& id) const {
    int64_t c = 0;
    for (const InvSlot& s : slots) if (!s.empty() && s.id == id) c += s.centi;
    return c / 100.0;
}

int Inventory::count(const std::string& id) const {
    int64_t c = 0;
    for (const InvSlot& s : slots) if (!s.empty() && s.id == id) c += s.centi;
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
    if (a == b || a < 0 || b < 0 || a >= kSlots || b >= kSlots) return;
    InvSlot& from = slots[a];
    InvSlot& to = slots[b];
    if (!from.empty() && !to.empty() && from.id == to.id) {
        const int64_t cap = (int64_t)StackLimit(from.id) * 100;
        const int64_t put = std::min(from.centi, cap - to.centi);
        to.centi += put; from.centi -= put;
        if (from.centi <= 0) from = InvSlot{};
        return;
    }
    std::swap(from, to);
}

void Inventory::sort() {
    std::unordered_map<std::string, int64_t> total;
    for (const InvSlot& s : slots) if (!s.empty()) total[s.id] += s.centi;
    std::vector<std::pair<std::string, int64_t>> items(total.begin(), total.end());
    std::sort(items.begin(), items.end(), [](const auto& x, const auto& y) {
        const ItemKind kx = KindOf(x.first), ky = KindOf(y.first);
        if (kx != ky) return (int)kx > (int)ky; // tools, then blocks, then resources
        const ItemDef* dx = itemDB.find(x.first);
        const ItemDef* dy = itemDB.find(y.first);
        return (dx ? dx->name : x.first) < (dy ? dy->name : y.first);
    });
    for (auto& s : slots) s = InvSlot{};
    size_t i = 0;
    for (auto& it : items) {
        const int64_t cap = (int64_t)StackLimit(it.first) * 100;
        for (int64_t left = it.second; left > 0 && i < slots.size(); ++i) {
            const int64_t put = std::min(left, cap);
            slots[i].id = it.first; slots[i].centi = put; left -= put;
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
    return OreKindAt(local, body.seed);
}
