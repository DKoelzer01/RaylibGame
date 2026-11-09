#include "inventory.h"
#include "utils/logger.h"

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
        d.stackLimit = std::max(1, atoi(f[2].c_str()));
        const auto rgb = Split(f[3], ',');
        if (rgb.size() == 3) d.color = { (unsigned char)atoi(rgb[0].c_str()), (unsigned char)atoi(rgb[1].c_str()), (unsigned char)atoi(rgb[2].c_str()), 255 };
        d.description = f[4];
        byId[d.id] = itemList.size();
        itemList.push_back(d);
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
    logger.logf("[Items] %zu items, %zu recipes\n", itemList.size(), recipeList.size());
    return true;
}

const ItemDef* ItemDB::find(const std::string& id) const {
    auto it = byId.find(id);
    return it == byId.end() ? nullptr : &itemList[it->second];
}

const char* ItemForMaterial(VoxelMaterial m) {
    switch (m) {
        case VoxelMaterial::Regolith: return "regolith";
        case VoxelMaterial::Ice: return "ice";
        case VoxelMaterial::Ore: return "iron_ore";
        default: return "stone";
    }
}

void Inventory::addFractional(const std::string& id, float amount) {
    float& p = partial[id];
    p += amount;
    const int whole = (int)floorf(p);
    if (whole > 0) { p -= whole; add(id, whole); }
}

void Inventory::add(const std::string& id, int n) {
    if (n <= 0) return;
    const ItemDef* d = itemDB.find(id);
    const int limit = d ? d->stackLimit : 999;
    int& c = counts[id];
    c = std::min(limit, c + n);
}

bool Inventory::remove(const std::string& id, int n) {
    auto it = counts.find(id);
    if (it == counts.end() || it->second < n) return false;
    it->second -= n;
    if (it->second == 0) counts.erase(it);
    return true;
}

int Inventory::count(const std::string& id) const {
    auto it = counts.find(id);
    return it == counts.end() ? 0 : it->second;
}

bool Inventory::canCraft(const Recipe& r) const {
    for (auto& in : r.inputs) if (count(in.first) < in.second) return false;
    // The whole output must fit under the stack limit (after its own inputs are consumed, in case
    // a recipe uses its output item), or the craft would silently destroy ingredients.
    const ItemDef* d = itemDB.find(r.output);
    const int limit = d ? d->stackLimit : 999;
    int have = count(r.output);
    for (auto& in : r.inputs) if (in.first == r.output) have -= in.second;
    return have + r.outputCount <= limit;
}

bool Inventory::craft(const Recipe& r) {
    if (!canCraft(r)) return false;
    for (auto& in : r.inputs) remove(in.first, in.second);
    add(r.output, r.outputCount);
    return true;
}
