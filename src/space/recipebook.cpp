#include "recipebook.h"
#include "gameui.h"
#include "inventory.h"
#include "survival.h"
#include "terraintool.h"

#include <raylib.h>
#include <cctype>
#include <cstdio>
#include <map>
#include <set>

namespace {
    struct Unlock { std::string kind, target; double n = 1.0; };
    std::map<std::string, Unlock> unlocks; // output -> its extra condition
    bool loaded = false;
    std::set<std::string> known;           // outputs known at the last check
    bool primed = false;                   // the first check after a reset is silent
    double lastCheck = -10.0;

    void Load() {
        if (loaded) return;
        loaded = true;
        FILE* f = fopen("resources/data/recipe_unlocks.txt", "r");
        if (!f) return;
        char buf[512];
        while (fgets(buf, sizeof buf, f)) {
            if (buf[0] == '#' || buf[0] == '\n' || buf[0] == '\r') continue;
            char out[64] = {}, kind[32] = {}, target[64] = {};
            double n = 1.0;
            const int got = sscanf(buf, " %63[^ |] | %31s %63s %lf", out, kind, target, &n);
            if (got >= 2) unlocks[out] = { kind, got >= 3 ? target : "", got >= 4 ? n : 1.0 };
        }
        fclose(f);
    }
    std::string Lower(std::string s) { for (char& c : s) c = (char)tolower((unsigned char)c); return s; }
}

bool RecipeBook_Discovered(const std::string& item) { return inventory.obtainedTotal(item) > 0.0 || inventory.count(item) > 0; }

bool RecipeBook_Known(const Recipe& r) {
    if (terrainTool.creative) return true;
    Load();
    for (const auto& in : r.inputs) if (!RecipeBook_Discovered(in.first)) return false;
    auto it = unlocks.find(r.output);
    return it == unlocks.end() || QuestBook_Met(it->second.kind, it->second.target, it->second.n);
}

bool RecipeBook_Matches(const Recipe& r, const std::string& q) {
    if (q.empty()) return true;
    auto has = [&](const std::string& id) {
        if (Lower(id).find(q) != std::string::npos) return true;
        const ItemDef* d = itemDB.find(id);
        return d && Lower(d->name).find(q) != std::string::npos;
    };
    if (has(r.output)) return true;
    for (const auto& in : r.inputs) if (has(in.first)) return true;
    return false;
}

void RecipeBook_Reset() { known.clear(); primed = false; lastCheck = -10.0; }

int RecipeBook_KnownCount() { int n = 0; for (const Recipe& r : itemDB.recipes()) n += RecipeBook_Known(r) ? 1 : 0; return n; }

void RecipeBook_Update() {
    if (GetTime() - lastCheck < 1.0 || terrainTool.creative) return;
    lastCheck = GetTime();
    std::set<std::string> now;
    for (const Recipe& r : itemDB.recipes()) if (RecipeBook_Known(r)) now.insert(r.output);
    if (primed) {
        std::string list;
        int fresh = 0;
        for (const std::string& o : now)
            if (!known.count(o)) {
                if (fresh < 3) { const ItemDef* d = itemDB.find(o); list += (list.empty() ? "" : ", ") + (d ? d->name : o); }
                ++fresh;
            }
        if (fresh > 0) GameUI_Toast(TextFormat("New recipe%s: %s%s (Tab)", fresh > 1 ? "s" : "", list.c_str(), fresh > 3 ? TextFormat(" and %d more", fresh - 3) : ""));
    }
    known.swap(now);
    primed = true;
}
