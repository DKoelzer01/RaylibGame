#include "recipebook.h"
#include "researchgate.h"
#include "playerinput.h"
#include "playersession.h"
#include "simclock.h"
#include "gameui.h"
#include "inventory.h"
#include "survival.h"
#include "terraintool.h"
#include "audio.h"
#include <algorithm>

#include <raylib.h>
#include <cctype>
#include <cstdio>
#include <map>
#include <set>

namespace {
    struct Unlock { std::string kind, target; double n = 1.0; };
    std::map<std::string, Unlock> unlocks; // output -> its extra condition
    bool loaded = false;
    RecipeBookState& R() { return LocalSession().recipes; } // the local player's (#313)

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
    return RecipeBook_Unlocked(r); // its recipe_unlocks condition, and (#1271, off until L6) its research node
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

void RecipeBook_Reset() { R().known.clear(); R().primed = false; R().lastCheck = -10.0; }

int RecipeBook_KnownCount() { int n = 0; for (const Recipe& r : itemDB.recipes()) n += RecipeBook_Known(r) ? 1 : 0; return n; }

void RecipeBook_Update() {
    if (SimTime() - R().lastCheck < 1.0 || terrainTool.creative) return;
    R().lastCheck = SimTime();
    std::set<std::string> now;
    for (const Recipe& r : itemDB.recipes()) if (RecipeBook_Known(r)) now.insert(r.output);
    if (R().primed) {
        std::string list;
        int fresh = 0;
        for (const std::string& o : now)
            if (!R().known.count(o)) {
                if (fresh < 3) { const ItemDef* d = itemDB.find(o); list += (list.empty() ? "" : ", ") + (d ? d->name : o); }
                ++fresh;
            }
        if (fresh > 0) { const std::string key = PlayerInput_KeyName(PlayerInput_Binding(Act::Crafting)); GameUI_Toast(TextFormat("New recipe%s: %s%s (%s: crafting)", fresh > 1 ? "s" : "", list.c_str(), fresh > 3 ? TextFormat(" and %d more", fresh - 3) : "", key.c_str())); } // (#827: the key as bound)
    }
    R().known.swap(now);
    R().primed = true;
}

bool RecipeBook_Unlocked(const Recipe& r) {
    Load();
    if (ResearchGate_On() && !ResearchGate_RecipeAllowed(ResearchGate_Team(), r.output)) return false; // (#1271) its research node
    auto it = unlocks.find(r.output);
    return it == unlocks.end() || QuestBook_Met(it->second.kind, it->second.target, it->second.n);
}

bool RecipeBook_LockedByResearch(const Recipe& r) {
    if (terrainTool.creative || !ResearchGate_On() || ResearchGate_RecipeAllowed(ResearchGate_Team(), r.output)) return false;
    Load();
    for (const auto& in : r.inputs) if (!RecipeBook_Discovered(in.first)) return false; // (not found yet: still to discover)
    auto it = unlocks.find(r.output);
    return it == unlocks.end() || QuestBook_Met(it->second.kind, it->second.target, it->second.n);
}

std::string RecipeBook_LockReason(const Recipe& r) {
    if (!RecipeBook_LockedByResearch(r)) return "";
    const ResearchNode* n = ResearchGate_NodeForRecipe(r.output);
    if (!n) n = ResearchGate_NodeForMachineItem(r.output);
    return n ? "Locked: research " + n->name + " (tier " + std::to_string(n->tier) + ", Research tab)" : "Locked: research it first";
}

bool Crafting_Request(int index, int count) {
    const auto& recipes = itemDB.recipes();
    if (index < 0 || index >= (int)recipes.size() || count < 1) return false;
    const Recipe& r = recipes[index];
    if (gCraftRequest) { gCraftRequest(index, r.output, std::min(count, 64), true, false); return true; } // (#854: from nearby chests too, as the screen counts them) // the server crafts and says so
    int made = 0;
    while (made < count && inventory.canCraft(r) && inventory.craft(r)) ++made;
    if (made > 0) Audio_Play(Sfx::Craft);
    return made > 0;
}

bool Crafting_RequestRecursive(int index, int count) {
    const auto& recipes = itemDB.recipes();
    if (!gCraftRequest || index < 0 || index >= (int)recipes.size() || count < 1) return false;
    gCraftRequest(index, recipes[index].output, std::min(count, 64), false, true); // (#946) one request: the server plans it whole
    return true;
}
