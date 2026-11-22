#include "machinetiers.h"
#include "inventory.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <mutex>
#include <sstream>

namespace {
struct Tier { std::string name; int volts = 0; Color color = GRAY; };
struct Family { PartType type = PartType::Count; int first = 0; std::vector<int> tiers; std::map<int, std::string> names; };
std::vector<Tier> tiers;
std::vector<Family> families;
MachineOverclock overclock;
int yields[6] = { 0, 1, 10, 50, 75, 100 };
bool loaded = false;
std::once_flag once;

std::string Lower(std::string s) { for (char& c : s) c = (char)std::tolower((unsigned char)c); return s; }
std::vector<std::string> Cols(const std::string& line) {
    std::vector<std::string> cols;
    std::stringstream ss(line);
    std::string c;
    while (std::getline(ss, c, '|')) {
        const size_t a = c.find_first_not_of(" \t"), b = c.find_last_not_of(" \t\r");
        cols.push_back(a == std::string::npos ? "" : c.substr(a, b - a + 1));
    }
    return cols;
}
int TierByName(const std::string& name) {
    const std::string n = Lower(name);
    for (int t = 0; t < (int)tiers.size(); ++t) if (Lower(tiers[t].name) == n) return t;
    return -1;
}
const Family* FamilyOf(PartType t) { MachineTiers_Load(); for (const Family& f : families) if (f.type == t) return &f; return nullptr; }

void Parse(std::istream& in) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> c = Cols(line);
        if (c.size() < 2) continue;
        if (c[0] == "tier" && c.size() >= 3) {
            Tier t{ c[1], atoi(c[2].c_str()) };
            int r = 0, g = 0, b = 0;
            if (c.size() >= 4 && sscanf(c[3].c_str(), "%d,%d,%d", &r, &g, &b) == 3) t.color = { (unsigned char)r, (unsigned char)g, (unsigned char)b, 255 }; // (L7d)
            tiers.push_back(t);
            continue;
        }
        if (c[0] == "overclock" && c.size() >= 4) { overclock = { std::max(1, atoi(c[1].c_str())), std::max(1, atoi(c[2].c_str())), std::max(0, atoi(c[3].c_str())) }; continue; }
        if (c[0] == "yield" && c.size() >= 7) { for (int k = 0; k < 6; ++k) yields[k] = std::clamp(atoi(c[k + 1].c_str()), 0, 100); continue; }
        if (c.size() < 3) continue;
        Family f;
        for (int t = 0; t < (int)PartType::Count; ++t) { const char* m = PartMachineName((PartType)t); if (m && c[0] == m) { f.type = (PartType)t; break; } }
        f.first = TierByName(c[1]);
        if (f.type == PartType::Count || f.first < 0) { printf("machine_tiers: skipped \"%s\"\n", line.c_str()); continue; }
        std::stringstream ts(c[2]);
        std::string tn;
        while (ts >> tn) { const int t = TierByName(tn); if (t >= f.first && std::find(f.tiers.begin(), f.tiers.end(), t) == f.tiers.end()) f.tiers.push_back(t); }
        std::sort(f.tiers.begin(), f.tiers.end());
        if (f.tiers.empty() || f.tiers[0] != f.first) f.tiers.insert(f.tiers.begin(), f.first);
        if (c.size() >= 4) {
            std::stringstream ns(c[3]);
            std::string kv;
            while (ns >> kv) { const size_t colon = kv.find(':'); if (colon != std::string::npos) { const int t = TierByName(kv.substr(0, colon)); if (t >= 0) f.names[t] = kv.substr(colon + 1); } }
        }
        families.push_back(f);
    }
}
} // namespace

bool MachineTiers_Load() {
    std::call_once(once, [] {
        std::ifstream f("resources/data/machine_tiers.txt");
        if (f) { Parse(f); loaded = true; }
        if (tiers.empty()) tiers = { { "LV", 32 }, { "MV", 128 }, { "HV", 512 }, { "EV", 2048 }, { "IV", 8192 } };
    });
    return loaded;
}
int MachineTiers_Count() { MachineTiers_Load(); return (int)tiers.size(); }
const char* MachineTiers_Name(int t) { MachineTiers_Load(); return t >= 0 && t < (int)tiers.size() ? tiers[t].name.c_str() : ""; }
int MachineTiers_Volts(int t) { MachineTiers_Load(); return t >= 0 && t < (int)tiers.size() ? tiers[t].volts : 0; }
int MachineTiers_ByName(const std::string& name) { MachineTiers_Load(); return TierByName(name); }
int MachineTiers_ForVolts(int volts) {
    MachineTiers_Load();
    int best = -1;
    for (int t = 0; t < (int)tiers.size(); ++t) if (tiers[t].volts <= volts) best = t;
    return best;
}
Color MachineTiers_Color(int t) { MachineTiers_Load(); return t >= 0 && t < (int)tiers.size() ? tiers[t].color : GRAY; }
const MachineOverclock& MachineTiers_Overclock() { MachineTiers_Load(); return overclock; }

bool MachineTiers_Tiered(PartType t) { return FamilyOf(t) != nullptr; }
int MachineTiers_First(PartType t) { const Family* f = FamilyOf(t); return f ? f->first : -1; }
int MachineTiers_MaxStep(PartType t) { const Family* f = FamilyOf(t); return f ? f->tiers.back() - f->first : 0; }
int MachineTiers_PartTier(PartType t, int step, int baseVolts) {
    const Family* f = FamilyOf(t);
    if (f) return f->first + std::clamp(step, 0, f->tiers.back() - f->first);
    return MachineTiers_ForVolts(baseVolts);
}
std::string MachineTiers_VariantItem(PartType t, int step) {
    const Family* f = FamilyOf(t);
    const char* base = PartCostItem(t);
    if (step == 0) return base ? base : "";
    if (!f || !base) return "";
    const int tier = f->first + step;
    if (std::find(f->tiers.begin(), f->tiers.end(), tier) == f->tiers.end()) return "";
    auto it = f->names.find(tier);
    return it != f->names.end() ? it->second : std::string(base) + "_" + Lower(tiers[tier].name);
}
bool MachineTiers_ItemPart(const std::string& item, PartType& type, int& step) {
    MachineTiers_Load();
    for (const Family& f : families)
        for (int tier : f.tiers)
            if (MachineTiers_VariantItem(f.type, tier - f.first) == item) { type = f.type; step = tier - f.first; return true; }
    return false;
}
std::string MachineTiers_ItemsText() {
    MachineTiers_Load();
    std::string out;
    for (const Family& f : families) {
        const ItemDef* base = itemDB.find(PartCostItem(f.type));
        for (int tier : f.tiers) {
            if (tier == f.first) continue;
            const std::string id = MachineTiers_VariantItem(f.type, tier - f.first);
            std::string name = base ? base->name : id;
            const std::string firstName = tiers[f.first].name + " ";
            if (name.compare(0, firstName.size(), firstName) == 0) name = name.substr(firstName.size()); // "LV assembler" -> "assembler"
            name = tiers[tier].name + " " + name;
            const Color b0 = base ? base->color : GRAY, tc = tiers[tier].color; // (L7d) its base's colour toward its tier's
            const Color c = { (unsigned char)((b0.r * 3 + tc.r * 2) / 5), (unsigned char)((b0.g * 3 + tc.g * 2) / 5), (unsigned char)((b0.b * 3 + tc.b * 2) / 5), 255 };
            char buf[512];
            snprintf(buf, sizeof buf, "%s | %s | block | %d,%d,%d | %s (%d V): its tier's recipes and every one below, those overclocked.\n",
                     id.c_str(), name.c_str(), c.r, c.g, c.b, name.c_str(), tiers[tier].volts);
            out += buf;
        }
    }
    return out;
}

int MachineTiers_Steps(int machineTier, int recipeTier) {
    if (recipeTier > machineTier) return -1;
    const int n = machineTier - recipeTier;
    const int cap = MachineTiers_Overclock().cap;
    return cap > 0 ? std::min(n, cap) : n;
}
int MachineTiers_Ticks(int ticks, int steps) {
    const int s = MachineTiers_Overclock().speed;
    long long d = 1;
    for (int k = 0; k < steps && d < 1000000; ++k) d *= s;
    return (int)std::max<long long>(1, ticks / d);
}
int MachineTiers_Eu(int eu, int steps) {
    const int m = MachineTiers_Overclock().eu;
    long long v = eu;
    for (int k = 0; k < steps && v < 1000000000; ++k) v *= m;
    return (int)std::min<long long>(v, 1000000000);
}
int MachineTiers_YieldPercent(int asmTier, int itemTier) {
    MachineTiers_Load();
    const int above = itemTier - asmTier; // the item's tiers over the assembler's
    if (above >= 3) return yields[0];
    if (above == 2) return yields[1];
    if (above == 1) return yields[2];
    if (above == 0) return yields[3];
    if (above == -1) return yields[4];
    return yields[5];
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>
TEST(MachineTiers, TheDataFileGivesTiersOverclockFamiliesAndYields) {
    ASSERT_TRUE(MachineTiers_Load());
    EXPECT_EQ(MachineTiers_ByName("lv"), 0);
    EXPECT_EQ(MachineTiers_ByName("MV"), 1);
    EXPECT_EQ(MachineTiers_Volts(MachineTiers_ByName("HV")), 512);
    EXPECT_EQ(MachineTiers_ForVolts(128), 1);
    EXPECT_EQ(MachineTiers_ForVolts(8), -1) << "steam machines are under LV";
    const MachineOverclock& oc = MachineTiers_Overclock();
    EXPECT_EQ(oc.speed, 2); EXPECT_EQ(oc.eu, 3); EXPECT_EQ(oc.cap, 0) << "the user's numbers: x2 speed for x3 EU/t, no cap";
    EXPECT_TRUE(MachineTiers_Tiered(PartType::Macerator));
    EXPECT_FALSE(MachineTiers_Tiered(PartType::Campfire));
    EXPECT_EQ(MachineTiers_VariantItem(PartType::Macerator, 0), "macerator");
    EXPECT_EQ(MachineTiers_VariantItem(PartType::Macerator, 1), "macerator_mv");
    EXPECT_EQ(MachineTiers_VariantItem(PartType::LvAssembler, 2), "hv_assembler");
    EXPECT_EQ(MachineTiers_VariantItem(PartType::Macerator, 9), "");
    PartType t = PartType::Count; int step = -1;
    ASSERT_TRUE(MachineTiers_ItemPart("mv_assembler", t, step));
    EXPECT_EQ(t, PartType::LvAssembler); EXPECT_EQ(step, 1);
    EXPECT_EQ(MachineTiers_PartTier(PartType::Electrolyzer, 1, 128), 2) << "an electrolyzer starts at MV: its first variant is HV";
    EXPECT_EQ(MachineTiers_PartTier(PartType::Campfire, 3, 0), -1) << "untiered: its voltage's tier, the step ignored";
    EXPECT_EQ(MachineTiers_YieldPercent(0, 0), 50);
    EXPECT_EQ(MachineTiers_YieldPercent(0, 1), 10);
    EXPECT_EQ(MachineTiers_YieldPercent(0, 2), 1);
    EXPECT_EQ(MachineTiers_YieldPercent(0, 3), 0);
    EXPECT_EQ(MachineTiers_YieldPercent(1, 0), 75);
    EXPECT_EQ(MachineTiers_YieldPercent(3, 0), 100);
}
TEST(MachineTiers, OverclockingHalvesTheTimeAndTriplesTheEuPerStep) {
    MachineTiers_Load();
    EXPECT_EQ(MachineTiers_Steps(1, 0), 1);
    EXPECT_EQ(MachineTiers_Steps(2, 0), 2);
    EXPECT_EQ(MachineTiers_Steps(0, 1), -1) << "an MV recipe never runs in an LV machine";
    EXPECT_EQ(MachineTiers_Ticks(200, 0), 200);
    EXPECT_EQ(MachineTiers_Ticks(200, 1), 100);
    EXPECT_EQ(MachineTiers_Ticks(200, 2), 50);
    EXPECT_EQ(MachineTiers_Ticks(3, 5), 1) << "never under a tick";
    EXPECT_EQ(MachineTiers_Eu(16, 1), 48);
    EXPECT_EQ(MachineTiers_Eu(16, 2), 144);
}
// (#1427 L7b) Every variant item exists (made from the table) and has its own recipe (machine_tier_recipes.txt), which
// never takes a machine of its family: each tier has its own recipe (the user's decision).
TEST(MachineTiers, EveryVariantHasItsOwnRecipeAndNoneTakesTheTierBelow) {
    ItemDB db;
    ASSERT_TRUE(db.load("resources/data/items.txt", "resources/data/recipes.txt"));
    ASSERT_TRUE(db.loadExtraText(MachineTiers_ItemsText()));
    ASSERT_TRUE(db.loadExtra("resources/data/machine_tier_recipes.txt"));
    int variants = 0;
    for (int t = 0; t < (int)PartType::Count; ++t) {
        const PartType type = (PartType)t;
        for (int step = 1; step <= MachineTiers_MaxStep(type); ++step) {
            const std::string id = MachineTiers_VariantItem(type, step);
            ASSERT_FALSE(id.empty());
            ++variants;
            EXPECT_NE(db.find(id), nullptr) << id << ": its item";
            const Recipe* r = nullptr;
            for (const Recipe& rc : db.recipes()) if (rc.output == id) { r = &rc; break; }
            ASSERT_NE(r, nullptr) << id << ": its recipe";
            for (const auto& in : r->inputs) {
                PartType it = PartType::Count; int is = 0;
                const bool aMachine = MachineTiers_ItemPart(in.first, it, is) && it == type;
                EXPECT_FALSE(aMachine) << id << " takes " << in.first << ": a tier's own recipe, never the tier below's machine";
            }
        }
    }
    EXPECT_EQ(variants, 12) << "MV and HV of the five LV families, HV of the two MV ones";
}
#endif
