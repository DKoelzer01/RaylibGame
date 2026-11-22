#include "orefamilies.h"
#include "ores.h"
#include "utils/logger.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <set>
#include <sstream>

namespace {
std::string Trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r')) ++a;
    while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r')) --b;
    return s.substr(a, b - a);
}
std::vector<std::string> Split(const std::string& s, char c) {
    std::vector<std::string> out;
    size_t start = 0;
    for (size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == c) { out.push_back(Trim(s.substr(start, i - start))); start = i + 1; }
    return out;
}
std::string Lower(std::string s) { for (char& c : s) c = (char)std::tolower((unsigned char)c); return s; }
const char* kShapes[] = { "blob", "sheet", "pipe", "seam", "stockwork" };
// The rock strata a family can sit in (VoxelMaterial, voxel.h), the surface biomes (Biome / BiomeName, voxel.h) and the
// underground regions (#1547), lower case
const std::set<std::string> kHosts = { "rock", "granite", "basalt", "sandstone", "limestone", "marble" };
const std::set<std::string> kBiomes = { "ocean", "jungle", "savanna", "desert", "badlands", "forest", "meadow", "swamp",
                                        "taiga", "tundra", "alpine", "polar" };
const std::set<std::string> kRegions = { "none", "karst", "volcanic", "frozen", "crystal", "flooded", "fungal" };

// "name mult, name mult" or "any"
std::vector<OreFamilyBias> Biases(const std::string& s) {
    std::vector<OreFamilyBias> out;
    if (s == "any") return out;
    for (const std::string& e : Split(s, ',')) {
        std::istringstream bs(e);
        OreFamilyBias b;
        if (!(bs >> b.name)) continue;
        if (!(bs >> b.factor)) b.factor = 1.0f;
        b.name = Lower(b.name);
        out.push_back(b);
    }
    return out;
}
float Factor(const std::vector<OreFamilyBias>& list, const std::string& name) {
    const std::string n = Lower(name);
    for (const OreFamilyBias& b : list) if (b.name == n) return b.factor;
    return 1.0f;
}

std::vector<OreFamily> Load() {
    std::ifstream f("resources/data/ore_families.txt");
    if (!f) { LOG_ERROR("OreFamilies", "resources/data/ore_families.txt: not found: no families\n"); return {}; }
    std::stringstream ss;
    ss << f.rdbuf();
    return OreFamilies_Parse(ss.str());
}
const std::vector<OreFamily>& T() {
    static const std::vector<OreFamily> t = Load(); // (thread-safe: the mesher's workers ask too)
    return t;
}
} // namespace

std::vector<OreFamily> OreFamilies_Parse(const std::string& text) {
    std::vector<OreFamily> out;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;
        // id | members | hosts | depth | biomes | regions | rarity | shape | size
        const auto p = Split(line, '|');
        if (p.size() < 9) { LOG_WARN("OreFamilies", "short line '%s'\n", line.c_str()); continue; }
        OreFamily fam;
        fam.id = p[0];
        for (const std::string& m : Split(p[1], ',')) {
            std::istringstream ms(m);
            OreFamilyMember mem;
            if (!(ms >> mem.ore)) continue;
            if (!(ms >> mem.weight)) mem.weight = 1.0f;
            for (int k = 0; k < Ores_Count(); ++k) if (Ores_Get(k)->id == mem.ore) { mem.kind = k; break; }
            fam.totalWeight += std::max(0.0f, mem.weight);
            fam.members.push_back(mem);
        }
        if (p[2] != "any") { std::istringstream hs(p[2]); for (std::string h; hs >> h;) fam.hosts.push_back(Lower(h)); }
        {
            std::istringstream ds(p[3]);
            std::string a, b;
            ds >> a >> b;
            fam.minDepth = strtof(a.c_str(), nullptr);
            if (!b.empty() && b != "-") fam.maxDepth = strtof(b.c_str(), nullptr);
        }
        fam.biomes = Biases(p[4]);
        fam.regions = Biases(p[5]);
        fam.rarity = strtof(p[6].c_str(), nullptr);
        bool shapeOk = false;
        for (int s = 0; s < (int)(sizeof(kShapes) / sizeof(kShapes[0])); ++s) if (p[7] == kShapes[s]) { fam.shape = (VeinShape)s; shapeOk = true; }
        if (!shapeOk) { LOG_WARN("OreFamilies", "%s: unknown shape '%s': skipped\n", fam.id.c_str(), p[7].c_str()); continue; }
        {
            std::istringstream zs(p[8]);
            float a = 0.0f, b = 0.0f;
            zs >> a;
            if (!(zs >> b)) b = a;
            fam.minSize = a; fam.maxSize = b;
        }
        out.push_back(fam);
    }
    return out;
}

const std::vector<OreFamily>& OreFamilies_All() { return T(); }
const OreFamily* OreFamilies_Find(const std::string& id) {
    for (const OreFamily& f : T()) if (f.id == id) return &f;
    return nullptr;
}
const OreFamilyMember* OreFamilies_Member(const OreFamily& f, float r) {
    if (f.members.empty() || f.totalWeight <= 0.0f) return nullptr;
    const float pick = std::clamp(r, 0.0f, 1.0f) * f.totalWeight;
    float acc = 0.0f;
    for (const OreFamilyMember& m : f.members) { acc += std::max(0.0f, m.weight); if (pick < acc) return &m; }
    return &f.members.back();
}
float OreFamilies_Rarity(const OreFamily& f, const std::string& biome, const std::string& region) {
    return f.rarity * Factor(f.biomes, biome) * Factor(f.regions, region);
}
bool OreFamilies_Hosts(const OreFamily& f, const std::string& rock) {
    if (f.hosts.empty()) return true;
    const std::string r = Lower(rock);
    return std::find(f.hosts.begin(), f.hosts.end(), r) != f.hosts.end();
}
const char* OreFamilies_ShapeName(VeinShape s) { return kShapes[(int)s]; }

std::string OreFamilies_Validate(const std::vector<OreFamily>& fams) {
    std::ostringstream e;
    std::set<std::string> seen;
    for (const OreFamily& f : fams) {
        if (!seen.insert(f.id).second) e << f.id << ": a second family of this id\n";
        if (f.members.empty()) e << f.id << ": no members\n";
        for (const OreFamilyMember& m : f.members) {
            if (m.kind < 0) e << f.id << ": " << m.ore << " is not in ores.txt\n";
            if (m.weight <= 0.0f) e << f.id << ": " << m.ore << " weight " << m.weight << " (over 0)\n";
        }
        for (const std::string& h : f.hosts) if (!kHosts.count(h)) e << f.id << ": unknown host rock '" << h << "'\n";
        for (const OreFamilyBias& b : f.biomes) {
            if (!kBiomes.count(b.name)) e << f.id << ": unknown biome '" << b.name << "'\n";
            if (b.factor < 0.0f) e << f.id << ": " << b.name << " factor under 0\n";
        }
        for (const OreFamilyBias& b : f.regions) {
            if (!kRegions.count(b.name)) e << f.id << ": unknown region '" << b.name << "'\n";
            if (b.factor < 0.0f) e << f.id << ": " << b.name << " factor under 0\n";
        }
        if (f.minDepth < 0.0f || f.minDepth >= f.maxDepth) e << f.id << ": depth " << f.minDepth << " to " << f.maxDepth << "\n";
        if (f.rarity <= 0.0f) e << f.id << ": rarity not over 0\n";
        if (f.minSize <= 0.0f || f.minSize > f.maxSize) e << f.id << ": size " << f.minSize << " to " << f.maxSize << "\n";
    }
    return e.str();
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"

TEST(OreFamilies, ParsesMembersHostsDepthBiomesRegionsShapeAndSize) {
    const auto f = OreFamilies_Parse(
        "# comment\n"
        "lead_zinc  | iron 5, copper 3, tin 0.5 | limestone marble | 30 140 | desert 2, badlands 1.5 | volcanic 3 | 0.8 | stockwork | 12 24\n"
        "coal_beds  | coal 1                    | any              | 4 -    | any                     | any        | 1   | seam      | 3\n"
        "broken     | iron 1                    | any              | 4 10   | any                     | any        | 1   | spiral    | 3\n");
    ASSERT_EQ(f.size(), 2u) << "an unknown shape is skipped";
    const OreFamily& g = f[0];
    ASSERT_EQ(g.members.size(), 3u);
    EXPECT_EQ(g.members[0].ore, "iron"); EXPECT_FLOAT_EQ(g.members[0].weight, 5.0f); EXPECT_EQ(g.members[0].kind, 0);
    EXPECT_FLOAT_EQ(g.members[2].weight, 0.5f) << "float weights";
    EXPECT_FLOAT_EQ(g.totalWeight, 8.5f);
    EXPECT_EQ(g.hosts, (std::vector<std::string>{ "limestone", "marble" }));
    EXPECT_TRUE(OreFamilies_Hosts(g, "Marble")) << "case-insensitive";
    EXPECT_FALSE(OreFamilies_Hosts(g, "granite")) << "hosts are hard eligibility";
    EXPECT_TRUE(OreFamilies_Hosts(f[1], "basalt")) << "any";
    EXPECT_FLOAT_EQ(g.minDepth, 30.0f); EXPECT_FLOAT_EQ(g.maxDepth, 140.0f);
    EXPECT_FLOAT_EQ(OreFamilies_Rarity(g, "Desert"), 1.6f) << "as BiomeName prints it";
    EXPECT_FLOAT_EQ(OreFamilies_Rarity(g, "forest"), 0.8f) << "unlisted biomes: x1";
    EXPECT_FLOAT_EQ(OreFamilies_Rarity(g, "desert", "volcanic"), 4.8f) << "biome x region";
    EXPECT_EQ(g.shape, VeinShape::Stockwork);
    EXPECT_FLOAT_EQ(g.minSize, 12.0f); EXPECT_FLOAT_EQ(g.maxSize, 24.0f);
    EXPECT_FLOAT_EQ(f[1].minSize, 3.0f); EXPECT_FLOAT_EQ(f[1].maxSize, 3.0f) << "one number: a fixed size";
    EXPECT_GT(f[1].maxDepth, 1e8f) << "'-': no lower limit";
    EXPECT_EQ(OreFamilies_Validate(f), "");
    // Members by weight: 5 / 8.5 iron, 3 / 8.5 copper, 0.5 / 8.5 tin
    int n[3] = {};
    for (int i = 0; i < 1700; ++i) { const OreFamilyMember* m = OreFamilies_Member(g, (i + 0.5f) / 1700.0f); n[m == &g.members[0] ? 0 : m == &g.members[1] ? 1 : 2]++; }
    EXPECT_EQ(n[0], 1000); EXPECT_EQ(n[1], 600); EXPECT_EQ(n[2], 100);
}

TEST(OreFamilies, ValidateNamesEachProblem) {
    const auto f = OreFamilies_Parse(
        "a | unobtainium 1, iron 0 | lava | 50 20 | moon 2 | sky 2 | 0 | blob | 9 3\n"
        "a | iron 1                | any  | 0 10  | any    | any   | 1 | blob | 4\n");
    const std::string e = OreFamilies_Validate(f);
    for (const char* want : { "unobtainium is not in ores.txt", "weight 0", "unknown host rock 'lava'", "unknown biome 'moon'",
                              "unknown region 'sky'", "depth 50 to 20", "rarity not over 0", "size 9 to 3", "a second family" })
        EXPECT_NE(e.find(want), std::string::npos) << want << " in:\n" << e;
}

// The game's file: every family checks out, and every ore in ores.txt with no deposit or pocket share of its own is in
// at least one family (else nothing would ever place it).
TEST(OreFamilies, TheGamesFamiliesAreValidAndPlaceEveryNewOre) {
    const auto& f = OreFamilies_All();
    ASSERT_FALSE(f.empty());
    EXPECT_EQ(OreFamilies_Validate(f), "");
    for (int k = 0; k < Ores_Count(); ++k) {
        const OreDef* o = Ores_Get(k);
        if (o->deposit > 0.0f || o->pocket > 0.0f) continue;
        bool found = false;
        for (const OreFamily& fam : f) for (const OreFamilyMember& m : fam.members) found = found || m.kind == k;
        EXPECT_TRUE(found) << o->id << ": in no family";
    }
}
#endif
