#include "biomevariants.h"
#include "envflag.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace {
    std::string Lower(std::string s) { for (char& c : s) c = (char)std::tolower((unsigned char)c); return s; }

    // A gate bound: a number, or "-" for open.
    bool Bound(const std::string& tok, float open, float& out) {
        if (tok == "-") { out = open; return true; }
        char* end = nullptr;
        out = strtof(tok.c_str(), &end);
        return end && *end == 0 && end != tok.c_str();
    }

    // (#2032) Each id's tints: its own, else its type's base's, else 1.
    void Tints(BiomeVariantTable& t) {
        t.grassTint.assign(t.defs.size(), Vector3{ 1.0f, 1.0f, 1.0f });
        t.leafTint.assign(t.defs.size(), Vector3{ 1.0f, 1.0f, 1.0f });
        for (size_t i = 0; i < t.defs.size(); ++i) {
            const int b = t.base(t.defs[i].type);
            const BiomeVariantDef& own = t.defs[i];
            const BiomeVariantDef* g = own.grassTint[0] >= 0.0f ? &own : b >= 0 && t.defs[b].grassTint[0] >= 0.0f ? &t.defs[b] : nullptr;
            const BiomeVariantDef* l = own.leafTint[0] >= 0.0f ? &own : b >= 0 && t.defs[b].leafTint[0] >= 0.0f ? &t.defs[b] : nullptr;
            if (g) t.grassTint[i] = { g->grassTint[0], g->grassTint[1], g->grassTint[2] };
            if (l) t.leafTint[i] = { l->leafTint[0], l->leafTint[1], l->leafTint[2] };
        }
    }

    void Bands(BiomeVariantTable& t) {
        t.bandLo.assign(t.defs.size(), 0.0f);
        t.bandHi.assign(t.defs.size(), 1.0f);
        for (int b = 0; b < (int)Biome::Count; ++b) {
            float total = 0.0f;
            for (int id : t.byType[b]) total += t.defs[id].weight;
            float at = 0.0f;
            for (int id : t.byType[b]) { t.bandLo[id] = at / total; at += t.defs[id].weight; t.bandHi[id] = at / total; }
            if (!t.byType[b].empty()) t.bandHi[t.byType[b].back()] = 1.0f; // no gap at the top from rounding
        }
    }
}

int BiomeVariantTable::find(const std::string& key) const {
    for (int i = 0; i < (int)defs.size(); ++i) if (defs[i].key == key) return i;
    return -1;
}

bool BiomeVariants_Parse(const std::string& text, BiomeVariantTable& out, std::string* err) {
    BiomeVariantTable t;
    auto fail = [&](int line, const std::string& why) { if (err) *err = "biomes.txt:" + std::to_string(line) + ": " + why; return false; };
    std::istringstream in(text);
    std::string raw;
    int lineNo = 0;
    while (std::getline(in, raw)) {
        ++lineNo;
        if (const size_t hash = raw.find('#'); hash != std::string::npos) raw.erase(hash);
        while (!raw.empty() && (raw.back() == '\r' || raw.back() == ' ' || raw.back() == '\t')) raw.pop_back();
        if (raw.empty()) continue;
        const bool indented = raw[0] == ' ' || raw[0] == '\t';
        std::istringstream ls(raw);
        std::string key;
        ls >> key;
        if (!indented) {
            if (key != "variant") return fail(lineNo, "expected 'variant <type>.<name>', got '" + key + "'");
            std::string id;
            ls >> id;
            const size_t dot = id.find('.');
            if (dot == std::string::npos || dot == 0 || dot + 1 == id.size()) return fail(lineNo, "variant '" + id + "' is not <type>.<name>");
            int type = -1;
            for (int b = 0; b < (int)Biome::Count; ++b) if (Lower(BiomeName((Biome)b)) == id.substr(0, dot)) type = b;
            if (type < 0) return fail(lineNo, "unknown biome type '" + id.substr(0, dot) + "'");
            if (t.find(id) >= 0) return fail(lineNo, "variant '" + id + "' listed twice");
            BiomeVariantDef d;
            d.key = id; d.name = id.substr(dot + 1); d.type = (Biome)type;
            t.byType[type].push_back((int)t.defs.size());
            t.defs.push_back(d);
            continue;
        }
        if (t.defs.empty()) return fail(lineNo, "'" + key + "' before any variant");
        BiomeVariantDef& d = t.defs.back();
        if (key == "name") {
            std::string rest;
            std::getline(ls, rest);
            const size_t s = rest.find_first_not_of(" \t");
            if (s == std::string::npos) return fail(lineNo, "empty name");
            d.name = rest.substr(s);
        } else if (key == "weight") {
            std::string v;
            ls >> v;
            if (!Bound(v, 0.0f, d.weight) || v == "-" || !(d.weight > 0.0f)) return fail(lineNo, "weight must be a number > 0");
        } else if (key == "elev" || key == "temp" || key == "moist") {
            std::string lo, hi, extra;
            ls >> lo >> hi >> extra;
            float* mn = key == "elev" ? &d.elevMin : key == "temp" ? &d.tempMin : &d.moistMin;
            float* mx = key == "elev" ? &d.elevMax : key == "temp" ? &d.tempMax : &d.moistMax;
            if (!extra.empty() || !Bound(lo, -1e9f, *mn) || !Bound(hi, 1e9f, *mx) || *mn > *mx) return fail(lineNo, key + " takes '<min> <max>' (numbers or -)");
        } else if (key == "flora") {
            std::string v[4], extra;
            ls >> v[0] >> v[1] >> v[2] >> v[3] >> extra;
            float* dst[4] = { &d.trees, &d.bushes, &d.stones, &d.size };
            for (int i = 0; i < 4; ++i)
                if (!Bound(v[i], 1.0f, *dst[i]) || *dst[i] < 0.0f) return fail(lineNo, "flora takes '<trees> <bushes> <stones> <size>' (numbers >= 0 or - for 1)");
            if (!extra.empty() || !(d.size > 0.0f)) return fail(lineNo, "flora takes '<trees> <bushes> <stones> <size>' (size > 0)");
        } else if (key == "sizes") { // (#1995 phase 3)
            std::string v[3], extra;
            ls >> v[0] >> v[1] >> v[2] >> extra;
            float* dst[3] = { &d.bushSize, &d.stoneSize, &d.grassSize };
            for (int i = 0; i < 3; ++i)
                if (!Bound(v[i], 0.0f, *dst[i]) || *dst[i] < -2.0f || *dst[i] > 2.0f) return fail(lineNo, "sizes takes '<bushes> <stones> <grass>' (numbers -2..2, or - for 0)");
            if (!extra.empty()) return fail(lineNo, "sizes takes '<bushes> <stones> <grass>' (three values)");
        } else if (key == "shrubs") { // (#1995 phase 3)
            std::string v[BiomeVariantDef::kShrubKinds], extra;
            ls >> v[0] >> v[1] >> v[2] >> extra;
            float total = 0.0f;
            for (int i = 0; i < BiomeVariantDef::kShrubKinds; ++i) {
                if (!Bound(v[i], 0.0f, d.shrubs[i]) || v[i] == "-" || d.shrubs[i] < 0.0f) return fail(lineNo, "shrubs takes '<round> <mat> <upright>' weights (numbers >= 0)");
                total += d.shrubs[i];
            }
            if (!extra.empty() || !(total > 0.0f)) return fail(lineNo, "shrubs takes '<round> <mat> <upright>' weights, at least one > 0");
        } else if (key == "rubber") { // (#2012)
            std::string extra;
            if (!(ls >> d.rubber) || (ls >> extra) || d.rubber < 0 || d.rubber > 100) return fail(lineNo, "rubber takes one whole percent 0..100");
        } else if (key == "tint") { // (#2032)
            std::string v[6], extra;
            ls >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5] >> extra;
            float* dst[6] = { &d.grassTint[0], &d.grassTint[1], &d.grassTint[2], &d.leafTint[0], &d.leafTint[1], &d.leafTint[2] };
            // Each group is three numbers 0..4, or "- - -" for the type's base's (a base: 1).
            const char* kUse = "tint takes '<grass r g b> <leaf r g b>' (numbers 0..4, or - - - to keep the base's)";
            for (int g = 0; g < 2; ++g) {
                const bool keep = v[3 * g] == "-";
                for (int i = 3 * g; i < 3 * g + 3; ++i) {
                    if (keep) { if (v[i] != "-") return fail(lineNo, kUse); continue; }
                    if (v[i] == "-" || !Bound(v[i], 1.0f, *dst[i]) || *dst[i] < 0.0f || *dst[i] > 4.0f) return fail(lineNo, kUse);
                }
            }
            if (!extra.empty()) return fail(lineNo, kUse);
        } else if (key == "cover") {
            std::string v, extra;
            ls >> v >> extra;
            if (!extra.empty() || v == "-" || !Bound(v, 0.0f, d.cover) || d.cover < 0.0f || d.cover > 1.0f) return fail(lineNo, "cover takes one number 0..1");
        } else if (key == "species") {
            std::string sp, w;
            int n = 0;
            while (ls >> sp) {
                const int k = WoodSpeciesByName(sp.c_str());
                if (k < 0) {
                    std::string known;
                    for (const char* name : kWoodSpeciesNames) known += std::string(known.empty() ? "" : " ") + name;
                    return fail(lineNo, "unknown species '" + sp + "' (" + known + ")");
                }
                if (!(ls >> w) || !Bound(w, 0.0f, d.species[k]) || w == "-" || d.species[k] < 0.0f) return fail(lineNo, "species takes '<name> <weight>' pairs");
                ++n;
            }
            if (n == 0 || !d.hasSpecies()) return fail(lineNo, "species needs at least one '<name> <weight>' with weight > 0");
        } else {
            return fail(lineNo, "unknown key '" + key + "'");
        }
    }
    for (int b = 0; b < (int)Biome::Count; ++b) {
        if (t.byType[b].empty()) return fail(lineNo, "no variants for " + Lower(BiomeName((Biome)b)));
        if (t.defs[t.byType[b][0]].gated()) return fail(lineNo, "the base (first) variant of " + Lower(BiomeName((Biome)b)) + " has a gate: it is the fallback");
    }
    if ((int)t.defs.size() > kMaxBiomeVariants) return fail(lineNo, "more than " + std::to_string(kMaxBiomeVariants) + " variants (biometint.glsl's arrays)");
    Bands(t);
    Tints(t);
    out = std::move(t);
    return true;
}

void BiomeVariants_Default(BiomeVariantTable& out) {
    out = BiomeVariantTable{};
    for (int b = 0; b < (int)Biome::Count; ++b) {
        BiomeVariantDef d;
        d.key = Lower(BiomeName((Biome)b)) + ".base"; d.name = BiomeName((Biome)b); d.type = (Biome)b;
        out.byType[b].push_back((int)out.defs.size());
        out.defs.push_back(d);
    }
    Bands(out);
    Tints(out);
}

const BiomeVariantTable& BiomeVariants() {
    static const BiomeVariantTable table = [] {
        BiomeVariantTable t;
        std::ifstream f("resources/data/biomes.txt", std::ios::binary);
        std::string err;
        if (!f) { fprintf(stderr, "biomes.txt: missing; one variant per biome\n"); BiomeVariants_Default(t); return t; }
        std::stringstream ss;
        ss << f.rdbuf();
        if (!BiomeVariants_Parse(ss.str(), t, &err)) { fprintf(stderr, "%s; one variant per biome\n", err.c_str()); BiomeVariants_Default(t); }
        return t;
    }();
    return table;
}

const char* BiomeVariantName(int id) {
    const BiomeVariantTable& t = BiomeVariants();
    return id >= 0 && id < (int)t.defs.size() ? t.defs[id].name.c_str() : "";
}

bool BiomeVariants_RubberRoll(int percent, uint64_t plantId) {
    if (percent <= 0) return false;
    if (percent >= 100) return true;
    uint64_t h = plantId * 0x9E3779B97F4A7C15ull; // a mix of the id: plant ids are spatial hashes, neighbours alike in their low bits
    h ^= h >> 31;
    return (int)(h % 100) < percent;
}

SmallFloraSizes SmallFloraSizesOf(int variant, Biome fallback) {
    const BiomeVariantTable& t = BiomeVariants();
    const BiomeVariantDef* d = variant >= 0 && variant < (int)t.defs.size() ? &t.defs[variant] : nullptr;
    // The types' own (GUESS, #1995 phase 3): the jungle's undergrowth is rank, the savanna's and meadow's grass tall, the
    // tundra's scrub and grass stunted, the mountains and badlands strewn with rocks and boulders, the desert with stones.
    SmallFloraSizes s;
    switch (d ? d->type : fallback) {
        case Biome::Jungle:   s = { 0.6f, 0.0f, 0.4f }; break;
        case Biome::Savanna:  s = { 0.0f, 0.2f, 0.9f }; break;
        case Biome::Desert:   s = { -0.5f, 0.3f, -0.5f }; break;
        case Biome::Badlands: s = { -0.3f, 0.6f, -0.3f }; break;
        case Biome::Meadow:   s = { 0.0f, -0.3f, 0.4f }; break;
        case Biome::Swamp:    s = { 0.3f, -0.5f, 0.6f }; break;
        case Biome::Taiga:    s = { -0.2f, 0.3f, -0.2f }; break;
        case Biome::Tundra:   s = { -0.7f, 0.3f, -0.7f }; break;
        case Biome::Alpine:   s = { -0.5f, 0.7f, -0.5f }; break;
        default: break; // forest: the mid mix
    }
    if (d) { s.bush += d->bushSize; s.stone += d->stoneSize; s.grass += d->grassSize; }
    return s;
}

int SmallFloraSizeClass(float shift, float patch, float u) {
    // Base shares 30 / 42 / 21 / 7 %, tilted by exp(t * (class - 1.5)): t > 0 moves weight to the big classes.
    static const float kBase[4] = { 0.30f, 0.42f, 0.21f, 0.07f };
    const float t = std::clamp(shift + 1.4f * (patch - 0.5f), -2.5f, 2.5f);
    float w[4], total = 0.0f;
    for (int k = 0; k < 4; ++k) total += (w[k] = kBase[k] * expf(t * ((float)k - 1.5f)));
    float at = u * total;
    for (int k = 0; k < 3; ++k) if ((at -= w[k]) < 0.0f) return k;
    return 3;
}

float BiomeClimateCode(float tempC, float moist, int variant) {
    const float step = floorf(std::clamp((tempC + 30.0f) / 70.0f, 0.0f, 1.0f) * 62.0f);
    const float id = variant >= 0 && variant < kMaxBiomeVariants && BiomeTintOn() ? (float)variant : 0.0f;
    return 2.0f + step + 64.0f * id + std::clamp(moist, 0.0f, 0.999f);
}

bool BiomeTintOn() { static const bool on = !(getenv("RAYLIBGAME_BIOME_TINT") && !EnvFlag("RAYLIBGAME_BIOME_TINT")); return on; }

// GUESS (#2032): the climate gradient keeps 40% of its pull inside a biome; the biome's tint carries the rest.
float BiomeClimateWeight() { return BiomeTintOn() ? 0.4f : 1.0f; }

void BiomeTint_Bind(const Shader& s) {
    // Locations per program (shaders hot-reload); the arrays are the same every frame, but a reload loses them.
    struct Locs { unsigned int prog = 0; int grass = -1, leaf = -1, weight = -1; };
    static Locs cache[8];
    Locs* L = nullptr;
    for (Locs& c : cache) if (c.prog == s.id) { L = &c; break; }
    if (!L) {
        static int next = 0;
        L = &cache[next++ % 8];
        L->prog = s.id;
        L->grass = GetShaderLocation(s, "biomeGrassTint");
        L->leaf = GetShaderLocation(s, "biomeLeafTint");
        L->weight = GetShaderLocation(s, "biomeClimateWeight");
    }
    static Vector3 grass[kMaxBiomeVariants], leaf[kMaxBiomeVariants];
    static const bool filled = [] {
        const BiomeVariantTable& t = BiomeVariants();
        for (int i = 0; i < kMaxBiomeVariants; ++i) {
            const bool on = BiomeTintOn() && i < (int)t.defs.size();
            grass[i] = on ? t.grassTint[i] : Vector3{ 1.0f, 1.0f, 1.0f };
            leaf[i] = on ? t.leafTint[i] : Vector3{ 1.0f, 1.0f, 1.0f };
        }
        return true;
    }();
    (void)filled;
    const float w = BiomeClimateWeight();
    if (L->grass >= 0) SetShaderValueV(s, L->grass, grass, SHADER_UNIFORM_VEC3, kMaxBiomeVariants);
    if (L->leaf >= 0) SetShaderValueV(s, L->leaf, leaf, SHADER_UNIFORM_VEC3, kMaxBiomeVariants);
    if (L->weight >= 0) SetShaderValue(s, L->weight, &w, SHADER_UNIFORM_FLOAT);
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>
namespace {
    std::string AllBases(const std::string& except = "") {
        std::string s;
        for (int b = 0; b < (int)Biome::Count; ++b) {
            const std::string k = Lower(BiomeName((Biome)b));
            if (k != except) s += "variant " + k + ".base\n";
        }
        return s;
    }
}

TEST(BiomeVariants, ParsesBlocksBandsAndGates) {
    BiomeVariantTable t;
    std::string err;
    const std::string text = AllBases("forest") + "# a comment\nvariant forest.mixed\n  name  Mixed woodland  # trailing\n  weight 3\n"
                             "variant forest.birch\r\n\tweight 1\n  elev - 400\n  temp 5 -\n";
    ASSERT_TRUE(BiomeVariants_Parse(text, t, &err)) << err;
    const int mixed = t.find("forest.mixed"), birch = t.find("forest.birch");
    ASSERT_GE(mixed, 0); ASSERT_GE(birch, 0);
    EXPECT_EQ(t.base(Biome::Forest), mixed);
    EXPECT_EQ(t.defs[mixed].name, "Mixed woodland");
    EXPECT_EQ(t.defs[birch].name, "birch"); // no name line: the id
    EXPECT_FLOAT_EQ(t.bandLo[mixed], 0.0f); EXPECT_FLOAT_EQ(t.bandHi[mixed], 0.75f);
    EXPECT_FLOAT_EQ(t.bandLo[birch], 0.75f); EXPECT_FLOAT_EQ(t.bandHi[birch], 1.0f);
    EXPECT_TRUE(t.defs[birch].gated());
    EXPECT_TRUE(t.defs[birch].admits(100.0f, 10.0f, 0.5f));
    EXPECT_FALSE(t.defs[birch].admits(500.0f, 10.0f, 0.5f));
    EXPECT_FALSE(t.defs[birch].admits(100.0f, 0.0f, 0.5f));
    EXPECT_FALSE(t.defs[mixed].gated());
}

// Flora keys (#1978): multipliers default to 1, "-" keeps 1, species by name; bad values are errors.
TEST(BiomeVariants, ParsesFlora) {
    BiomeVariantTable t;
    std::string err;
    ASSERT_TRUE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  flora 1.5 - 0.5 1.2\n  cover 0.3\n  species birch 3 pine 1\n", t, &err)) << err;
    const BiomeVariantDef& d = t.defs[t.find("forest.x")];
    EXPECT_FLOAT_EQ(d.trees, 1.5f); EXPECT_FLOAT_EQ(d.bushes, 1.0f); EXPECT_FLOAT_EQ(d.stones, 0.5f); EXPECT_FLOAT_EQ(d.size, 1.2f);
    EXPECT_FLOAT_EQ(d.cover, 0.3f);
    EXPECT_FLOAT_EQ(d.species[1], 3.0f); EXPECT_FLOAT_EQ(d.species[3], 1.0f); EXPECT_FLOAT_EQ(d.species[0], 0.0f);
    const BiomeVariantDef& base = t.defs[t.base(Biome::Forest)];
    EXPECT_FLOAT_EQ(base.trees, 1.0f); EXPECT_LT(base.cover, 0.0f); EXPECT_FALSE(base.hasSpecies());
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  species fungus 2\n", t, &err));
    EXPECT_NE(err.find("fungus"), std::string::npos) << err;
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  species oak\n", t, &err));
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  flora 1 1 1\n", t, &err));
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  flora 1 1 1 0\n", t, &err));
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  cover 2\n", t, &err));
}

TEST(BiomeVariants, RefusesBadFiles) {
    BiomeVariantTable t;
    std::string err;
    EXPECT_FALSE(BiomeVariants_Parse(AllBases("swamp"), t, &err)); // a type without variants
    EXPECT_NE(err.find("swamp"), std::string::npos) << err;
    EXPECT_FALSE(BiomeVariants_Parse(AllBases("swamp") + "variant swamp.mangrove\n  elev - 20\n", t, &err)); // gated base
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  wieght 2\n", t, &err)); // a typo is an error, not a default
    EXPECT_NE(err.find("wieght"), std::string::npos) << err;
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant jungel.x\n", t, &err));
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  weight 0\n", t, &err));
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  temp 10 5\n", t, &err));
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant forest.base\n", t, &err)); // twice
    EXPECT_FALSE(BiomeVariants_Parse("  weight 2\n" + AllBases(), t, &err));
}

// The shipped file parses, every type has its base first, and it holds the plan's 34 land variants + the ocean (#1964).
TEST(BiomeVariants, ShippedFileLoads) {
    std::ifstream f("resources/data/biomes.txt", std::ios::binary);
    ASSERT_TRUE(f) << "resources/data/biomes.txt (tests run from the repo root)";
    std::stringstream ss;
    ss << f.rdbuf();
    BiomeVariantTable t;
    std::string err;
    ASSERT_TRUE(BiomeVariants_Parse(ss.str(), t, &err)) << err;
    EXPECT_EQ((int)t.defs.size(), 35);
    EXPECT_EQ((int)t.byType[(int)Biome::Ocean].size(), 1);
    for (const char* k : { "forest.old_growth", "swamp.mangrove", "savanna.acacia" }) { // (#2012) the user's ruling on #1964
        ASSERT_GE(t.find(k), 0) << k;
        EXPECT_EQ(t.defs[t.find(k)].rubber, 33) << k;
    }
    EXPECT_EQ(t.defs[t.find("forest.mixed")].rubber, 0);
}

// (#2012) The rubber key: a whole percent 0..100, default 0; the roll gives about that share of trees, the same tree
// always the same answer.
TEST(BiomeVariants, RubberPercentAndRoll) {
    BiomeVariantTable t;
    std::string err;
    ASSERT_TRUE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  rubber 33\n", t, &err)) << err;
    EXPECT_EQ(t.defs[t.find("forest.x")].rubber, 33);
    EXPECT_EQ(t.defs[t.base(Biome::Forest)].rubber, 0);
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  rubber 120\n", t, &err));
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  rubber 33 5\n", t, &err));
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  rubber lots\n", t, &err));
    int yes = 0;
    for (uint64_t id = 0; id < 3000; ++id) {
        const uint64_t plant = id * 4096 + 7; // spaced ids, as spatial hashes are
        const bool r = BiomeVariants_RubberRoll(33, plant);
        EXPECT_EQ(r, BiomeVariants_RubberRoll(33, plant));
        yes += r;
    }
    EXPECT_NEAR(yes, 990, 90) << "about a third of 3000 trees";
    EXPECT_FALSE(BiomeVariants_RubberRoll(0, 5));
    EXPECT_TRUE(BiomeVariants_RubberRoll(100, 5));
}

// (#2032) Tints: a base's is its type's; a variant takes its own group, or its base's for "- - -", or 1; bad lines fail.
TEST(BiomeVariants, Tints) {
    BiomeVariantTable t;
    std::string err;
    ASSERT_TRUE(BiomeVariants_Parse(AllBases("savanna") + "variant savanna.base\n  tint 2 1.5 0.5  1.5 1.2 0.6\n"
                                    "variant savanna.x\n  tint - - -  0.8 0.9 1\nvariant savanna.y\n", t, &err)) << err;
    const int b = t.find("savanna.base"), x = t.find("savanna.x"), y = t.find("savanna.y"), f = t.find("forest.base");
    EXPECT_FLOAT_EQ(t.grassTint[b].x, 2.0f); EXPECT_FLOAT_EQ(t.leafTint[b].z, 0.6f);
    EXPECT_FLOAT_EQ(t.grassTint[x].y, 1.5f) << "grass from the base"; EXPECT_FLOAT_EQ(t.leafTint[x].x, 0.8f) << "its own leaves";
    EXPECT_FLOAT_EQ(t.grassTint[y].z, 0.5f); EXPECT_FLOAT_EQ(t.leafTint[y].y, 1.2f);
    EXPECT_FLOAT_EQ(t.grassTint[f].x, 1.0f); EXPECT_FLOAT_EQ(t.leafTint[f].y, 1.0f);
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  tint 1 1 1\n", t, &err));
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  tint 1 - 1  1 1 1\n", t, &err));
    EXPECT_FALSE(BiomeVariants_Parse(AllBases() + "variant forest.x\n  tint 1 1 9  1 1 1\n", t, &err));
}

// (#2032) The climate code the shaders decode (biometint.glsl): the step, the variant id and the moisture survive.
TEST(BiomeVariants, ClimateCodeRoundTrip) {
    if (!BiomeTintOn()) GTEST_SKIP() << "RAYLIBGAME_BIOME_TINT=0";
    for (int id : { -1, 0, 7, 63 })
        for (float temp : { -30.0f, 4.3f, 39.9f })
            for (float moist : { 0.0f, 0.37f, 0.999f }) {
                const float code = BiomeClimateCode(temp, moist, id);
                const float n = floorf(code) - 2.0f, gotId = floorf(n / 64.0f), step = n - 64.0f * gotId;
                EXPECT_EQ((int)gotId, std::max(id, 0)) << code;
                EXPECT_EQ(step, floorf(std::clamp((temp + 30.0f) / 70.0f, 0.0f, 1.0f) * 62.0f)) << code;
                EXPECT_NEAR(code - floorf(code), moist, 1e-3f) << code; // float precision at codes up to ~4100
            }
}
#endif
