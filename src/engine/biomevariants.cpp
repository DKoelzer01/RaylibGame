#include "biomevariants.h"
#include <cctype>
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
        } else {
            return fail(lineNo, "unknown key '" + key + "'");
        }
    }
    for (int b = 0; b < (int)Biome::Count; ++b) {
        if (t.byType[b].empty()) return fail(lineNo, "no variants for " + Lower(BiomeName((Biome)b)));
        if (t.defs[t.byType[b][0]].gated()) return fail(lineNo, "the base (first) variant of " + Lower(BiomeName((Biome)b)) + " has a gate: it is the fallback");
    }
    Bands(t);
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
}
#endif
