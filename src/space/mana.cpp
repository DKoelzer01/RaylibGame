#include "mana.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

namespace {
std::vector<ManaKind> kinds;
std::map<std::string, float> rates;
bool loaded = false;

std::string Trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}

// A kind line keeps the names of its parts until every kind is read (a compound may be listed before its bases).
struct Pending { ManaKind k; std::string a, b; int line = 0; };
}

bool Mana_LoadText(const std::string& text, std::string& error) {
    std::vector<Pending> pending;
    std::map<std::string, float> newRates;
    std::istringstream in(text);
    int lineNo = 0;
    auto fail = [&](const std::string& why) { error = "mana.txt line " + std::to_string(lineNo) + ": " + why; return false; };
    for (std::string line; std::getline(in, line);) {
        ++lineNo;
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("rate ", 0) == 0) {
            std::istringstream ls(line.substr(5));
            std::string name; float v = 0.0f;
            if (!(ls >> name >> v) || v <= 0.0f || v > 1.0f) return fail("a rate is 'rate <name> <0..1>'");
            newRates[name] = v;
            continue;
        }
        std::vector<std::string> fields;
        { std::istringstream fs(line); for (std::string f; std::getline(fs, f, '|');) fields.push_back(Trim(f)); }
        if (fields.size() < 2 || fields[0].rfind("kind ", 0) != 0) return fail("not a kind line");
        Pending p;
        p.line = lineNo;
        p.k.id = Trim(fields[0].substr(5));
        p.k.name = fields[1];
        if (p.k.id.empty()) return fail("a kind without an id");
        bool hasOp = false;
        for (size_t f = 2; f < fields.size(); ++f) {
            std::istringstream fs(fields[f]);
            std::string key; fs >> key;
            std::string rest; std::getline(fs, rest); rest = Trim(rest);
            auto op = [&](ManaOp o) { if (hasOp) return false; hasOp = true; p.k.op = o; return true; };
            if (key == "latin") p.k.latin = rest;
            else if (key == "base") { if (!op(ManaOp::Base)) return fail("more than one of base/combine/invert/cancel"); }
            else if (key == "combine") {
                if (!op(ManaOp::Combine)) return fail("more than one of base/combine/invert/cancel");
                std::istringstream rs(rest);
                if (!(rs >> p.a >> p.b)) return fail("combine needs two kinds");
            } else if (key == "invert") {
                if (!op(ManaOp::Invert)) return fail("more than one of base/combine/invert/cancel");
                std::istringstream rs(rest);
                if (!(rs >> p.a)) return fail("invert needs a kind");
            } else if (key == "cancel") { if (!op(ManaOp::Cancel)) return fail("more than one of base/combine/invert/cancel"); }
            else if (key == "status") p.k.status = rest;
            else if (key == "colour") p.k.colour = (uint32_t)strtoul(rest.c_str(), nullptr, 16);
            else if (key == "slice") p.k.slice = atoi(rest.c_str()) == 1;
            else return fail("unknown field '" + key + "'");
        }
        if (!hasOp) return fail("a kind needs one of base, combine, invert, cancel");
        pending.push_back(std::move(p));
    }
    std::vector<ManaKind> out;
    for (const Pending& p : pending) out.push_back(p.k);
    auto find = [&](const std::string& id) { for (size_t i = 0; i < out.size(); ++i) if (out[i].id == id) return (int)i; return -1; };
    for (size_t i = 0; i < pending.size(); ++i) {
        lineNo = pending[i].line;
        if (!pending[i].a.empty() && (out[i].a = find(pending[i].a)) < 0) return fail("no kind '" + pending[i].a + "'");
        if (!pending[i].b.empty() && (out[i].b = find(pending[i].b)) < 0) return fail("no kind '" + pending[i].b + "'");
    }
    kinds = std::move(out);
    rates = std::move(newRates);
    loaded = true;
    return true;
}

bool Mana_Load(const char* path) {
    std::ifstream f(path);
    if (!f) { printf("mana: cannot read %s\n", path); return false; }
    std::stringstream ss; ss << f.rdbuf();
    std::string error;
    if (!Mana_LoadText(ss.str(), error)) { printf("mana: %s\n", error.c_str()); return false; }
    const std::string problems = Mana_Validate();
    if (!problems.empty()) printf("mana: %s", problems.c_str());
    return true;
}

bool Mana_Loaded() { return loaded; }
int Mana_Count() { return (int)kinds.size(); }
const ManaKind& Mana_Kind(int k) { static const ManaKind none; return k >= 0 && k < (int)kinds.size() ? kinds[k] : none; }
int Mana_Find(const std::string& id) {
    for (size_t i = 0; i < kinds.size(); ++i) if (kinds[i].id == id) return (int)i;
    return -1;
}

std::string Mana_Validate() {
    std::string out;
    std::set<std::string> ids;
    int bases = 0, cancels = 0;
    std::set<std::pair<int, int>> pairs;
    std::map<int, int> inverted; // kind -> how many kinds invert it
    for (size_t i = 0; i < kinds.size(); ++i) {
        const ManaKind& k = kinds[i];
        if (!ids.insert(k.id).second) out += "kind " + k.id + ": the id is used twice\n";
        switch (k.op) {
            case ManaOp::Base: ++bases; break;
            case ManaOp::Cancel: ++cancels; break;
            case ManaOp::Combine: {
                const ManaKind& a = Mana_Kind(k.a); const ManaKind& b = Mana_Kind(k.b);
                if (k.a == k.b) out += "kind " + k.id + ": combines a kind with itself\n";
                else if (a.op != ManaOp::Base || b.op != ManaOp::Base) out += "kind " + k.id + ": combines a kind that is not a base (no compound of compounds in v1)\n";
                else if (!pairs.insert({ std::min(k.a, k.b), std::max(k.a, k.b) }).second) out += "kind " + k.id + ": that pair already makes another kind\n";
                break;
            }
            case ManaOp::Invert: {
                const ManaKind& a = Mana_Kind(k.a);
                if (a.op != ManaOp::Base && a.op != ManaOp::Combine) out += "kind " + k.id + ": inverts " + a.id + ", neither a base nor a compound\n";
                if (++inverted[k.a] > 1) out += "kind " + k.id + ": " + a.id + " already has an inverse\n";
                break;
            }
        }
    }
    if (bases != 4) out += "the table has " + std::to_string(bases) + " base kinds (it needs 4)\n";
    if (cancels != 1) out += "the table has " + std::to_string(cancels) + " cancel kinds (it needs exactly 1: Aether)\n";
    if ((int)kinds.size() > kManaMaxKinds) out += "the table has " + std::to_string(kinds.size()) + " kinds (at most " + std::to_string(kManaMaxKinds) + " in v1)\n";
    return out;
}

int Mana_Inverse(int k) {
    const ManaKind& a = Mana_Kind(k);
    if (a.op == ManaOp::Invert) return a.a;
    for (size_t i = 0; i < kinds.size(); ++i) if (kinds[i].op == ManaOp::Invert && kinds[i].a == k) return (int)i;
    return -1;
}

bool Mana_Cancels(int a, int b) { return a >= 0 && b >= 0 && Mana_Inverse(a) == b; }

int Mana_Aether() {
    for (size_t i = 0; i < kinds.size(); ++i) if (kinds[i].op == ManaOp::Cancel) return (int)i;
    return -1;
}

int Mana_Combine(int a, int b) {
    if (a < 0 || b < 0 || a == b) return -1;
    if (Mana_Cancels(a, b)) return Mana_Aether();
    for (size_t i = 0; i < kinds.size(); ++i)
        if (kinds[i].op == ManaOp::Combine && ((kinds[i].a == a && kinds[i].b == b) || (kinds[i].a == b && kinds[i].b == a))) return (int)i;
    return -1;
}

float Mana_Rate(const std::string& name) { auto it = rates.find(name); return it == rates.end() ? 0.0f : it->second; }

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1248): the shipped kinds table and its algebra ----
#include <gtest/gtest.h>

TEST(Mana, TheShippedTableLoadsAndValidates) {
    ASSERT_TRUE(Mana_Load()) << "resources/data/mana.txt (tests run from the repo root)";
    EXPECT_EQ(Mana_Validate(), "");
    EXPECT_EQ(Mana_Count(), 16) << "the doc's 16 kinds";
    int slice = 0;
    for (int k = 0; k < Mana_Count(); ++k) slice += Mana_Kind(k).slice;
    EXPECT_EQ(slice, 9) << "the first slice: 4 bases, Aether, Frost, Storm, Shard, Rot (decided 2026-11-21)";
    for (const char* id : { "earth", "water", "fire", "air", "aether", "frost", "storm", "shard", "rot" }) {
        const int k = Mana_Find(id);
        ASSERT_GE(k, 0) << id;
        EXPECT_TRUE(Mana_Kind(k).slice) << id;
    }
    for (const char* r : { "combine", "split", "invert", "cancel", "alone", "from_aether" }) EXPECT_GT(Mana_Rate(r), 0.0f) << r;
}

TEST(Mana, TheAlgebraFollowsTheDoc) {
    ASSERT_TRUE(Mana_Load());
    auto K = [](const char* id) { return Mana_Find(id); };
    EXPECT_EQ(Mana_Combine(K("fire"), K("water")), K("steam"));
    EXPECT_EQ(Mana_Combine(K("water"), K("fire")), K("steam")) << "combine is symmetric";
    EXPECT_EQ(Mana_Combine(K("earth"), K("fire")), K("magma"));
    EXPECT_EQ(Mana_Combine(K("earth"), K("water")), K("life"));
    EXPECT_EQ(Mana_Combine(K("fire"), K("air")), K("storm"));
    EXPECT_EQ(Mana_Combine(K("water"), K("air")), K("mist"));
    EXPECT_EQ(Mana_Combine(K("earth"), K("air")), K("shard"));
    EXPECT_EQ(Mana_Inverse(K("fire")), K("frost"));
    EXPECT_EQ(Mana_Inverse(K("frost")), K("fire")) << "inversion is symmetric";
    EXPECT_EQ(Mana_Inverse(K("life")), K("rot"));
    EXPECT_EQ(Mana_Inverse(K("storm")), -1) << "no inverses of other compounds in v1";
    EXPECT_TRUE(Mana_Cancels(K("air"), K("still")));
    EXPECT_EQ(Mana_Combine(K("fire"), K("frost")), K("aether")) << "A + ~A cancels into Aether";
    EXPECT_EQ(Mana_Combine(K("life"), K("rot")), K("aether"));
    EXPECT_EQ(Mana_Combine(K("storm"), K("shard")), -1) << "no compound of compounds";
    EXPECT_EQ(Mana_Combine(K("fire"), K("fire")), -1);
    EXPECT_EQ(Mana_Kind(K("frost")).status, "chill");
    EXPECT_EQ(Mana_Kind(K("earth")).latin, "Terra") << "Latin is flavour beside the plain name";
}

TEST(Mana, ValidationCatchesBadTables) {
    std::string error;
    ASSERT_TRUE(Mana_LoadText(
        "kind a | A | base\nkind b | B | base\nkind c | C | base\n"          // three bases
        "kind ab | AB | combine a b\nkind ba | BA | combine b a\n"         // the same pair twice
        "kind abc | ABC | combine ab c\n"                                   // a compound of a compound
        "kind na | ~A | invert a\nkind na2 | ~A2 | invert a\n"              // a inverted twice
        "kind aa | AA | combine a a\n"                                      // a kind with itself
        "kind x | X | cancel\nkind y | Y | cancel\n",                       // two cancel kinds
        error)) << error;
    const std::string v = Mana_Validate();
    EXPECT_NE(v.find("3 base kinds"), std::string::npos) << v;
    EXPECT_NE(v.find("already makes another kind"), std::string::npos) << v;
    EXPECT_NE(v.find("not a base"), std::string::npos) << v;
    EXPECT_NE(v.find("already has an inverse"), std::string::npos) << v;
    EXPECT_NE(v.find("with itself"), std::string::npos) << v;
    EXPECT_NE(v.find("2 cancel kinds"), std::string::npos) << v;
    std::string many;
    for (int i = 0; i < 17; ++i) many += "kind k" + std::to_string(i) + " | K | " + (i < 4 ? "base" : i == 4 ? "cancel" : "combine k0 k1") + "\n";
    ASSERT_TRUE(Mana_LoadText(many, error)) << error;
    EXPECT_NE(Mana_Validate().find("at most 16"), std::string::npos);
    EXPECT_FALSE(Mana_LoadText("kind x | X | colour red | flavour\n", error)) << "an unknown field is refused";
    EXPECT_FALSE(Mana_LoadText("kind x | X | latin Ex\n", error)) << "a kind needs an operation";
    EXPECT_FALSE(Mana_LoadText("kind x | X | combine a b\n", error)) << "parts must exist";
    EXPECT_FALSE(Mana_LoadText("kind x | X | base | cancel\n", error)) << "one operation";
    EXPECT_FALSE(Mana_LoadText("rate combine 1.5\n", error)) << "a rate is 0..1";
    ASSERT_TRUE(Mana_Load()); // (the shipped table back for the other tests)
}
#endif
