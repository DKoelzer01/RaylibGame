#include "cacheopeners.h"
#include "inventory.h"
#include "plots.h"

#include <cstdio>
#include <fstream>
#include <sstream>

namespace {
struct Catalogue {
    bool loaded = false;
    std::vector<CacheOpener> openers;
    std::string defaults[3];          // by CacheKind: the opener id of a cache whose part names none
} C;
CachePowerHooks power;

int KindNamed(const std::string& k) { return k == "small" ? 0 : k == "vault" ? 1 : k == "deep_vault" ? 2 : -1; }

void LoadOnce() {
    if (C.loaded) return;
    C.loaded = true;
    std::ifstream in("resources/data/cache_openers.txt");
    if (!in) { fprintf(stderr, "cache openers: resources/data/cache_openers.txt not found: every cache stays sealed\n"); return; }
    std::stringstream ss; ss << in.rdbuf();
    std::string error;
    if (!CacheOpeners_Parse(ss.str(), &error)) fprintf(stderr, "cache openers: %s\n", error.c_str());
}

// 2000 -> "2,000"
std::string Thousands(long long v) {
    std::string s = std::to_string(v);
    for (int i = (int)s.size() - 3; i > 0; i -= 3) s.insert((size_t)i, ",");
    return s;
}
std::string ItemName(const std::string& id) { const ItemDef* d = itemDB.find(id); return d ? d->name : id; }
// "a Crowbar", "an Iron ingot"
std::string WithArticle(const std::string& name) {
    const char c = name.empty() ? 'x' : (char)tolower((unsigned char)name[0]);
    return (std::string("aeiou").find(c) != std::string::npos ? "an " : "a ") + name;
}
} // namespace

bool CacheOpeners_Parse(const std::string& text, std::string* error) {
    C = Catalogue{};
    C.loaded = true;
    std::istringstream in(text);
    int lineNo = 0;
    auto bad = [&](const std::string& why) { if (error) *error = "line " + std::to_string(lineNo) + ": " + why; return false; };
    for (std::string line; std::getline(in, line);) {
        ++lineNo;
        if (const size_t hash = line.find_first_of("#\r"); hash != std::string::npos) line.resize(hash); // a comment; a CRLF checkout's CR on Linux
        std::istringstream ls(line);
        std::string cmd;
        if (!(ls >> cmd)) continue;
        if (cmd == "opener") {
            CacheOpener o;
            std::string kind;
            if (!(ls >> o.index >> o.id >> kind)) return bad("opener <n> <id> <power|tool|part> ...");
            if (o.index < 1 || o.index > 7) return bad("an opener's number is 1-7 (the part's material keeps it in 3 bits)");
            if (CacheOpener_ByIndex(o.index) || CacheOpener_Find(o.id)) return bad("opener " + o.id + ": its number or id is used already");
            if (kind == "power") {
                o.kind = OpenerKind::Power;
                if (!(ls >> o.eu >> o.capEUt >> o.minVolts) || o.eu <= 0 || o.capEUt <= 0 || o.minVolts < 0) return bad("power <EU> <cap EU/t> <min V>");
            } else if (kind == "tool" || kind == "part") {
                o.kind = kind == "tool" ? OpenerKind::Tool : OpenerKind::Part;
                if (!(ls >> o.item)) return bad(kind + " <item>");
            } else return bad("an opener is power, tool or part, not " + kind);
            C.openers.push_back(o);
        } else if (cmd == "default") {
            std::string kind, id;
            if (!(ls >> kind >> id) || KindNamed(kind) < 0) return bad("default <small|vault|deep_vault> <opener id>");
            if (!CacheOpener_Find(id)) return bad("default " + kind + ": no opener " + id + " above");
            C.defaults[KindNamed(kind)] = id;
        } else return bad("unknown line " + cmd);
    }
    return true;
}

void CacheOpeners_Reset() { C = Catalogue{}; }

const CacheOpener* CacheOpener_Find(const std::string& id) {
    LoadOnce();
    for (const CacheOpener& o : C.openers) if (o.id == id) return &o;
    return nullptr;
}

const CacheOpener* CacheOpener_ByIndex(int index) {
    LoadOnce();
    for (const CacheOpener& o : C.openers) if (o.index == index) return &o;
    return nullptr;
}

const CacheOpener* CacheOpenerOf(const Part& p) {
    LoadOnce(); // (before C.defaults is read below: the first call could find "" otherwise)
    if (const int n = (p.material & kCacheOpenerMask) >> kCacheOpenerShift; n != 0)
        if (const CacheOpener* o = CacheOpener_ByIndex(n)) return o;
    const int k = p.material & 3; // the kind (ruincache.h's CacheKindOf)
    return CacheOpener_Find(C.defaults[k > 2 ? 0 : k]);
}

std::string CacheOpener_NeedText(const CacheOpener& o) {
    switch (o.kind) {
    case OpenerKind::Power: return "needs " + Thousands(o.eu) + " EU at " + std::to_string(o.minVolts) + " V or more (a battery bank or a cable)";
    case OpenerKind::Tool: return "needs " + WithArticle(ItemName(o.item)) + " (kept)";
    case OpenerKind::Part: return "needs " + WithArticle(ItemName(o.item)) + " (used up)";
    }
    return "";
}

void CacheOpeners_SetPower(CachePowerHooks hooks) { power = std::move(hooks); }
float CacheOpener_Charge(const Planetoid& body, const Plot& plot, int part) { return power.charged ? power.charged(body, plot, part) : 0.0f; }
void CacheOpener_Spend(Planetoid& body, Plot& plot, int part) { if (power.spend) power.spend(body, plot, part); }

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>

namespace {
const char* kTestText =
    "opener 1 small_power power 2000 8 8\n"
    "opener 2 crowbar tool crowbar   # kept\n"
    "opener 3 ruin_key part ruin_key\n"
    "default small small_power\n"
    "default vault crowbar\n";
Part CachePart(int kind, int opener) { Part p; p.type = PartType::RuinCache; p.material = (uint8_t)(kind | CacheOpenerBits(opener)); return p; }
// Each test parses its own catalogue; afterwards the game's file is read again on next use.
class CacheOpenersTest : public ::testing::Test { protected: void TearDown() override { CacheOpeners_Reset(); } };
}

TEST_F(CacheOpenersTest, ParsesTheThreeKinds) {
    std::string err;
    ASSERT_TRUE(CacheOpeners_Parse(kTestText, &err)) << err;
    const CacheOpener* p = CacheOpener_Find("small_power");
    ASSERT_NE(p, nullptr);
    EXPECT_EQ(p->kind, OpenerKind::Power);
    EXPECT_EQ(p->eu, 2000); EXPECT_EQ(p->capEUt, 8); EXPECT_EQ(p->minVolts, 8);
    ASSERT_NE(CacheOpener_ByIndex(2), nullptr);
    EXPECT_EQ(CacheOpener_ByIndex(2)->kind, OpenerKind::Tool);
    EXPECT_EQ(CacheOpener_ByIndex(3)->kind, OpenerKind::Part);
    EXPECT_EQ(CacheOpener_ByIndex(3)->item, "ruin_key");
    EXPECT_EQ(CacheOpener_NeedText(*p), "needs 2,000 EU at 8 V or more (a battery bank or a cable)");
    EXPECT_EQ(CacheOpener_NeedText(*CacheOpener_ByIndex(3)), "needs a ruin_key (used up)"); // (no item data here: the id)
    ASSERT_TRUE(CacheOpeners_Parse("opener 1 ingot part iron_ingot\n"));
    EXPECT_EQ(CacheOpener_NeedText(*CacheOpener_ByIndex(1)), "needs an iron_ingot (used up)");
}

TEST_F(CacheOpenersTest, APartKeepsItsOpenerElseItsKindsDefault) {
    ASSERT_TRUE(CacheOpeners_Parse(kTestText));
    EXPECT_EQ(CacheOpenerOf(CachePart(0, 0))->id, "small_power"); // small, none named: the default
    EXPECT_EQ(CacheOpenerOf(CachePart(1, 0))->id, "crowbar");     // vault: its default
    EXPECT_EQ(CacheOpenerOf(CachePart(0, 3))->id, "ruin_key");    // named in the part
    EXPECT_EQ(CacheOpenerOf(CachePart(0 | 0x10, 2))->id, "crowbar"); // the opened bit doesn't disturb it
    EXPECT_EQ(CacheOpenerOf(CachePart(2, 0)), nullptr);           // deep vault: no default in this catalogue
    EXPECT_EQ(CacheOpenerOf(CachePart(0, 6))->id, "small_power"); // an unknown number: the kind's default
}

TEST_F(CacheOpenersTest, RefusesBadLines) {
    std::string err;
    EXPECT_FALSE(CacheOpeners_Parse("opener 8 x tool crowbar\n", &err)); // 3 bits
    EXPECT_FALSE(CacheOpeners_Parse("opener 1 a tool crowbar\nopener 1 b tool crowbar\n", &err)); // number reused
    EXPECT_FALSE(CacheOpeners_Parse("opener 1 a power 0 8 8\n", &err));
    EXPECT_FALSE(CacheOpeners_Parse("default small nothing\n", &err));
    EXPECT_FALSE(CacheOpeners_Parse("opener 1 a key ruin_key\n", &err));
}

TEST_F(CacheOpenersTest, TheGamesFileLoads) {
    std::ifstream in("resources/data/cache_openers.txt");
    ASSERT_TRUE(in.good());
    std::stringstream ss; ss << in.rdbuf();
    std::string err;
    ASSERT_TRUE(CacheOpeners_Parse(ss.str(), &err)) << err;
    for (int k = 0; k < 3; ++k) EXPECT_NE(CacheOpenerOf(CachePart(k, 0)), nullptr) << "kind " << k << " has no default";
}
#endif
