#include "ores.h"
#include "utils/logger.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

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
bool PatternOf(const std::string& s, OrePattern& p) {
    static const char* names[] = { "veins", "veins_crust", "speckle", "seams", "nodules", "crystals" };
    for (int i = 0; i < (int)(sizeof(names) / sizeof(names[0])); ++i) if (s == names[i]) { p = (OrePattern)i; return true; }
    return false;
}

struct Table {
    std::vector<OreDef> ores;
    // Cumulative shares as floats, with their kinds.
    std::vector<std::pair<float, int>> deposit, pocket;
};

Table Load(const char* path) {
    Table t;
    std::ifstream f(path);
    if (!f) { LOG_ERROR("Ores", "%s: not found: no ores\n", path); return t; }
    std::string line;
    std::vector<double> depShare, pocShare; // (parsed in double: summed so, 0.22 + 0.22 lands on the same float as 0.44f)
    while (std::getline(f, line)) {
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;
        const auto p = Split(line, '|');
        if (p.size() < 8) { LOG_WARN("Ores", "%s: short line '%s'\n", path, line.c_str()); continue; }
        if ((int)t.ores.size() >= kOresMax - 1) { LOG_WARN("Ores", "%s: more than %d kinds: the rest ignored\n", path, kOresMax - 1); break; }
        OreDef o;
        o.id = p[0]; o.item = p[1];
        const auto rgb = Split(p[2], ',');
        if (rgb.size() == 3) for (int k = 0; k < 3; ++k) o.colour[k] = strtof(rgb[k].c_str(), nullptr);
        if (!PatternOf(p[3], o.pattern)) LOG_WARN("Ores", "%s: %s: unknown pattern '%s' (veins)\n", path, o.id.c_str(), p[3].c_str());
        o.minDepth = strtof(p[4].c_str(), nullptr);
        o.deposit = strtof(p[5].c_str(), nullptr);
        o.roll = atoi(p[6].c_str());
        o.pocket = strtof(p[7].c_str(), nullptr);
        depShare.push_back(strtod(p[5].c_str(), nullptr));
        pocShare.push_back(strtod(p[7].c_str(), nullptr));
        t.ores.push_back(o);
    }
    auto cumulative = [&](const std::vector<int>& order, const std::vector<double>& shares, std::vector<std::pair<float, int>>& out) {
        double sum = 0.0;
        for (int k : order) {
            const double share = shares[k];
            if (share <= 0.0) continue;
            sum += share;
            out.push_back({ (float)sum, k });
        }
    };
    std::vector<int> rows(t.ores.size());
    for (size_t i = 0; i < rows.size(); ++i) rows[i] = (int)i;
    std::vector<int> rolled = rows;
    std::stable_sort(rolled.begin(), rolled.end(), [&](int a, int b) { return t.ores[a].roll < t.ores[b].roll; });
    cumulative(rolled, depShare, t.deposit);
    cumulative(rows, pocShare, t.pocket);
    return t;
}

const Table& T() {
    static const Table t = Load("resources/data/ores.txt"); // (thread-safe: the mesher's workers ask too)
    return t;
}

int Pick(const std::vector<std::pair<float, int>>& cum, float r) {
    if (cum.empty()) return -1;
    for (size_t i = 0; i + 1 < cum.size(); ++i) if (r < cum[i].first) return cum[i].second;
    return cum.back().second; // the last takes the rest
}
} // namespace

const std::vector<OreDef>& Ores_All() { return T().ores; }
int Ores_Count() { return (int)T().ores.size(); }
const OreDef* Ores_Get(int kind) { return kind >= 0 && kind < Ores_Count() ? &T().ores[kind] : nullptr; }
int Ores_DepositKind(float r) { return Pick(T().deposit, r); }
int Ores_PocketKind(float r) { return Pick(T().pocket, r); }

void Ores_ShaderTable(float* colours, int* patterns) {
    memset(colours, 0, sizeof(float) * 3 * kOresMax);
    memset(patterns, 0, sizeof(int) * kOresMax);
    const auto& o = T().ores;
    for (size_t k = 0; k < o.size() && (int)k < kOresMax; ++k) {
        for (int c = 0; c < 3; ++c) colours[k * 3 + c] = o[k].colour[c];
        patterns[k] = (int)o[k].pattern;
    }
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"
// Today's rules before the table (#1466): the table must reproduce them exactly.
static int OldDepositKind(float r) { return r < 0.22f ? 3 : r < 0.44f ? 0 : r < 0.64f ? 1 : r < 0.76f ? 2 : r < 0.90f ? 4 : 5; }
static int OldPocketKind(float r) { return r < 0.35f ? 0 : r < 0.60f ? 1 : r < 0.75f ? 2 : 3; }

TEST(Ores, TheTableHoldsTodaysSixKinds) {
    ASSERT_EQ(Ores_Count(), 6);
    const char* items[] = { "iron_ore", "copper_ore", "tin_ore", "coal", "bauxite_ore", "ilmenite_ore" };
    const float minDepth[] = { 4, 4, 4, 4, 20, 60 };
    for (int k = 0; k < 6; ++k) { EXPECT_EQ(Ores_Get(k)->item, items[k]); EXPECT_EQ(Ores_Get(k)->minDepth, minDepth[k]) << k; }
    // the shader's colours and patterns, as lighting.fs had them per kind
    const float col[6][3] = { { 0.58f, 0.27f, 0.08f }, { 0.10f, 0.42f, 0.30f }, { 0.62f, 0.62f, 0.64f }, { 0.03f, 0.03f, 0.035f }, { 0.62f, 0.22f, 0.16f }, { 0.13f, 0.08f, 0.16f } };
    float c[3 * kOresMax]; int pat[kOresMax];
    Ores_ShaderTable(c, pat);
    for (int k = 0; k < 6; ++k) {
        for (int i = 0; i < 3; ++i) EXPECT_EQ(c[k * 3 + i], col[k][i]) << k;
        EXPECT_EQ(pat[k], k) << "kind k drew pattern k";
    }
}

TEST(Ores, DepositAndPocketRollsMatchTheOldRules) {
    // every float r in [0, 1) near the thresholds, and a sweep
    std::vector<float> rs;
    for (int i = 0; i < 100000; ++i) rs.push_back(i / 100000.0f);
    for (float th : { 0.22f, 0.44f, 0.64f, 0.76f, 0.90f, 0.35f, 0.60f, 0.75f }) {
        float a = th, b = th;
        for (int s = 0; s < 4; ++s) { a = nextafterf(a, 0.0f); b = nextafterf(b, 1.0f); rs.push_back(a); rs.push_back(b); }
        rs.push_back(th);
    }
    for (float r : rs) {
        ASSERT_EQ(Ores_DepositKind(r), OldDepositKind(r)) << r;
        ASSERT_EQ(Ores_PocketKind(r), OldPocketKind(r)) << r;
    }
}
#endif
