#include "deposits.h"
#include "utils/logger.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

namespace {
int BandIndex(DepthBand b) { return b == DepthBand::Deep ? 2 : b == DepthBand::Middle ? 1 : 0; } // the surface counts as shallow
int BandByName(const std::string& s) { return s == "shallow" ? 0 : s == "middle" ? 1 : s == "deep" ? 2 : -1; }

DepositTable Load() {
    std::ifstream f("resources/data/deposits.txt");
    if (!f) { LOG_ERROR("Deposits", "resources/data/deposits.txt: not found: deposits hold nothing\n"); return {}; }
    std::stringstream ss;
    ss << f.rdbuf();
    return Deposits_Parse(ss.str());
}
} // namespace

DepositTable Deposits_Parse(const std::string& text) {
    DepositTable t;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == ' ')) line.pop_back();
        std::istringstream ls(line);
        std::string tag;
        if (!(ls >> tag) || tag[0] == '#') continue;
        if (tag == "pool") {
            std::string band; long long s = 0, m = 0, l = 0; double grade = 1.0;
            const int b = (ls >> band >> s >> m >> l >> grade) ? BandByName(band) : -1;
            if (b < 0 || s < 0 || m < 0 || l < 0 || grade < 1.0) { LOG_WARN("Deposits", "bad pool line '%s': skipped\n", line.c_str()); continue; }
            t.pool[b][0] = s; t.pool[b][1] = m; t.pool[b][2] = l;
            t.gradeCenti[b] = (int)std::lround(grade * 100.0);
        } else if (tag == "rare") {
            if (!(ls >> t.rareFactor) || t.rareFactor <= 0.0f) { LOG_WARN("Deposits", "bad rare line '%s': skipped\n", line.c_str()); t.rareFactor = 1.0f; continue; }
            for (std::string ore; ls >> ore;) t.rare.push_back(ore);
        } else if (tag == "drill") {
            DepositDrill d;
            if (!(ls >> d.part >> d.reach) || d.reach <= 0.0f) { LOG_WARN("Deposits", "bad drill line '%s': skipped\n", line.c_str()); continue; }
            if (!(ls >> d.spacing)) d.spacing = 8.0f;
            t.drills.push_back(d);
        } else {
            LOG_WARN("Deposits", "unknown line '%s': skipped\n", line.c_str());
        }
    }
    return t;
}

const DepositTable& Deposits_All() {
    static const DepositTable t = Load();
    return t;
}

long long Deposits_Pool(const DepositTable& t, DepthBand band, float sizeRoll, const std::string& mainOre) {
    const int size = sizeRoll < 1.0f / 3.0f ? 0 : sizeRoll < 2.0f / 3.0f ? 1 : 2;
    long long units = t.pool[BandIndex(band)][size];
    if (std::find(t.rare.begin(), t.rare.end(), mainOre) != t.rare.end()) units = (long long)std::llround(units * (double)t.rareFactor);
    return units;
}

int Deposits_GradeCenti(const DepositTable& t, DepthBand band) { return t.gradeCenti[BandIndex(band)]; }

const DepositDrill* Deposits_Drill(const DepositTable& t, const std::string& part) {
    for (const DepositDrill& d : t.drills) if (d.part == part) return &d;
    return nullptr;
}

uint64_t Deposits_Id(int cx, int cy, int cz) {
    auto axis = [](int v) { return (uint64_t)((uint32_t)v & 0x1fffffu); };
    return (1ull << 63) | (axis(cx) << 42) | (axis(cy) << 21) | axis(cz);
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>

TEST(Deposits, ParsesPoolsGradesRareOresAndDrills) {
    const DepositTable t = Deposits_Parse(
        "# comment\r\n"
        "pool shallow 2000 5000 10000 1.0\n"
        "pool middle 10000 25000 50000 1.25\n"
        "pool deep 50000 150000 400000 1.5\n"
        "rare 0.25 gold_ore\n"
        "drill miner 128 8\n"
        "pool nowhere 1 2 3 1.0\n");
    EXPECT_EQ(Deposits_Pool(t, DepthBand::Surface, 0.1f, "iron_ore"), 2000);
    EXPECT_EQ(Deposits_Pool(t, DepthBand::Shallow, 0.5f, "iron_ore"), 5000);
    EXPECT_EQ(Deposits_Pool(t, DepthBand::Middle, 0.9f, "iron_ore"), 50000);
    EXPECT_EQ(Deposits_Pool(t, DepthBand::Deep, 0.5f, "iron_ore"), 150000);
    EXPECT_EQ(Deposits_Pool(t, DepthBand::Deep, 0.99f, "gold_ore"), 100000); // a quarter
    EXPECT_EQ(Deposits_GradeCenti(t, DepthBand::Shallow), 100);
    EXPECT_EQ(Deposits_GradeCenti(t, DepthBand::Middle), 125);
    EXPECT_EQ(Deposits_GradeCenti(t, DepthBand::Deep), 150);
    ASSERT_NE(Deposits_Drill(t, "miner"), nullptr);
    EXPECT_FLOAT_EQ(Deposits_Drill(t, "miner")->reach, 128.0f);
    EXPECT_FLOAT_EQ(Deposits_Drill(t, "miner")->spacing, 8.0f);
    EXPECT_EQ(Deposits_Drill(t, "furnace"), nullptr);
}

TEST(Deposits, IdsAreDistinctPerCellAndNeverZero) {
    EXPECT_NE(Deposits_Id(0, 0, 0), 0u);
    EXPECT_NE(Deposits_Id(-1, 0, 0), Deposits_Id(1, 0, 0));
    EXPECT_NE(Deposits_Id(0, -1, 0), Deposits_Id(0, 0, -1));
    EXPECT_NE(Deposits_Id(600, -600, 3), Deposits_Id(-600, 600, 3));
}

TEST(Deposits, TheShippedTableHasEveryBandAndTheMiner) {
    const DepositTable& t = Deposits_All();
    for (int b = 0; b < 3; ++b) for (int s = 0; s < 3; ++s) EXPECT_GT(t.pool[b][s], 0) << "band " << b << " size " << s;
    EXPECT_NE(Deposits_Drill(t, "miner"), nullptr);
}
#endif
