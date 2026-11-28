#ifndef DEPOSITS_H
#define DEPOSITS_H

// (#2071, the underground plan #1465) Ore deposits as drills see them: finite pools of ore units by depth band and size,
// their grade, and each drill's reach (resources/data/deposits.txt). The factory keeps what is left of each deposit a drill
// has touched; the generator (oreveins.h) says where deposits are. Pure functions of the table and their inputs.

#include "depthtable.h"
#include <cstdint>
#include <string>
#include <vector>

struct DepositDrill { std::string part; float reach = 0.0f, spacing = 8.0f; }; // m below its base; m between two heads
struct DepositTable {
    long long pool[3][3] = {};        // [band: shallow, middle, deep][size: small, medium, large] units
    int gradeCenti[3] = { 100, 100, 100 }; // ore out per unit mined, in hundredths
    float rareFactor = 1.0f;
    std::vector<std::string> rare;    // main ore items whose deposits hold rareFactor x the pool
    std::vector<DepositDrill> drills;
};

// Parse deposits.txt (comments and blank lines skipped; bad lines skipped with a warning).
DepositTable Deposits_Parse(const std::string& text);
const DepositTable& Deposits_All();
// A deposit's pool: its band (the surface counts as shallow), its size roll in [0, 1) (thirds: small, medium, large), its
// main ore item.
long long Deposits_Pool(const DepositTable& t, DepthBand band, float sizeRoll, const std::string& mainOre);
int Deposits_GradeCenti(const DepositTable& t, DepthBand band);
const DepositDrill* Deposits_Drill(const DepositTable& t, const std::string& part); // nullptr: not a drill
// A deposit's id: its 96 m cell, packed (21 bits an axis, signed). Never 0 (0 means "no deposit").
uint64_t Deposits_Id(int cx, int cy, int cz);

#endif
