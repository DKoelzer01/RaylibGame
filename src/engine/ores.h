#ifndef ORES_H
#define ORES_H

// (#1466, the underground plan #1465) Ore kinds from resources/data/ores.txt: what PlanetField::oreKind places (deposit
// and pocket shares, depth) and how the terrain shader draws each (colour, pattern). A kind is its row index; the mesher
// sends it to the shader as (kind + 1) / 256, so up to kOresMax kinds. Loaded once, on first use (thread-safe).

#include <string>
#include <vector>

enum class OrePattern { Veins = 0, VeinsCrust, Speckle, Seams, Nodules, Crystals }; // lighting.fs, its ore block

struct OreDef {
    std::string id, item;
    float colour[3] = { 0.5f, 0.5f, 0.5f };
    OrePattern pattern = OrePattern::Veins;
    float minDepth = 4.0f;
    float deposit = 0.0f, pocket = 0.0f; // shares
    int roll = 0;                        // the deposit roll's order
};

constexpr int kOresMax = 64; // the shader's ore table (lighting.fs kOres)

const std::vector<OreDef>& Ores_All();
int Ores_Count();
const OreDef* Ores_Get(int kind);
// The deposit and pocket rolls: a uniform r in [0, 1) -> a kind (cumulative shares, in roll order / row order; the last
// with a share takes the rest). -1 when no ore has a share.
int Ores_DepositKind(float r);
int Ores_PocketKind(float r);
// The shader's table: colours (3 per kind) and pattern ids, kOresMax entries (unused ones zero).
void Ores_ShaderTable(float* colours, int* patterns);

#endif
