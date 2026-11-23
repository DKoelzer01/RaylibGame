#ifndef FORGE_H
#define FORGE_H

// (#1613, anvil PR 1) The forging minigame core: a 12x8 height grid (0-3 per cell), hammer strokes that move metal
// between cells, item heat as the clock, scoring against a pattern, and the replay the server runs on a client's stroke
// list. Pure and deterministic, integer math only: the same strokes give the same result everywhere. No UI and no game
// code calls it yet. Patterns come from resources/data/forge_patterns.txt.
//
// Grid index: y * kForgeW + x. y = 0 is the north row.

#include <array>
#include <cstdint>
#include <string>
#include <vector>

constexpr int kForgeW = 12, kForgeH = 8, kForgeMaxHeight = 3;
constexpr int kForgeHeatFull = 1000, kForgeWorkingHeat = 400; // permille; under working heat strokes do nothing
constexpr int kForgeHeatPerStroke = 25;                        // 1/40 a stroke
constexpr int kForgeUnitsPerIngot = 36;                        // a 4x3 block, 3 high

enum class ForgeStroke : uint8_t { Flatten, Draw, Upset, Punch, Reheat };

struct ForgeMove {
    ForgeStroke type;
    int8_t x = 0, y = 0;
    int8_t dx = 0, dy = 0;
    uint32_t tick = 0; // 20 Hz time since start
};

struct ForgePattern {
    std::string id;
    int ingots = 1;
    int par = 0;
    std::array<uint8_t, kForgeW * kForgeH> h{};
    int units() const;
};

struct ForgeState {
    std::array<uint8_t, kForgeW * kForgeH> h{};
    int heat = 0;      // permille
    uint32_t tick = 0; // the last stroke's tick (time cools between strokes)
    int strokes = 0;   // strokes that changed something (no-ops don't count)
    int waste = 0;     // units punched out
    int reheats = 0;
    int inputUnits = 0;
};

// 1 ingot: a 4x3 block 3 high at columns 4..7, rows 2..4 (36 units); 2 ingots: 8x3 at columns 2..9, rows 2..4 (72).
// Other counts are clamped to 1..2.
ForgeState Forge_Start(int ingots, int heatPermille);

// Cools the workpiece for the time since the last move, then does the stroke. Returns whether the grid changed.
// Metal is conserved: sum(h) + waste == inputUnits.
bool Forge_Apply(ForgeState& s, const ForgeMove& m, int coolPerSecondPermille, int hammerTier, int metalTier);

// Permille per second: tier 1 (copper, bronze) 20, tier 2 (iron) 17, tier 3+ (steel and up) 14.
int Forge_CoolRate(int metalTier);

// Scoring: 75..130, capped by the anvil and the player's skill (the caps come in later PRs; pass 130 for none).
int Forge_Quality(const ForgeState& s, const ForgePattern& p, int anvilCap, int skillCap);
const char* Forge_Grade(int q); // "Crude" < 90, "Common" 90-109, "Fine" 110-124, "Masterwork" >= 125
int Forge_QuickQuality();       // the Quick forge: 95

struct ForgeResult {
    int q = 75;
    int waste = 0;
    int strokes = 0;
    int reheats = 0;
    bool valid = true;
};
// What the server runs: Forge_Start with the pattern's ingot count, every move in order, then score. valid = false
// (and q = 75) for more than 4096 moves or ticks past 20 minutes. outState (optional) receives the final state.
ForgeResult Forge_Replay(const ForgePattern& p, int heatPermille, int metalTier, int hammerTier, int anvilCap,
                         int skillCap, const std::vector<ForgeMove>& moves, ForgeState* outState = nullptr);

// Patterns: parse (errors appended as "line N: message"; bad blocks are skipped) and the game file, loaded once.
std::vector<ForgePattern> Forge_ParsePatterns(const std::string& text, std::vector<std::string>* errors);
const std::vector<ForgePattern>& Forge_AllPatterns(); // resources/data/forge_patterns.txt
const ForgePattern* Forge_Pattern(const std::string& id);

// A greedy solver (tests now, hints later): Draw strokes that carry surplus units to deficit cells along paths of
// material, Reheats when the heat would drop under working heat, ticks +10 per move. Never punches. Gives up after
// 2000 moves and returns what it has. Deterministic.
std::vector<ForgeMove> Forge_Solve(const ForgePattern& p, int metalTier);

#endif
