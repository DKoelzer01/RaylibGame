#include "forge.h"
#include "utils/logger.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace {
constexpr int kCells = kForgeW * kForgeH;
constexpr int kMaxReplayMoves = 4096;
constexpr uint32_t kMaxReplayTicks = 20u * 60u * 20u; // 20 minutes at 20 Hz
constexpr int kSolveMaxMoves = 2000;
constexpr uint32_t kSolveTickStep = 10;

// Neighbour order everywhere: -y, +x, +y, -x.
constexpr int kDx[4] = { 0, 1, 0, -1 };
constexpr int kDy[4] = { -1, 0, 1, 0 };

bool InGrid(int x, int y) { return x >= 0 && y >= 0 && x < kForgeW && y < kForgeH; }
int Idx(int x, int y) { return y * kForgeW + x; }

std::string Trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r' || s[a] == '\n')) ++a;
    while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r' || s[b - 1] == '\n')) --b;
    return s.substr(a, b - a);
}

std::string ReadFile(const char* path) {
    std::ifstream f(path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

void FillBlock(ForgeState& s, int x0, int x1, int y0, int y1) {
    for (int y = y0; y <= y1; ++y)
        for (int x = x0; x <= x1; ++x) s.h[Idx(x, y)] = kForgeMaxHeight;
}
} // namespace

int ForgePattern::units() const {
    int n = 0;
    for (uint8_t v : h) n += v;
    return n;
}

ForgeState Forge_Start(int ingots, int heatPermille) {
    ForgeState s;
    ingots = std::max(1, std::min(2, ingots));
    if (ingots == 1) FillBlock(s, 4, 7, 2, 4);
    else FillBlock(s, 2, 9, 2, 4);
    s.heat = heatPermille;
    s.inputUnits = ingots * kForgeUnitsPerIngot;
    return s;
}

int Forge_CoolRate(int metalTier) {
    if (metalTier <= 1) return 20;
    if (metalTier == 2) return 17;
    return 14;
}

int Forge_MetalTier(const std::string& itemOrMaterial) {
    std::string s = itemOrMaterial;
    if (s.compare(0, 6, "rough_") == 0) s = s.substr(6);
    const std::string m = s.substr(0, s.find('_'));
    if (m == "copper" || m == "tin" || m == "bronze") return 1;
    if (m == "iron") return 2;
    if (m == "steel") return 3;
    if (m == "titanium" || m == "tungsten") return 4;
    return 0;
}

int Forge_JobTier(const std::string& jobItem) {
    const int t = Forge_MetalTier(jobItem);
    const std::string suffix = "_hammer_head";
    const bool hammerHead = jobItem.size() > suffix.size() && jobItem.compare(jobItem.size() - suffix.size(), suffix.size(), suffix) == 0;
    return std::max(0, hammerHead ? t - 1 : t);
}

bool Forge_Apply(ForgeState& s, const ForgeMove& m, int coolPerSecondPermille, int hammerTier, int metalTier) {
    // Time first. Ticks never go backwards.
    const uint32_t t = m.tick < s.tick ? s.tick : m.tick;
    const int64_t loss = (int64_t)coolPerSecondPermille * (int64_t)(t - s.tick) / 20;
    s.heat = (int)std::max<int64_t>(0, (int64_t)s.heat - loss);
    s.tick = t;

    if (m.type == ForgeStroke::Reheat) {
        s.heat = kForgeHeatFull;
        ++s.reheats;
        return false;
    }
    if (hammerTier < metalTier || s.heat < kForgeWorkingHeat) return false;
    const int x = m.x, y = m.y;
    if (!InGrid(x, y)) return false;

    bool changed = false;
    uint8_t& here = s.h[Idx(x, y)];
    switch (m.type) {
    case ForgeStroke::Flatten: {
        if (here == 0) break;
        int best = -1, bestH = kForgeMaxHeight; // the lowest neighbour that is < max; ties: first in order
        for (int d = 0; d < 4; ++d) {
            const int nx = x + kDx[d], ny = y + kDy[d];
            if (!InGrid(nx, ny)) continue;
            const int nh = s.h[Idx(nx, ny)];
            if (nh < bestH) { bestH = nh; best = Idx(nx, ny); }
        }
        if (best < 0) break;
        --here;
        ++s.h[best];
        changed = true;
        break;
    }
    case ForgeStroke::Draw: {
        const int adx = std::abs((int)m.dx), ady = std::abs((int)m.dy);
        if (adx + ady != 1) break;
        const int nx = x + m.dx, ny = y + m.dy;
        if (here == 0 || !InGrid(nx, ny) || s.h[Idx(nx, ny)] >= kForgeMaxHeight) break;
        --here;
        ++s.h[Idx(nx, ny)];
        changed = true;
        break;
    }
    case ForgeStroke::Upset: {
        if (here >= kForgeMaxHeight) break;
        int best = -1, bestH = 0; // the highest neighbour that is > 0; ties: first in order
        for (int d = 0; d < 4; ++d) {
            const int nx = x + kDx[d], ny = y + kDy[d];
            if (!InGrid(nx, ny)) continue;
            const int nh = s.h[Idx(nx, ny)];
            if (nh > bestH) { bestH = nh; best = Idx(nx, ny); }
        }
        if (best < 0) break;
        --s.h[best];
        ++here;
        changed = true;
        break;
    }
    case ForgeStroke::Punch:
        if (here == 0) break;
        --here;
        ++s.waste;
        changed = true;
        break;
    case ForgeStroke::Reheat: break;
    }
    if (changed) {
        ++s.strokes;
        s.heat = std::max(0, s.heat - kForgeHeatPerStroke);
    }
    return changed;
}

int Forge_Quality(const ForgeState& s, const ForgePattern& p, int anvilCap, int skillCap) {
    int covered = 0, match = 0;
    for (int i = 0; i < kCells; ++i) {
        if (p.h[i] > 0 || s.h[i] > 0) {
            ++covered;
            if (p.h[i] == s.h[i]) ++match;
        }
    }
    const int extra = std::max(0, s.strokes - p.par);
    const int wasteOver = std::max(0, s.waste - (s.inputUnits - p.units()));
    int raw = match * 1000 / std::max(1, covered) - 20 * extra - 30 * wasteOver - 20 * s.reheats;
    raw = std::max(0, std::min(1000, raw));
    return std::min(std::min(75 + 55 * raw / 1000, anvilCap), skillCap);
}

const char* Forge_Grade(int q) {
    if (q < 90) return "Crude";
    if (q < 110) return "Common";
    if (q < 125) return "Fine";
    return "Masterwork";
}

int Forge_QuickQuality() { return 95; }

ForgeResult Forge_Replay(const ForgePattern& p, int heatPermille, int metalTier, int hammerTier, int anvilCap,
                         int skillCap, const std::vector<ForgeMove>& moves, ForgeState* outState) {
    ForgeResult r;
    ForgeState s = Forge_Start(p.ingots, heatPermille);
    bool valid = moves.size() <= (size_t)kMaxReplayMoves;
    for (const ForgeMove& m : moves)
        if (m.tick > kMaxReplayTicks) valid = false;
    if (!valid) {
        r.valid = false;
        r.q = 75;
        if (outState) *outState = s;
        return r;
    }
    const int cool = Forge_CoolRate(metalTier);
    for (const ForgeMove& m : moves) Forge_Apply(s, m, cool, hammerTier, metalTier);
    r.q = Forge_Quality(s, p, anvilCap, skillCap);
    r.waste = s.waste;
    r.strokes = s.strokes;
    r.reheats = s.reheats;
    r.valid = true;
    if (outState) *outState = s;
    return r;
}

// ---- Patterns ----

std::vector<ForgePattern> Forge_ParsePatterns(const std::string& text, std::vector<std::string>* errors) {
    std::vector<ForgePattern> out;
    auto err = [&](int ln, const std::string& msg) {
        if (errors) errors->push_back("line " + std::to_string(ln) + ": " + msg);
    };
    std::istringstream in(text);
    std::string raw;
    int ln = 0;
    bool inBlock = false, bad = false, hasPar = false;
    int rows = 0, startLine = 0;
    ForgePattern cur;

    auto finish = [&](int endLine) {
        if (rows != kForgeH) {
            err(endLine, "pattern " + cur.id + " has " + std::to_string(rows) + " rows, expected " + std::to_string(kForgeH));
            bad = true;
        }
        if (!hasPar) {
            err(startLine, "pattern " + cur.id + " is missing par");
            bad = true;
        }
        if (!bad && cur.units() != cur.ingots * kForgeUnitsPerIngot) {
            err(startLine, "pattern " + cur.id + " has " + std::to_string(cur.units()) + " units, " +
                               std::to_string(cur.ingots) + " ingot(s) give " +
                               std::to_string(cur.ingots * kForgeUnitsPerIngot));
            bad = true;
        }
        if (!bad) out.push_back(cur);
        inBlock = false;
    };

    while (std::getline(in, raw)) {
        ++ln;
        const std::string line = Trim(raw);
        if (line.compare(0, 8, "pattern ") == 0) {
            if (inBlock) {
                err(ln, "pattern " + cur.id + " is missing 'end'");
                inBlock = false;
            }
            inBlock = true;
            bad = false;
            hasPar = false;
            rows = 0;
            startLine = ln;
            cur = ForgePattern();
            std::vector<std::string> parts;
            size_t a = 0;
            for (size_t i = 0; i <= line.size(); ++i)
                if (i == line.size() || line[i] == '|') { parts.push_back(Trim(line.substr(a, i - a))); a = i + 1; }
            cur.id = Trim(parts[0].substr(8));
            if (cur.id.empty()) { err(ln, "pattern without an id"); bad = true; }
            for (size_t i = 1; i < parts.size(); ++i) {
                std::istringstream ts(parts[i]);
                std::string key;
                int val = 0;
                if (!(ts >> key >> val)) { err(ln, "bad field '" + parts[i] + "'"); bad = true; continue; }
                if (key == "ingots") {
                    if (val < 1 || val > 2) { err(ln, "ingots must be 1 or 2"); bad = true; }
                    else cur.ingots = val;
                } else if (key == "par") {
                    if (val < 0) { err(ln, "par must not be negative"); bad = true; }
                    else { cur.par = val; hasPar = true; }
                } else {
                    err(ln, "unknown field '" + key + "'");
                    bad = true;
                }
            }
            continue;
        }
        if (!inBlock) {
            if (!line.empty() && line[0] != '#') err(ln, "expected 'pattern <id> | ingots <n> | par <strokes>'");
            continue;
        }
        if (line.empty() || line[0] == '#') continue;
        if (line == "end") { finish(ln); continue; }
        // A row of 12 cells.
        if (rows >= kForgeH) {
            if (rows == kForgeH) { err(ln, "too many rows (expected " + std::to_string(kForgeH) + ")"); bad = true; }
            ++rows;
            continue;
        }
        if ((int)line.size() != kForgeW) {
            err(ln, "row has " + std::to_string(line.size()) + " cells, expected " + std::to_string(kForgeW));
            bad = true;
        } else {
            for (int x = 0; x < kForgeW; ++x) {
                const char c = line[x];
                if (c == '.') cur.h[rows * kForgeW + x] = 0;
                else if (c >= '0' && c <= '3') cur.h[rows * kForgeW + x] = (uint8_t)(c - '0');
                else {
                    err(ln, std::string("bad height '") + c + "' (use . or 0-3)");
                    bad = true;
                    break;
                }
            }
        }
        ++rows;
    }
    if (inBlock) err(ln, "pattern " + cur.id + " is missing 'end'");
    return out;
}

const std::vector<ForgePattern>& Forge_AllPatterns() {
    static const std::vector<ForgePattern> all = [] {
        std::vector<std::string> errors;
        std::vector<ForgePattern> v = Forge_ParsePatterns(ReadFile("resources/data/forge_patterns.txt"), &errors);
        for (const std::string& e : errors) LOG_WARN("Forge", "%s\n", e.c_str());
        return v;
    }();
    return all;
}

const ForgePattern* Forge_Pattern(const std::string& id) {
    for (const ForgePattern& p : Forge_AllPatterns())
        if (p.id == id) return &p;
    return nullptr;
}

// ---- Solver ----

std::vector<ForgeMove> Forge_Solve(const ForgePattern& p, int metalTier) {
    std::vector<ForgeMove> out;
    ForgeState s = Forge_Start(p.ingots, kForgeHeatFull);
    const int cool = Forge_CoolRate(metalTier);
    uint32_t tick = 0;

    // Appends one move (a Reheat first if the stroke would meet a too-cold piece) and applies it.
    auto emit = [&](ForgeStroke type, int x, int y, int dx, int dy) -> bool {
        if (type != ForgeStroke::Reheat) {
            const int predicted = s.heat - cool * (int)kSolveTickStep / 20;
            if (predicted < kForgeWorkingHeat) {
                tick += kSolveTickStep;
                ForgeMove r;
                r.type = ForgeStroke::Reheat;
                r.tick = tick;
                out.push_back(r);
                Forge_Apply(s, r, cool, metalTier, metalTier);
            }
        }
        tick += kSolveTickStep;
        ForgeMove m;
        m.type = type;
        m.x = (int8_t)x;
        m.y = (int8_t)y;
        m.dx = (int8_t)dx;
        m.dy = (int8_t)dy;
        m.tick = tick;
        out.push_back(m);
        return Forge_Apply(s, m, cool, metalTier, metalTier);
    };

    while ((int)out.size() < kSolveMaxMoves) {
        // The shortest path over material from a surplus cell to a deficit cell (the first found wins ties).
        std::vector<int> bestPath;
        for (int src = 0; src < kCells; ++src) {
            if (s.h[src] <= p.h[src]) continue;
            int parent[kCells];
            bool seen[kCells];
            std::fill(parent, parent + kCells, -1);
            std::fill(seen, seen + kCells, false);
            std::vector<int> queue;
            queue.push_back(src);
            seen[src] = true;
            int found = -1;
            for (size_t qi = 0; qi < queue.size() && found < 0; ++qi) {
                const int cur = queue[qi];
                const int cx = cur % kForgeW, cy = cur / kForgeW;
                for (int d = 0; d < 4; ++d) {
                    const int nx = cx + kDx[d], ny = cy + kDy[d];
                    if (!InGrid(nx, ny)) continue;
                    const int n = Idx(nx, ny);
                    if (seen[n]) continue;
                    seen[n] = true;
                    parent[n] = cur;
                    if (s.h[n] < p.h[n]) { found = n; break; }
                    if (s.h[n] >= 1) queue.push_back(n);
                }
            }
            if (found < 0) continue;
            std::vector<int> path;
            for (int c = found; c >= 0; c = parent[c]) path.push_back(c);
            std::reverse(path.begin(), path.end());
            if (bestPath.empty() || path.size() < bestPath.size()) bestPath = path;
        }
        if (bestPath.empty()) break;
        // Draw from the far end back, so every step finds room and the cells between end up unchanged.
        bool ok = true;
        for (int i = (int)bestPath.size() - 2; i >= 0 && ok; --i) {
            const int a = bestPath[i], b = bestPath[i + 1];
            if ((int)out.size() >= kSolveMaxMoves) return out;
            ok = emit(ForgeStroke::Draw, a % kForgeW, a / kForgeW, b % kForgeW - a % kForgeW, b / kForgeW - a / kForgeW);
        }
        if (!ok) break; // cannot happen; never loop forever
    }
    return out;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"

namespace {
int Sum(const ForgeState& s) {
    int n = 0;
    for (uint8_t v : s.h) n += v;
    return n;
}
ForgeMove Mv(ForgeStroke t, int x, int y, int dx = 0, int dy = 0, uint32_t tick = 0) {
    ForgeMove m;
    m.type = t;
    m.x = (int8_t)x;
    m.y = (int8_t)y;
    m.dx = (int8_t)dx;
    m.dy = (int8_t)dy;
    m.tick = tick;
    return m;
}
// A small hand-built state: all zero, full heat, the given sum counted as input.
ForgeState Empty() {
    ForgeState s;
    s.heat = kForgeHeatFull;
    return s;
}
void SetInput(ForgeState& s) { s.inputUnits = Sum(s); }
} // namespace

TEST(Forge, StartBlocks) {
    const ForgeState a = Forge_Start(1, kForgeHeatFull);
    EXPECT_EQ(Sum(a), 36);
    EXPECT_EQ(a.inputUnits, 36);
    for (int y = 0; y < kForgeH; ++y)
        for (int x = 0; x < kForgeW; ++x)
            EXPECT_EQ(a.h[y * kForgeW + x], (x >= 4 && x <= 7 && y >= 2 && y <= 4) ? 3 : 0);
    const ForgeState b = Forge_Start(2, 700);
    EXPECT_EQ(Sum(b), 72);
    EXPECT_EQ(b.inputUnits, 72);
    EXPECT_EQ(b.heat, 700);
    for (int y = 0; y < kForgeH; ++y)
        for (int x = 0; x < kForgeW; ++x)
            EXPECT_EQ(b.h[y * kForgeW + x], (x >= 2 && x <= 9 && y >= 2 && y <= 4) ? 3 : 0);
    EXPECT_EQ(Sum(Forge_Start(0, 1000)), 36); // clamped
    EXPECT_EQ(Sum(Forge_Start(9, 1000)), 72);
}

TEST(Forge, MetalIsConserved) {
    uint32_t rng = 12345;
    auto next = [&]() { rng = rng * 1664525u + 1013904223u; return rng >> 8; };
    ForgeState s = Forge_Start(2, kForgeHeatFull);
    uint32_t tick = 0;
    for (int i = 0; i < 5000; ++i) {
        ForgeMove m;
        m.type = (ForgeStroke)(next() % 5);
        m.x = (int8_t)((int)(next() % 14) - 1); // sometimes off the grid
        m.y = (int8_t)((int)(next() % 10) - 1);
        m.dx = (int8_t)((int)(next() % 3) - 1);
        m.dy = (int8_t)((int)(next() % 3) - 1);
        const uint32_t r = next() % 8;
        tick = r == 0 && tick > 5 ? tick - 5 : tick + r; // sometimes backwards
        m.tick = tick;
        Forge_Apply(s, m, Forge_CoolRate(2), 2, 2);
        ASSERT_EQ(Sum(s) + s.waste, s.inputUnits) << "move " << i;
        for (uint8_t v : s.h) ASSERT_LE(v, kForgeMaxHeight);
    }
}

TEST(Forge, StrokesFollowTheRules) {
    // Flatten: tie order -y, +x, +y, -x.
    {
        ForgeState s = Empty();
        s.h[4 * kForgeW + 5] = 2;
        SetInput(s);
        EXPECT_TRUE(Forge_Apply(s, Mv(ForgeStroke::Flatten, 5, 4), 0, 1, 1));
        EXPECT_EQ(s.h[3 * kForgeW + 5], 1); // north first
        EXPECT_EQ(s.h[4 * kForgeW + 5], 1);
        EXPECT_EQ(s.strokes, 1);
    }
    {
        ForgeState s = Empty();
        s.h[4 * kForgeW + 5] = 3;
        s.h[3 * kForgeW + 5] = 2; // N
        s.h[4 * kForgeW + 6] = 1; // E (tied lowest with S)
        s.h[5 * kForgeW + 5] = 1; // S
        s.h[4 * kForgeW + 4] = 3; // W
        SetInput(s);
        EXPECT_TRUE(Forge_Apply(s, Mv(ForgeStroke::Flatten, 5, 4), 0, 1, 1));
        EXPECT_EQ(s.h[4 * kForgeW + 6], 2); // east beats south
        EXPECT_EQ(s.h[5 * kForgeW + 5], 1);
        EXPECT_EQ(s.h[4 * kForgeW + 5], 2);
    }
    {   // Flatten with nowhere to go, or an empty cell: no-op.
        ForgeState s = Empty();
        s.h[4 * kForgeW + 5] = 2;
        s.h[3 * kForgeW + 5] = 3;
        s.h[4 * kForgeW + 6] = 3;
        s.h[5 * kForgeW + 5] = 3;
        s.h[4 * kForgeW + 4] = 3;
        SetInput(s);
        EXPECT_FALSE(Forge_Apply(s, Mv(ForgeStroke::Flatten, 5, 4), 0, 1, 1));
        EXPECT_FALSE(Forge_Apply(s, Mv(ForgeStroke::Flatten, 9, 7), 0, 1, 1));
        EXPECT_EQ(s.strokes, 0);
    }
    // Draw: moves one unit; blocked at max; not a unit step.
    {
        ForgeState s = Empty();
        s.h[4 * kForgeW + 5] = 2;
        s.h[4 * kForgeW + 6] = 2;
        SetInput(s);
        EXPECT_TRUE(Forge_Apply(s, Mv(ForgeStroke::Draw, 5, 4, 1, 0), 0, 1, 1));
        EXPECT_EQ(s.h[4 * kForgeW + 5], 1);
        EXPECT_EQ(s.h[4 * kForgeW + 6], 3);
        EXPECT_FALSE(Forge_Apply(s, Mv(ForgeStroke::Draw, 5, 4, 1, 0), 0, 1, 1)); // target at max
        EXPECT_FALSE(Forge_Apply(s, Mv(ForgeStroke::Draw, 5, 4, 1, 1), 0, 1, 1));  // diagonal
        EXPECT_FALSE(Forge_Apply(s, Mv(ForgeStroke::Draw, 5, 4, 0, 0), 0, 1, 1));  // no step
        EXPECT_FALSE(Forge_Apply(s, Mv(ForgeStroke::Draw, 0, 0, -1, 0), 0, 1, 1)); // off the grid
        EXPECT_FALSE(Forge_Apply(s, Mv(ForgeStroke::Draw, 5, 5, 0, 1), 0, 1, 1));  // empty source
        EXPECT_EQ(s.strokes, 1);
    }
    // Upset: from the highest neighbour (ties in order).
    {
        ForgeState s = Empty();
        s.h[3 * kForgeW + 5] = 1; // N
        s.h[4 * kForgeW + 6] = 2; // E
        s.h[5 * kForgeW + 5] = 2; // S (ties with E)
        s.h[4 * kForgeW + 4] = 1; // W
        SetInput(s);
        EXPECT_TRUE(Forge_Apply(s, Mv(ForgeStroke::Upset, 5, 4), 0, 1, 1));
        EXPECT_EQ(s.h[4 * kForgeW + 5], 1);
        EXPECT_EQ(s.h[4 * kForgeW + 6], 1); // east gave it
        EXPECT_EQ(s.h[5 * kForgeW + 5], 2);
        ForgeState full = Empty();
        full.h[4 * kForgeW + 5] = 3;
        full.h[4 * kForgeW + 6] = 3;
        SetInput(full);
        EXPECT_FALSE(Forge_Apply(full, Mv(ForgeStroke::Upset, 5, 4), 0, 1, 1)); // already max
        ForgeState lone = Empty();
        lone.h[4 * kForgeW + 5] = 1;
        SetInput(lone);
        EXPECT_FALSE(Forge_Apply(lone, Mv(ForgeStroke::Upset, 8, 4), 0, 1, 1)); // no neighbour with metal
    }
    // Punch: counts waste.
    {
        ForgeState s = Empty();
        s.h[4 * kForgeW + 5] = 2;
        SetInput(s);
        EXPECT_TRUE(Forge_Apply(s, Mv(ForgeStroke::Punch, 5, 4), 0, 1, 1));
        EXPECT_EQ(s.h[4 * kForgeW + 5], 1);
        EXPECT_EQ(s.waste, 1);
        EXPECT_EQ(Sum(s) + s.waste, s.inputUnits);
        EXPECT_FALSE(Forge_Apply(s, Mv(ForgeStroke::Punch, 0, 0), 0, 1, 1));
    }
}

TEST(Forge, HeatIsTheClock) {
    // Under working heat a stroke does nothing.
    ForgeState cold = Forge_Start(1, kForgeWorkingHeat - 1);
    EXPECT_FALSE(Forge_Apply(cold, Mv(ForgeStroke::Punch, 4, 2), 0, 1, 1));
    EXPECT_EQ(cold.strokes, 0);
    // Exactly working heat works.
    ForgeState warm = Forge_Start(1, kForgeWorkingHeat);
    EXPECT_TRUE(Forge_Apply(warm, Mv(ForgeStroke::Punch, 4, 2), 0, 1, 1));
    EXPECT_EQ(warm.heat, kForgeWorkingHeat - kForgeHeatPerStroke);
    // Time cools at the metal's rate: 1 s at tier 1 = 20 permille; the stroke then costs 25.
    ForgeState s = Forge_Start(1, kForgeHeatFull);
    EXPECT_TRUE(Forge_Apply(s, Mv(ForgeStroke::Punch, 4, 2, 0, 0, 20), Forge_CoolRate(1), 1, 1));
    EXPECT_EQ(s.heat, 1000 - 20 - kForgeHeatPerStroke);
    EXPECT_EQ(s.tick, 20u);
    // A move with an earlier tick counts as the same tick (no negative cooling).
    const int before = s.heat;
    EXPECT_TRUE(Forge_Apply(s, Mv(ForgeStroke::Punch, 5, 2, 0, 0, 5), Forge_CoolRate(1), 1, 1));
    EXPECT_EQ(s.heat, before - kForgeHeatPerStroke);
    EXPECT_EQ(s.tick, 20u);
    // Cold by the clock: 500 permille, 101 ticks at 20/s cool 101 -> 399.
    ForgeState c = Forge_Start(1, 500);
    EXPECT_FALSE(Forge_Apply(c, Mv(ForgeStroke::Punch, 4, 2, 0, 0, 101), 20, 1, 1));
    EXPECT_EQ(c.heat, 399);
    // Reheat restores the heat and counts, but it is not a stroke.
    EXPECT_FALSE(Forge_Apply(c, Mv(ForgeStroke::Reheat, 0, 0, 0, 0, 101), 20, 1, 1));
    EXPECT_EQ(c.heat, kForgeHeatFull);
    EXPECT_EQ(c.reheats, 1);
    EXPECT_EQ(c.strokes, 0);
    EXPECT_TRUE(Forge_Apply(c, Mv(ForgeStroke::Punch, 4, 2, 0, 0, 102), 20, 1, 1));
    // A hammer below the metal's tier does nothing.
    ForgeState weak = Forge_Start(1, kForgeHeatFull);
    EXPECT_FALSE(Forge_Apply(weak, Mv(ForgeStroke::Punch, 4, 2), 0, 1, 2));
    EXPECT_EQ(weak.strokes, 0);
    EXPECT_EQ(Sum(weak), 36);
    EXPECT_TRUE(Forge_Apply(weak, Mv(ForgeStroke::Punch, 4, 2), 0, 2, 2));
    // Rates by tier.
    EXPECT_EQ(Forge_CoolRate(1), 20);
    EXPECT_EQ(Forge_CoolRate(2), 17);
    EXPECT_EQ(Forge_CoolRate(3), 14);
    EXPECT_EQ(Forge_CoolRate(5), 14);
}

TEST(Forge, Scoring) {
    // A pattern equal to the start block with no moves: a perfect forge.
    ForgePattern block;
    block.id = "block";
    block.ingots = 1;
    block.par = 10;
    block.h = Forge_Start(1, kForgeHeatFull).h;
    const ForgeState fresh = Forge_Start(1, kForgeHeatFull);
    EXPECT_EQ(Forge_Quality(fresh, block, 130, 130), 130);
    EXPECT_STREQ(Forge_Grade(130), "Masterwork");
    EXPECT_EQ(Forge_Quality(fresh, block, 118, 130), 118);
    EXPECT_EQ(Forge_Quality(fresh, block, 130, 100), 100);
    // Punching 3 units (3 cells 3 -> 2): match 9 of 12 = 750, minus 30 * 3 waste over, no extra strokes (par 10): 660.
    ForgeState s = fresh;
    EXPECT_TRUE(Forge_Apply(s, Mv(ForgeStroke::Punch, 4, 2), 0, 1, 1));
    EXPECT_TRUE(Forge_Apply(s, Mv(ForgeStroke::Punch, 5, 2), 0, 1, 1));
    EXPECT_TRUE(Forge_Apply(s, Mv(ForgeStroke::Punch, 6, 2), 0, 1, 1));
    EXPECT_EQ(Forge_Quality(s, block, 130, 130), 75 + 55 * 660 / 1000);
    // Extra strokes and reheats cost 20 each: par 1 with 3 strokes = 2 extra = 40, one reheat 20.
    block.par = 1;
    s.reheats = 1;
    EXPECT_EQ(Forge_Quality(s, block, 130, 130), 75 + 55 * (750 - 90 - 40 - 20) / 1000);
    // The worst case floors at 75.
    ForgeState ruined = fresh;
    for (uint8_t& v : ruined.h) v = 0;
    ruined.waste = 36;
    ruined.strokes = 50;
    EXPECT_EQ(Forge_Quality(ruined, block, 130, 130), 75);
    // Waste the pattern needs is free: a pattern of 30 units from 36 may punch 6.
    ForgePattern small = block;
    small.par = 100;
    for (int x = 4; x <= 7; ++x) small.h[2 * kForgeW + x] = 0; // drops the top row: 24 units left
    ForgeState t = fresh;
    for (int x = 4; x <= 7; ++x)
        for (int i = 0; i < 3; ++i) Forge_Apply(t, Mv(ForgeStroke::Punch, x, 2), 0, 1, 1);
    EXPECT_EQ(t.waste, 12);
    EXPECT_EQ(small.units(), 24);
    EXPECT_EQ(Forge_Quality(t, small, 130, 130), 130); // 12 punched = 36 - 24: no waste over; strokes 12 <= par
    // Grades.
    EXPECT_STREQ(Forge_Grade(75), "Crude");
    EXPECT_STREQ(Forge_Grade(89), "Crude");
    EXPECT_STREQ(Forge_Grade(90), "Common");
    EXPECT_STREQ(Forge_Grade(109), "Common");
    EXPECT_STREQ(Forge_Grade(110), "Fine");
    EXPECT_STREQ(Forge_Grade(124), "Fine");
    EXPECT_STREQ(Forge_Grade(125), "Masterwork");
    EXPECT_EQ(Forge_MetalTier("bronze"), 1);
    EXPECT_EQ(Forge_MetalTier("rough_bronze_mace_head"), 1);
    EXPECT_EQ(Forge_MetalTier("iron_bloom"), 2);
    EXPECT_EQ(Forge_MetalTier("iron_ingot"), 2);
    EXPECT_EQ(Forge_MetalTier("steel_hammer_head"), 3);
    EXPECT_EQ(Forge_MetalTier("tungsten_steel_ingot"), 4);
    EXPECT_EQ(Forge_MetalTier("weapon_sword"), 0);
    EXPECT_EQ(Forge_MetalTier(""), 0);
    EXPECT_EQ(Forge_JobTier("steel_pick_head"), 3);
    EXPECT_EQ(Forge_JobTier("steel_hammer_head"), 2) << "a bronze or iron hammer makes the steel one";
    EXPECT_EQ(Forge_JobTier("rough_bronze_hammer_head"), 0);
    EXPECT_EQ(Forge_JobTier("iron_ingot"), 2);
    EXPECT_EQ(Forge_QuickQuality(), 95);
    EXPECT_STREQ(Forge_Grade(Forge_QuickQuality()), "Common");
}

TEST(Forge, PatternsParse) {
    const std::vector<ForgePattern>& all = Forge_AllPatterns();
    EXPECT_EQ(all.size(), 8u);
    for (const char* id : { "pick_head", "axe_head", "shovel_head", "hoe_head", "blade", "spearhead", "mace_head", "hammer_head" }) {
        const ForgePattern* p = Forge_Pattern(id);
        ASSERT_NE(p, nullptr) << id;
        EXPECT_EQ(p->units(), p->ingots * kForgeUnitsPerIngot) << id;
        EXPECT_GT(p->par, 0) << id;
    }
    EXPECT_EQ(Forge_Pattern("blade")->ingots, 2);
    EXPECT_EQ(Forge_Pattern("mace_head")->ingots, 2);
    EXPECT_EQ(Forge_Pattern("spearhead")->ingots, 1);
    EXPECT_EQ(Forge_Pattern("nope"), nullptr);

    // Bad input errors with line numbers and drops the block.
    const std::string row36 = "....3333....\n....3333....\n....3333....\n";
    const std::string empty = "............\n";
    auto rows = [&](int nonEmptyBefore) { // 3 full rows (36 units) padded to 8 rows
        std::string s;
        for (int i = 0; i < nonEmptyBefore; ++i) s += empty;
        s += row36;
        for (int i = 0; i < 5 - nonEmptyBefore; ++i) s += empty;
        return s;
    };
    {
        std::vector<std::string> errors;
        auto v = Forge_ParsePatterns("pattern ok | ingots 1 | par 5\n" + rows(2) + "end\n", &errors);
        EXPECT_TRUE(errors.empty());
        ASSERT_EQ(v.size(), 1u);
        EXPECT_EQ(v[0].par, 5);
    }
    struct Bad { const char* what; std::string text; };
    const Bad bads[] = {
        { "too few rows", "pattern a | ingots 1 | par 5\n" + row36 + "end\n" },
        { "too many rows", "pattern a | ingots 1 | par 5\n" + rows(2) + empty + "end\n" },
        { "short row", "pattern a | ingots 1 | par 5\n...\n" + rows(2) + "end\n" },
        { "bad digit", "pattern a | ingots 1 | par 5\n....3393....\n" + rows(1) + "end\n" },
        { "units mismatch", "pattern a | ingots 2 | par 5\n" + rows(2) + "end\n" },
        { "missing par", "pattern a | ingots 1\n" + rows(2) + "end\n" },
        { "missing end", "pattern a | ingots 1 | par 5\n" + rows(2) },
        { "stray line", "hello\n" },
    };
    for (const Bad& b : bads) {
        std::vector<std::string> errors;
        auto v = Forge_ParsePatterns(b.text, &errors);
        EXPECT_TRUE(v.empty()) << b.what;
        ASSERT_FALSE(errors.empty()) << b.what;
        EXPECT_EQ(errors[0].rfind("line ", 0), 0u) << b.what << ": " << errors[0];
    }
    // A bad block does not take its neighbours down.
    std::vector<std::string> errors;
    auto v = Forge_ParsePatterns("pattern bad | ingots 1\n" + rows(2) + "end\npattern good | ingots 1 | par 3\n" + rows(2) + "end\n", &errors);
    ASSERT_EQ(v.size(), 1u);
    EXPECT_EQ(v[0].id, "good");
    EXPECT_FALSE(errors.empty());
}

TEST(Forge, SolverReachesEveryPattern) {
    ASSERT_EQ(Forge_AllPatterns().size(), 8u);
    for (const ForgePattern& p : Forge_AllPatterns()) {
        const int tier = 2;
        const std::vector<ForgeMove> moves = Forge_Solve(p, tier);
        EXPECT_LT(moves.size(), 2000u) << p.id;
        ForgeState s;
        const ForgeResult r = Forge_Replay(p, kForgeHeatFull, tier, tier, 130, 130, moves, &s);
        EXPECT_TRUE(r.valid) << p.id;
        EXPECT_EQ(s.h, p.h) << p.id << " not reached";
        EXPECT_EQ(r.waste, 0) << p.id;
        EXPECT_GE(r.q, 120) << p.id << " q " << r.q << " strokes " << r.strokes << " par " << p.par << " reheats "
                            << r.reheats;
        EXPECT_LE(r.strokes, p.par) << p.id;
    }
    // Another metal tier also gets there (the heat plan differs).
    for (const ForgePattern& p : Forge_AllPatterns()) {
        ForgeState s;
        Forge_Replay(p, kForgeHeatFull, 1, 1, 130, 130, Forge_Solve(p, 1), &s);
        EXPECT_EQ(s.h, p.h) << p.id << " tier 1";
    }
}

TEST(Forge, ReplayIsDeterministic) {
    const ForgePattern* p = Forge_Pattern("axe_head");
    ASSERT_NE(p, nullptr);
    const std::vector<ForgeMove> moves = Forge_Solve(*p, 2);
    ForgeState a, b;
    const ForgeResult ra = Forge_Replay(*p, kForgeHeatFull, 2, 2, 130, 130, moves, &a);
    const ForgeResult rb = Forge_Replay(*p, kForgeHeatFull, 2, 2, 130, 130, moves, &b);
    EXPECT_EQ(a.h, b.h);
    EXPECT_EQ(a.heat, b.heat);
    EXPECT_EQ(a.tick, b.tick);
    EXPECT_EQ(ra.q, rb.q);
    EXPECT_EQ(ra.strokes, rb.strokes);
    EXPECT_EQ(ra.reheats, rb.reheats);
    EXPECT_EQ(ra.waste, rb.waste);
    EXPECT_TRUE(ra.valid);
    // A hostile client: too many moves, or ticks past 20 minutes.
    std::vector<ForgeMove> many(4097, Mv(ForgeStroke::Reheat, 0, 0));
    const ForgeResult rm = Forge_Replay(*p, kForgeHeatFull, 2, 2, 130, 130, many);
    EXPECT_FALSE(rm.valid);
    EXPECT_EQ(rm.q, 75);
    std::vector<ForgeMove> okMoves(4096, Mv(ForgeStroke::Reheat, 0, 0));
    EXPECT_TRUE(Forge_Replay(*p, kForgeHeatFull, 2, 2, 130, 130, okMoves).valid);
    std::vector<ForgeMove> late = { Mv(ForgeStroke::Punch, 4, 2, 0, 0, 20u * 60u * 20u + 1u) };
    const ForgeResult rl = Forge_Replay(*p, kForgeHeatFull, 2, 2, 130, 130, late);
    EXPECT_FALSE(rl.valid);
    EXPECT_EQ(rl.q, 75);
}
#endif
