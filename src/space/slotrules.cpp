#include "slotrules.h"
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace {
    enum class Drag { None, Split, One, Gather, Shift };
    InvSlot* held = nullptr;
    std::vector<SlotGrid*> grids;
    Drag drag = Drag::None;
    InvSlot dragStack;                                  // Split: the held stack when the drag began
    std::vector<std::pair<SlotGrid*, int>> visited;     // the slots dragged across, in order
    std::string gatherItem;                             // Gather: the item picked up by the press

    bool Visited(SlotGrid* g, int i) { return std::find(visited.begin(), visited.end(), std::make_pair(g, i)) != visited.end(); }
    void Settle() { if (held->centi <= 0) *held = InvSlot{}; }
    int64_t One() { return std::min<int64_t>(100, held->centi); } // one unit (or the fraction left)

    // Up to centi of the held stack into slot i (what fits); the hand loses what went.
    int64_t PutHeld(SlotGrid* g, int i, int64_t centi) {
        if (held->empty()) return 0;
        const int64_t n = g->infinite() ? held->centi : std::min(centi, g->room(i, held->id));
        if (n <= 0) return 0;
        const int64_t went = g->put(i, *held, n);
        held->centi -= went;
        Settle();
        return went;
    }
    // The held stack and slot i's change places (a click with another item). A container answers a tick later: its stack
    // comes onto the emptied hand through the cursor sink, after the put that follows it in the same queue.
    void Swap(SlotGrid* g, int i) {
        const InvSlot at = g->at(i);
        if (g->infinite()) { *held = InvSlot{}; return; } // the Catalog: a held stack put on it is deleted
        const InvSlot mine = *held;
        *held = InvSlot{};
        g->take(i, at.centi, *held);
        const int64_t went = g->put(i, mine, mine.centi);
        if (went < mine.centi && held->empty()) { *held = mine; held->centi -= went; } // (didn't fit: the rest stays held)
    }
    // A click with a held stack: put it all down / merge (what doesn't fit stays held), or swap with another item.
    void ClickHeld(SlotGrid* g, int i) {
        const InvSlot at = g->at(i);
        if (at.empty() || at.id == held->id || g->infinite()) PutHeld(g, i, held->centi);
        else Swap(g, i);
    }
    // Double-click: every stack of the held item onto the hand, up to a full stack; partial stacks first, as Minecraft does.
    void Gather() {
        if (held->empty()) return;
        const std::string id = held->id;
        int64_t want = 0;
        for (SlotGrid* g : grids) if (g && !g->infinite()) { want = g->capOf(id) - held->centi; break; }
        for (int pass = 0; pass < 2 && want > 0; ++pass)
            for (SlotGrid* g : grids) {
                if (!g || g->infinite()) continue;
                for (int j = 0; j < g->size() && want > 0; ++j) {
                    const InvSlot s = g->at(j);
                    if (s.empty() || s.id != id || (pass == 0) == (s.centi >= g->capOf(id))) continue;
                    const int64_t k = std::min(want, s.centi);
                    g->take(j, k, *held);
                    want -= k; // (a container's comes next tick: counted now, so the hand never overfills)
                }
            }
    }
}

void SlotRules_SetHeld(InvSlot* h) { held = h; }
void SlotRules_SetGrids(const std::vector<SlotGrid*>& g) { grids = g; }
void SlotRules_Cancel() { drag = Drag::None; visited.clear(); }

bool SlotRules_Frame(const SlotInput& in, SlotGrid* g, int i) {
    if (!held) return false;
    if (g && (i < 0 || i >= g->size())) g = nullptr;
    bool did = false;
    const InvSlot at = g ? g->at(i) : InvSlot{};
    // Keys over a slot (an empty hand): 1-0 and the off hand swap with that cell; Q drops one, Ctrl+Q the stack.
    if (g && held->empty() && drag == Drag::None) {
        if (in.hotbarKey >= 0 && in.hotbarKey < Inventory::kHotbar) { g->swapCell(i, Inventory::kGrid + in.hotbarKey); did = true; }
        if (in.offHandKey) { g->swapCell(i, Inventory::kOffHand); did = true; }
        if (in.dropKey && !at.empty() && !g->infinite()) { g->drop(i, in.ctrl ? at.centi : std::min<int64_t>(100, at.centi)); did = true; }
        // The wheel (Mouse Tweaks): down sends one (Shift: the stack) to the other side, up pulls one in.
        const int notches = (int)std::lround(std::fabs(in.wheel));
        if (notches > 0 && !g->infinite()) {
            const bool out = (in.wheel < 0.0f) != in.wheelInvert;
            for (int k = 0; k < notches; ++k) {
                const InvSlot now = g->at(i);
                if (out) { if (!now.empty()) g->moveOut(i, in.shift ? now.centi : std::min<int64_t>(100, now.centi)); }
                else g->pullIn(i, 100);
            }
            did = true;
        }
    }
    // Middle-click in creative: a full stack of that item.
    if (in.mPress && g && in.creative && held->empty() && !at.empty()) { *held = at; held->centi = g->capOf(at.id); did = true; }
    // Left press
    if (in.lPress && g) {
        did = true;
        if (in.lDouble && !held->empty()) { Gather(); drag = Drag::None; visited.clear(); }
        else if (in.shift && held->empty()) {
            if (!at.empty() && !g->infinite()) g->moveOut(i, at.centi);
            if (!g->infinite()) { drag = Drag::Shift; visited = { { g, i } }; }
        } else if (held->empty()) {
            if (!at.empty()) {
                g->take(i, g->infinite() ? g->capOf(at.id) : at.centi, *held);
                if (!g->infinite()) { drag = Drag::Gather; gatherItem = at.id; visited = { { g, i } }; }
            }
        } else { drag = Drag::Split; dragStack = *held; visited = { { g, i } }; } // a click or a split: decided on release
    }
    // Right press
    if (in.rPress && g) {
        did = true;
        if (held->empty()) {
            if (!at.empty()) {
                const int64_t units = at.centi / 100;
                const int64_t half = g->infinite() ? 100 : units >= 1 ? ((units + 1) / 2) * 100 : at.centi; // the larger half
                g->take(i, half, *held);
            }
        } else if (!at.empty() && at.id != held->id && !g->infinite()) Swap(g, i);
        else { PutHeld(g, i, One()); drag = Drag::One; visited = { { g, i } }; }
    }
    // Dragging across new slots
    if (g && drag != Drag::None && !Visited(g, i)) {
        switch (drag) {
            case Drag::Split:
                if (in.lDown && !g->infinite() && g->room(i, dragStack.id) > 0) visited.push_back({ g, i });
                break;
            case Drag::One:
                if (in.rDown && !held->empty() && g->room(i, held->id) > 0 && !g->infinite()) { PutHeld(g, i, One()); visited.push_back({ g, i }); did = true; }
                break;
            case Drag::Gather:
                if (in.lDown && !at.empty() && at.id == gatherItem && !g->infinite() && (held->empty() || held->id == gatherItem)) {
                    const int64_t want = g->capOf(gatherItem) - (held->empty() ? 0 : held->centi);
                    if (want > 0) { g->take(i, std::min(want, at.centi), *held); did = true; }
                    visited.push_back({ g, i });
                }
                break;
            case Drag::Shift:
                if (in.lDown && !at.empty() && !g->infinite()) { g->moveOut(i, at.centi); visited.push_back({ g, i }); did = true; }
                break;
            default: break;
        }
    }
    // Releases
    if (in.lRelease && drag != Drag::None && drag != Drag::One) {
        if (drag == Drag::Split && !held->empty()) {
            if (visited.size() <= 1) { if (!visited.empty()) ClickHeld(visited[0].first, visited[0].second); }
            else { // split evenly in whole units (fewer units than slots: one each to the first)
                const int64_t units = dragStack.centi / 100;
                const int64_t share = std::max<int64_t>(1, units / (int64_t)visited.size()) * 100;
                for (const auto& [vg, vi] : visited) { if (held->empty()) break; PutHeld(vg, vi, std::min(share, held->centi)); }
            }
            did = true;
        }
        drag = Drag::None; visited.clear();
    }
    if (in.rRelease && drag == Drag::One) { drag = Drag::None; visited.clear(); }
    if (!in.lDown && !in.rDown && !in.lRelease && !in.rRelease && drag != Drag::None && !in.lPress && !in.rPress) { drag = Drag::None; visited.clear(); } // (a release missed)
    return did;
}

std::vector<SlotDragPreview> SlotRules_DragPreview() {
    std::vector<SlotDragPreview> out;
    if (drag != Drag::Split || visited.size() < 2) return out;
    const int64_t share = std::max<int64_t>(1, (dragStack.centi / 100) / (int64_t)visited.size()) * 100;
    int64_t left = dragStack.centi;
    for (const auto& [g, i] : visited) { const int64_t k = std::min(share, left); if (k <= 0) break; out.push_back({ g, i, k }); left -= k; }
    return out;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#862): every rule on two in-memory grids (A the "inventory", B the "chest"; each one's other side is the
// other) and a Catalog. Stacks of 64 (the grids' own cap: no item tables needed).
#include "unit_helpers.h"

namespace {
struct TestGrid : SlotGrid {
    std::vector<InvSlot> cells;
    TestGrid* other = nullptr;
    bool catalog = false;
    std::vector<std::pair<int, int64_t>> dropped;
    std::vector<std::pair<int, int>> swapped;
    explicit TestGrid(int n) : cells(n) {}
    int size() const override { return (int)cells.size(); }
    InvSlot at(int i) const override { return cells[i]; }
    int64_t capOf(const std::string&) const override { return 6400; }
    int64_t room(int i, const std::string& id) const override {
        if (catalog) return 0;
        return cells[i].empty() ? 6400 : cells[i].id == id ? 6400 - cells[i].centi : 0;
    }
    int64_t put(int i, const InvSlot& h, int64_t centi) override {
        if (catalog) return h.centi;
        const int64_t k = std::min(centi, room(i, h.id));
        if (k <= 0) return 0;
        if (cells[i].empty()) { cells[i] = h; cells[i].centi = 0; }
        cells[i].centi += k;
        return k;
    }
    int64_t take(int i, int64_t centi, InvSlot& h) override {
        const InvSlot s = cells[i];
        if (s.empty() || (!h.empty() && h.id != s.id)) return 0;
        const int64_t k = catalog ? centi : std::min(centi, s.centi);
        if (h.empty()) { h = s; h.centi = 0; }
        h.centi += k;
        if (!catalog) { cells[i].centi -= k; if (cells[i].centi <= 0) cells[i] = InvSlot{}; }
        return k;
    }
    void moveOut(int i, int64_t centi) override { // onto the same item first, then the first empty slot
        const InvSlot s = cells[i];
        int64_t left = std::min(centi, s.centi);
        for (int pass = 0; pass < 2 && left > 0; ++pass)
            for (int j = 0; j < other->size() && left > 0; ++j) {
                InvSlot& c = other->cells[j];
                if (pass == 0 ? (c.empty() || c.id != s.id) : !c.empty()) continue;
                const int64_t k = std::min(left, other->room(j, s.id));
                if (c.empty()) { c = s; c.centi = 0; }
                c.centi += k; left -= k; cells[i].centi -= k;
            }
        if (cells[i].centi <= 0) cells[i] = InvSlot{};
    }
    void pullIn(int i, int64_t centi) override {
        if (cells[i].empty()) return;
        for (InvSlot& o : other->cells)
            if (!o.empty() && o.id == cells[i].id) {
                const int64_t k = std::min({ centi, o.centi, room(i, cells[i].id) });
                o.centi -= k; cells[i].centi += k;
                if (o.centi <= 0) o = InvSlot{};
                return;
            }
    }
    void swapCell(int i, int cell) override { swapped.push_back({ i, cell }); }
    void drop(int i, int64_t centi) override { dropped.push_back({ i, centi }); cells[i].centi -= centi; if (cells[i].centi <= 0) cells[i] = InvSlot{}; }
    bool infinite() const override { return catalog; }
};

InvSlot S(const char* id, int n) { InvSlot s; s.id = id; s.centi = (int64_t)n * 100; return s; }
int N(const InvSlot& s, const char* id) { return !s.empty() && s.id == id ? (int)(s.centi / 100) : 0; }

struct SlotRulesTest : ::testing::Test {
    TestGrid A{ 9 }, B{ 9 }, C{ 3 };
    InvSlot hand;
    void SetUp() override {
        A.other = &B; B.other = &A;
        C.catalog = true; C.cells[0] = S("diamond", 1);
        SlotRules_SetHeld(&hand);
        SlotRules_SetGrids({ &A, &B, &C });
        SlotRules_Cancel();
    }
    void TearDown() override { SlotRules_SetHeld(nullptr); SlotRules_SetGrids({}); }
    // A press and its release on one slot.
    void Click(SlotGrid* g, int i, bool right = false, bool shift = false, bool dbl = false) {
        SlotInput p; p.shift = shift; (right ? p.rPress : p.lPress) = true; (right ? p.rDown : p.lDown) = true; p.lDouble = dbl && !right;
        SlotRules_Frame(p, g, i);
        SlotInput r; r.shift = shift; (right ? r.rRelease : r.lRelease) = true;
        SlotRules_Frame(r, g, i);
    }
    // A press on the first slot, held across the rest, released on the last.
    void Drag(SlotGrid* g, const std::vector<int>& slots, bool right = false, bool shift = false) {
        SlotInput p; p.shift = shift; (right ? p.rPress : p.lPress) = true; (right ? p.rDown : p.lDown) = true;
        SlotRules_Frame(p, g, slots[0]);
        for (size_t k = 1; k < slots.size(); ++k) { SlotInput m; m.shift = shift; (right ? m.rDown : m.lDown) = true; SlotRules_Frame(m, g, slots[k]); }
        SlotInput r; r.shift = shift; (right ? r.rRelease : r.lRelease) = true;
        SlotRules_Frame(r, g, slots.back());
    }
    void Key(SlotGrid* g, int i, const SlotInput& in) { SlotRules_Frame(in, g, i); }
};
}

// -- Minecraft

TEST_F(SlotRulesTest, LeftClickPicksUpAndPutsDownTheWholeStack) {
    A.cells[0] = S("stone", 10);
    Click(&A, 0);
    EXPECT_EQ(N(hand, "stone"), 10); EXPECT_TRUE(A.cells[0].empty());
    Click(&A, 1);
    EXPECT_EQ(N(A.cells[1], "stone"), 10); EXPECT_TRUE(hand.empty());
}

TEST_F(SlotRulesTest, LeftClickMergesKeepingWhatDoesNotFit) {
    hand = S("stone", 10); A.cells[2] = S("stone", 60);
    Click(&A, 2);
    EXPECT_EQ(N(A.cells[2], "stone"), 64); EXPECT_EQ(N(hand, "stone"), 6);
}

TEST_F(SlotRulesTest, LeftClickOnAnotherItemSwaps) {
    hand = S("stone", 6); A.cells[3] = S("oak_log", 3);
    Click(&A, 3);
    EXPECT_EQ(N(A.cells[3], "stone"), 6); EXPECT_EQ(N(hand, "oak_log"), 3);
}

TEST_F(SlotRulesTest, RightClickTakesTheLargerHalf) {
    A.cells[0] = S("stone", 9);
    Click(&A, 0, true);
    EXPECT_EQ(N(hand, "stone"), 5); EXPECT_EQ(N(A.cells[0], "stone"), 4);
}

TEST_F(SlotRulesTest, RightClickHoldingAStackPutsOneDown) {
    hand = S("stone", 5);
    Click(&A, 1, true);
    EXPECT_EQ(N(A.cells[1], "stone"), 1); EXPECT_EQ(N(hand, "stone"), 4);
}

TEST_F(SlotRulesTest, RightClickOnAnotherItemSwaps) {
    hand = S("stone", 4); A.cells[2] = S("oak_log", 2);
    Click(&A, 2, true);
    EXPECT_EQ(N(hand, "oak_log"), 2); EXPECT_EQ(N(A.cells[2], "stone"), 4);
}

TEST_F(SlotRulesTest, ShiftClickQuickMovesOntoTheSameItemFirst) {
    A.cells[0] = S("stone", 10); B.cells[4] = S("stone", 60);
    Click(&A, 0, false, true);
    EXPECT_TRUE(A.cells[0].empty());
    EXPECT_EQ(N(B.cells[4], "stone"), 64);
    EXPECT_EQ(N(B.cells[0], "stone"), 6);
}

TEST_F(SlotRulesTest, DoubleClickGathersPartialStacksFirstUpToAStack) {
    A.cells[0] = S("stone", 10); A.cells[3] = S("stone", 64); B.cells[1] = S("stone", 20);
    Click(&A, 0);
    Click(&A, 0, false, false, true);
    EXPECT_EQ(N(hand, "stone"), 64);
    EXPECT_TRUE(B.cells[1].empty()) << "the partial stack goes first";
    EXPECT_EQ(N(A.cells[3], "stone"), 30) << "then the full one, only what fits";
}

TEST_F(SlotRulesTest, LeftDragSplitsEvenlyKeepingTheRest) {
    hand = S("stone", 10);
    Drag(&A, { 0, 1, 2 });
    EXPECT_EQ(N(A.cells[0], "stone"), 3); EXPECT_EQ(N(A.cells[1], "stone"), 3); EXPECT_EQ(N(A.cells[2], "stone"), 3);
    EXPECT_EQ(N(hand, "stone"), 1);
}

TEST_F(SlotRulesTest, LeftDragSplitSkipsSlotsWithAnotherItem) {
    hand = S("stone", 10); A.cells[1] = S("oak_log", 1);
    Drag(&A, { 0, 1, 2 });
    EXPECT_EQ(N(A.cells[0], "stone"), 5); EXPECT_EQ(N(A.cells[1], "oak_log"), 1); EXPECT_EQ(N(A.cells[2], "stone"), 5);
}

TEST_F(SlotRulesTest, LeftDragPreviewShowsEachShare) {
    hand = S("stone", 7);
    SlotInput p; p.lPress = p.lDown = true; SlotRules_Frame(p, &A, 0);
    SlotInput m; m.lDown = true; SlotRules_Frame(m, &A, 1);
    const auto preview = SlotRules_DragPreview();
    ASSERT_EQ(preview.size(), 2u);
    EXPECT_EQ(preview[0].centi, 300); EXPECT_EQ(preview[1].centi, 300);
    EXPECT_TRUE(A.cells[0].empty()) << "nothing moves before the release";
}

TEST_F(SlotRulesTest, RightDragPutsOneIntoEachSlot) {
    hand = S("stone", 10);
    Drag(&A, { 0, 1, 2, 3 }, true);
    for (int k = 0; k < 4; ++k) EXPECT_EQ(N(A.cells[k], "stone"), 1);
    EXPECT_EQ(N(hand, "stone"), 6);
}

TEST_F(SlotRulesTest, NumberAndOffHandKeysSwapWithThoseCells) {
    A.cells[0] = S("stone", 5);
    SlotInput k; k.hotbarKey = 3; Key(&A, 0, k);
    ASSERT_FALSE(A.swapped.empty());
    EXPECT_EQ(A.swapped.back(), std::make_pair(0, Inventory::kGrid + 3));
    SlotInput o; o.offHandKey = true; Key(&A, 0, o);
    EXPECT_EQ(A.swapped.back(), std::make_pair(0, Inventory::kOffHand));
}

TEST_F(SlotRulesTest, QDropsOneAndCtrlQTheStack) {
    A.cells[0] = S("stone", 5);
    SlotInput q; q.dropKey = true; Key(&A, 0, q);
    EXPECT_EQ(N(A.cells[0], "stone"), 4); EXPECT_EQ(A.dropped.back().second, 100);
    SlotInput cq; cq.dropKey = true; cq.ctrl = true; Key(&A, 0, cq);
    EXPECT_TRUE(A.cells[0].empty()); EXPECT_EQ(A.dropped.back().second, 400);
}

TEST_F(SlotRulesTest, MiddleClickGivesAFullStackInCreativeOnly) {
    A.cells[0] = S("stone", 5);
    SlotInput m; m.mPress = true; Key(&A, 0, m);
    EXPECT_TRUE(hand.empty());
    m.creative = true; Key(&A, 0, m);
    EXPECT_EQ(N(hand, "stone"), 64); EXPECT_EQ(N(A.cells[0], "stone"), 5);
}

// -- the Catalog

TEST_F(SlotRulesTest, CatalogGivesFreeStacksAndDeletesWhatIsPutOnIt) {
    Click(&C, 0);
    EXPECT_EQ(N(hand, "diamond"), 64); EXPECT_EQ(N(C.cells[0], "diamond"), 1);
    Click(&C, 0);
    EXPECT_TRUE(hand.empty());
    Click(&C, 0, true);
    EXPECT_EQ(N(hand, "diamond"), 1);
}

// -- Mouse Tweaks

TEST_F(SlotRulesTest, LeftDragAfterAPickUpGathersThatItem) {
    A.cells[0] = S("stone", 10); A.cells[1] = S("stone", 20); A.cells[2] = S("oak_log", 4); A.cells[3] = S("stone", 50);
    Drag(&A, { 0, 1, 2, 3 });
    EXPECT_EQ(N(hand, "stone"), 64);
    EXPECT_EQ(N(A.cells[2], "oak_log"), 4);
    EXPECT_EQ(N(A.cells[3], "stone"), 16);
}

TEST_F(SlotRulesTest, ShiftDragQuickMovesEverySlotCrossed) {
    A.cells[0] = S("stone", 10); A.cells[1] = S("oak_log", 3); A.cells[2] = S("sand", 2);
    Drag(&A, { 0, 1, 2 }, false, true);
    EXPECT_TRUE(A.cells[0].empty()); EXPECT_TRUE(A.cells[1].empty()); EXPECT_TRUE(A.cells[2].empty());
    EXPECT_EQ(N(B.cells[1], "oak_log"), 3);
}

TEST_F(SlotRulesTest, WheelDownMovesOneOrTheStackWithShift) {
    A.cells[0] = S("stone", 10);
    SlotInput w; w.wheel = -1.0f; Key(&A, 0, w);
    EXPECT_EQ(N(A.cells[0], "stone"), 9); EXPECT_EQ(N(B.cells[0], "stone"), 1);
    w.shift = true; Key(&A, 0, w);
    EXPECT_TRUE(A.cells[0].empty()); EXPECT_EQ(N(B.cells[0], "stone"), 10);
}

TEST_F(SlotRulesTest, WheelUpPullsOnePerNotchAndTheDirectionCanBeSwapped) {
    A.cells[5] = S("stone", 1); B.cells[0] = S("stone", 10);
    SlotInput w; w.wheel = 2.0f; Key(&A, 5, w);
    EXPECT_EQ(N(A.cells[5], "stone"), 3); EXPECT_EQ(N(B.cells[0], "stone"), 8);
    SlotInput inv; inv.wheel = 1.0f; inv.wheelInvert = true; Key(&A, 5, inv);
    EXPECT_EQ(N(A.cells[5], "stone"), 2);
}

TEST_F(SlotRulesTest, NoHeldStackNoGridDoesNothing) {
    SlotRules_SetHeld(nullptr);
    SlotInput p; p.lPress = true;
    EXPECT_FALSE(SlotRules_Frame(p, &A, 0));
}
#endif
