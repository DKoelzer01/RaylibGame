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

// ---- --test-slot-rules: every rule on two in-memory grids (A the "inventory", B the "chest"; each one's other side is
// the other), plus a Catalog. Stacks of 64.
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
        void moveOut(int i, int64_t centi) override {
            InvSlot s = cells[i];
            int64_t left = std::min(centi, s.centi);
            for (int pass = 0; pass < 2 && left > 0; ++pass)
                for (int j = 0; j < other->size() && left > 0; ++j) {
                    InvSlot& c = other->cells[j];
                    if (pass == 0 ? (c.empty() || c.id != s.id) : !c.empty()) continue;
                    const int64_t k = std::min(left, other->room(j, s.id));
                    if (c.empty()) { c = s; c.centi = 0; }
                    c.centi += k; left -= k;
                    cells[i].centi -= k;
                }
            if (cells[i].centi <= 0) cells[i] = InvSlot{};
        }
        void pullIn(int i, int64_t centi) override {
            if (cells[i].empty()) return;
            for (int j = 0; j < other->size(); ++j)
                if (other->cells[j].id == cells[i].id && !other->cells[j].empty()) {
                    const int64_t k = std::min({ centi, other->cells[j].centi, room(i, cells[i].id) });
                    other->cells[j].centi -= k; cells[i].centi += k;
                    if (other->cells[j].centi <= 0) other->cells[j] = InvSlot{};
                    return;
                }
        }
        void swapCell(int i, int cell) override { swapped.push_back({ i, cell }); }
        void drop(int i, int64_t centi) override { dropped.push_back({ i, centi }); cells[i].centi -= centi; if (cells[i].centi <= 0) cells[i] = InvSlot{}; }
        bool infinite() const override { return catalog; }
    };
    InvSlot S(const char* id, int n) { InvSlot s; s.id = id; s.centi = (int64_t)n * 100; return s; }
    int N(const InvSlot& s, const char* id) { return !s.empty() && s.id == id ? (int)(s.centi / 100) : 0; }
}

bool SlotRules_SelfTest() {
    bool ok = true;
    auto check = [&](bool c, const char* what) { printf("bench: slot rules: %s: %s\n", what, c ? "yes" : "NO"); ok = ok && c; };
    TestGrid A(9), B(9), C(3);
    A.other = &B; B.other = &A;
    C.catalog = true; C.cells[0] = S("diamond", 1);
    InvSlot hand;
    SlotRules_SetHeld(&hand);
    SlotRules_SetGrids({ &A, &B });
    auto reset = [&]() { for (auto* g : { &A, &B }) for (auto& c : g->cells) c = InvSlot{}; hand = InvSlot{}; SlotRules_Cancel(); };
    // input helpers: a press and release on one slot, a drag across several
    auto click = [&](SlotGrid* g, int i, bool right = false, bool shift = false, bool dbl = false) {
        SlotInput p; p.shift = shift; (right ? p.rPress : p.lPress) = true; (right ? p.rDown : p.lDown) = true; p.lDouble = dbl && !right;
        SlotRules_Frame(p, g, i);
        SlotInput r; r.shift = shift; (right ? r.rRelease : r.lRelease) = true;
        SlotRules_Frame(r, g, i);
    };
    auto dragAcross = [&](SlotGrid* g, std::vector<int> slots, bool right = false, bool shift = false) {
        SlotInput p; p.shift = shift; (right ? p.rPress : p.lPress) = true; (right ? p.rDown : p.lDown) = true;
        SlotRules_Frame(p, g, slots[0]);
        for (size_t k = 1; k < slots.size(); ++k) { SlotInput m; m.shift = shift; (right ? m.rDown : m.lDown) = true; SlotRules_Frame(m, g, slots[k]); }
        SlotInput r; r.shift = shift; (right ? r.rRelease : r.lRelease) = true;
        SlotRules_Frame(r, g, slots.back());
    };
    auto key = [&](SlotGrid* g, int i, SlotInput in) { SlotRules_Frame(in, g, i); };

    // -- Minecraft
    reset(); A.cells[0] = S("stone", 10);
    click(&A, 0);                 check(N(hand, "stone") == 10 && A.cells[0].empty(), "LMB on a stack picks it all up");
    click(&A, 1);                 check(N(A.cells[1], "stone") == 10 && hand.empty(), "LMB holding a stack puts it all down");
    A.cells[2] = S("stone", 60); click(&A, 1); click(&A, 2);
    check(N(A.cells[2], "stone") == 64 && N(hand, "stone") == 6, "LMB on the same item merges; what doesn't fit stays held");
    A.cells[3] = S("oak_log", 3); click(&A, 3);
    check(N(A.cells[3], "stone") == 6 && N(hand, "oak_log") == 3, "LMB on another item swaps");
    reset(); A.cells[0] = S("stone", 9);
    click(&A, 0, true);           check(N(hand, "stone") == 5 && N(A.cells[0], "stone") == 4, "RMB on a stack picks up half, rounded up");
    click(&A, 1, true);           check(N(A.cells[1], "stone") == 1 && N(hand, "stone") == 4, "RMB holding a stack puts one down");
    A.cells[2] = S("oak_log", 2); click(&A, 2, true);
    check(N(hand, "oak_log") == 2 && N(A.cells[2], "stone") == 4, "RMB on another item swaps");
    reset(); A.cells[0] = S("stone", 10); B.cells[4] = S("stone", 60);
    click(&A, 0, false, true);    check(A.cells[0].empty() && N(B.cells[4], "stone") == 64 && N(B.cells[0], "stone") == 6, "Shift+LMB quick-moves to the other side (onto the same item first)");
    reset(); A.cells[0] = S("stone", 10); A.cells[3] = S("stone", 64); B.cells[1] = S("stone", 20);
    click(&A, 0); click(&A, 0, false, false, true);
    check(N(hand, "stone") == 64 && N(B.cells[1], "stone") == 0 && N(A.cells[3], "stone") == 30, "a double-click gathers that item onto the hand (partial stacks first), up to a stack");
    reset(); hand = S("stone", 10);
    dragAcross(&A, { 0, 1, 2 });  check(N(A.cells[0], "stone") == 3 && N(A.cells[1], "stone") == 3 && N(A.cells[2], "stone") == 3 && N(hand, "stone") == 1, "LMB drag splits the stack evenly (the rest stays held)");
    reset(); hand = S("stone", 10); A.cells[1] = S("oak_log", 1);
    dragAcross(&A, { 0, 1, 2 });  check(N(A.cells[0], "stone") == 5 && N(A.cells[1], "oak_log") == 1 && N(A.cells[2], "stone") == 5, "...skipping slots holding another item");
    reset(); hand = S("stone", 10);
    dragAcross(&A, { 0, 1, 2, 3 }, true); check(N(A.cells[0], "stone") == 1 && N(A.cells[3], "stone") == 1 && N(hand, "stone") == 6, "RMB drag puts one into each slot");
    reset(); A.cells[0] = S("stone", 5);
    { SlotInput k; k.hotbarKey = 3; key(&A, 0, k); } check(!A.swapped.empty() && A.swapped.back() == std::make_pair(0, Inventory::kGrid + 3), "4 over a slot swaps it with hotbar cell 4");
    { SlotInput k; k.offHandKey = true; key(&A, 0, k); } check(A.swapped.back() == std::make_pair(0, Inventory::kOffHand), "the off-hand key swaps it with the off hand");
    { SlotInput k; k.dropKey = true; key(&A, 0, k); } check(N(A.cells[0], "stone") == 4 && A.dropped.back().second == 100, "Q drops one");
    { SlotInput k; k.dropKey = true; k.ctrl = true; key(&A, 0, k); } check(A.cells[0].empty() && A.dropped.back().second == 400, "Ctrl+Q drops the stack");
    reset(); A.cells[0] = S("stone", 5);
    { SlotInput k; k.mPress = true; key(&A, 0, k); } check(hand.empty(), "middle-click in survival takes nothing");
    { SlotInput k; k.mPress = true; k.creative = true; key(&A, 0, k); } check(N(hand, "stone") == 64 && N(A.cells[0], "stone") == 5, "middle-click in creative: a full stack of it");
    // -- the Catalog
    reset(); SlotRules_SetGrids({ &A, &B, &C });
    click(&C, 0);                 check(N(hand, "diamond") == 64 && N(C.cells[0], "diamond") == 1, "Catalog: LMB takes a full stack for free");
    click(&C, 0);                 check(hand.empty(), "Catalog: a held stack put on it is deleted");
    click(&C, 0, true);           check(N(hand, "diamond") == 1, "Catalog: RMB takes one");
    SlotRules_SetGrids({ &A, &B });
    // -- Mouse Tweaks
    reset(); A.cells[0] = S("stone", 10); A.cells[1] = S("stone", 20); A.cells[2] = S("oak_log", 4); A.cells[3] = S("stone", 50);
    dragAcross(&A, { 0, 1, 2, 3 }); check(N(hand, "stone") == 64 && N(A.cells[2], "oak_log") == 4 && N(A.cells[3], "stone") == 16, "LMB drag after a pick-up gathers that item from the slots crossed (up to a stack)");
    reset(); A.cells[0] = S("stone", 10); A.cells[1] = S("oak_log", 3); A.cells[2] = S("sand", 2);
    dragAcross(&A, { 0, 1, 2 }, false, true); check(A.cells[0].empty() && A.cells[1].empty() && A.cells[2].empty() && N(B.cells[1], "oak_log") == 3, "Shift+LMB drag quick-moves every slot crossed");
    reset(); A.cells[0] = S("stone", 10);
    { SlotInput k; k.wheel = -1.0f; key(&A, 0, k); } check(N(A.cells[0], "stone") == 9 && N(B.cells[0], "stone") == 1, "the wheel down moves one to the other side");
    { SlotInput k; k.wheel = -1.0f; k.shift = true; key(&A, 0, k); } check(A.cells[0].empty() && N(B.cells[0], "stone") == 10, "...Shift+wheel the stack");
    A.cells[5] = S("stone", 1);
    { SlotInput k; k.wheel = 2.0f; key(&A, 5, k); } check(N(A.cells[5], "stone") == 3 && N(B.cells[0], "stone") == 8, "the wheel up pulls one in per notch");
    { SlotInput k; k.wheel = 1.0f; k.wheelInvert = true; key(&A, 5, k); } check(N(A.cells[5], "stone") == 2, "the direction can be swapped");
    SlotRules_SetHeld(nullptr);
    printf("bench: slot rules (#862): %s\n", ok ? "PASS" : "FAIL");
    return ok;
}
