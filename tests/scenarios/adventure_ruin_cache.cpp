// adventure/ruin-cache (#2072): a sealed ruin cache and a vault, opened through the real input path (E held, on foot).
// Never broken; E without what it needs opens nothing; a small cache spends a power cell, a vault keeps the crowbar; the
// hold takes kCacheHoldSeconds (half of it opens nothing); the rewards land in the inventory and the result shows; opened,
// it stays open (once per world) and E again gives nothing. The rewards are caches.txt's roll (#2072 slice 2).
#include "scenario.h"
#include "buildtool.h"
#include "gameui.h"
#include "inventory.h"
#include "player.h"
#include "playerinput.h"
#include "plots.h"
#include "ruincache.h"
#include "world.h"
#include <raymath.h>

namespace {
uint32_t plotId = 0;
int small = -1, vault = -1;
constexpr uint64_t kE = 1ull << (int)Act::Use;

Plot* ThePlot(ScenarioCtx& c) { return c.body->plots ? c.body->plots->byId(plotId) : nullptr; }
// The part at a cell (indices move when a part is placed again: opening re-places the cache).
int At(Plot& p, Int3 cell) { return p.partAt(cell); }
const Int3 kSmall{ 0, 0, 0 }, kVault{ 3, 0, 0 };

Vector3 CellCentre(Plot& p, Int3 cell) { return p.toBody({ (cell.x + 0.5f) * kCellSize, (cell.y + 0.4f) * kCellSize, (cell.z + 0.5f) * kCellSize }); }
// Each frame: the player's view on the cache at `cell` (the player stands; only the view turns).
void Aim(ScenarioCtx& c, Int3 cell) {
    Plot* p = ThePlot(c);
    if (!p) return;
    const Vector3 at = CellCentre(*p, cell);
    const Vector3 feet = c.body->toLocalD(player.position).toFloat();
    const Vector3 eye = Vector3Add(feet, Vector3Scale(Vector3Normalize(feet), 1.6f));
    player.forward = c.body->dirToWorld(Vector3Normalize(Vector3Subtract(at, eye)));
}
bool Opened(ScenarioCtx& c, Int3 cell) { Plot* p = ThePlot(c); const int i = p ? At(*p, cell) : -1; return i >= 0 && CacheOpened(p->parts[i]); }

// Holds E for `seconds` (sim time) from the step's start, then lets go; true once `after` seconds more have passed.
bool HoldE(ScenarioCtx& c, double seconds, double after) {
    if (c.t < seconds) { gFakeHeld |= kE; return false; }
    gFakeHeld &= ~kE;
    return c.t >= seconds + after;
}
} // namespace

SCENARIO("adventure/ruin-cache") {
    s.tags = { "adventure", "ui", "world" }; // (world: the runner's part 2)
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;
    s.step("a cache and a vault on the ground, the player on foot before them", 20.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const Vector3 dir = c.spawnDirLocal();
        const float ground = c.field->radius + c.field->elevation(dir);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        Plot* plot = c.body->plots->plotFor(Vector3Scale(dir, ground + 0.05f), 0.0f, true);
        plotId = plot->id;
        small = plot->place({ PartType::RuinCache, kSmall, 0, (uint8_t)CacheKind::Small });
        vault = plot->place({ PartType::RuinCache, kVault, 0, (uint8_t)CacheKind::Vault });
        if (!c.check("both placed", small >= 0 && vault >= 0)) return true;
        c.check("a cache can't be broken", PartUnbreakable(plot->parts[At(*plot, kSmall)]));
        const uint64_t idA = CacheId(*c.body, *plot, plot->parts[At(*plot, kSmall)]), idB = CacheId(*c.body, *plot, plot->parts[At(*plot, kVault)]);
        c.check("each cache has its own id", idA != idB);
        c.check("an id is the same asked twice", idA == CacheId(*c.body, *plot, plot->parts[At(*plot, kSmall)]));
        // The player 2.5 m out from between them, on foot.
        const Vector3 mid = CellCentre(*plot, { 2, 0, 0 });
        const Vector3 up = Vector3Normalize(mid), out = Vector3Normalize(Vector3CrossProduct(up, plot->east));
        player.flyMode = false;
        player.placeAt(c.world, *c.body, c.body->toWorld(Vector3Add(Vector3Scale(up, Vector3Length(mid) + 0.3f), Vector3Scale(out, 2.5f))));
        inventory.remove("power_cell", inventory.count("power_cell"));
        inventory.remove("crowbar", inventory.count("crowbar"));
        gScenarioInput = true; // (the scripted E is used: scenarios are benches, which take no input otherwise)
        return true;
    });
    s.step("aimed at the small cache, on the ground", 10.0, [](ScenarioCtx& c) {
        Aim(c, kSmall);
        c.waitingFor("grounded, the cache aimed");
        if (c.t < 1.0 || !player.grounded || BuildTool_AimedType() != PartType::RuinCache) return false;
        c.check("E on it is the cache's, not the game window's", BuildTool_AimingUsable());
        c.check("it says what it needs", !CacheMissing(inventory, CacheKind::Small).empty(), CacheMissing(inventory, CacheKind::Small).c_str());
        return true;
    });
    s.step("E held without a power cell: still sealed", 10.0, [](ScenarioCtx& c) {
        Aim(c, kSmall);
        if (!HoldE(c, kCacheHoldSeconds + 1.0, 0.5)) return false;
        c.check("sealed", !Opened(c, kSmall));
        c.check("the inventory never opened", !GameUI_InventoryOpen());
        return true;
    });
    s.step("with a power cell: half the hold opens nothing", 10.0, [](ScenarioCtx& c) {
        if (inventory.count("power_cell") == 0) inventory.add("power_cell", 1);
        Aim(c, kSmall);
        if (!HoldE(c, kCacheHoldSeconds * 0.5, 0.5)) return false;
        c.check("still sealed after half a hold", !Opened(c, kSmall));
        c.eq("the power cell not spent", inventory.count("power_cell"), 1);
        return true;
    });
    s.step("the full hold opens it: the cell spent, the rewards in, the result shown", 12.0, [](ScenarioCtx& c) {
        Aim(c, kSmall);
        if (!HoldE(c, kCacheHoldSeconds + 0.4, 1.0)) return false;
        c.check("opened", Opened(c, kSmall));
        c.eq("the power cell spent", inventory.count("power_cell"), 0);
        // (#2072 slice 2) What it holds is caches.txt's roll (unit-tested: Caches.*); a small cache holds 2 rewards.
        c.eq("a small cache's 2 result lines", (int)RuinCache_LastResult().size(), 2);
        for (const std::string& l : RuinCache_LastResult()) c.note("reward: " + l);
        c.metric("result lines", (double)RuinCache_LastResult().size(), "");
        return true;
    });
    s.step("E on it again: empty, nothing more", 10.0, [](ScenarioCtx& c) {
        static double shownBefore = -1.0; // the result screen's lines at the step's first frame
        static std::vector<std::string> before;
        if (shownBefore < 0.0) { before = RuinCache_LastResult(); shownBefore = 1.0; }
        Aim(c, kSmall);
        if (!HoldE(c, kCacheHoldSeconds + 0.5, 0.5)) return false;
        c.check("still opened", Opened(c, kSmall));
        c.check("no new rewards", RuinCache_LastResult() == before);
        return true;
    });
    s.step("the vault: no crowbar, still sealed", 10.0, [](ScenarioCtx& c) {
        Aim(c, kVault);
        if (!HoldE(c, kCacheHoldSeconds + 0.5, 0.5)) return false;
        c.check("aimed at the vault", BuildTool_AimedType() == PartType::RuinCache);
        c.check("sealed", !Opened(c, kVault));
        return true;
    });
    s.step("with a crowbar it opens, and the crowbar is kept", 12.0, [](ScenarioCtx& c) {
        static bool added = false; // (once, at the step's first frame: the crowbar is kept, so its count proves it)
        if (!added) { inventory.add("crowbar", 1); added = true; }
        Aim(c, kVault);
        if (!HoldE(c, kCacheHoldSeconds + 0.4, 1.5)) return false;
        c.check("opened", Opened(c, kVault));
        c.eq("the crowbar kept", inventory.count("crowbar"), 1);
        c.keepCapture("the vault's result on screen");
        gScenarioInput = false;
        return true;
    });
}
