// adventure/ruin-cache (#2072): sealed ruin caches opened through the real input path (E held, on foot), one per opener kind
// (cache_openers.txt): a small cache by power (its default), a vault by a crowbar (a tool, kept) and a third by a key (a
// part, spent). Never broken; E without what it needs opens nothing; a lock part-charged says how far; the hold takes
// kCacheHoldSeconds (half of it opens nothing); the rewards land in the inventory and the result shows; opened, it stays
// open (once per world) and E again gives nothing. The factory's charge is Progression's: here a stand-in reports it.
#include "scenario.h"
#include "buildtool.h"
#include "cacheopeners.h"
#include "factory.h"
#include "gameui.h"
#include "inventory.h"
#include "player.h"
#include "playerinput.h"
#include "plots.h"
#include "ruincache.h"
#include "world.h"
#include <raymath.h>
#include <fstream>
#include <sstream>

namespace {
uint32_t plotId = 0;
int small = -1, vault = -1, keyed = -1;
constexpr uint64_t kE = 1ull << (int)Act::Use;
// The stand-in for the factory's lock charge (Progression's side): the small cache's charge, and whether it was spent.
float charge = 0.0f;
int spends = 0;
constexpr int kKeyOpener = 7; // the scenario's own opener: a key (a part, spent), appended to the game's catalogue
constexpr const char* kKey = "iron_ingot";

Plot* ThePlot(ScenarioCtx& c) { return c.body->plots ? c.body->plots->byId(plotId) : nullptr; }
// The part at a cell (indices move when a part is placed again: opening re-places the cache).
int At(Plot& p, Int3 cell) { return p.partAt(cell); }
const Int3 kSmall{ 0, 0, 0 }, kVault{ 3, 0, 0 }, kKeyed{ 6, 0, 0 };

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
        // The game's openers, and a key of the scenario's own (no key item exists yet: an ingot stands in).
        std::ifstream in("resources/data/cache_openers.txt");
        std::stringstream text; text << in.rdbuf() << "\nopener " << kKeyOpener << " test_key part " << kKey << "\n";
        std::string err;
        if (!c.check("the openers load", CacheOpeners_Parse(text.str(), &err), err)) return true;
        const CacheOpener* crowbar = CacheOpener_Find("crowbar");
        if (!c.check("a crowbar opener", crowbar && crowbar->kind == OpenerKind::Tool)) return true;
        auto isSmall = [](const Plot& p, int part) {
            if (p.id != plotId || part < 0 || part >= (int)p.parts.size()) return false;
            const Int3 o = p.parts[part].origin;
            return o.x == kSmall.x && o.y == kSmall.y && o.z == kSmall.z;
        };
        CacheOpeners_SetPower({ [isSmall](const Planetoid&, const Plot& p, int part) { return isSmall(p, part) ? charge : 0.0f; },
                                [](Planetoid&, Plot&, int) { charge = 0.0f; ++spends; } });
        small = plot->place({ PartType::RuinCache, kSmall, 0, (uint8_t)CacheKind::Small }); // its default: power
        vault = plot->place({ PartType::RuinCache, kVault, 0, (uint8_t)((int)CacheKind::Vault | CacheOpenerBits(crowbar->index)) });
        keyed = plot->place({ PartType::RuinCache, kKeyed, 0, (uint8_t)((int)CacheKind::Small | CacheOpenerBits(kKeyOpener)) });
        if (!c.check("three placed", small >= 0 && vault >= 0 && keyed >= 0)) return true;
        const CacheOpener* def = CacheOpenerOf(plot->parts[At(*plot, kSmall)]);
        c.check("a small cache opens by power by default", def && def->kind == OpenerKind::Power);
        c.check("the vault keeps its crowbar opener", CacheOpenerOf(plot->parts[At(*plot, kVault)]) == crowbar);
        c.check("a cache can't be broken", PartUnbreakable(plot->parts[At(*plot, kSmall)]));
        const uint64_t idA = CacheId(*c.body, *plot, plot->parts[At(*plot, kSmall)]), idB = CacheId(*c.body, *plot, plot->parts[At(*plot, kVault)]);
        c.check("each cache has its own id", idA != idB);
        c.check("an id is the same asked twice", idA == CacheId(*c.body, *plot, plot->parts[At(*plot, kSmall)]));
        // The player 2.5 m out from between them, on foot.
        const Vector3 mid = CellCentre(*plot, { 2, 0, 0 });
        const Vector3 up = Vector3Normalize(mid), out = Vector3Normalize(Vector3CrossProduct(up, plot->east));
        player.flyMode = false;
        player.placeAt(c.world, *c.body, c.body->toWorld(Vector3Add(Vector3Scale(up, Vector3Length(mid) + 0.3f), Vector3Scale(out, 2.5f))));
        inventory.remove("crowbar", inventory.count("crowbar"));
        inventory.remove(kKey, inventory.count(kKey));
        gScenarioInput = true; // (the scripted E is used: scenarios are benches, which take no input otherwise)
        return true;
    });
    s.step("aimed at the small cache, on the ground", 10.0, [](ScenarioCtx& c) {
        Aim(c, kSmall);
        c.waitingFor("grounded, the cache aimed");
        if (c.t < 1.0 || !player.grounded || BuildTool_AimedType() != PartType::RuinCache) return false;
        c.check("E on it is the cache's, not the game window's", BuildTool_AimingUsable());
        Plot* p = ThePlot(c);
        const std::string need = p ? CacheMissing(*c.body, *p, At(*p, kSmall), inventory) : "";
        c.check("it says what it needs: power", need.find(" EU at ") != std::string::npos, need);
        return true;
    });
    s.step("E held with no charge in its lock: still sealed", 10.0, [](ScenarioCtx& c) {
        Aim(c, kSmall);
        if (!HoldE(c, kCacheHoldSeconds + 1.0, 0.5)) return false;
        c.check("sealed", !Opened(c, kSmall));
        c.check("the inventory never opened", !GameUI_InventoryOpen());
        return true;
    });
    s.step("its lock 40% charged: sealed, and the prompt says how far", 10.0, [](ScenarioCtx& c) {
        charge = 0.4f;
        Aim(c, kSmall);
        if (!HoldE(c, kCacheHoldSeconds + 0.5, 0.5)) return false;
        Plot* p = ThePlot(c);
        const std::string need = p ? CacheMissing(*c.body, *p, At(*p, kSmall), inventory) : "";
        c.check("the prompt reads \"Charging 40%\"", need.rfind("Charging 40%", 0) == 0, need);
        c.check("still sealed", !Opened(c, kSmall));
        return true;
    });
    s.step("charged: half the hold opens nothing", 10.0, [](ScenarioCtx& c) {
        charge = 1.0f;
        Aim(c, kSmall);
        if (!HoldE(c, kCacheHoldSeconds * 0.5, 0.5)) return false;
        c.check("still sealed after half a hold", !Opened(c, kSmall));
        c.eq("the charge not spent", spends, 0);
        return true;
    });
    s.step("the full hold opens it: the charge spent, the rewards in, the result shown", 12.0, [](ScenarioCtx& c) {
        static int stoneBefore = -1; // (at the step's first frame)
        if (stoneBefore < 0) stoneBefore = inventory.count("stone");
        Aim(c, kSmall);
        if (!HoldE(c, kCacheHoldSeconds + 0.4, 1.0)) return false;
        c.check("opened", Opened(c, kSmall));
        c.eq("the charge spent once", spends, 1);
        c.eq("the stub's 4 stone in the inventory", inventory.count("stone") - stoneBefore, 4);
        c.eq("one result line", (int)RuinCache_LastResult().size(), 1);
        c.metric("result lines", (double)RuinCache_LastResult().size(), "");
        return true;
    });
    s.step("E on it again: empty, nothing more", 10.0, [](ScenarioCtx& c) {
        static int stoneBefore = -1;
        if (stoneBefore < 0) stoneBefore = inventory.count("stone");
        Aim(c, kSmall);
        if (!HoldE(c, kCacheHoldSeconds + 0.5, 0.5)) return false;
        c.check("still opened", Opened(c, kSmall));
        c.eq("no more stone", inventory.count("stone") - stoneBefore, 0);
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
        return true;
    });
    s.step("the keyed cache: no key, still sealed", 10.0, [](ScenarioCtx& c) {
        Aim(c, kKeyed);
        if (!HoldE(c, kCacheHoldSeconds + 0.5, 0.5)) return false;
        Plot* p = ThePlot(c);
        const std::string need = p ? CacheMissing(*c.body, *p, At(*p, kKeyed), inventory) : "";
        c.check("it names the key, used up", need.find("(used up)") != std::string::npos, need);
        c.check("sealed", !Opened(c, kKeyed));
        return true;
    });
    s.step("with the key it opens, and the key is spent", 12.0, [](ScenarioCtx& c) {
        static bool added = false;
        if (!added) { inventory.add(kKey, 2); added = true; }
        Aim(c, kKeyed);
        if (!HoldE(c, kCacheHoldSeconds + 0.4, 1.0)) return false;
        c.check("opened", Opened(c, kKeyed));
        c.eq("one key spent", inventory.count(kKey), 1);
        gScenarioInput = false;
        Factory::InstallCachePower(); // the game's own hooks again (the factory's lock charge)
        CacheOpeners_Reset(); // the game's own catalogue again, without the test key
        return true;
    });
}
