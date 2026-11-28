// factory/cache-power (#2072): a sealed cache's power lock (cache_openers.txt) charged by the factory, in game. A generator
// burning coal on a cable that runs past a small cache (its default opener: 2000 EU at 8 EU a tick, 8 V or more): the
// lock charges (the prompt says how far), fills, and the cache then opens through the real input path (E held, on foot);
// the charge is spent with it. The game's own hooks (Factory::InstallCachePower), no stand-in.
#include "scenario.h"
#include "cacheopeners.h"
#include "factory.h"
#include "inventory.h"
#include "player.h"
#include "playerinput.h"
#include "plots.h"
#include "ruincache.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>

namespace {
uint32_t plotId = 0;
int gen = -1;
constexpr uint64_t kE = 1ull << (int)Act::Use;
const Int3 kGen{ 0, 0, 0 }, kCache{ 3, 0, 0 };

Plot* ThePlot(ScenarioCtx& c) { return c.body->plots ? c.body->plots->byId(plotId) : nullptr; }
Vector3 CellCentre(Plot& p, Int3 cell) { return p.toBody({ (cell.x + 0.5f) * kCellSize, (cell.y + 0.4f) * kCellSize, (cell.z + 0.5f) * kCellSize }); }
void Aim(ScenarioCtx& c) {
    Plot* p = ThePlot(c);
    if (!p) return;
    const Vector3 at = CellCentre(*p, kCache);
    const Vector3 feet = c.body->toLocalD(player.position).toFloat();
    const Vector3 eye = Vector3Add(feet, Vector3Scale(Vector3Normalize(feet), 1.6f));
    player.forward = c.body->dirToWorld(Vector3Normalize(Vector3Subtract(at, eye)));
}
int CachePart(ScenarioCtx& c) { Plot* p = ThePlot(c); return p ? p->partAt(kCache) : -1; }
float Charge(ScenarioCtx& c) { Plot* p = ThePlot(c); const int i = CachePart(c); return p && i >= 0 ? CacheOpener_Charge(*c.body, *p, i) : -1.0f; }
std::string Missing(ScenarioCtx& c) { Plot* p = ThePlot(c); const int i = CachePart(c); return p && i >= 0 ? CacheMissing(*c.body, *p, i, inventory) : "?"; }
bool Opened(ScenarioCtx& c) { Plot* p = ThePlot(c); const int i = CachePart(c); return p && i >= 0 && CacheOpened(p->parts[i]); }
} // namespace

SCENARIO("factory/cache-power") {
    s.tags = { "factory", "adventure", "power" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day", "--record-hud", "1" }; // the recording shows the cache prompt's charge
    s.timeout = 90.0;
    s.step("a generator, a cable and a sealed cache; the player on foot before it", 20.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        terrainTool.creative = true;
        const Vector3 dir = c.spawnDirLocal();
        const float ground = c.field->radius + c.field->elevation(dir);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        Plot* plot = c.body->plots->plotFor(Vector3Scale(dir, ground + 0.05f), 0.0f, true);
        plotId = plot->id;
        gen = plot->place({ PartType::Generator, kGen, 0 });
        for (int x = 0; x <= 4; ++x) plot->place({ PartType::Cable, { x, 0, -1 }, 0 });
        const int cache = plot->place({ PartType::RuinCache, kCache, 0, (uint8_t)CacheKind::Small }); // its default: power
        if (!c.check("placed", gen >= 0 && cache >= 0)) return true;
        const CacheOpener* o = CacheOpenerOf(plot->parts[cache]);
        c.check("a small cache opens by power", o && o->kind == OpenerKind::Power);
        const Vector3 mid = CellCentre(*plot, kCache);
        const Vector3 up = Vector3Normalize(mid), out = Vector3Normalize(Vector3CrossProduct(up, plot->east));
        player.flyMode = false;
        player.placeAt(c.world, *c.body, c.body->toWorld(Vector3Add(Vector3Scale(up, Vector3Length(mid) + 0.3f), Vector3Scale(out, 3.0f))));
        gScenarioInput = true;
        Factory::Sync(c.world);
        Factory::WaitTicks(3);
        c.check("sealed, no charge yet", Charge(c) == 0.0f && Missing(c).rfind("Sealed", 0) == 0, Missing(c).c_str());
        Factory::Insert(c.body->id, plotId, gen, "coal", 16);
        return true;
    });
    s.step("the generator charges the lock: the prompt says how far", 15.0, [](ScenarioCtx& c) {
        Aim(c);
        c.waitingFor(TextFormat("charge %.3f", Charge(c)));
        if (Charge(c) <= 0.0f) return false;
        c.check("charging, not open yet", Charge(c) < 1.0f && Missing(c).rfind("Charging", 0) == 0, Missing(c).c_str());
        return c.t >= 1.5; // (a moment on the prompt for the recording)
    });
    s.step("full: nothing missing", 30.0, [](ScenarioCtx& c) {
        Aim(c);
        Factory::SetSpeed(20);
        c.waitingFor(TextFormat("charge %.3f", Charge(c)));
        if (Charge(c) < 1.0f) return false;
        Factory::SetSpeed(1);
        c.check("nothing missing", Missing(c).empty(), Missing(c).c_str());
        return true;
    });
    s.step("held E opens it, and the charge is spent", 12.0, [](ScenarioCtx& c) {
        Aim(c);
        if (c.t < kCacheHoldSeconds + 0.4) { gFakeHeld |= kE; return false; }
        gFakeHeld &= ~kE;
        if (c.t < kCacheHoldSeconds + 1.4) return false;
        c.check("opened", Opened(c));
        Plot* p = ThePlot(c);
        c.check("the charge is spent", p && Factory::CacheCharged(c.body->id, p->id, CachePart(c)) == 0.0f);
        gScenarioInput = false;
        return true;
    });
}
