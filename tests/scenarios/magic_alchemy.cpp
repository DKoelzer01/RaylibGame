// magic/alchemy (#1254): with RAYLIBGAME_MAGIC on, in the sandbox: a distillery powered by a generator boils eggs into
// Life essence, a pipe run carries it to an essence tank, a second tank holding Fire on the same run refuses the run (aspects
// never mix), a second run feeds the tank's Life to a cauldron that brews the cleansing draught (Order poured in as from
// phials), the draught lifts 15 burden once a day, and essence spilled from an outlet fizzles away. The camera watches the
// row of machines (the factory at 20x).
#include "scenario.h"
#include "alchemy.h"
#include "buildtool.h"
#include "factory.h"
#include "fluid.h"
#include "inventory.h"
#include "magicsight.h"
#include "playersession.h"
#include "plots.h"
#include "potions.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>
#include <vector>

namespace {
Plot* plot = nullptr;
int distillery = -1, generator = -1, tank = -1, tank2 = -1, cauldron = -1, outlet = -1;
Vector3 ground{}, up{}, east{};
int fizzledBefore = 0;
int PartAt(Int3 c) { return plot ? plot->partAt(c) : -1; }
int Count(const std::vector<std::pair<std::string, int>>& slots, const std::string& item) { int n = 0; for (const auto& [id, k] : slots) if (id == item) n += k; return n; }
Factory::MachineKey Key(ScenarioCtx& c, int part) { return std::make_tuple(c.body->id, plot->id, part); }
void Watch(ScenarioCtx& c) { // 4 m in front of the row, 3 m up, looking at its middle
    const Vector3 mid = Vector3Add(ground, Vector3Scale(east, 2.0f));
    c.camera(Vector3Add(mid, Vector3Add(Vector3Scale(Vector3CrossProduct(east, up), 5.0f), Vector3Scale(up, 3.0f))), Vector3Add(mid, Vector3Scale(up, 0.5f)), up);
}
int TankLife(ScenarioCtx& c) { return Count(Factory::View(c.body->id, plot->id, tank).inSlots, "essence_life"); }
}

SCENARIO("magic/alchemy") {
    s.tags = { "magic", "alchemy", "machines", "factory", "fluid" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.env = { { "RAYLIBGAME_MAGIC", "1" } };
    s.timeout = 150.0;
    s.step("a distillery on a generator, a pipe to a tank, a pipe to a cauldron", 15.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        up = c.spawnDirLocal();
        ground = Vector3Scale(up, c.field->radius + c.field->elevation(up) + 0.05f);
        east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, up));
        terrainTool.creative = true;
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        bool placed = BuildTool_Place(*c.body, *plot, { PartType::Distillery, { 0, 0, 0 }, 0 }, false)
                   && BuildTool_Place(*c.body, *plot, { PartType::Cable, { 0, 0, 2 }, 0 }, false)
                   && BuildTool_Place(*c.body, *plot, { PartType::Generator, { 0, 0, 3 }, 0 }, false)
                   && BuildTool_Place(*c.body, *plot, { PartType::Pipe, { 2, 0, 0 }, 0 }, false)
                   && BuildTool_Place(*c.body, *plot, { PartType::Pipe, { 3, 0, 0 }, 0 }, false)
                   && BuildTool_Place(*c.body, *plot, { PartType::EssenceTank, { 4, 0, 0 }, 0 }, false)
                   && BuildTool_Place(*c.body, *plot, { PartType::Pipe, { 6, 0, 0 }, 0 }, false)
                   && BuildTool_Place(*c.body, *plot, { PartType::Cauldron, { 7, 0, 0 }, 0 }, false);
        if (!c.check("the alchemy row placed", placed)) return true;
        distillery = PartAt({ 0, 0, 0 }); generator = PartAt({ 0, 0, 3 }); tank = PartAt({ 4, 0, 0 }); cauldron = PartAt({ 7, 0, 0 });
        Factory::Sync(c.world);
        Factory::WaitTicks(3);
        Factory::SetSpeed(20);
        Factory::InsertFuel(c.body->id, plot->id, generator, "coal", 8);
        Factory::Insert(c.body->id, plot->id, distillery, "egg", 12); // egg: life 2 (aspects.txt)
        return true;
    });
    s.step("the distillery's Life essence runs down the pipe into the tank", 40.0, [](ScenarioCtx& c) {
        Watch(c);
        Factory::Sync(c.world);
        Alchemy_Update(c.world, c.T);
        const int life = TankLife(c);
        c.waitingFor(TextFormat("tank: %d Life", life));
        if (life < 8) return false;
        c.check("the run carries Life", Alchemy_PanelLine(c.body->id, plot->id, tank).find("Life") != std::string::npos, Alchemy_PanelLine(c.body->id, plot->id, tank));
        c.check("essence never leaves as an item", !Factory::Accepts(PartType::Chest, "essence_life"));
        c.metric("tank Life after the run filled it", life, "measures", false);
        // A second tank holding Fire on the same run: the run is refused.
        if (!c.check("a second tank on the first run", BuildTool_Place(*c.body, *plot, { PartType::EssenceTank, { 2, 0, -2 }, 0 }, false))) return true;
        tank2 = PartAt({ 2, 0, -2 });
        Factory::Sync(c.world);
        Factory::WaitTicks(2);
        Factory::FluidMove fm; fm.hasTo = true; fm.to = Key(c, tank2); fm.toWhich = 0; fm.item = "essence_fire"; fm.count = 5;
        c.check("fire poured into the second tank", Factory::MoveFluids({ fm })[0] == 5);
        return true;
    });
    s.step("two aspects on one run: refused", 10.0, [](ScenarioCtx& c) {
        Watch(c);
        Factory::Sync(c.world);
        Alchemy_Update(c.world, c.T);
        const std::string line = Alchemy_PanelLine(c.body->id, plot->id, tank);
        c.waitingFor(line);
        if (line.find("refused") == std::string::npos) return false;
        const int before = TankLife(c);
        c.check("Fire never reached the Life tank", Count(Factory::View(c.body->id, plot->id, tank).inSlots, "essence_fire") == 0);
        c.note("refused: " + line + TextFormat(" (tank Life %d)", before));
        // Empty the Fire again (as with phials): the run carries Life again.
        Factory::FluidMove fm; fm.hasFrom = true; fm.from = Key(c, tank2); fm.fromWhich = 0; fm.item = "essence_fire"; fm.count = 5;
        Factory::MoveFluids({ fm });
        // The cauldron: water, empty phials, fuel, and Order poured in (as from phials); Life comes from the tank's run.
        Factory::Insert(c.body->id, plot->id, cauldron, "water", 2);
        Factory::Insert(c.body->id, plot->id, cauldron, "phial", 4);
        Factory::InsertFuel(c.body->id, plot->id, cauldron, "coal", 4);
        Factory::FluidMove order; order.hasTo = true; order.to = Key(c, cauldron); order.toWhich = 0; order.item = "essence_order"; order.count = 4; order.cap = 16;
        c.check("Order poured into the cauldron", Factory::MoveFluids({ order })[0] == 4);
        return true;
    });
    s.step("the cauldron brews the cleansing draught", 40.0, [](ScenarioCtx& c) {
        Watch(c);
        Factory::Sync(c.world);
        Alchemy_Update(c.world, c.T);
        const MachineView v = Factory::View(c.body->id, plot->id, cauldron);
        const int made = Count(v.outSlots, "potion_cleansing");
        c.waitingFor(TextFormat("cauldron: %s -> %s", v.input.c_str(), v.output.c_str()));
        if (made < 1) return false;
        c.check("the Life tank's run is not refused any more", Alchemy_PanelLine(c.body->id, plot->id, tank).find("refused") == std::string::npos);
        // Drink it: 15 burden off, and only once a day.
        MagicPlayerState& m = LocalSession().magic;
        m.burden = 30.0f;
        Potion_Local() = PotionState{};
        inventory.add("potion_cleansing", 2);
        c.check("the draught is drunk", Potion_UseLocal("potion_cleansing"));
        c.check("15 burden lifted", fabsf(m.burden - 15.0f) < 0.01f, TextFormat("%.1f", m.burden));
        c.check("not twice in a day", !Potion_UseLocal("potion_cleansing") && fabsf(m.burden - 15.0f) < 0.01f);
        // An outlet on the cauldron's run: what the tank holds beyond the cauldron's feed spills and fizzles.
        if (!c.check("an outlet on the second run", BuildTool_Place(*c.body, *plot, { PartType::Outlet, { 6, 0, 2 }, 0 }, false)
                     && BuildTool_Place(*c.body, *plot, { PartType::Pipe, { 6, 0, 1 }, 0 }, false))) return true;
        outlet = PartAt({ 6, 0, 2 });
        fizzledBefore = Alchemy_Stats().fizzled;
        return true;
    });
    s.step("spilled essence fizzles away", 25.0, [](ScenarioCtx& c) {
        Watch(c);
        Factory::Sync(c.world);
        Alchemy_Update(c.world, c.T);
        const AlchemyStats st = Alchemy_Stats();
        c.waitingFor(TextFormat("spilled %d, fizzled %d", st.spilled, st.fizzled - fizzledBefore));
        static double fizzledAt = -1.0;
        if (st.spilled <= 0 || st.fizzled - fizzledBefore < 1) return false;
        if (fizzledAt < 0.0) fizzledAt = c.T;
        if (c.T < fizzledAt + 2.0) return false; // the takes reach the world over the next fluid steps
        const Vector3 at = plot->toBody({ 6.5f * kCellSize, 0.5f * kCellSize, 2.5f * kCellSize });
        const float left = Fluid::FluidAt(*c.body, at, 10.0f, Alchemy_EssenceKind());
        const float spilledM3 = st.spilled * Alchemy_Settings().spillM3;
        c.metric("essence left in the world after the fizzle", left, "m^3");
        c.check("the spill fizzled away (under a tenth left)", left < 0.1f * spilledM3, TextFormat("%.3f of %.2f m^3", left, spilledM3));
        Factory::SetSpeed(1);
        return true;
    });
}
