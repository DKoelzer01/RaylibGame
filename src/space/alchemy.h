#pragma once

// Alchemy (#1254, magic step 7; the magic doc section 4), behind RAYLIBGAME_MAGIC. Decided (user, 2026-11-21): ONE generic
// essence fluid; its aspect is a data-driven tag on what holds it. A tank, a phial or one connected pipe run holds one aspect
// at a time; a run that would mix two moves nothing (as fluid kinds never mix). Spilled into the world it loses its aspect
// and fizzles away within seconds. Look and use per aspect come from data (aspects.txt colours, alchemy.txt).
//
// Where the tag lives:
// - machines (distillery outputs, an essence tank, a cauldron's inputs): a factory stack named "essence_<aspect>" whose
//   count is in measures. The stack's name is the tag; it is saved and sent to clients with the machine like any slot.
//   "essence_" items are fluid to the factory (Factory::AddFluidPrefix): no conveyor, chest or hand moves them.
// - phials: items "phial_<aspect>" (one measure each), made per aspect at load; "phial" is the empty one.
// - pipe runs: the existing Pipe and Outlet parts. A run's aspect is worked out each pass from what is on it (its tanks,
//   else its distilleries' first essence), so it has nothing to save. Two aspects in its tanks, or a water pump, and the run
//   is refused.
// - spills: an outlet on a run pours "essence" fluid (fluids.txt), with no aspect; the pass takes it back after a few seconds.

#include "factory.h"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

class Scene;
enum class PartType : uint8_t;

struct AlchemySettings {
    int tankMeasures = 256, pipeMeasures = 8, feedMeasures = 16;
    float spillM3 = 0.05f, fizzleSeconds = 3.0f;
    int distillTicks = 100, distillEu = 4, brewTicks = 300;
};
struct EssenceLook { std::string aspect, text; float glow = 0.5f; };
enum class PotionEffect : uint8_t { Heal, NightVision, WaterBreathing, ColdResist, MiningHaste, CleanseBurden };
struct PotionDef {
    std::string id, name, text;                          // id without the "potion_" of its item
    std::vector<std::pair<std::string, int>> brew;       // aspect, measures (one or two)
    int makes = 2;
    PotionEffect effect = PotionEffect::Heal;
    float amount = 0.0f, seconds = 0.0f;
    bool onceADay = false;
    unsigned colour = 0xFFFFFF;
};

// ---- Data (pure) ----
// alchemy.txt; aspect names are checked against aspects.txt (Tome_Load first). Replaces the tables only on success.
bool Alchemy_LoadText(const std::string& text, std::string& error);
bool Alchemy_Load(const char* path = "resources/data/alchemy.txt");
const AlchemySettings& Alchemy_Settings();
const std::vector<EssenceLook>& Alchemy_Looks();
const EssenceLook* Alchemy_Look(const std::string& aspect);
const std::vector<PotionDef>& Alchemy_Potions();
const PotionDef* Alchemy_PotionOfItem(const std::string& item);  // "potion_<id>", else nullptr
inline std::string Alchemy_EssenceItem(const std::string& aspect) { return "essence_" + aspect; }
inline std::string Alchemy_PhialItem(const std::string& aspect) { return "phial_" + aspect; }
inline std::string Alchemy_PotionItem(const PotionDef& p) { return "potion_" + p.id; }
std::string Alchemy_AspectOfItem(const std::string& item);       // essence_<a> or phial_<a> -> a; else ""
// The items made from the data, in items.txt's format: an essence and a full phial per aspect, a potion per brew.
std::string Alchemy_ItemsText();
// The machine recipes made from the data: the distillery (every "of item:" line of aspects.txt: one item into a measure of
// essence per point of each aspect, at most four) and the cauldron (its brews: the essences, water 1 and an empty phial
// per potion made).
void Alchemy_Recipes(std::vector<Factory::RecipeInfo>& out);

// ---- Pipe runs (pure) ----
enum class RunRole : uint8_t { Distillery, Tank, Cauldron, Outlet, WaterPump };
struct RunMember {
    RunRole role = RunRole::Tank;
    std::vector<std::pair<std::string, int>> essence;   // its essence stacks (distillery: outputs; tank, cauldron: inputs)
    std::vector<std::string> wants;                      // a cauldron: the essence items its brews take
};
struct RunMove { int from = -1, to = -1; std::string item; int count = 0; }; // indices into the members; to an outlet: a spill
struct RunPlan {
    std::string aspect;          // what the run carries this pass ("" idle)
    bool refused = false;
    std::string why;             // refused: why ("Fire and Life would mix")
    std::vector<RunMove> moves;
};
// One pass of a run: refused with a water pump on it (kinds never mix) or two aspects in its tanks; else its aspect is its
// tanks' (or its first distillery's first essence) and up to pipeMeasures move: distilleries to cauldrons that want it, then
// tanks with room, then outlets; tanks to cauldrons, then outlets. Other aspects wait where they are.
RunPlan Alchemy_PlanRun(const std::vector<RunMember>& members, const AlchemySettings& s);

// ---- The game ----
bool Alchemy_IsVessel(PartType t);                  // distillery, essence tank, cauldron
void Alchemy_Install();                             // with RAYLIBGAME_MAGIC, after Tome_Load and before the factory starts
void Alchemy_Update(Scene& world, double now);      // pipe runs and spills, twice a second (not on a joined client)
// RMB on a vessel with phials in hand: empty ones fill from it, full ones pour in. False: not phials (open its screen).
bool Alchemy_UsePhials(uint64_t body, uint32_t plot, int part, PartType type, const std::string& held);
// The run state shown on a vessel's look-at panel ("" none): "Run: Fire essence", "Run refused: Fire and Life would mix".
std::string Alchemy_PanelLine(uint64_t body, uint32_t plot, int part);
int Alchemy_EssenceKind();                          // the fluid kind of spills (-1 without one)
struct AlchemyStats { int runs = 0, refused = 0, moved = 0, spilled = 0, fizzled = 0; };
AlchemyStats Alchemy_Stats();                       // totals since the start (runs, refused: the last pass)
