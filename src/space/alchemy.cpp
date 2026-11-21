#include "alchemy.h"

#include "plots.h"
#include "tome.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

namespace {
AlchemySettings settings;
std::vector<EssenceLook> looks;
std::vector<PotionDef> potions;

std::string Trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}
std::vector<std::string> Fields(const std::string& line) {
    std::vector<std::string> out;
    std::istringstream fs(line);
    for (std::string f; std::getline(fs, f, '|');) out.push_back(Trim(f));
    return out;
}
bool StartsWith(const std::string& s, const char* p) { return s.rfind(p, 0) == 0; }
}

bool Alchemy_LoadText(const std::string& text, std::string& error) {
    AlchemySettings outS;
    std::vector<EssenceLook> outL;
    std::vector<PotionDef> outP;
    std::istringstream in(text);
    int lineNo = 0;
    auto fail = [&](const std::string& why) { error = "alchemy.txt line " + std::to_string(lineNo) + ": " + why; return false; };
    auto knownAspect = [](const std::string& a) { return Aspects_Find(a) != nullptr; };
    for (std::string line; std::getline(in, line);) {
        ++lineNo;
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;
        if (StartsWith(line, "set ")) {
            std::istringstream ss(line.substr(4));
            std::string name; double v = 0.0;
            if (!(ss >> name >> v)) return fail("'set <name> <value>'");
            if (name == "tank_measures") outS.tankMeasures = (int)v;
            else if (name == "pipe_measures") outS.pipeMeasures = (int)v;
            else if (name == "feed_measures") outS.feedMeasures = (int)v;
            else if (name == "spill_m3") outS.spillM3 = (float)v;
            else if (name == "fizzle_seconds") outS.fizzleSeconds = (float)v;
            else if (name == "distill_ticks") outS.distillTicks = (int)v;
            else if (name == "distill_eu") outS.distillEu = (int)v;
            else if (name == "brew_ticks") outS.brewTicks = (int)v;
            else return fail("unknown setting '" + name + "'");
            if (v <= 0.0) return fail("'" + name + "' must be above 0");
            continue;
        }
        const auto f = Fields(line);
        if (StartsWith(f[0], "essence ")) {
            EssenceLook l;
            l.aspect = Trim(f[0].substr(8));
            if (!knownAspect(l.aspect)) return fail("unknown aspect '" + l.aspect + "'");
            for (const EssenceLook& o : outL) if (o.aspect == l.aspect) return fail("essence '" + l.aspect + "' twice");
            for (size_t k = 1; k < f.size(); ++k) {
                if (StartsWith(f[k], "glow ")) l.glow = std::clamp((float)atof(f[k].c_str() + 5), 0.0f, 1.0f);
                else if (StartsWith(f[k], "text ")) l.text = Trim(f[k].substr(5));
                else return fail("unknown field '" + f[k] + "'");
            }
            outL.push_back(l);
        } else if (StartsWith(f[0], "potion ")) {
            if (f.size() < 2) return fail("a potion needs a name");
            PotionDef p;
            p.id = Trim(f[0].substr(7));
            p.name = f[1];
            if (p.id.empty() || p.id.find_first_of(" \t") != std::string::npos) return fail("bad potion id '" + p.id + "'");
            for (const PotionDef& o : outP) if (o.id == p.id) return fail("potion '" + p.id + "' twice");
            bool hasEffect = false;
            for (size_t k = 2; k < f.size(); ++k) {
                const std::string& fk = f[k];
                if (StartsWith(fk, "brew ")) {
                    std::istringstream ps(fk.substr(5));
                    for (std::string part; std::getline(ps, part, ',');) {
                        std::istringstream ws(part);
                        std::string a; int n = 0;
                        if (!(ws >> a >> n) || n <= 0) return fail("'" + Trim(part) + "' is not '<aspect> <n>'");
                        if (!knownAspect(a)) return fail("unknown aspect '" + a + "'");
                        for (const auto& b : p.brew) if (b.first == a) return fail("aspect '" + a + "' twice in a brew");
                        p.brew.push_back({ a, n });
                    }
                } else if (StartsWith(fk, "makes ")) p.makes = atoi(fk.c_str() + 6);
                else if (StartsWith(fk, "effect ")) {
                    std::istringstream es(fk.substr(7));
                    std::string kind; float amount = 0.0f;
                    es >> kind; es >> amount;
                    if (kind == "heal") p.effect = PotionEffect::Heal;
                    else if (kind == "night_vision") p.effect = PotionEffect::NightVision;
                    else if (kind == "water_breathing") p.effect = PotionEffect::WaterBreathing;
                    else if (kind == "cold_resist") p.effect = PotionEffect::ColdResist;
                    else if (kind == "mining_haste") p.effect = PotionEffect::MiningHaste;
                    else if (kind == "cleanse_burden") p.effect = PotionEffect::CleanseBurden;
                    else return fail("unknown effect '" + kind + "'");
                    p.amount = amount;
                    hasEffect = true;
                } else if (StartsWith(fk, "seconds ")) p.seconds = (float)atof(fk.c_str() + 8);
                else if (StartsWith(fk, "colour ")) p.colour = (unsigned)strtoul(fk.c_str() + 7, nullptr, 16);
                else if (fk == "once_a_day") p.onceADay = true;
                else if (StartsWith(fk, "text ")) p.text = Trim(fk.substr(5));
                else return fail("unknown field '" + fk + "'");
            }
            if (p.brew.empty() || p.brew.size() > 2) return fail("potion '" + p.id + "' needs one or two essences");
            if (!hasEffect) return fail("potion '" + p.id + "' has no effect");
            if (p.makes < 1 || p.makes > 8) return fail("potion '" + p.id + "' makes 1-8");
            const bool timed = p.effect == PotionEffect::NightVision || p.effect == PotionEffect::WaterBreathing ||
                               p.effect == PotionEffect::ColdResist || p.effect == PotionEffect::MiningHaste;
            if (timed && p.seconds <= 0.0f) return fail("potion '" + p.id + "' lasts no time");
            if (!timed && p.amount <= 0.0f) return fail("potion '" + p.id + "' needs an amount");
            outP.push_back(p);
        } else return fail("not a set, essence or potion line");
    }
    if (outS.feedMeasures > outS.tankMeasures) { error = "alchemy.txt: feed_measures above tank_measures"; return false; }
    settings = outS;
    looks = std::move(outL);
    potions = std::move(outP);
    return true;
}

bool Alchemy_Load(const char* path) {
    std::ifstream f(path);
    if (!f) { printf("alchemy: cannot read %s\n", path); return false; }
    std::stringstream ss; ss << f.rdbuf();
    std::string error;
    if (!Alchemy_LoadText(ss.str(), error)) { printf("alchemy: %s\n", error.c_str()); return false; }
    return true;
}

const AlchemySettings& Alchemy_Settings() { return settings; }
const std::vector<EssenceLook>& Alchemy_Looks() { return looks; }
const EssenceLook* Alchemy_Look(const std::string& aspect) { for (const EssenceLook& l : looks) if (l.aspect == aspect) return &l; return nullptr; }
const std::vector<PotionDef>& Alchemy_Potions() { return potions; }
const PotionDef* Alchemy_PotionOfItem(const std::string& item) {
    if (!StartsWith(item, "potion_")) return nullptr;
    for (const PotionDef& p : potions) if (item.compare(7, std::string::npos, p.id) == 0) return &p;
    return nullptr;
}
std::string Alchemy_AspectOfItem(const std::string& item) {
    std::string a;
    if (StartsWith(item, "essence_")) a = item.substr(8);
    else if (StartsWith(item, "phial_")) a = item.substr(6);
    return Aspects_Find(a) ? a : std::string();
}

std::string Alchemy_ItemsText() {
    std::string out;
    char buf[512];
    for (const AspectDef& a : Aspects_All()) {
        const int r = (a.colour >> 16) & 255, g = (a.colour >> 8) & 255, b = a.colour & 255;
        const EssenceLook* l = Alchemy_Look(a.id);
        const std::string use = l && !l->text.empty() ? " " + l->text : std::string();
        snprintf(buf, sizeof buf, "%s | %s essence | 999 | %d,%d,%d | %s essence (%s), held in a tank, a phial or a pipe run. Spilled, it fizzles.%s\n",
                 Alchemy_EssenceItem(a.id).c_str(), a.name.c_str(), r, g, b, a.name.c_str(), a.latin.c_str(), use.c_str());
        out += buf;
        snprintf(buf, sizeof buf, "%s | Phial of %s essence | 64 | %d,%d,%d | One measure of %s essence. Right-click an essence tank or a cauldron with it in hand to pour it in.%s\n",
                 Alchemy_PhialItem(a.id).c_str(), a.name.c_str(), r, g, b, a.name.c_str(), use.c_str());
        out += buf;
    }
    for (const PotionDef& p : potions) {
        const int r = (p.colour >> 16) & 255, g = (p.colour >> 8) & 255, b = p.colour & 255;
        snprintf(buf, sizeof buf, "%s | %s | 16 | %d,%d,%d | %s Drink it: bind it to the hotbar and press its key.\n", Alchemy_PotionItem(p).c_str(), p.name.c_str(), r, g, b, p.text.c_str());
        out += buf;
    }
    return out;
}

void Alchemy_Recipes(std::vector<Factory::RecipeInfo>& out) {
    for (const auto& [type, aspects] : Aspects_Sources()) {
        if (!StartsWith(type, "item:") || aspects.empty()) continue;
        Factory::RecipeInfo r;
        r.machine = PartType::Distillery;
        r.in.push_back({ type.substr(5), 1 });
        // At most four outputs (a machine's output slots): the strongest aspects (ties: by name).
        std::vector<std::pair<std::string, int>> sorted(aspects.begin(), aspects.end());
        std::stable_sort(sorted.begin(), sorted.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
        for (size_t k = 0; k < sorted.size() && k < 4; ++k) if (Aspects_Find(sorted[k].first)) r.out.push_back({ Alchemy_EssenceItem(sorted[k].first), sorted[k].second });
        if (r.out.empty()) continue;
        r.ticks = settings.distillTicks;
        r.eu = settings.distillEu;
        out.push_back(r);
    }
    for (const PotionDef& p : potions) {
        Factory::RecipeInfo r;
        r.machine = PartType::Cauldron;
        for (const auto& [a, n] : p.brew) r.in.push_back({ Alchemy_EssenceItem(a), n });
        r.in.push_back({ "water", 1 });
        r.in.push_back({ "phial", p.makes });
        r.out.push_back({ Alchemy_PotionItem(p), p.makes });
        r.ticks = settings.brewTicks;
        r.eu = 0;
        out.push_back(r);
    }
}

RunPlan Alchemy_PlanRun(const std::vector<RunMember>& members, const AlchemySettings& s) {
    RunPlan plan;
    auto aspectOf = [](const std::string& item) { return item.rfind("essence_", 0) == 0 ? item.substr(8) : std::string(); };
    for (const RunMember& m : members)
        if (m.role == RunRole::WaterPump) { plan.refused = true; plan.why = "a water pump is on it: water and essence never mix"; return plan; }
    // The tanks' aspects: one at most.
    std::vector<std::string> tankAspects;
    for (const RunMember& m : members) {
        if (m.role != RunRole::Tank) continue;
        for (const auto& [item, n] : m.essence)
            if (n > 0 && std::find(tankAspects.begin(), tankAspects.end(), aspectOf(item)) == tankAspects.end()) tankAspects.push_back(aspectOf(item));
    }
    if (tankAspects.size() > 1) {
        plan.refused = true;
        plan.why = tankAspects[0] + " and " + tankAspects[1] + " would mix";
        return plan;
    }
    if (!tankAspects.empty()) plan.aspect = tankAspects[0];
    else
        for (const RunMember& m : members)
            if (m.role == RunRole::Distillery && plan.aspect.empty())
                for (const auto& [item, n] : m.essence) if (n > 0) { plan.aspect = aspectOf(item); break; }
    if (plan.aspect.empty()) return plan; // nothing to carry
    const std::string item = "essence_" + plan.aspect;
    auto held = [&](const RunMember& m) { int c = 0; for (const auto& [it, n] : m.essence) if (it == item) c += n; return c; };
    std::vector<int> have(members.size()), room(members.size(), 0);
    for (size_t i = 0; i < members.size(); ++i) {
        have[i] = held(members[i]);
        const RunMember& m = members[i];
        if (m.role == RunRole::Tank) room[i] = std::max(0, s.tankMeasures - have[i]);
        else if (m.role == RunRole::Cauldron && std::find(m.wants.begin(), m.wants.end(), item) != m.wants.end()) room[i] = std::max(0, s.feedMeasures - have[i]);
        else if (m.role == RunRole::Outlet) room[i] = 1 << 30;
    }
    int budget = s.pipeMeasures;
    auto send = [&](RunRole fromRole, RunRole toRole) {
        for (size_t f = 0; f < members.size() && budget > 0; ++f) {
            if (members[f].role != fromRole) continue;
            for (size_t t = 0; t < members.size() && budget > 0 && have[f] > 0; ++t) {
                if (t == f || members[t].role != toRole || room[t] <= 0) continue;
                const int n = std::min({ have[f], room[t], budget });
                plan.moves.push_back({ (int)f, (int)t, item, n });
                have[f] -= n; room[t] -= n; budget -= n;
            }
        }
    };
    send(RunRole::Distillery, RunRole::Cauldron);
    send(RunRole::Distillery, RunRole::Tank);
    send(RunRole::Distillery, RunRole::Outlet);
    send(RunRole::Tank, RunRole::Cauldron);
    send(RunRole::Tank, RunRole::Outlet);
    return plan;
}

// ---------------------------------------------------------------------------------------------
// The game
// ---------------------------------------------------------------------------------------------
#include "fluid.h"
#include "gameui.h"
#include "inventory.h"
#include "magic.h"
#include "net.h"
#include "particles.h"
#include "plots.h"
#include "world.h"
#include "worldtypes.h"

#include <deque>
#include <unordered_map>
#include <unordered_set>

namespace {
double nextAt = 0.0;
AlchemyStats stats;
int essenceKind = -1;
constexpr float kFizzleReach = 8.0f; // m: a spill fizzles wherever it ran within this of its outlet
struct Spill { uint64_t body; Vector3 local; float m3; double at; int tries = 0; };
std::vector<Spill> spills;
std::map<Factory::MachineKey, std::string> panel; // the last pass's run state per vessel
inline int64_t CellKey(Int3 c) { return ((int64_t)(c.x + 100000) << 40) | ((int64_t)(c.y + 100000) << 20) | (int64_t)(c.z + 100000); }
// Where a vessel keeps its essence: a distillery in its outputs, a tank and a cauldron in their inputs.
int WhichOf(PartType t) { return t == PartType::Distillery ? 1 : 0; }
std::string AspectName(const std::string& aspect) { const AspectDef* a = Aspects_Find(aspect); return a ? a->name : aspect; }
}

bool Alchemy_IsVessel(PartType t) { return t == PartType::Distillery || t == PartType::EssenceTank || t == PartType::Cauldron; }
int Alchemy_EssenceKind() { return essenceKind; }
AlchemyStats Alchemy_Stats() { return stats; }

void Alchemy_Install() {
    static bool done = false;
    if (done) return;
    done = true;
    Factory::AddFluidPrefix("essence_");
    if (!Alchemy_Load()) return;
    itemDB.loadExtraText(Alchemy_ItemsText());
    gFactoryExtraRecipes = [](std::vector<Factory::RecipeInfo>& out) { Alchemy_Recipes(out); };
}

std::string Alchemy_PanelLine(uint64_t body, uint32_t plot, int part) {
    auto it = panel.find(std::make_tuple(body, plot, part));
    return it == panel.end() ? std::string() : it->second;
}

void Alchemy_Update(Scene& world, double now) {
    if (!Magic_Enabled() || gNetRole == NetRole::Client) return;
    if (nextAt - now > 1.0) { nextAt = now; spills.clear(); } // (the sim clock went back: a new world or a load)
    if (now < nextAt) return;
    nextAt = now + 0.5;
    if (essenceKind < 0) essenceKind = Fluid::KindByName("essence");
    const AlchemySettings& s = Alchemy_Settings();
    // Spills fizzle: their fluid is taken back out of the world (a few tries while it still runs off).
    for (size_t i = 0; i < spills.size();) {
        Spill& sp = spills[i];
        if (now - sp.at < s.fizzleSeconds) { ++i; continue; }
        Planetoid* b = Planetoid::find(sp.body);
        if (b && essenceKind >= 0) {
            Fluid::Take(*b, sp.local, std::max(sp.m3 * 1.5f, Fluid::FluidAt(*b, sp.local, kFizzleReach, essenceKind)), kFizzleReach, essenceKind); // (the spill runs out several metres)
            ParticleBurst pb; pb.kind = ParticleKind::Spark; pb.color = { 215, 195, 255, 255 }; pb.count = 10; pb.size = 0.06f; pb.life = 0.6f; pb.speed = 1.2f; pb.spread = 0.9f;
            Particles_Emit(*b, sp.local, pb);
        }
        if (b && essenceKind >= 0 && ++sp.tries < 12 && Fluid::FluidAt(*b, sp.local, kFizzleReach, essenceKind) > 0.005f) { sp.at = now - s.fizzleSeconds + 0.5; ++i; continue; }
        ++stats.fizzled;
        spills.erase(spills.begin() + i);
    }
    // The runs of every plot: pipe cells joined face to face, and the vessels, outlets and pumps they touch.
    struct Run { Planetoid* body; Plot* plot; std::vector<int> parts; std::vector<RunRole> roles; };
    std::vector<Run> runs;
    std::vector<std::pair<Int3, CellShape>> fp;
    for (auto& obj : world.objects) {
        Planetoid* body = dynamic_cast<Planetoid*>(obj.get());
        if (!body || !body->plots) continue;
        for (auto& pl : body->plots->plots) {
            std::unordered_map<int64_t, int> pipes, ends; // ends: vessels, outlets and pumps by cell
            bool anyVessel = false;
            for (int i = 0; i < (int)pl->parts.size(); ++i) {
                const PartType t = pl->parts[i].type;
                const bool vessel = Alchemy_IsVessel(t);
                anyVessel = anyVessel || vessel;
                if (t != PartType::Pipe && !vessel && t != PartType::Outlet && t != PartType::WaterPump && t != PartType::HandPump) continue;
                Plot::Footprint(pl->parts[i], fp);
                for (auto& c : fp) (t == PartType::Pipe ? pipes : ends)[CellKey(c.first)] = i;
            }
            if (!anyVessel || pipes.empty()) continue;
            std::unordered_set<int64_t> seen;
            const int d[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
            for (int i = 0; i < (int)pl->parts.size(); ++i) {
                if (pl->parts[i].type != PartType::Pipe) continue;
                Plot::Footprint(pl->parts[i], fp);
                if (fp.empty() || seen.count(CellKey(fp[0].first))) continue;
                Run run{ body, pl.get(), {}, {} };
                std::deque<Int3> q;
                for (auto& c : fp) if (seen.insert(CellKey(c.first)).second) q.push_back(c.first);
                std::set<int> members;
                while (!q.empty() && seen.size() < 20000) {
                    const Int3 c = q.front(); q.pop_front();
                    for (auto& o : d) {
                        const Int3 n = { c.x + o[0], c.y + o[1], c.z + o[2] };
                        const int64_t k = CellKey(n);
                        if (auto e = ends.find(k); e != ends.end()) { members.insert(e->second); continue; }
                        if (pipes.count(k) && seen.insert(k).second) q.push_back(n);
                    }
                }
                bool hasVessel = false;
                for (int m : members) {
                    const PartType t = pl->parts[m].type;
                    run.parts.push_back(m);
                    run.roles.push_back(t == PartType::Distillery ? RunRole::Distillery : t == PartType::EssenceTank ? RunRole::Tank : t == PartType::Cauldron ? RunRole::Cauldron
                                      : t == PartType::Outlet ? RunRole::Outlet : RunRole::WaterPump);
                    hasVessel = hasVessel || Alchemy_IsVessel(t);
                }
                if (hasVessel) runs.push_back(std::move(run));
            }
        }
    }
    // The vessels' essence, in one read.
    std::vector<Factory::MachineKey> keys;
    std::vector<std::vector<int>> keyOf(runs.size());
    for (size_t r = 0; r < runs.size(); ++r)
        for (size_t m = 0; m < runs[r].parts.size(); ++m) {
            const bool vessel = runs[r].roles[m] != RunRole::Outlet && runs[r].roles[m] != RunRole::WaterPump;
            keyOf[r].push_back(vessel ? (int)keys.size() : -1);
            if (vessel) keys.push_back(std::make_tuple(runs[r].body->id, runs[r].plot->id, runs[r].parts[m]));
        }
    const std::vector<Factory::FluidSlots> slots = Factory::ReadFluids(keys);
    // The cauldron's wants: the essences its brews take.
    std::vector<std::string> cauldronWants;
    for (const PotionDef& p : Alchemy_Potions()) for (const auto& [a, n] : p.brew) cauldronWants.push_back(Alchemy_EssenceItem(a));
    std::vector<Factory::FluidMove> moves;
    struct SpillAt { size_t move; Planetoid* body; Plot* plot; int part; };
    std::vector<SpillAt> spillMoves;
    panel.clear();
    stats.runs = (int)runs.size(); stats.refused = 0;
    for (size_t r = 0; r < runs.size(); ++r) {
        Run& run = runs[r];
        std::vector<RunMember> members(run.parts.size());
        for (size_t m = 0; m < run.parts.size(); ++m) {
            members[m].role = run.roles[m];
            const int k = keyOf[r][m];
            if (k >= 0 && slots[k].valid) members[m].essence = WhichOf(run.plot->parts[run.parts[m]].type) ? slots[k].out : slots[k].in;
            if (run.roles[m] == RunRole::Cauldron) members[m].wants = cauldronWants;
        }
        const RunPlan plan = Alchemy_PlanRun(members, s);
        if (plan.refused) ++stats.refused;
        const std::string line = plan.refused ? "Pipe run refused: " + plan.why : plan.aspect.empty() ? std::string("Pipe run: empty") : "Pipe run: " + AspectName(plan.aspect) + " essence";
        for (size_t m = 0; m < run.parts.size(); ++m) { // a vessel on several runs (a tank between two) shows a refused one first
            if (keyOf[r][m] < 0) continue;
            std::string& shown = panel[keys[keyOf[r][m]]];
            if (shown.find("refused") == std::string::npos) shown = line;
        }
        for (const RunMove& mv : plan.moves) {
            Factory::FluidMove fm;
            fm.item = mv.item; fm.count = mv.count;
            fm.hasFrom = true; fm.from = keys[keyOf[r][mv.from]]; fm.fromWhich = WhichOf(run.plot->parts[run.parts[mv.from]].type);
            if (members[mv.to].role == RunRole::Outlet) spillMoves.push_back({ moves.size(), run.body, run.plot, run.parts[mv.to] });
            else {
                fm.hasTo = true; fm.to = keys[keyOf[r][mv.to]]; fm.toWhich = 0;
                fm.cap = members[mv.to].role == RunRole::Tank ? s.tankMeasures : s.feedMeasures;
            }
            moves.push_back(fm);
        }
    }
    if (moves.empty()) return;
    const std::vector<int> done = Factory::MoveFluids(moves);
    for (int n : done) stats.moved += n;
    for (const SpillAt& sp : spillMoves) {
        const int n = done[sp.move];
        if (n <= 0 || essenceKind < 0) continue;
        stats.spilled += n;
        Plot::Footprint(sp.plot->parts[sp.part], fp);
        const Int3 c = fp.empty() ? sp.plot->parts[sp.part].origin : fp[0].first;
        const Vector3 at = sp.plot->toBody({ (c.x + 0.5f) * kCellSize, (c.y + 0.5f) * kCellSize, (c.z + 0.5f) * kCellSize });
        Fluid::Pour(*sp.body, at, n * s.spillM3, essenceKind);
        spills.push_back({ sp.body->id, at, n * s.spillM3, now });
    }
}

bool Alchemy_UsePhials(uint64_t body, uint32_t plot, int part, PartType type, const std::string& held) {
    if (!Alchemy_IsVessel(type)) return false;
    const std::string heldAspect = Alchemy_AspectOfItem(held);
    const bool empty = held == "phial", full = !heldAspect.empty() && held.rfind("phial_", 0) == 0;
    if (!empty && !full) return false;
    if (gNetRole == NetRole::Client) { GameUI_Toast("Phials work only on the host for now"); return true; }
    const int have = inventory.count(held);
    if (have <= 0) return true;
    const Factory::MachineKey key = std::make_tuple(body, plot, part);
    const std::vector<Factory::FluidSlots> slots = Factory::ReadFluids({ key });
    if (slots.empty() || !slots[0].valid) return true;
    const int which = WhichOf(type);
    const auto& stacks = which ? slots[0].out : slots[0].in;
    const AlchemySettings& s = Alchemy_Settings();
    if (empty) { // fill: from its biggest essence
        const std::pair<std::string, int>* best = nullptr;
        for (const auto& st : stacks) if (!best || st.second > best->second) best = &st;
        if (!best || best->second <= 0) { GameUI_Toast("Nothing to fill a phial with"); return true; }
        Factory::FluidMove fm; fm.hasFrom = true; fm.from = key; fm.fromWhich = which; fm.item = best->first; fm.count = std::min(have, best->second);
        const int n = Factory::MoveFluids({ fm })[0];
        if (n > 0 && inventory.remove("phial", n)) inventory.add(Alchemy_PhialItem(Alchemy_AspectOfItem(fm.item)), n);
        return true;
    }
    // Pour: a tank of one aspect, or a cauldron's brew; never into a distillery.
    if (type == PartType::Distillery) { GameUI_Toast("A distillery only gives essence"); return true; }
    const std::string item = Alchemy_EssenceItem(heldAspect);
    int cap = s.feedMeasures;
    if (type == PartType::EssenceTank) {
        for (const auto& st : stacks)
            if (st.second > 0 && st.first != item) { GameUI_Toast(TextFormat("This tank holds %s essence: aspects never mix", AspectName(Alchemy_AspectOfItem(st.first)).c_str())); return true; }
        cap = s.tankMeasures;
    } else {
        bool wanted = false;
        for (const PotionDef& p : Alchemy_Potions()) for (const auto& [a, n] : p.brew) wanted = wanted || a == heldAspect;
        if (!wanted) { GameUI_Toast(TextFormat("No brew uses %s essence", AspectName(heldAspect).c_str())); return true; }
    }
    Factory::FluidMove fm; fm.hasTo = true; fm.to = key; fm.toWhich = 0; fm.item = item; fm.count = have; fm.cap = cap;
    const int n = Factory::MoveFluids({ fm })[0];
    if (n <= 0) { GameUI_Toast("It's full"); return true; }
    if (inventory.remove(held, n)) inventory.add("phial", n);
    return true;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1254) ----
#include <gtest/gtest.h>

namespace {
void LoadShipped() {
    ASSERT_TRUE(Tome_Load());
    ASSERT_TRUE(Alchemy_Load()) << "resources/data/alchemy.txt (tests run from the repo root)";
}
RunMember Member(RunRole role, std::vector<std::pair<std::string, int>> essence = {}, std::vector<std::string> wants = {}) {
    RunMember m; m.role = role; m.essence = std::move(essence); m.wants = std::move(wants); return m;
}
}

TEST(Alchemy, ShippedDataLoadsAndEveryAspectHasALook) {
    LoadShipped();
    EXPECT_EQ(Alchemy_Looks().size(), Aspects_All().size()) << "a look line per aspect";
    EXPECT_EQ(Alchemy_Potions().size(), 6u);
    const PotionDef* c = Alchemy_PotionOfItem("potion_cleansing");
    ASSERT_NE(c, nullptr);
    EXPECT_EQ(c->effect, PotionEffect::CleanseBurden);
    EXPECT_FLOAT_EQ(c->amount, 15.0f);
    EXPECT_TRUE(c->onceADay);
    EXPECT_EQ(Alchemy_PotionOfItem("potion_nothing"), nullptr);
    EXPECT_EQ(Alchemy_AspectOfItem("phial_fire"), "fire");
    EXPECT_EQ(Alchemy_AspectOfItem("essence_life"), "life");
    EXPECT_EQ(Alchemy_AspectOfItem("phial_nothing"), "");
    EXPECT_EQ(Alchemy_AspectOfItem("phial"), "");
}

TEST(Alchemy, BadDataIsRefused) {
    LoadShipped();
    std::string error;
    EXPECT_FALSE(Alchemy_LoadText("essence nothing | glow 1\n", error)) << "unknown aspect";
    EXPECT_FALSE(Alchemy_LoadText("set nothing 3\n", error));
    EXPECT_FALSE(Alchemy_LoadText("set tank_measures 0\n", error));
    EXPECT_FALSE(Alchemy_LoadText("potion a | A | brew life 2 | effect fly\n", error)) << "unknown effect";
    EXPECT_FALSE(Alchemy_LoadText("potion a | A | brew life 2, fire 1, water 1 | effect heal 5\n", error)) << "three essences";
    EXPECT_FALSE(Alchemy_LoadText("potion a | A | brew life 2 | effect night_vision\n", error)) << "a timed effect needs seconds";
    EXPECT_FALSE(Alchemy_LoadText("potion a | A | brew life 2, life 1 | effect heal 5\n", error)) << "an aspect twice";
    EXPECT_TRUE(Alchemy_LoadText("potion a | A | brew life 2 | effect heal 5\n", error)) << error;
    EXPECT_EQ(Alchemy_Potions().size(), 1u);
    LoadShipped(); // back to the shipped tables for the other tests
}

TEST(Alchemy, RecipesComeFromTheAspectsAndTheBrews) {
    LoadShipped();
    std::vector<Factory::RecipeInfo> rs;
    Alchemy_Recipes(rs);
    int distill = 0, brew = 0;
    const Factory::RecipeInfo* oak = nullptr; const Factory::RecipeInfo* cleanse = nullptr;
    for (const auto& r : rs) {
        if (r.machine == PartType::Distillery) { ++distill; if (r.in[0].first == "oak_log") oak = &r; EXPECT_LE(r.out.size(), 4u); }
        if (r.machine == PartType::Cauldron) { ++brew; if (r.out[0].first == "potion_cleansing") cleanse = &r; }
    }
    int items = 0;
    for (const auto& [type, a] : Aspects_Sources()) items += type.rfind("item:", 0) == 0;
    EXPECT_EQ(distill, items) << "every item with aspects distils";
    EXPECT_EQ(brew, (int)Alchemy_Potions().size());
    ASSERT_NE(oak, nullptr);
    ASSERT_EQ(oak->out.size(), 1u);
    EXPECT_EQ(oak->out[0], (std::pair<std::string, int>{ "essence_plant", 2 }));
    ASSERT_NE(cleanse, nullptr);
    std::map<std::string, int> in(cleanse->in.begin(), cleanse->in.end());
    EXPECT_EQ(in["essence_life"], 2);
    EXPECT_EQ(in["essence_order"], 2);
    EXPECT_EQ(in["water"], 1);
    EXPECT_EQ(in["phial"], 1) << "an empty phial per potion";
    EXPECT_LE(cleanse->in.size(), 4u) << "fits a machine's four input slots";
}

TEST(Alchemy, ItemsLoadIntoTheItemTable) {
    LoadShipped();
    ItemDB db;
    ASSERT_TRUE(db.load("resources/data/items.txt", "resources/data/recipes.txt"));
    ASSERT_TRUE(db.loadExtra("resources/data/magic_items.txt"));
    ASSERT_TRUE(db.loadExtraText(Alchemy_ItemsText()));
    ASSERT_NE(db.find("phial"), nullptr);
    ASSERT_NE(db.find("phial_fire"), nullptr);
    ASSERT_NE(db.find("essence_order"), nullptr);
    ASSERT_NE(db.find("potion_healing"), nullptr);
    EXPECT_EQ(db.find("potion_cleansing")->name, "Cleansing draught");
    EXPECT_EQ(db.find("distillery")->kind, ItemKind::Block);
    EXPECT_EQ(db.find("cauldron")->kind, ItemKind::Block);
    EXPECT_EQ(db.find("essence_tank")->kind, ItemKind::Block);
    const Color fire = db.find("phial_fire")->color;
    EXPECT_EQ(fire.r, 0xe8); EXPECT_EQ(fire.g, 0x64); EXPECT_EQ(fire.b, 0x2c) << "the aspect's colour";
}

TEST(Alchemy, EssenceIsFluidToTheFactory) {
    Factory::AddFluidPrefix("essence_");
    EXPECT_TRUE(Factory::IsFluidItem("essence_fire"));
    EXPECT_FALSE(Factory::IsFluidItem("phial_fire"));
    EXPECT_FALSE(Factory::Accepts(PartType::Chest, "essence_fire")) << "no chest takes it from a conveyor or a hand";
    EXPECT_TRUE(Factory::Accepts(PartType::Chest, "stone"));
}

TEST(Alchemy, ARunCarriesOneAspect) {
    AlchemySettings s; s.pipeMeasures = 8; s.tankMeasures = 20; s.feedMeasures = 4;
    // A distillery, a tank, a cauldron that wants Life: Life goes to the cauldron first, then the tank.
    std::vector<RunMember> run = { Member(RunRole::Distillery, { { "essence_life", 10 }, { "essence_plant", 3 } }), Member(RunRole::Tank),
                                   Member(RunRole::Cauldron, {}, { "essence_life", "essence_order" }) };
    RunPlan p = Alchemy_PlanRun(run, s);
    EXPECT_FALSE(p.refused);
    EXPECT_EQ(p.aspect, "life");
    int toCauldron = 0, toTank = 0;
    for (const RunMove& m : p.moves) { EXPECT_EQ(m.item, "essence_life") << "plant waits"; (m.to == 2 ? toCauldron : toTank) += m.count; }
    EXPECT_EQ(toCauldron, 4) << "up to its feed";
    EXPECT_EQ(toTank, 4) << "the rest of the pass's budget";
    // The tank's aspect wins over the distillery's first essence: Plant waits in the distillery.
    run[1].essence = { { "essence_plant", 5 } };
    p = Alchemy_PlanRun(run, s);
    EXPECT_EQ(p.aspect, "plant");
    for (const RunMove& m : p.moves) EXPECT_EQ(m.item, "essence_plant");
    // Nothing on it: idle.
    p = Alchemy_PlanRun({ Member(RunRole::Tank), Member(RunRole::Outlet) }, s);
    EXPECT_TRUE(p.aspect.empty());
    EXPECT_TRUE(p.moves.empty());
}

TEST(Alchemy, MixingIsRefused) {
    AlchemySettings s;
    RunPlan p = Alchemy_PlanRun({ Member(RunRole::Tank, { { "essence_fire", 5 } }), Member(RunRole::Tank, { { "essence_life", 5 } }), Member(RunRole::Outlet) }, s);
    EXPECT_TRUE(p.refused) << "two aspects in its tanks";
    EXPECT_TRUE(p.moves.empty());
    EXPECT_NE(p.why.find("would mix"), std::string::npos) << p.why;
    p = Alchemy_PlanRun({ Member(RunRole::Distillery, { { "essence_fire", 5 } }), Member(RunRole::Tank), Member(RunRole::WaterPump) }, s);
    EXPECT_TRUE(p.refused) << "water and essence never mix";
    EXPECT_TRUE(p.moves.empty());
    // A full tank of the run's aspect takes no more; an outlet spills what is left.
    s.tankMeasures = 5; s.pipeMeasures = 8;
    p = Alchemy_PlanRun({ Member(RunRole::Distillery, { { "essence_fire", 6 } }), Member(RunRole::Tank, { { "essence_fire", 5 } }), Member(RunRole::Outlet) }, s);
    ASSERT_FALSE(p.refused);
    int fromDistillery = 0, fromTank = 0;
    for (const RunMove& m : p.moves) { EXPECT_EQ(m.to, 2) << "the full tank takes none"; (m.from == 0 ? fromDistillery : fromTank) += m.count; }
    EXPECT_EQ(fromDistillery, 6) << "the distillery's essence spills";
    EXPECT_EQ(fromTank, 2) << "the tank drains through the outlet with the rest of the pass's 8";
    // The budget caps a pass: a tank emptying through an outlet.
    s.pipeMeasures = 3;
    p = Alchemy_PlanRun({ Member(RunRole::Tank, { { "essence_fire", 50 } }), Member(RunRole::Outlet) }, s);
    ASSERT_EQ(p.moves.size(), 1u);
    EXPECT_EQ(p.moves[0].count, 3);
}
#endif
