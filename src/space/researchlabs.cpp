#include "researchlabs.h"
#include "factory.h"
#include "gameui.h"
#include "net.h"
#include "plots.h"
#include "research.h"
#include "researchtab.h"
#include "researchgate.h"
#include "world.h"

#include <raylib.h>
#include <cstdlib>
#include <string>

namespace {
double nextAt = 0.0;
int unitsCounted = 0;
int LabTier(PartType t) { return t == PartType::SteamLab ? 3 : t == PartType::LvLab ? 4 : t == PartType::MvLab ? 5 : t == PartType::HvLab ? 7 : 0; }
}

int ResearchLabs_Update(Scene& world, double now, bool force) {
    if (gNetRole == NetRole::Client) return 0;            // (the server's: synced in a later step)
    if (nextAt - now > 1.0) nextAt = now; // (the sim clock went back: a new world or a load)
    if (!force && now < nextAt) return 0;
    nextAt = now + 0.5;
    ResearchState* st = nullptr;                           // the team's (made on first need: a team of one in single player)
    int added = 0;
    for (auto& o : world.objects) {
        Planetoid* body = dynamic_cast<Planetoid*>(o.get());
        if (!body || !body->plots) continue;
        for (auto& plot : body->plots->plots) {
            for (size_t i = 0; i < plot->parts.size(); ++i) {
                const int maxTier = LabTier(plot->parts[i].type);
                if (!maxTier) continue;
                if (!st) st = &Research_Team(ResearchTab_Team());
                // Steer it: the recipe of the first queued node it can run, else an output nothing makes (it idles).
                const ResearchNode* head = Research_HeadFor(*st, maxTier);
                Factory::SetVariant(body->id, plot->id, (int)i, head ? "research_t" + std::to_string(head->tier) : std::string("research_none"));
                // Collect: its units, onto the queued node of their tier.
                for (const auto& [item, count] : Factory::DrainOutputs(body->id, plot->id, (int)i, "research_t")) {
                    const int tier = atoi(item.c_str() + 10);
                    added += count;
                    for (const std::string& name : Research_AddUnits(*st, tier, count)) {
                        GameUI_Toast(TextFormat("Research done: %s", name.c_str()));
                        // (#1271 L6b-2) tier 3's working build: a steam lab, on a plot whose steam engine turns, completed a node
                        if (plot->parts[i].type == PartType::SteamLab)
                            for (int e = 0; e < (int)plot->parts.size(); ++e)
                                if (plot->parts[e].type == PartType::SteamEngine && Factory::RotationAt(body->id, plot->id, e).rpm != 0.0f) { ResearchGate_MarkBuild("lab_node"); break; }
                    }
                }
            }
        }
    }
    unitsCounted += added;
    return added;
}

int ResearchLabs_UnitsCounted() { return unitsCounted; }
void ResearchLabs_Reset() { nextAt = 0.0; unitsCounted = 0; }
