#include "skillplay.h"
#include "skills.h"
#include "playersession.h"
#include "terraintool.h"
#include "tools.h"
#include "gameui.h"
#include "simclock.h"
#include "net.h"
#include "world.h"
#include "player.h"
#include <algorithm>
#include <raymath.h>
#include <cmath>
#include <cstring>
#include <vector>

namespace {
    // A stroke's or fill's cell for the anti-farming share: the body-local point on a 1 m grid (a fill sphere's size).
    uint64_t CellKey(Vector3 local) {
        const auto q = [](float v) { return (uint64_t)(uint32_t)(int32_t)std::floor(v) & 0x1FFFFFull; };
        return q(local.x) | (q(local.y) << 21) | (q(local.z) << 42);
    }
    bool IsStroke(const char* e) { return !strcmp(e, "rock_stroke") || !strcmp(e, "ore_stroke") || !strcmp(e, "soil_stroke"); }
    int levelUps = 0;

    void OnEvent(const char* event, double n, Vector3 local) {
        if (gNetRole == NetRole::Client) return; // (#920) joined: the server counts our XP (its SkillState replaces ours)
        SkillState& s = LocalSession().skills;
        if (!strcmp(event, "fill")) { // (no XP: it marks the ground as ours: the cells round its centre, where a dig would hit it)
            for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx)
                Skills_NoteFill(s, CellKey(Vector3Add(local, Vector3{ (float)dx, (float)dy, (float)dz })), SimTime());
            return;
        }
        if (IsStroke(event)) n *= Skills_DigShare(s, CellKey(local), SimTime()); // ground we filled, dug out again, pays little
        std::vector<std::string> up;
        const int k = Skills_Event(s, event, n, worldXpRate, &up);
        if (k <= 0) return;
        levelUps += k;
        for (const std::string& id : up) {
            const SkillDef* d = Skills_Find(id);
            const int lv = Skills_Level(s, id);
            GameUI_Toast(TextFormat("%s %d: a point to spend (Skills)", d ? d->name.c_str() : id.c_str(), lv));
        }
    }
    float OnStat(const char* stat) { return (float)Skills_Stat(LocalSession().skills, stat); }
    double wearCarry = 0.0;
}

void SkillPlay_Init() {
    Skills_Load();
    gSkillEvent = OnEvent;
    gSkillStat = OnStat;
}

int SkillPlay_WearStrokes(const std::string& toolId, int strokes) {
    if (strokes <= 0) return 0;
    // Each tool's wear stat: Mining's Light Hands (picks), (#919) Digging's for shovels and hoes, Woodcutting's for axes.
    const ToolType t = ToolOf(toolId).type;
    const char* stat = t == ToolType::Pickaxe ? "pick_wear" : t == ToolType::Shovel || t == ToolType::Hoe ? "shovel_wear" : t == ToolType::Axe ? "axe_wear" : nullptr;
    if (!stat) return strokes;
    wearCarry += strokes * Skills_Stat(LocalSession().skills, stat); // fewer strokes wear the tool
    const int n = (int)wearCarry;
    wearCarry -= n;
    return n;
}

int SkillPlay_LevelUps() { return levelUps; }

// (#919) Athletics and Survival XP from the player's own step: distance walked on the ground (not flying, not in a
// vehicle: this runs on foot only), and a night lived through.
void SkillPlay_Step(Scene& world, float dt) {
    static double walked = 0.0, nightFor = 0.0;
    static bool night = false;
    Player& p = LocalSession().player;
    if (LocalSession().survival.dead) { nightFor = 0.0; night = false; return; } // (a death breaks the night)
    if (p.grounded && !p.flyMode) {
        const Vector3 flat = Vector3Subtract(p.velocity, Vector3Scale(p.up, Vector3DotProduct(p.velocity, p.up)));
        walked += std::min(Vector3Length(flat), 20.0f) * dt; // (capped: a fling across the ground is not a walk)
        while (walked >= 10.0) { walked -= 10.0; if (gSkillEvent) gSkillEvent("walk_10m", 1.0, Vector3{}); }
    }
    // The sun's height over the player's horizon: down below -0.05, up over +0.05 (the band between: dusk and dawn).
    const float sunUp = Vector3DotProduct(Vector3Normalize(p.up), Vector3Normalize(world.sunDirection));
    if (sunUp < -0.05f) { night = true; nightFor += dt; }
    else if (sunUp > 0.05f && night) {
        if (nightFor >= 60.0 && gSkillEvent) gSkillEvent("night", 1.0, Vector3{});
        night = false; nightFor = 0.0;
    }
}
