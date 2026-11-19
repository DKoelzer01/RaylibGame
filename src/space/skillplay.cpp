#include "skillplay.h"
#include "skills.h"
#include "playersession.h"
#include "terraintool.h"
#include "tools.h"
#include "gameui.h"
#include "simclock.h"
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
    if (ToolOf(toolId).type != ToolType::Pickaxe) return strokes;
    wearCarry += strokes * Skills_Stat(LocalSession().skills, "pick_wear"); // Light Hands: fewer strokes wear the pick
    const int n = (int)wearCarry;
    wearCarry -= n;
    return n;
}

int SkillPlay_LevelUps() { return levelUps; }
