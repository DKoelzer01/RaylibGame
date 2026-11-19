#include "skillplay.h"
#include "skills.h"
#include "playersession.h"
#include "terraintool.h"
#include "tools.h"
#include "gameui.h"
#include "simclock.h"
#include <cmath>
#include <cstring>
#include <vector>

namespace {
    // A stroke's cell for the anti-farming share: the body-local point on a 0.5 m grid.
    uint64_t CellKey(Vector3 local) {
        const auto q = [](float v) { return (uint64_t)(uint32_t)(int32_t)std::floor(v * 2.0f) & 0x1FFFFFull; };
        return q(local.x) | (q(local.y) << 21) | (q(local.z) << 42);
    }
    bool IsStroke(const char* e) { return !strcmp(e, "rock_stroke") || !strcmp(e, "ore_stroke") || !strcmp(e, "soil_stroke"); }

    void OnEvent(const char* event, double n, Vector3 local) {
        SkillState& s = LocalSession().skills;
        if (IsStroke(event)) n *= Skills_CellShare(s, CellKey(local), SimTime()); // the same cell again within a minute pays little
        std::vector<std::string> up;
        if (Skills_Event(s, event, n, worldXpRate, &up) <= 0) return;
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
