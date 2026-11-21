#include "laddertasks.h"
#include "factory.h"
#include "plots.h"
#include "scene.h"
#include "simclock.h"
#include "survival.h"
#include "worldtypes.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>

// Tier tasks that watch the factory (#1271 L5). See laddertasks.h.

bool LadderTasks_WheelPressSecond(bool pressWorking, float pressRpm, float sourceRpm) {
    return pressWorking && fabsf(pressRpm) > 0.0f && fabsf(sourceRpm) > 0.0f;
}

void LadderTasks_Update(Scene& world) {
    static double next = 0.0;
    const double now = SimTime();
    if (next - now > 2.0) next = now; // (the clock went back: a load)
    if (now < next) return;
    next = now + 1.0;
    // The seconds of factory time since the last look (a fast factory clock, tests at 60x, counts its own time).
    static long long lastTicks = -1;
    const long long ticks = Factory::Stats().ticks;
    const double seconds = lastTicks < 0 || ticks < lastTicks ? 0.0 : (double)(ticks - lastTicks) / kTicksPerSecond;
    lastTicks = ticks;
    if (seconds <= 0.0) return;
    for (const auto& obj : world.objects) {
        auto* body = dynamic_cast<Planetoid*>(obj.get());
        if (!body || !body->plots) continue;
        const Vector3 cam = body->toLocal(world.camera.position);
        for (const auto& pl : body->plots->plots) {
            if (Vector3Distance(pl->toBody({ 0, 0, 0 }), cam) > 1000.0f) continue;
            float source = 0.0f, pressRpm = 0.0f;
            bool working = false;
            for (int i = 0; i < (int)pl->parts.size(); ++i) {
                const PartType t = pl->parts[i].type;
                if (t == PartType::WaterWheel || t == PartType::Windmill) source = std::max(source, fabsf(Factory::RotationAt(body->id, pl->id, i).rpm));
                else if (t == PartType::Press && Factory::Working(body->id, pl->id, i)) { working = true; pressRpm = std::max(pressRpm, fabsf(Factory::RotationAt(body->id, pl->id, i).rpm)); }
            }
            if (LadderTasks_WheelPressSecond(working, pressRpm, source)) { QuestBook_Count("wheel_press", seconds); return; } // (one plot is enough)
        }
    }
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>
TEST(LadderTasks, AWheelPressSecondNeedsAWorkingTurningPressAndATurningSource) {
    EXPECT_TRUE(LadderTasks_WheelPressSecond(true, 16.0f, 16.0f));
    EXPECT_FALSE(LadderTasks_WheelPressSecond(false, 16.0f, 16.0f)) << "an idle press";
    EXPECT_FALSE(LadderTasks_WheelPressSecond(true, 0.0f, 16.0f)) << "a press that doesn't turn (another shaft, a crank let go)";
    EXPECT_FALSE(LadderTasks_WheelPressSecond(true, 16.0f, 0.0f)) << "no wheel or windmill turning (a crank or a steam engine)";
}
#endif
