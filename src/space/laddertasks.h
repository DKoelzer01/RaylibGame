#pragma once

// Tier tasks that watch the factory (#1271, ladder L5): "do wheel_press N" counts the seconds a mechanical press works while a
// water wheel or windmill turns on the same plot (tier 2's working build: a wheel or windmill turning a press under load).
// Plots near the player, once a second of sim time, into the selected quest book.

class Scene;
void LadderTasks_Update(Scene& world);
// The rule itself (tests): a second counts when a press works and turns, and a wheel or windmill turns, on the plot.
bool LadderTasks_WheelPressSecond(bool pressWorking, float pressRpm, float sourceRpm);
