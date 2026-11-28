#ifndef WATERCAVES_H
#define WATERCAVES_H

// Water and caves (#1242, audio plan #1003 step 6), positional sounds in the world:
//  - rivers: the nearest river near the listener sounds along its course, as loud as its flow (its width, depth and fall,
//    times the share of its natural flow still in it: a diverted river goes quiet, #656);
//  - simulated water (dug channels, pours, released reaches): its flow (Fluid::FlowAt) on an 8 m grid around the listener;
//  - waterfalls: their own emitters at the pool (waterfalls.cpp), by the fall's power;
//  - caves: drips from the ceiling underground, wind at the cave mouths near you, a deep rumble under volcanoes;
//  (the sea's surf is the shore bed, soundscape.cpp; reverb by enclosure is Steam Audio's, #1294.)
// Audio v2 only.

class Scene;

// Pure helpers (tested).
// 0..1: how loud a river reach is, by its half width and depth (m), its slope (fall per metre) and the share of its
// natural flow left (1 as generated, 0 dammed dry).
float WaterCaves_RiverLoudness(float halfWidth, float depth, float slope, float flowShare);
// 0..1: a simulated water flow's loudness from its discharge (m^3/s, as Fluid::FlowAt reads it).
float WaterCaves_FlowLoudness(float discharge);
// Seconds to the next cave drip at this depth below the surface (m; deeper caves drip more, up to a point).
float WaterCaves_DripInterval(float depth, float random01);

void WaterCaves_Update(Scene& world, float dt);   // once a frame (v2)
void WaterCaves_Reset();                          // a new or unloaded world: drop the emitters

#endif
