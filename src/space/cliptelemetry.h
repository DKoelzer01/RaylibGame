#ifndef CLIPTELEMETRY_H
#define CLIPTELEMETRY_H

#include <string>

class Scene;

// Clip telemetry (#1442): in benches, scenarios, recordings and profiled runs, every frame checks whether the camera eye
// (and a 0.15 m ring around it, the near plane's reach) or the player's capsule is inside solid ground or a built cube, on
// the CPU density field with the player's edits. A run then knows its view is not normal gameplay: #1428 chased "black
// tubes" that were a bench camera 26 m inside a hill. Changes print as
//   bench: camera inside solid (d=... m) at t=... s   /   bench: camera out of solid at t=... s (... s inside)
// (the player likewise); scenarios get the counters in results.json, and a recorded frame gets a "CLIPPED" tag.
// Off (no samples) outside those runs. RAYLIBGAME_CLIP_TELEMETRY=0 turns it off, =1 on in any run.
struct ClipStats {
    long cameraFrames = 0, playerFrames = 0;      // frames with the camera / the player inside solid
    double cameraSeconds = 0.0, playerSeconds = 0.0;
    double cameraFirst = -1.0, cameraLast = -1.0; // game seconds of the first and last such frame (-1: none)
    double playerFirst = -1.0, playerLast = -1.0;
    float cameraDepthMax = 0.0f, playerDepthMax = 0.0f; // the deepest seen (m under the open surface above it)
};

// Once a frame after the world is drawn (the camera that drew it). t: game seconds (for the messages).
void ClipTelemetry_Frame(Scene& world, double t, double dt);
bool ClipTelemetry_CameraInside();          // this frame
float ClipTelemetry_CameraDepth();          // m, this frame (0 when outside)
const ClipStats& ClipTelemetry_Stats();     // since the start (or the last reset)
void ClipTelemetry_Reset();                 // a scenario's own count starts here
void ClipTelemetry_DrawTag();               // "CLIPPED d=... m" in the frame's corner while the camera is inside (recordings)

#endif
