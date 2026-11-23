#ifndef FRAMEPACER_H
#define FRAMEPACER_H

// (#1583) The frame cap, paced by the game after EndDrawing instead of inside raylib's EndDrawing: raylib's wait (Sleep for
// 95% of the time left, 1 ms timer resolution, then a 0.35 ms spin at 144 FPS) overshot by up to ~1 ms (capped frames ran
// long and uneven), and it sat inside the GPU timer's "swap+poll" segment, so profiles showed the cap's padding as GPU
// swap time. Here: a high-resolution waitable timer to 0.5 ms before the deadline, then a spin; main.cpp wraps it in its own
// profiler zone ("frame cap"). raylib's own cap is set to 0 (its frame time still counts this wait: it measures from the
// end of EndDrawing to the next BeginDrawing).

void FramePacer_SetTarget(int fps); // 0: uncapped (also sets raylib's SetTargetFPS(0))
int FramePacer_Target();
void FramePacer_Wait();             // once a frame, after EndDrawing: until this frame's deadline

#endif
