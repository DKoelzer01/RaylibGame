#ifndef RECORDER_H
#define RECORDER_H

// In-game frame recording for visual debugging (#375). Agents read stills, not video, so a recording is PNG frames, a
// contact sheet (one PNG of up to 25 evenly spaced frames, each labelled with its frame number and time) and a GIF, all
// written by the game itself (no ffmpeg, no Python).
//
// Command line (works with every bench / test / shot mode):
//   --record <name>      frames to <root>/<YYYY-MM-DD>/<name>/frame_0000.png ... (an absolute path: there); implies --fixed-dt 1/60
//                        root = RAYLIBGAME_CAPTURE_ROOT, else //DEWEYSERVER/files/claude-media/raylib-captures when
//                        reachable, else %USERPROFILE%/RaylibCaptures; kept under 10 GB (RAYLIBGAME_CAPTURE_MAX_GB),
//                        oldest captures deleted first
//   --record-compare <n> diff each frame against an earlier capture (share of pixels that differ)
//   --frames N           frames to capture (default 120)
//   --record-every K     one frame every K rendered frames (default 1)
//   --record-start S     wait S seconds of play (and for the world to be built) before the first frame (default 0)
//   --record-hud 0|1     include the HUD and UI (default 1; 0 grabs the frame before them)
//   --record-continue    keep running after the last frame (default: exit when done)
//   --record-on-cue      hold the first frame until the scenario calls Recorder_Cue() (a test's event after a long setup,
//                        e.g. --test-dam-break cues it as the wall goes)
//   --gif                also <dir>/capture.gif (640 px wide)
//   --sheet              also <dir>/sheet.png
//   --fixed-dt <s>       simulate with this frame time whatever the real one (reproducible captures with --seed)
// Normal play: F10 starts and stops a recording into <root>/<YYYY-MM-DD>/<HHMMSS>/ (every 2nd frame until F10 again, at most 30 min; the GIF is its first 30 s; the sheet spans it all).
// When a recording ends it prints "record: <dir> frames=N gif=<path|-> sheet=<path|-> hash=<first>..<last>".

#include <functional>
#include <string>
#include <raylib.h>

bool Recorder_ParseArg(int argc, char** argv, int& i); // consumes one of the flags above (true when it did)
void Recorder_Start();                    // after the window opens: seeds, fixed dt
void Recorder_Update(bool worldReady);    // once a frame before drawing: start delay, F10
bool Recorder_WantsFrame(bool hud);       // this frame is captured at this point (hud = after the UI)
void Recorder_Capture();                  // grab the back buffer now (call where WantsFrame said so)
bool Recorder_Finished();                 // a command-line recording is done and the game should exit
void Recorder_Shutdown();                 // finish writing (frames still queued, GIF, sheet)
bool Recorder_Active();
void Recorder_Cue();                      // --record-on-cue: start capturing now
void Recorder_Screenshot();             // F2 (#579): the next finished frame (HUD included) to <root>/<date>/shot-<HHMMSS>/screenshot.png
// (#910) F2 also writes <shot dir>/lod.txt from this (called on the main thread with the screenshot image)
void Recorder_SetShotExtra(std::function<std::string(const Image&)> f);
std::string Recorder_LastDir();           // the last finished recording's folder ("" if none): bug reports attach it (#169)
void Recorder_FrameDone();               // end of every frame: counts rendered frames for --record-every
void Recorder_DrawIndicator();          // (#1119) F10: REC, a red dot and the time, after the capture (never in the recording)

#endif
