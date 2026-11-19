#ifndef BUGREPORT_H
#define BUGREPORT_H

// In-game bug reports (#169): a window from the pause menu ("Report a bug"): a title, what happened, tags (they become the
// issue's labels) and a severity, with attachments: the screenshot taken as the game paused, log.txt, the last F10
// capture (contact sheet and GIF) and, optionally, the save. Sent to the bug relay (tools/bugrelay: tailnet-only on
// deweyserver for now), which stores the files in the Claude media store and opens a GitHub issue. The relay URL is
// RAYLIBGAME_BUG_RELAY (default http://deweyserver:8765/report). Sending runs on a thread with the system's curl.exe;
// the result comes back as a toast. (#617) A report that fails waits in saves/bug-outbox (a notice says so until
// dismissed) and is sent again at start-up and every 5 minutes; the pause menu shows how many wait.

class Scene;

void BugReport_CaptureScene();        // the world as the game paused (call right after the scene is drawn, before the menu)
void BugReport_Open(Scene& world);    // the window (from the pause menu)
bool BugReport_IsOpen();
void BugReport_Draw();                // inside the ImGui frame
void BugReport_Update();              // once a frame: a finished send's toast
void BugReport_TestSend(Scene& world); // --test-bug-report: a scripted report (screenshot if one was taken, log.txt)
bool BugReport_Sending();
int BugReport_Unsent();             // (#617) reports waiting in the outbox (saves/bug-outbox: sent again at start-up and every 5 min)

#endif
