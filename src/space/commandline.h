#ifndef COMMANDLINE_H
#define COMMANDLINE_H

// The command line (#823): Enter (or /) opens a line at the bottom left; Enter sends it, Esc closes it. Lines starting
// with "/" are commands: in single player this game runs them (/gamemode creative|survival, /help); joined to a server,
// every line goes to the server (chat, and commands it allows: /gamemode for ops). A full chat box is #342.

#include <string>

class Scene;

bool CommandLine_Open();                    // the line is open (gameplay keys don't apply)
void CommandLine_Update(Scene& world);      // per frame (before ImGui): opens on its key
void CommandLine_Draw();                    // inside ImGui
// Runs a command line as typed (with or without the "/"): the reply ("" none). Single player's commands.
std::string CommandLine_Run(const std::string& line);
// Game mode (#823): creative (build, craft and place for free; the Catalog screen) or survival. Set by /gamemode.
void GameMode_SetCreative(bool on);
void CommandLine_TestOpen(const char* text); // tests: open the line with this typed
void CommandLine_Submit(const std::string& line); // what Enter does with the line (closes it; runs or sends it)
// (#1563) Single player's history (typed lines and their answers, the oldest first; up to 200), shown above the open line.
int CommandLine_HistoryCount();
std::string CommandLine_HistoryLine(int k);

#endif
