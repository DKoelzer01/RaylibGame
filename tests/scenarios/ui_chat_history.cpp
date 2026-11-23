// ui/chat-history (#1563): in single player the command line kept only its last answer, for a few seconds, so a long help
// text was gone before it was read. Lines are sent as Enter sends them (/help, an unknown command, plain text); the history
// must hold each typed line and its answer, and show them (a scrollable pane) when the line is opened again.
#include "scenario.h"
#include "commandline.h"
#include "world.h"
#include "imgui.h"
#include "imgui_internal.h"
#include <cstring>
#include <string>

SCENARIO("ui/chat-history") {
    s.tags = { "world", "ui" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.step("settle", 30.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T >= 1.0 && WorldIdle(); });
    s.step("three lines sent, then the line opened again", 10.0, [](ScenarioCtx& c) {
        static int n0 = 0;
        if (c.t < 0.05) {
            n0 = CommandLine_HistoryCount();
            CommandLine_Submit("/help");
            CommandLine_Submit("/nosuchcommand");
            CommandLine_Submit("hello");
            CommandLine_TestOpen("");
            return false;
        }
        if (c.t < 0.5) return false;
        const int n = CommandLine_HistoryCount() - n0;
        c.eq("lines kept (3 typed, 3 answers)", n, 6);
        c.check("the typed /help kept", CommandLine_HistoryLine(n0) == "/help");
        c.check("its answer kept (the command list)", CommandLine_HistoryLine(n0 + 1).find("/gamemode") != std::string::npos, CommandLine_HistoryLine(n0 + 1));
        c.check("an error kept", CommandLine_HistoryLine(n0 + 3).find("Unknown command") != std::string::npos, CommandLine_HistoryLine(n0 + 3));
        const ImGuiWindow* w = ImGui::FindWindowByName("##cmdline");
        bool pane = false;
        if (w) for (const ImGuiWindow* ch : w->DC.ChildWindows) pane = pane || strstr(ch->Name, "cmdhistory") != nullptr;
        c.check("the history pane shows with the line open", pane);
        if (w) c.check("the line's window on screen", w->Pos.y >= 0.0f && w->Pos.y + w->Size.y <= ImGui::GetIO().DisplaySize.y + 0.5f,
                       TextFormat("y %.0f, height %.0f", w->Pos.y, w->Size.y));
        c.grab("chat-history"); // (RAYLIBGAME_SLIDE_SHOTS: a picture to look at)
        CommandLine_Submit(""); // (closes it)
        return true;
    });
}
