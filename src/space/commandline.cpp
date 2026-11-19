#include "commandline.h"
#include "gameui.h"
#include "netclient.h"
#include "playerinput.h"
#include "terraintool.h"

#include "imgui.h"
#include <raylib.h>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <sstream>

namespace {
    bool open = false;
    bool focusNext = false;
    char text[256] = {};
    std::string lastReply;
    double replyUntil = 0.0;
    std::string Lower(std::string s) { for (char& c : s) c = (char)tolower((unsigned char)c); return s; }
}

bool CommandLine_Open() { return open; }

void GameMode_SetCreative(bool on) {
    terrainTool.creative = on; // (flying follows it, main.cpp; saved with the player's record)
    if (!on && GameUI_Screen() == GameScreen::Catalog) GameUI_ShowScreen(GameScreen::Inventory);
}

std::string CommandLine_Run(const std::string& typed) {
    std::string line = typed;
    while (!line.empty() && (line[0] == '/' || line[0] == ' ')) line.erase(0, 1);
    std::istringstream in(line);
    std::string cmd, arg;
    in >> cmd >> arg;
    cmd = Lower(cmd); arg = Lower(arg);
    if (cmd == "gamemode" || cmd == "gm") {
        if (arg == "creative" || arg == "c" || arg == "1") { GameMode_SetCreative(true); return "Game mode: creative (build, craft and place for free; F flies; the Catalog screen)"; }
        if (arg == "survival" || arg == "s" || arg == "0") { GameMode_SetCreative(false); return "Game mode: survival"; }
        return terrainTool.creative ? "Game mode: creative (/gamemode survival to leave it)" : "Game mode: survival (/gamemode creative for free building)";
    }
    if (cmd == "help" || cmd == "?") return "Commands: /gamemode creative|survival, /help";
    return "Unknown command: /" + cmd + " (/help lists them)";
}

void CommandLine_Update(Scene& world) {
    if (open || GameUI_InventoryOpen() || GameUI_TextInputActive()) return;
    // Enter opens it empty; / opens it with the slash typed.
    if (PlayerInput_BindingPressed(Act::Chat)) { open = true; focusNext = true; text[0] = 0; }
    else if (IsKeyPressed(KEY_SLASH)) { open = true; focusNext = true; strcpy(text, "/"); }
    if (open) { while (GetCharPressed() != 0) {} } // (the key that opened it types nothing)
}

void CommandLine_Draw() {
    const float s = GameUI_Scale();
    if (!open) {
        if (!lastReply.empty() && GetTime() < replyUntil) { // the last answer, for a while
            ImGui::SetNextWindowPos(ImVec2(20 * s, ImGui::GetIO().DisplaySize.y - 200 * s));
            ImGui::SetNextWindowBgAlpha(0.55f);
            ImGui::Begin("##cmdreply", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoInputs | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing);
            ImGui::TextUnformatted(lastReply.c_str());
            ImGui::End();
        }
        return;
    }
    const ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowPos(ImVec2(20 * s, io.DisplaySize.y - 160 * s));
    ImGui::SetNextWindowSize(ImVec2(std::min(760.0f * s, io.DisplaySize.x - 40 * s), 0.0f));
    ImGui::SetNextWindowBgAlpha(0.8f);
    ImGui::Begin("##cmdline", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoSavedSettings);
    const bool joined = gNetClient && gNetClient->state == NetClient::State::Joined;
    ImGui::TextDisabled(joined ? "Say something, or a /command (Enter sends, Esc closes)" : "A /command (Enter runs it, Esc closes; /help lists them)");
    if (focusNext) { ImGui::SetKeyboardFocusHere(); focusNext = false; }
    ImGui::SetNextItemWidth(-1.0f);
    const bool sent = ImGui::InputText("##cmd", text, sizeof(text), ImGuiInputTextFlags_EnterReturnsTrue);
    if (sent) {
        const std::string line = text;
        open = false;
        if (!line.empty()) {
            if (joined) NetClient_SendChat(*gNetClient, line); // (the server answers commands it allows, and says chat to everyone)
            else if (line[0] == '/') { lastReply = CommandLine_Run(line); replyUntil = GetTime() + 6.0; }
            else { lastReply = "Nobody else is here: lines starting with / are commands (/help)"; replyUntil = GetTime() + 5.0; }
        }
    }
    if (ImGui::IsKeyPressed(ImGuiKey_Escape)) open = false;
    ImGui::End();
}

void CommandLine_TestOpen(const char* t) { open = true; focusNext = true; snprintf(text, sizeof(text), "%s", t ? t : ""); }
