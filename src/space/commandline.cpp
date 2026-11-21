#include "commandline.h"
#include "gameui.h"
#include "netclient.h"
#include "playerinput.h"
#include "terraintool.h"
#include "locate.h"
#include "navigation.h"
#include "playersession.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <vector>

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

// (#1171) /locate <feature> [n]: the nearest n, distance and bearing, a waypoint on the nearest. /tp <lat> <lon> [height] or
// /tp <feature> [height]: there (on the ground, or height metres over it); the loading screen holds the player until the
// place has loaded. Creative only (finding things is the game in survival). Joined, the line goes to the server, which
// doesn't run them yet (it owns the player's position).
namespace {
Scene* gCmdWorld = nullptr; // (CommandLine_Update's: the commands need the player's body)
std::string LocateCommand(const std::string& cmd, const std::string& line) {
    std::istringstream in(line);
    std::string word; in >> word; // (the command)
    std::vector<std::string> args;
    while (in >> word) args.push_back(word);
    Player& p = LocalSession().player;
    if (!gCmdWorld || !p.nearest) return "No world here";
    if (!terrainTool.creative) return "/" + cmd + " is for creative mode (/gamemode creative)";
    Planetoid& body = *p.nearest;
    const Vector3 from = body.toLocalD(p.position).toFloat();
    if (cmd == "locate") {
        LocateHit nearest;
        nearest.distance = -1.0f;
        const std::string reply = Locate_LocateCommand(body, from, args, &nearest);
        if (nearest.distance < 0.0f) return reply;
        Nav_AddWaypoint(body.name, Vector3Normalize(nearest.at), "/locate " + Lower(args[0]), false); // (the compass and the map show it)
        return reply + "\n(a waypoint on the nearest; /tp " + Lower(args[0]) + " goes there)";
    }
    Vector3 feet{}; bool hover = false; std::string reply;
    if (!Locate_TpTarget(body, from, args, feet, hover, reply)) return reply;
    p.flyMode = hover;
    p.placeAt(*gCmdWorld, body, (body.toWorldD(DVec3(feet)) - gCmdWorld->renderOrigin).toFloat()); // (its loading screen holds the player)
    return reply;
}
}

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
    if (cmd == "locate" || cmd == "tp") return LocateCommand(cmd, line);
    if (cmd == "help" || cmd == "?") return "Commands: /gamemode creative|survival, /locate <feature> [n], /tp <lat> <lon> [height] | /tp <feature>, /help";
    return "Unknown command: /" + cmd + " (/help lists them)";
}

void CommandLine_Update(Scene& world) {
    gCmdWorld = &world; // (#1171: /locate and /tp)
    if (open || GameUI_InventoryOpen() || GameUI_TextInputActive()) return;
    // Enter opens it empty; / opens it with the slash typed.
    if (PlayerInput_BindingPressed(Act::Chat)) { open = true; focusNext = true; text[0] = 0; }
    else if (IsKeyPressed(KEY_SLASH)) { open = true; focusNext = true; strcpy(text, "/"); }
    if (open) { while (GetCharPressed() != 0) {} } // (the key that opened it types nothing)
}

void CommandLine_Draw() {
    const float s = GameUI_HudScale(); // (#925) the HUD's size, not the menus'
    if (!open) {
        if (!lastReply.empty() && GetTime() < replyUntil) { // the last answer, for a while
            ImGui::SetNextWindowPos(ImVec2(20 * s, ImGui::GetIO().DisplaySize.y - 200 * s));
            ImGui::SetNextWindowBgAlpha(0.55f);
            ImGui::Begin("##cmdreply", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoInputs | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing);
            ImGui::SetWindowFontScale(GameUI_HudScale() / GameUI_Scale()); // (#925: the style's fonts are at the menus' 2x)
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
    ImGui::SetWindowFontScale(GameUI_HudScale() / GameUI_Scale()); // (#925: the style's fonts are at the menus' 2x)
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
