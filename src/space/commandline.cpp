#include "commandline.h"
#include "gameui.h"
#include "netclient.h"
#include "playerinput.h"
#include "terraintool.h"
#include "player.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <cstdio>

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
    Scene* gCmdWorld = nullptr; // the scene CommandLine_Update last ran in (for /locate)

    // (#1181) /locate KIND: the nearest cave entrance of a kind within 5 km of the player - its distance and compass
    // bearing (from north along the ground). Kinds: cave (any entrance), hillside, sinkhole (and karst sinks, skylights:
    // every shaft), resurgence, chamber.
    std::string Locate(const std::string& what) {
        int kind = -1;
        const char* label = "cave entrance";
        if (what == "hillside" || what == "mouth") { kind = 0; label = "hillside cave mouth"; }
        else if (what == "sinkhole" || what == "sink" || what == "shaft") { kind = 1; label = "sinkhole or shaft"; }
        else if (what == "resurgence") { kind = 2; label = "resurgence"; }
        else if (what == "chamber") { kind = 3; label = "cave chamber"; }
        else if (!what.empty() && what != "cave" && what != "entrance") return "/locate cave|hillside|sinkhole|resurgence|chamber";
        Planetoid* b = player.nearest;
        const PlanetField* pf = b ? dynamic_cast<const PlanetField*>(b->field.get()) : nullptr;
        if (!pf || !gCmdWorld) return "Nothing to search here (stand on a planet)";
        const Vector3 local = b->toLocal((player.position - gCmdWorld->renderOrigin).toFloat());
        const float len = Vector3Length(local);
        if (len < 1.0f) return "Nothing to search here";
        const Vector3 up = Vector3Scale(local, 1.0f / len);
        Vector3 at{}, out{};
        if (!pf->nearestCaveMouth(up, 5000.0f, at, out, kind)) return std::string("No ") + label + " within 5 km";
        const Vector3 atS = Vector3Scale(Vector3Normalize(at), pf->radius), meS = Vector3Scale(up, pf->radius);
        const float dist = Vector3Distance(atS, meS);
        Vector3 to = Vector3Subtract(atS, meS);
        to = Vector3Subtract(to, Vector3Scale(up, Vector3DotProduct(to, up)));
        Vector3 north = Vector3Subtract({ 0, 1, 0 }, Vector3Scale(up, up.y));                // toward the body's +Y pole
        if (Vector3Length(north) < 1e-4f) north = Vector3Subtract({ 0, 0, 1 }, Vector3Scale(up, up.z));
        north = Vector3Normalize(north);
        const Vector3 east = Vector3CrossProduct(north, up);
        float bearing = atan2f(Vector3DotProduct(to, east), Vector3DotProduct(to, north)) * RAD2DEG;
        if (bearing < 0.0f) bearing += 360.0f;
        static const char* dirs[8] = { "N", "NE", "E", "SE", "S", "SW", "W", "NW" };
        const int oct = (int)floorf((bearing + 22.5f) / 45.0f) % 8;
        char buf[160];
        snprintf(buf, sizeof buf, "Nearest %s: %.0f m %s (bearing %.0f deg)", label, dist, dirs[oct], bearing);
        return buf;
    }
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
    if (cmd == "locate" || cmd == "loc") return Locate(arg);
    if (cmd == "help" || cmd == "?") return "Commands: /gamemode creative|survival, /locate cave|hillside|sinkhole|resurgence|chamber, /help";
    return "Unknown command: /" + cmd + " (/help lists them)";
}

void CommandLine_Update(Scene& world) {
    gCmdWorld = &world;
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
