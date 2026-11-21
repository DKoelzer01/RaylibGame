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
        if (args.empty()) return "/locate <feature> [n]: " + Locate_FeatureList();
        const int n = args.size() > 1 ? std::max(1, atoi(args[1].c_str())) : 1;
        std::vector<LocateHit> hits;
        std::string why;
        if (!Locate_Find(body, from, args[0], n, hits, why)) return why;
        std::string reply;
        for (size_t i = 0; i < hits.size(); ++i) {
            float la, lo; Locate_LatLonOf(hits[i].at, la, lo);
            const std::string dist = hits[i].distance < 1000.0f ? TextFormat("%.0f m", hits[i].distance) : TextFormat("%.1f km", hits[i].distance / 1000.0f);
            reply += TextFormat("%s%s: %s %s (%.0f deg), lat %.4f lon %.4f", i ? "\n" : "", hits[i].what.c_str(), dist.c_str(),
                                Locate_CompassPoint(hits[i].bearing), hits[i].bearing, la, lo);
        }
        Nav_AddWaypoint(body.name, Vector3Normalize(hits[0].at), "/locate " + Lower(args[0]), false); // (the compass and the map show it)
        if (!why.empty()) reply += "\n" + why;
        return reply + "\n(a waypoint on the nearest; /tp " + Lower(args[0]) + " goes there)";
    }
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    if (!pf) return "Nowhere to go here (not a planet)";
    Vector3 to{}; float height = 0.0f; std::string where;
    char* end = nullptr;
    const float lat = args.empty() ? 0.0f : strtof(args[0].c_str(), &end);
    if (args.size() >= 2 && end && *end == 0) { // lat lon [height]
        if (fabsf(lat) > 90.0f) return "Latitude is -90 to 90";
        const float lon = (float)atof(args[1].c_str());
        height = args.size() > 2 ? std::max(0.0f, (float)atof(args[2].c_str())) : 0.0f;
        const Vector3 d = Locate_DirFromLatLon(lat, lon);
        to = Vector3Scale(d, pf->radius + pf->elevation(d));
        where = TextFormat("lat %.4f lon %.4f", lat, lon);
    } else if (!args.empty()) { // a feature: its nearest
        std::vector<LocateHit> hits; std::string why;
        if (!Locate_Find(body, from, args[0], 1, hits, why)) return why;
        to = hits[0].at; where = hits[0].what;
        height = args.size() > 1 ? std::max(0.0f, (float)atof(args[1].c_str())) : 0.0f;
    } else return "/tp <lat> <lon> [height] or /tp <feature> [height] (" + Locate_FeatureList() + ")";
    // On the ground (a cave feature: at it, underground) 0.3 m up so the feet don't start in the soil, or height m over it.
    const Vector3 up = Vector3Normalize(to);
    const Vector3 feet = Vector3Add(to, Vector3Scale(up, 0.3f + height));
    p.flyMode = height > 2.0f; // (creative: hover there rather than fall)
    p.placeAt(*gCmdWorld, body, (body.toWorldD(DVec3(feet)) - gCmdWorld->renderOrigin).toFloat());
    return "Teleported to " + where + (height > 0.0f ? std::string(TextFormat(", %.0f m up", height)) : std::string());
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
