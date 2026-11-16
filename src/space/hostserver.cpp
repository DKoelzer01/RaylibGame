#include "hostserver.h"
#include "events.h"
#include "joininfo.h"
#include "net.h"
#include "scene.h"
#include "world.h"

#include <raylib.h>
#include <cstdio>
#include <string>

namespace {
    bool listening = false;
    int joins = 0;
    uint32_t localClient = 0;
    std::vector<NetMessage> localInbox;
    double autosaveSeconds = 0.0, lastSave = -1.0;
    int autosaves = 0;

    void SendText(uint32_t conn, const std::string& text) { Net_Send(conn, text.data(), (uint32_t)text.size(), true); }
}

bool Host_Start(uint16_t port, std::string& status) {
    gNetRole = NetRole::Host; // the authority either way: the world runs and rolls here
    lastSave = GetTime();
    if (!Net_Available()) { status = "built without ENABLE_GNS: the world runs, nobody can connect"; return false; }
    std::string err;
    if (!Net_Init(err)) { status = "network init failed: " + err; return false; }
    if (!Net_Host(port)) { status = "could not listen on port " + std::to_string(port); return false; }
    listening = true;
    status = "listening on port " + std::to_string(port);
    return true;
}

void Host_SetAutosaveMinutes(double minutes) { autosaveSeconds = minutes > 0.0 ? minutes * 60.0 : 0.0; }

void Host_Update(Scene& world) {
    // Autosave (user, 2026-11-15: every 30 minutes by default; --autosave-minutes N, 0 = off). Wall-clock time: a
    // server's saves protect against crashes, whatever the game clock does.
    if (autosaveSeconds > 0.0 && lastSave >= 0.0 && GetTime() - lastSave >= autosaveSeconds) {
        lastSave = GetTime();
        const bool ok = SaveWorld(world);
        if (ok) ++autosaves;
        printf("host: autosave %d to %s %s\n", autosaves, WorldSavePath(), ok ? "done" : "FAILED");
    }
    if (!listening) return;
    Net_Poll();
    std::vector<NetMessage> in;
    Net_Receive(in);
    for (NetMessage& m : in) {
        if (localClient && m.conn == localClient) { localInbox.push_back(std::move(m)); continue; }
        const std::string text(m.data.begin(), m.data.end());
        if (text.rfind("join", 0) == 0) {
            SendText(m.conn, "J" + JoinInfo_Write(JoinInfo_FromWorld(world)));
            const EventsShower s = Events_GetShower();
            char buf[96];
            snprintf(buf, sizeof buf, "S %.3f %.3f %u", s.start, s.end, s.id);
            SendText(m.conn, buf);
            ++joins;
            printf("host: answered a join (connection %u, %d so far)\n", m.conn, joins);
        }
    }
}

bool Host_SaveOnShutdown(Scene& world) {
    const bool ok = SaveWorld(world);
    printf("host: saved on shutdown to %s %s\n", WorldSavePath(), ok ? "done" : "FAILED");
    return ok;
}

void Host_Stop() {
    if (listening) Net_Shutdown();
    listening = false;
    localClient = 0;
    localInbox.clear();
    gNetRole = NetRole::Offline;
}

bool Host_Listening() { return listening; }
int Host_Joins() { return joins; }
int Host_Autosaves() { return autosaves; }
void Host_SetLocalClient(uint32_t conn) { localClient = conn; }
std::vector<NetMessage> Host_TakeLocalClientMessages() { std::vector<NetMessage> out; out.swap(localInbox); return out; }
