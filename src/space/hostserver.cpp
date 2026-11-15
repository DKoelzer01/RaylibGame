#include "hostserver.h"
#include "events.h"
#include "joininfo.h"
#include "net.h"
#include "scene.h"

#include <cstdio>
#include <string>

namespace {
    bool listening = false;
    int joins = 0;
    uint32_t localClient = 0;
    std::vector<NetMessage> localInbox;

    void SendText(uint32_t conn, const std::string& text) { Net_Send(conn, text.data(), (uint32_t)text.size(), true); }
}

bool Host_Start(uint16_t port, std::string& status) {
    gNetRole = NetRole::Host; // the authority either way: the world runs and rolls here
    if (!Net_Available()) { status = "built without ENABLE_GNS: the world runs, nobody can connect"; return false; }
    std::string err;
    if (!Net_Init(err)) { status = "network init failed: " + err; return false; }
    if (!Net_Host(port)) { status = "could not listen on port " + std::to_string(port); return false; }
    listening = true;
    status = "listening on port " + std::to_string(port);
    return true;
}

void Host_Update(Scene& world) {
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

void Host_Stop() {
    if (listening) Net_Shutdown();
    listening = false;
    localClient = 0;
    localInbox.clear();
    gNetRole = NetRole::Offline;
}

bool Host_Listening() { return listening; }
int Host_Joins() { return joins; }
void Host_SetLocalClient(uint32_t conn) { localClient = conn; }
std::vector<NetMessage> Host_TakeLocalClientMessages() { std::vector<NetMessage> out; out.swap(localInbox); return out; }
