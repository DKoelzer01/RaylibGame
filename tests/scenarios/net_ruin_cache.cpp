// net/ruin-cache (#2072): a ruin cache opened by a joined client goes through the host. The game hosts (--host; needs a
// GameNetworkingSockets build, else one SKIP check passes); bob and carol join. The host puts a small cache at bob's feet.
// bob asks to open it with no power cell: refused, still sealed. With one (the host gives it): the host opens it, spends
// the cell and gives the stub pool's 4 stone; bob hears the result lines, and his game's result panel takes them (fed
// here as main.cpp's chat loop does); carol hears the part edit that marks it opened, and no result. carol then stands at
// the cache with a power cell of her own: the host refuses her (opened already) and she keeps the cell.
// The clients run in this process, so they don't apply the world (the host's is the real one): the checks read the host's
// plot and inventories, and what each client received.
#include "scenario.h"
#include "hostserver.h"
#include "machineact.h"
#include "netclient.h"
#include "net.h"
#include "plots.h"
#include "ruincache.h"
#include "world.h"
#include "worldtypes.h"
#include <algorithm>
#include <string>
#include <thread>

namespace {
// A frame at most every 1/60 s of wall time (as net/team-store: a client polled thousands of times a second floods).
void RealTime() {
    static double next = 0.0;
    const double now = GetTime();
    if (next > now) std::this_thread::sleep_for(std::chrono::duration<double>(next - now));
    next = std::max(now, next) + 1.0 / 60.0;
}
constexpr uint16_t kPort = 27956;
struct State {
    NetClient bob, carol;
    uint64_t body = 0;
    uint32_t plot = 0;
    Int3 cell{ 0, 0, 0 }; // the cache's cell: the one bob stands in (the host's 4 m reach check)
    DVec3 feet;           // bob's spot on the host when the cache was placed (body-local)
    double t0 = 0.0;
    bool skip = false, asked = false;
} S;
bool Joined(const NetClient& n) { return n.state == NetClient::State::Joined; }
void UpdateAll() { RealTime(); for (NetClient* n : { &S.bob, &S.carol }) if (n->state != NetClient::State::Idle) NetClient_Update(*n); }
Plot* ThePlot() { Planetoid* b = Planetoid::findNet(S.body); return b && b->plots ? b->plots->byId(S.plot) : nullptr; }
int CacheIndex() { Plot* p = ThePlot(); return p ? p->partAt(S.cell) : -1; }
bool Opened() { Plot* p = ThePlot(); const int i = CacheIndex(); return i >= 0 && CacheOpened(p->parts[i]); }
// The request, from the cache's side. These clients don't play, so they send the host no player states: the host's copy
// of the player has no body and a stale place. A real client's states keep both current; here the host puts them there.
void Ask(ScenarioCtx& c, NetClient& n) {
    Host_TeleportPlayer(c.world, n.welcome.session, S.body, S.feet, false);
    NetClient_MachineAction(n, S.body, S.plot, CacheIndex(), kMachineActOpenCache, std::string(), 0, 0, 0); }
bool Heard(const NetClient& n, const std::string& text) { return std::any_of(n.chat.begin(), n.chat.end(), [&](const MsgChat& m) { return m.from == 0 && m.text == text; }); }
bool HeardStart(const NetClient& n, const std::string& head) { return std::any_of(n.chat.begin(), n.chat.end(), [&](const MsgChat& m) { return m.from == 0 && m.text.rfind(head, 0) == 0; }); }
// What the server said to this client (a failure message: why a wait timed out).
std::string Said(const NetClient& n) { std::string s; for (const MsgChat& m : n.chat) if (m.from == 0) s += "[" + m.text + "] "; return s.empty() ? "nothing" : s; }
// The opening's part edit: a cache placed into our plot with the opened bit.
bool SawOpened(const NetClient& n) {
    return std::any_of(n.plotLog.begin(), n.plotLog.end(), [](const MsgPlotEdit& e) {
        return e.plot == S.plot && e.op == NetPartOp::Place && e.part.type == (uint8_t)PartType::RuinCache && (e.part.material & kCacheOpened);
    });
}
}

SCENARIO("net/ruin-cache") {
    s.tags = { "multiplayer", "net", "adventure" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--host", "saves/scenario-net-ruin-cache.sav", "--port", std::to_string(kPort) };
    s.timeout = 90.0;
    s.step("bob and carol join; a sealed cache at bob's feet", 40.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        if (!Host_Listening()) { // (the harness has no SKIP verdict: one passing check, said plainly)
            c.check("SKIP: not listening (a build without GameNetworkingSockets)", true);
            S.skip = true; return true;
        }
        if (S.bob.state == NetClient::State::Idle) {
            S.bob.keepPlotLog = S.carol.keepPlotLog = true;
            const std::string addr = TextFormat("127.0.0.1:%d", kPort);
            NetClient_Connect(S.bob, addr.c_str(), "bob");
            NetClient_Connect(S.carol, addr.c_str(), "carol");
            S.t0 = c.t;
        }
        UpdateAll();
        c.waitingFor(TextFormat("bob %s, carol %s", NetClient_StateName(S.bob.state), NetClient_StateName(S.carol.state)));
        DVec3& feet = S.feet;
        if (!(Joined(S.bob) && Joined(S.carol) && Host_PlayerLocal(c.world, S.bob.welcome.session, S.body, feet))) {
            if (c.t - S.t0 > 30.0) c.fail("bob and carol not both joined after 30 s");
            return false;
        }
        Planetoid* b = Planetoid::findNet(S.body);
        if (!c.check("bob's body on the host", b != nullptr)) return true;
        if (!b->plots) b->plots = std::make_unique<PlotSet>();
        // (the host's save outlives a run on the runner: a fresh plot each time, at bob's feet now)
        Plot* plot = b->plots->plotFor(feet.toFloat(), 0.0f, true);
        S.plot = plot->id;
        S.cell = plot->cellAt(plot->toPlot(feet));
        const int i = plot->partAt(S.cell);
        if (i >= 0) plot->remove(i);
        c.check("the cache placed on the host", plot->place({ PartType::RuinCache, S.cell, 0, (uint8_t)CacheKind::Small }) >= 0);
        for (const char* item : { "power_cell", "stone" }) {
            const int have = Host_ItemCount(S.bob.welcome.session, item);
            if (have > 0) Host_TestGive(S.bob.welcome.session, item, -have);
            const int carolHas = Host_ItemCount(S.carol.welcome.session, item);
            if (carolHas > 0) Host_TestGive(S.carol.welcome.session, item, -carolHas);
        }
        S.bob.chat.clear(); S.carol.chat.clear();
        Ask(c, S.bob); // no power cell yet
        S.t0 = c.t;
        return true;
    });
    s.step("no power cell: refused, still sealed", 10.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        UpdateAll();
        c.waitingFor("the host's answer to bob");
        if (!Heard(S.bob, "Sealed: its lock needs a power cell")) {
            if (c.t - S.t0 > 5.0) c.fail("bob heard no refusal 5 s after asking without a power cell; heard: " + Said(S.bob));
            return false;
        }
        c.check("still sealed on the host", !Opened());
        c.check("no result lines", !HeardStart(S.bob, "Opened the "));
        c.check("bob gets a power cell", Host_TestGive(S.bob.welcome.session, "power_cell", 1));
        S.asked = false; S.t0 = c.t;
        return true;
    });
    s.step("with a power cell: the host opens it; bob gets the result, carol sees it opened", 15.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        UpdateAll();
        if (!S.asked && c.t - S.t0 > 0.5) { Ask(c, S.bob); S.asked = true; S.t0 = c.t; } // (after the give reached the host's copy)
        c.waitingFor("bob's result and carol's part edit");
        if (!S.asked || !(Heard(S.bob, "  4 stone") && SawOpened(S.carol))) {
            if (S.asked && c.t - S.t0 > 10.0) c.fail(std::string(TextFormat("10 s after asking: bob's result %d, carol saw it opened %d, opened on the host %d", (int)Heard(S.bob, "  4 stone"), (int)SawOpened(S.carol), (int)Opened())) + "; bob heard: " + Said(S.bob));
            return false;
        }
        c.check("opened on the host", Opened());
        c.eq("bob's power cell spent (the host's copy)", Host_ItemCount(S.bob.welcome.session, "power_cell"), 0);
        c.eq("bob's 4 stone (the host's copy)", Host_ItemCount(S.bob.welcome.session, "stone"), 4);
        c.check("bob heard the result's title", Heard(S.bob, "Opened the sealed cache:"));
        c.check("bob saw it opened too", SawOpened(S.bob));
        c.check("carol heard no result", !HeardStart(S.carol, "Opened the ") && !Heard(S.carol, "  4 stone"));
        // bob's game: its chat loop gives each line from the server to the result panel (main.cpp)
        int taken = 0;
        for (const MsgChat& m : S.bob.chat) if (m.from == 0 && RuinCache_TakeHostLine(m.text)) ++taken;
        c.eq("bob's result lines the panel took (title + 1)", taken, 2);
        c.check("bob's result panel shows the reward", RuinCache_LastResult() == std::vector<std::string>{ "4 stone" });
        c.check("a refusal line isn't taken", !RuinCache_TakeHostLine("Sealed: its lock needs a power cell"));
        // carol comes to the cache (Ask puts her at bob's spot) with a power cell of her own
        c.check("carol gets a power cell", Host_TestGive(S.carol.welcome.session, "power_cell", 1));
        S.carol.chat.clear();
        S.asked = false; S.t0 = c.t;
        return true;
    });
    s.step("carol can't open it again: refused, her cell kept", 15.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        UpdateAll();
        if (!S.asked && c.t - S.t0 > 0.5) { Ask(c, S.carol); S.asked = true; S.t0 = c.t; }
        c.waitingFor("the host's answer to carol");
        if (!S.asked || !Heard(S.carol, "Opened already: empty")) {
            if (S.asked && c.t - S.t0 > 5.0) c.fail("carol heard no refusal 5 s after asking at the opened cache; heard: " + Said(S.carol));
            return false;
        }
        c.check("still opened", Opened());
        c.eq("carol's power cell kept", Host_ItemCount(S.carol.welcome.session, "power_cell"), 1);
        c.eq("no stone for carol", Host_ItemCount(S.carol.welcome.session, "stone"), 0);
        c.check("no result for carol", !HeardStart(S.carol, "Opened the "));
        c.eq("sends the transport refused", Net_SendFailures(), 0);
        for (NetClient* n : { &S.bob, &S.carol }) NetClient_Disconnect(*n);
        return true;
    });
}
