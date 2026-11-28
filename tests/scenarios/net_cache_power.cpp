// net/cache-power (#2181): a joined client's sealed-cache prompt shows the host's lock charge. The game hosts (--host;
// needs a GameNetworkingSockets build, else one SKIP check passes); bob joins. The host puts a small cache (its default
// opener: a power lock) at bob's feet, a generator and a cable beside it, and burns coal. bob watches the cache, as his
// prompt does when he aims at it (Factory::View on a client asks the host to watch it): the host's views come back with
// the lock's size (energyMax) and its charge, as the host has it. Full, bob asks to open it: the host opens it and spends
// the charge.
// The client runs in this process (its factory is the host's), so the client side of the prompt (CacheCharged reading
// the host's view in remote mode) is the unit test Factory.AJoinedClientReadsACachesChargeFromTheHostsView; this checks
// what crosses the wire.
#include "scenario.h"
#include "cacheopeners.h"
#include "factory.h"
#include "hostserver.h"
#include "machineact.h"
#include "netclient.h"
#include "net.h"
#include "plots.h"
#include "ruincache.h"
#include "world.h"
#include "worldtypes.h"
#include <algorithm>
#include <cmath>
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
constexpr uint16_t kPort = 27957;
struct State {
    NetClient bob;
    uint64_t body = 0;
    uint32_t plot = 0;
    Int3 cell{ 0, 0, 0 }; // the cache's cell: the one bob stands in (the host's reach checks)
    DVec3 feet;           // bob's spot on the host when the cache was placed (body-local)
    int gen = -1;
    double t0 = 0.0;
    float seen = 0.0f;    // the highest charge bob's views showed
    bool skip = false, asked = false;
} S;
bool Joined(const NetClient& n) { return n.state == NetClient::State::Joined; }
void UpdateAll() { RealTime(); if (S.bob.state != NetClient::State::Idle) NetClient_Update(S.bob); }
Planetoid* TheBody() { return Planetoid::findNet(S.body); }
Plot* ThePlot() { Planetoid* b = TheBody(); return b && b->plots ? b->plots->byId(S.plot) : nullptr; }
int CacheIndex() { Plot* p = ThePlot(); return p ? p->partAt(S.cell) : -1; }
bool Opened() { Plot* p = ThePlot(); const int i = CacheIndex(); return i >= 0 && CacheOpened(p->parts[i]); }
float HostCharge() { Planetoid* b = TheBody(); Plot* p = ThePlot(); return b && p ? Factory::CacheCharged(b->id, p->id, CacheIndex()) : -1.0f; }
// bob's last view of the cache, as a charge (0..1); -1 if the last view isn't the cache's or carries no lock
float BobCharge() {
    const MsgMachineView& v = S.bob.lastMachineView;
    if (v.body != S.body || v.plot != S.plot || v.part != CacheIndex() || v.view.energyMax <= 0) return -1.0f;
    return std::min(1.0f, (float)v.view.energy / (float)v.view.energyMax);
}
// These clients don't play, so they send the host no player states: the host puts bob at the cache (its reach checks).
void BobAct(ScenarioCtx& c, int action) {
    Host_TeleportPlayer(c.world, S.bob.welcome.session, S.body, S.feet, false);
    NetClient_MachineAction(S.bob, S.body, S.plot, CacheIndex(), action, std::string(), 0, 0, 0);
}
std::string Said(const NetClient& n) { std::string s; for (const MsgChat& m : n.chat) if (m.from == 0) s += "[" + m.text + "] "; return s.empty() ? "nothing" : s; }
}

SCENARIO("net/cache-power") {
    s.tags = { "multiplayer", "net", "adventure", "power" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--host", "saves/scenario-net-cache-power.sav", "--port", std::to_string(kPort) };
    s.timeout = 120.0;
    s.step("bob joins; a sealed small cache at his feet, a generator burning coal beside it", 40.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        if (!Host_Listening()) { // (the harness has no SKIP verdict: one passing check, said plainly)
            c.check("SKIP: not listening (a build without GameNetworkingSockets)", true);
            S.skip = true; return true;
        }
        if (S.bob.state == NetClient::State::Idle) {
            NetClient_Connect(S.bob, TextFormat("127.0.0.1:%d", kPort), "bob");
            S.t0 = c.t;
        }
        UpdateAll();
        c.waitingFor(TextFormat("bob %s", NetClient_StateName(S.bob.state)));
        if (!(Joined(S.bob) && Host_PlayerLocal(c.world, S.bob.welcome.session, S.body, S.feet))) {
            if (c.t - S.t0 > 30.0) c.fail("bob not joined after 30 s");
            return false;
        }
        Planetoid* b = TheBody();
        if (!c.check("bob's body on the host", b != nullptr)) return true;
        if (!b->plots) b->plots = std::make_unique<PlotSet>();
        // (the host's save outlives a run on the runner: a fresh plot each time, at bob's feet now)
        Plot* plot = b->plots->plotFor(S.feet.toFloat(), 0.0f, true);
        S.plot = plot->id;
        S.cell = plot->cellAt(plot->toPlot(S.feet));
        const Int3 g{ S.cell.x - 3, S.cell.y, S.cell.z };
        for (int x = g.x; x <= S.cell.x + 1; ++x) for (const Int3 at : { Int3{ x, S.cell.y, S.cell.z }, Int3{ x, S.cell.y, S.cell.z - 1 } })
            if (const int i = plot->partAt(at); i >= 0) plot->remove(i);
        S.gen = plot->place({ PartType::Generator, g, 0 });
        for (int x = g.x; x <= S.cell.x + 1; ++x) plot->place({ PartType::Cable, { x, S.cell.y, S.cell.z - 1 }, 0 });
        const int cache = plot->place({ PartType::RuinCache, S.cell, 0, (uint8_t)CacheKind::Small }); // its default: power
        if (!c.check("placed on the host", S.gen >= 0 && cache >= 0)) return true;
        const CacheOpener* o = CacheOpenerOf(plot->parts[cache]);
        c.check("a small cache opens by power", o && o->kind == OpenerKind::Power);
        Factory::Sync(c.world);
        Factory::WaitTicks(3);
        Factory::Insert(b->id, S.plot, S.gen, "coal", 16);
        BobAct(c, 6); // bob aims at it: his prompt watches it
        S.bob.chat.clear();
        S.t0 = c.t;
        return true;
    });
    s.step("bob's views carry the lock and its charge, as the host has it", 20.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        UpdateAll();
        const float bob = BobCharge(), host = HostCharge();
        S.seen = std::max(S.seen, bob);
        c.waitingFor(TextFormat("bob's view %.3f (%d views), the host's %.3f", bob, S.bob.machineViews, host));
        if (!(bob > 0.0f && bob < 1.0f)) {
            if (c.t - S.t0 > 15.0) c.fail(std::string(TextFormat("no charging view after 15 s: bob's %.3f, the host's %.3f, %d views", bob, host, S.bob.machineViews)) + "; bob heard: " + Said(S.bob));
            return false;
        }
        const CacheOpener* o = CacheOpener_Find("small_power");
        c.eq("the view's lock size (EU)", S.bob.lastMachineView.view.energyMax, o ? o->eu : -1);
        c.check("bob's charge is the host's (within its 5-a-second views)", std::fabs(bob - host) < 0.1f, TextFormat("bob %.3f, host %.3f", bob, host));
        c.check("no refusal for the watch", Said(S.bob) == "nothing", Said(S.bob).c_str());
        Factory::SetSpeed(20);
        S.t0 = c.t;
        return true;
    });
    s.step("full on bob's view too", 40.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        UpdateAll();
        c.waitingFor(TextFormat("bob's view %.3f, the host's %.3f", BobCharge(), HostCharge()));
        if (BobCharge() < 1.0f) {
            if (c.t - S.t0 > 35.0) c.fail(TextFormat("bob's view not full after 35 s: %.3f (the host's %.3f)", BobCharge(), HostCharge()));
            return false;
        }
        Factory::SetSpeed(1);
        c.check("full on the host", HostCharge() >= 1.0f);
        S.t0 = c.t;
        return true;
    });
    s.step("bob opens it: the host opens it and spends the charge", 15.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        UpdateAll();
        if (!S.asked) { BobAct(c, kMachineActOpenCache); S.asked = true; S.t0 = c.t; }
        c.waitingFor("the cache opened on the host");
        if (!Opened()) {
            if (c.t - S.t0 > 10.0) c.fail("not opened 10 s after bob asked; bob heard: " + Said(S.bob));
            return false;
        }
        c.check("the charge is spent", HostCharge() == 0.0f);
        BobAct(c, 7); // bob looks away: no more views
        c.eq("sends the transport refused", Net_SendFailures(), 0);
        NetClient_Disconnect(S.bob);
        return true;
    });
}
