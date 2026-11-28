// net/alt-recipe (#2072 slice 2b): a joined player chooses a ruin cache's alternate recipe for a machine (MachineAction
// 14), and the host checks the requester's team has it. The game hosts (--host; needs a GameNetworkingSockets build, else
// one SKIP check passes); bob and carol join, each in a team of one, and stand by a coke oven on the host. bob's team has
// packed_coke, carol's doesn't: carol's request is refused (the oven keeps its base recipes, she hears why), bob's sets
// it, and bob's "" request goes back to the base recipes.
#include "scenario.h"
#include "buildtool.h"
#include "caches.h"
#include "factory.h"
#include "hostserver.h"
#include "machineact.h"
#include "netclient.h"
#include "net.h"
#include "playerrecord.h"
#include "plots.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>
#include <thread>

namespace {
// A frame at most every 1/60 s of wall time (as net/late-join-lag: a client polled thousands of times a second floods).
void RealTime() {
    static double next = 0.0;
    const double now = GetTime();
    if (next > now) std::this_thread::sleep_for(std::chrono::duration<double>(next - now));
    next = std::max(now, next) + 1.0 / 60.0;
}
constexpr uint16_t kPort = 27956;
struct State {
    NetClient bob, carol;
    uint64_t bodyNet = 0;
    uint32_t plotId = 0;
    int oven = -1;
    double t0 = 0.0;
    bool skip = false;
} S;
bool Joined(const NetClient& n) { return n.state == NetClient::State::Joined; }
void UpdateAll() { RealTime(); for (NetClient* n : { &S.bob, &S.carol }) if (n->state != NetClient::State::Idle) NetClient_Update(*n); }
bool Heard(const NetClient& n, const char* text) { for (const MsgChat& m : n.chat) if (m.text.find(text) != std::string::npos) return true; return false; }
std::string OvenAlt(ScenarioCtx& c) {
    Factory::View(c.body->id, S.plotId, S.oven); Factory::WaitTicks(2); // (View follows one machine: the first call after a switch is stale)
    return Factory::View(c.body->id, S.plotId, S.oven).alt;
}
void Send(NetClient& n, const std::string& id) { NetClient_MachineAction(n, S.bodyNet, S.plotId, S.oven, kMachineActChooseAlt, id, 0, 0, 0); }
}

SCENARIO("net/alt-recipe") {
    s.tags = { "multiplayer", "net" }; // (not "factory": the runner's part 1 takes that tag, part 2 takes "net")
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--host", "saves/scenario-net-alt-recipe.sav", "--port", std::to_string(kPort) };
    s.timeout = 120.0;
    s.realTime = true; // (#1972) joins over real sockets
    s.step("bob and carol join, each in a team of one", 40.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        if (!Host_Listening()) { // (the harness has no SKIP verdict: one passing check, said plainly)
            c.check("SKIP: not listening (a build without GameNetworkingSockets)", true);
            S.skip = true; return true;
        }
        if (S.bob.state == NetClient::State::Idle) {
            for (const char* n : { "bob", "carol" }) Teams_CreateSolo(PlayerId_FromName(n), n); // (fresh teams: the save outlives a run)
            const std::string addr = TextFormat("127.0.0.1:%d", kPort);
            NetClient_Connect(S.bob, addr.c_str(), "bob");
            NetClient_Connect(S.carol, addr.c_str(), "carol");
            S.t0 = c.clock;
        }
        UpdateAll();
        c.waitingFor(TextFormat("bob %s, carol %s", NetClient_StateName(S.bob.state), NetClient_StateName(S.carol.state)));
        if (!(Joined(S.bob) && Joined(S.carol))) {
            if (c.clock - S.t0 > 30.0) c.fail("bob and carol not both joined after 30 s");
            return false;
        }
        const int bobTeam = Teams_Of(PlayerId_FromName("bob")), carolTeam = Teams_Of(PlayerId_FromName("carol"));
        c.check("two teams", bobTeam && carolTeam && bobTeam != carolTeam, TextFormat("bob %d, carol %d", bobTeam, carolTeam));
        Caches_ResetTeams();
        Caches_AddTeamAlt(bobTeam, "packed_coke"); // as a cache's alt reward does (Caches_ApplyTeam)
        return true;
    });
    s.step("a coke oven by them on the host", 10.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        UpdateAll();
        Host_TestPutPlayerAt(c.world, S.carol.welcome.session, S.bob.welcome.session);
        uint64_t body = 0; DVec3 local;
        if (!c.check("bob is in the world on the host", Host_PlayerLocal(c.world, S.bob.welcome.session, body, local))) return true;
        Planetoid* b = Planetoid::findNet(body);
        if (!c.check("on the scenario's body", b == c.body)) return true;
        S.bodyNet = body;
        terrainTool.creative = true;
        const Vector3 up = Vector3Normalize(local.toFloat());
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, up));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(up, Vector3Length(local.toFloat())), Vector3Scale(east, 3.0f)));
        const Vector3 ground = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        Plot* plot = c.body->plots->plotFor(ground, 0.0f, true);
        if (!c.check("coke oven placed", BuildTool_Place(*c.body, *plot, { PartType::CokeOven, { 0, 0, 0 }, 0 }, false))) return true;
        S.plotId = plot->id;
        S.oven = plot->partAt({ 0, 0, 0 });
        Factory::Sync(c.world);
        Factory::WaitTicks(3);
        c.check("no alternate chosen yet", OvenAlt(c).empty());
        Send(S.carol, "packed_coke");
        S.t0 = c.t;
        return true;
    });
    s.step("carol's team hasn't found it: refused", 10.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        UpdateAll();
        c.waitingFor("carol hears the refusal");
        if (!Heard(S.carol, "hasn't found that recipe")) {
            if (c.t - S.t0 > 5.0) c.fail(TextFormat("no refusal for carol after 5 s (%zu chat lines; last: %s)", S.carol.chat.size(), S.carol.chat.empty() ? "-" : S.carol.chat.back().text.c_str()));
            return false;
        }
        c.check("the oven keeps its base recipes", OvenAlt(c).empty(), OvenAlt(c).c_str());
        Send(S.bob, "packed_coke");
        S.t0 = c.t;
        return true;
    });
    s.step("bob's team has it: the oven runs packed coke", 10.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        UpdateAll();
        const std::string alt = OvenAlt(c);
        c.waitingFor(TextFormat("the oven's alternate '%s'", alt.c_str()));
        if (alt != "packed_coke") {
            if (c.t - S.t0 > 5.0) c.fail(TextFormat("bob's choice not applied after 5 s (bob's last chat: %s)", S.bob.chat.empty() ? "-" : S.bob.chat.back().text.c_str()));
            return false;
        }
        c.check("bob heard no refusal", !Heard(S.bob, "hasn't found") && !Heard(S.bob, "out of reach") && !Heard(S.bob, "no such recipe"));
        Send(S.bob, "");
        S.t0 = c.t;
        return true;
    });
    s.step("bob's \"\" goes back to the base recipes", 10.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        UpdateAll();
        const std::string alt = OvenAlt(c);
        c.waitingFor(TextFormat("the oven's alternate '%s'", alt.c_str()));
        if (!alt.empty()) {
            if (c.t - S.t0 > 5.0) c.fail("bob's \"\" not applied after 5 s");
            return false;
        }
        c.eq("sends the transport refused", Net_SendFailures(), 0);
        for (NetClient* n : { &S.bob, &S.carol }) NetClient_Disconnect(*n);
        return true;
    });
}
