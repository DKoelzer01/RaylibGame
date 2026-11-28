// net/team-store (#1886): the host sends each player their team's progression store (research, explored map, blueprint
// library) and keeps two teams apart. The game hosts (--host; needs a GameNetworkingSockets build, else one SKIP check
// passes); bob and carol join, each in a team of one. The host's teams then research different nodes: each client hears
// only its own team's within the resend interval (5 s). dave joins late into bob's team (RAYLIBGAME_HOST_TEAM): he gets
// bob's research with his join. carol joins bob's team by the real commands (/team invite, /team accept): she gets the
// new team's store at once, and her old team's research doesn't come with her (the user, on #1886).
// The clients run in this process, so they don't apply the text (the host's state is the real one here): the checks read
// what each client received. Map pins (#1886): bob pins a spot; dave's late join brings it, carol's team doesn't see it
// until she joins it, and her removing it takes it from bob's and dave's stores at the next tick. Research from a client
// (#1886): bob queues a node his team sees and the host keeps it; carol's ask for a node her team can't see is refused.
#include "scenario.h"
#include "hostserver.h"
#include "netclient.h"
#include "net.h"
#include "playerrecord.h"
#include "research.h"
#include "survival.h"
#include "teamquests.h"
#include "world.h"
#include <algorithm>
#include <cstring>
#include <sstream>
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
constexpr uint16_t kPort = 27955;
struct State {
    NetClient bob, carol, dave;
    int bobTeam = 0, carolTeam = 0;
    int carolStoresBefore = 0;
    double t0 = 0.0, daveJoinedAt = -1.0;
    bool skip = false, invited = false, accepted = false;
} S;
bool Joined(const NetClient& n) { return n.state == NetClient::State::Joined; }
// The research section's "done a b c" line names the node.
bool Has(const NetClient& n, const char* node) {
    std::istringstream in(n.teamStoreText);
    bool research = false;
    for (std::string line; std::getline(in, line);) {
        if (line.rfind("section ", 0) == 0) research = line == "section research";
        else if (research && line.rfind("done ", 0) == 0) {
            std::istringstream words(line.substr(5));
            for (std::string w; words >> w;) if (w == node) return true;
        }
    }
    return false;
}
// (#1886) The pins section has the pin of this name.
bool HasPin(const NetClient& n, const char* name) {
    const size_t a = n.teamStoreText.find("section pins\n");
    if (a == std::string::npos) return false;
    std::istringstream in(n.teamStoreText.substr(a + 13, n.teamStoreText.find("section_end\n", a) - a - 13));
    for (std::string line; std::getline(in, line);) if (line.rfind("pin ", 0) == 0 && line.size() > strlen(name) && line.compare(line.size() - strlen(name), std::string::npos, name) == 0) return true;
    return false;
}
const Vector3 kPinDir = { 0.0f, 0.6f, 0.8f }; // body-local
void UpdateAll() { RealTime(); for (NetClient* n : { &S.bob, &S.carol, &S.dave }) if (n->state != NetClient::State::Idle) NetClient_Update(*n); }
}

SCENARIO("net/team-store") {
    s.tags = { "multiplayer", "net" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--host", "saves/scenario-net-team-store.sav", "--port", std::to_string(kPort) };
    s.env = { { "RAYLIBGAME_HOST_TEAM", "bob,dave" } };
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
            // The host's save outlives a run where saves/ persists (the runner): a previous run left carol in bob's team.
            // A fresh team of one each (new ids, so no research carried over either); dave joins bob's again by the env.
            for (const char* n : { "bob", "carol", "dave" }) Teams_CreateSolo(PlayerId_FromName(n), n);
            const std::string addr = TextFormat("127.0.0.1:%d", kPort);
            NetClient_Connect(S.bob, addr.c_str(), "bob");
            NetClient_Connect(S.carol, addr.c_str(), "carol");
            S.t0 = c.clock;
        }
        UpdateAll();
        c.waitingFor(TextFormat("bob %s, carol %s", NetClient_StateName(S.bob.state), NetClient_StateName(S.carol.state)));
        if (!(Joined(S.bob) && Joined(S.carol) && S.bob.teamStores > 0 && S.carol.teamStores > 0)) {
            if (c.clock - S.t0 > 30.0) c.fail("bob and carol not both joined with a team store after 30 s");
            return false;
        }
        S.bobTeam = Teams_Of(PlayerId_FromName("bob"));
        S.carolTeam = Teams_Of(PlayerId_FromName("carol"));
        c.check("the host made two teams", S.bobTeam && S.carolTeam && S.bobTeam != S.carolTeam, TextFormat("bob %d, carol %d", S.bobTeam, S.carolTeam));
        c.eq("bob's store is his team's", S.bob.teamStoreTeam, S.bobTeam);
        c.eq("carol's store is her team's", S.carol.teamStoreTeam, S.carolTeam);
        // The host's teams research different nodes (as labs would finish them).
        Research_Team(S.bobTeam).done.insert("smelting");
        Research_Team(S.carolTeam).done.insert("basic_tools");
        S.t0 = c.t;
        return true;
    });
    s.step("each hears only its own team's research", 15.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        UpdateAll();
        c.waitingFor(TextFormat("bob has smelting %d, carol has basic_tools %d", (int)Has(S.bob, "smelting"), (int)Has(S.carol, "basic_tools")));
        if (!(Has(S.bob, "smelting") && Has(S.carol, "basic_tools"))) {
            if (c.t - S.t0 > 10.0) c.fail("a team's research didn't reach its member within 10 s");
            return false;
        }
        c.note(TextFormat("heard after %.1f s (the host resends a changed store at most every 5 s)", c.t - S.t0));
        c.check("bob doesn't hear carol's team's research", !Has(S.bob, "basic_tools"));
        c.check("carol doesn't hear bob's team's research", !Has(S.carol, "smelting"));
        // (#1886) bob pins a spot on his map: a pin for his team (as a click on the map does in the game)
        c.check("bob's map pin is sent", NetClient_SendMapPin(S.bob, true, "terra", kPinDir, "Copper hill"));
        // (#1886) research from a client: bob queues a tier 3 node for his team (at tier 3: it shows in full); carol asks for
        // a tier 5 node her team can't see, which the host refuses
        TeamQuests_Book(S.bobTeam).tier = 3; TeamQuests_Book(S.carolTeam).tier = 3;
        c.check("bob's research queue change is sent", NetClient_SendResearchQueue(S.bob, (uint8_t)ResearchQueueOp::Add, "steelmaking"));
        c.check("carol's is sent", NetClient_SendResearchQueue(S.carol, (uint8_t)ResearchQueueOp::Add, "aluminium"));
        NetClient_Connect(S.dave, TextFormat("127.0.0.1:%d", kPort), "dave");
        S.t0 = c.clock;
        return true;
    });
    s.step("dave joins bob's team late and gets its research with his join", 30.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        UpdateAll();
        if (Joined(S.dave) && S.daveJoinedAt < 0.0) S.daveJoinedAt = c.t;
        c.waitingFor(TextFormat("dave %s, %d stores", NetClient_StateName(S.dave.state), S.dave.teamStores));
        if (!(Joined(S.dave) && S.dave.teamStores > 0)) {
            if (c.clock - S.t0 > 25.0) c.fail("dave not joined with a team store after 25 s");
            return false;
        }
        c.eq("dave's store is bob's team's", S.dave.teamStoreTeam, S.bobTeam);
        c.check("dave has bob's team's research", Has(S.dave, "smelting"));
        c.check("...and not carol's", !Has(S.dave, "basic_tools"));
        c.check("dave's join brings bob's map pin (the team's pins)", HasPin(S.dave, "Copper_hill"), S.dave.teamStoreText.substr(0, 300));
        c.check("bob has his pin back from the host", HasPin(S.bob, "Copper_hill"));
        c.check("carol's team doesn't see it", !HasPin(S.carol, "Copper_hill"));
        const auto& bq = Research_Team(S.bobTeam).queue;
        c.check("the host queued bob's research for his team", std::find(bq.begin(), bq.end(), "steelmaking") != bq.end());
        c.check("...and dave's store shows it", S.dave.teamStoreText.find("queue steelmaking") != std::string::npos);
        c.check("the host refused carol's hidden node", Research_Team(S.carolTeam).queue.empty());
        c.lt("s from dave's join to his store (it comes with the join, not the 5 s resend)", c.t - S.daveJoinedAt, 1.0);
        S.carolStoresBefore = S.carol.teamStores;
        NetClient_SendChat(S.bob, "/team invite carol");
        S.invited = true; S.t0 = c.t;
        return true;
    });
    s.step("carol joins bob's team: the new team's store, without her old research", 15.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        UpdateAll();
        if (!S.accepted && c.t - S.t0 > 1.0) { NetClient_SendChat(S.carol, "/team accept"); S.accepted = true; S.t0 = c.t; }
        c.waitingFor(TextFormat("carol's store team %d (bob's %d)", S.carol.teamStoreTeam, S.bobTeam));
        if (!S.accepted || S.carol.teamStoreTeam != S.bobTeam) {
            if (S.accepted && c.t - S.t0 > 10.0) c.fail(TextFormat("carol's store still from team %d 10 s after /team accept (host: carol in %d)", S.carol.teamStoreTeam, Teams_Of(PlayerId_FromName("carol"))));
            return false;
        }
        c.lt("s from /team accept to carol's new store (sent with the team change)", c.t - S.t0, 2.0);
        c.eq("the host has carol in bob's team", Teams_Of(PlayerId_FromName("carol")), S.bobTeam);
        c.check("carol has her new team's research", Has(S.carol, "smelting"));
        c.check("carol's old team's research stayed behind", !Has(S.carol, "basic_tools"));
        c.check("carol sees her new team's map pin", HasPin(S.carol, "Copper_hill"));
        c.check("carol removes the pin", NetClient_SendMapPin(S.carol, false, "terra", kPinDir, "Copper_hill"));
        S.t0 = c.t;
        return true;
    });
    s.step("a pin removed by one member goes from the others' maps", 10.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        UpdateAll();
        c.waitingFor(TextFormat("bob has the pin %d, dave %d", (int)HasPin(S.bob, "Copper_hill"), (int)HasPin(S.dave, "Copper_hill")));
        if (HasPin(S.bob, "Copper_hill") || HasPin(S.dave, "Copper_hill")) {
            if (c.t - S.t0 > 5.0) c.fail("the removed pin is still in bob's or dave's store 5 s later");
            return false;
        }
        c.lt("s from carol's removal to bob's and dave's stores (sent at the next tick, not the 5 s resend)", c.t - S.t0, 2.0);
        c.eq("sends the transport refused", Net_SendFailures(), 0);
        for (NetClient* n : { &S.bob, &S.carol, &S.dave }) NetClient_Disconnect(*n);
        return true;
    });
}
