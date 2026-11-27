// net/team-store (#1886): the host sends each player their team's progression store (research, explored map, blueprint
// library) and keeps two teams apart. The game hosts (--host; needs a GameNetworkingSockets build, else one SKIP check
// passes); bob and carol join, each in a team of one. The host's teams then research different nodes: each client hears
// only its own team's within the resend interval (5 s). dave joins late into bob's team (RAYLIBGAME_HOST_TEAM): he gets
// bob's research with his join. carol joins bob's team by the real commands (/team invite, /team accept): she gets the
// new team's store at once, and her old team's research doesn't come with her (the user, on #1886).
// The clients run in this process, so they don't apply the text (the host's state is the real one here): the checks read
// what each client received.
#include "scenario.h"
#include "hostserver.h"
#include "netclient.h"
#include "net.h"
#include "playerrecord.h"
#include "research.h"
#include "world.h"
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
void UpdateAll() { RealTime(); for (NetClient* n : { &S.bob, &S.carol, &S.dave }) if (n->state != NetClient::State::Idle) NetClient_Update(*n); }
}

SCENARIO("net/team-store") {
    s.tags = { "multiplayer", "net" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--host", "saves/scenario-net-team-store.sav", "--port", std::to_string(kPort) };
    s.env = { { "RAYLIBGAME_HOST_TEAM", "bob,dave" } };
    s.timeout = 120.0;
    s.step("bob and carol join, each in a team of one", 40.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        if (!Host_Listening()) { // (the harness has no SKIP verdict: one passing check, said plainly)
            c.check("SKIP: not listening (a build without GameNetworkingSockets)", true);
            S.skip = true; return true;
        }
        if (S.bob.state == NetClient::State::Idle) {
            const std::string addr = TextFormat("127.0.0.1:%d", kPort);
            NetClient_Connect(S.bob, addr.c_str(), "bob");
            NetClient_Connect(S.carol, addr.c_str(), "carol");
            S.t0 = c.t;
        }
        UpdateAll();
        c.waitingFor(TextFormat("bob %s, carol %s", NetClient_StateName(S.bob.state), NetClient_StateName(S.carol.state)));
        if (!(Joined(S.bob) && Joined(S.carol) && S.bob.teamStores > 0 && S.carol.teamStores > 0)) {
            if (c.t - S.t0 > 30.0) c.fail("bob and carol not both joined with a team store after 30 s");
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
        NetClient_Connect(S.dave, TextFormat("127.0.0.1:%d", kPort), "dave");
        S.t0 = c.t;
        return true;
    });
    s.step("dave joins bob's team late and gets its research with his join", 30.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        UpdateAll();
        if (Joined(S.dave) && S.daveJoinedAt < 0.0) S.daveJoinedAt = c.t;
        c.waitingFor(TextFormat("dave %s, %d stores", NetClient_StateName(S.dave.state), S.dave.teamStores));
        if (!(Joined(S.dave) && S.dave.teamStores > 0)) {
            if (c.t - S.t0 > 25.0) c.fail("dave not joined with a team store after 25 s");
            return false;
        }
        c.eq("dave's store is bob's team's", S.dave.teamStoreTeam, S.bobTeam);
        c.check("dave has bob's team's research", Has(S.dave, "smelting"));
        c.check("...and not carol's", !Has(S.dave, "basic_tools"));
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
        c.eq("sends the transport refused", Net_SendFailures(), 0);
        for (NetClient* n : { &S.bob, &S.carol, &S.dave }) NetClient_Disconnect(*n);
        return true;
    });
}
