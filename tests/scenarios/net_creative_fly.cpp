// net/creative-fly (#1922): a creative client flies and the server's copy of it follows. Player::update's fly path moves
// world.camera (its own player's view), so a server stepping a flying remote flew the server's camera instead and the
// remote stayed behind; the client's states from the server then pulled it back (rubber-banding, worse on a slow link).
// Now each flying input step carries the client's pose and the server takes it for creative players.
//
// The game hosts (--host, needs a GameNetworkingSockets build; without one the scenario notes SKIP and passes) and one
// in-process client, "flyer", joins, is made creative and flies 75 m along the ground 30 m up. The client's fly camera is
// world.camera only while it steps; between its steps the camera is parked 60 m off, looking at the flight, so a server
// that flew the remote with world.camera would put it at the parked eye. Checked: (1) the client is barely corrected
// during the flight, (2) the server's copy is where the client is (and moved), (3) the server reports it flying.
#include "scenario.h"
#include "hostserver.h"
#include "netclient.h"
#include "commandline.h"
#include "player.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <memory>
#include <string>
#include <thread>

namespace {
constexpr uint16_t kPort = 27953;
// A frame at most every 1/60 s of wall time: scenarios step a fixed 1/60 s a frame as fast as the machine goes (thousands
// a second headless), and a client sending an input a frame that fast is a flood the server kicks. A game is real time.
void RealTime() {
    static double next = 0.0;
    const double now = GetTime();
    if (next > now) std::this_thread::sleep_for(std::chrono::duration<double>(next - now));
    next = std::max(now, next) + 1.0 / 60.0;
}
constexpr double kFlySeconds = 5.0, kFlySpeed = 15.0, kHeight = 30.0; // s, m/s, m above the ground
struct State {
    NetClient cl;
    std::unique_ptr<Player> p;
    long long stepped = 0;
    double t0 = 0.0, tFly = -1.0;
    Vector3 up{}, east{};                 // body-local
    DVec3 start{};                        // body-local: the flight's first point (feet)
    int correctionsAtStart = 0;
    double worstError = 0.0;              // m: the largest prediction error the client saw while flying
    bool skip = false;
} S;

// The flight's feet at time t (body-local): along `east`, kHeight over the start; then held.
DVec3 FlightAt(double t) {
    const double d = kFlySpeed * std::clamp(t, 0.0, kFlySeconds);
    return S.start + DVec3(S.east.x * d, S.east.y * d, S.east.z * d);
}
// The host's own camera, parked 60 m to the side of the flight's middle, looking at it.
void Park(ScenarioCtx& c) {
    const DVec3 mid = FlightAt(kFlySeconds * 0.5);
    const Vector3 side = Vector3Normalize(Vector3CrossProduct(S.up, S.east));
    const Vector3 m = mid.toFloat();
    c.camera(Vector3Add(m, Vector3Add(Vector3Scale(side, 60.0f), Vector3Scale(S.up, 10.0f))), m, S.up);
}
}

static void BuildFly(ScenarioDef& s) {
    s.tags = { "multiplayer", "net" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--host", "saves/scenario-net-creative-fly.sav", "--port", std::to_string(kPort) };
    s.timeout = 120.0;
    s.realTime = true; // (#1972) the join is over real sockets
    s.step("join as flyer", 40.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        if (!Host_Listening()) { // (the harness has no SKIP verdict: one passing check, said plainly)
            c.check("SKIP: not listening (a build without GameNetworkingSockets)", true);
            S.skip = true; return true;
        }
        if (S.cl.state == NetClient::State::Idle) {
            NetClient_Connect(S.cl, TextFormat("127.0.0.1:%d", kPort), "flyer");
            S.t0 = c.clock;
        }
        RealTime(); NetClient_Update(S.cl);
        c.waitingFor("the client to join and get its state");
        if (S.cl.state != NetClient::State::Joined || !S.cl.haveOwn) {
            if (c.clock - S.t0 > 30.0) c.fail(TextFormat("not joined after 30 s (state %d)", (int)S.cl.state));
            return false;
        }
        Planetoid* b = Planetoid::findNet(S.cl.own.body);
        if (!b) { c.fail("the client's body isn't here"); return true; }
        S.p = std::make_unique<Player>();
        S.p->spawnOn(c.world, *b, b->dirToWorld(Vector3Normalize(S.cl.own.local.toFloat())));
        S.p->applyNetState(c.world, NetClient_FromMsg(S.cl.own));
        c.check("the server makes flyer creative (/gamemode)", Host_SetCreative("flyer", true));
        GameMode_SetCreative(true); // (the process's mode: the client's code reads it; the server's message says the same)
        S.up = Vector3Normalize(S.cl.own.local.toFloat());
        const Vector3 pole = fabsf(S.up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 };
        S.east = Vector3Normalize(Vector3CrossProduct(pole, S.up));
        S.start = S.cl.own.local + DVec3(S.up.x * kHeight, S.up.y * kHeight, S.up.z * kHeight);
        S.stepped = SimTick();
        return true;
    });
    s.step("fly 75 m, then hover", 30.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        RealTime(); NetClient_Update(S.cl);
        if (S.cl.state != NetClient::State::Joined) { c.fail("the client lost the server"); return true; }
        Planetoid* b = Planetoid::findNet(S.cl.own.body);
        if (!b) { c.fail("the client's body is gone"); return true; }
        if (c.t < 0.5) { Park(c); return false; } // (the creative flag reaches the client)
        if (S.tFly < 0.0) {
            S.tFly = c.t; S.correctionsAtStart = S.cl.corrections;
            S.p->flyAllowed = true; S.p->flyMode = true; S.p->flyFollow = false;
        }
        const double ft = c.t - S.tFly;
        // The client's steps, each with its fly camera where the flight is at that step (as its own game would have it).
        for (; S.stepped < SimTick(); ++S.stepped) {
            const double st = ft - (double)(SimTick() - S.stepped - 1) * SimDt();
            const DVec3 feet = FlightAt(st);
            const DVec3 eye = feet + DVec3(S.up.x * S.p->eyeHeight, S.up.y * S.p->eyeHeight, S.up.z * S.p->eyeHeight);
            c.world.camera.position = (b->toWorldD(eye) - c.world.renderOrigin).toFloat();
            c.world.camera.target = Vector3Add(c.world.camera.position, b->dirToWorld(S.east));
            c.world.camera.up = b->dirToWorld(S.up);
            PlayerInput in; in.enabled = true;
            NetClient_Step(S.cl, c.world, *S.p, in);
            if (ft > 0.3) S.worstError = std::max(S.worstError, S.cl.lastError); // (past the first states, sent before it flew)
        }
        Park(c); // the host's own view between the client's steps: the server must not fly the remote with it
        c.waitingFor("the flight and a second of hover");
        return ft > kFlySeconds + 1.5;
    });
    s.step("the server's flyer is where the client is", 10.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        Planetoid* b = Planetoid::findNet(S.cl.own.body);
        uint64_t sb = 0; DVec3 sl{};
        const bool have = b && Host_PlayerLocal(c.world, S.cl.welcome.session, sb, sl);
        c.check("the server has the flyer", have && sb == b->netId);
        if (!have) return true;
        // The client's last stepped position, body-local (not its world position now: no step since, and the planet's spin
        // has carried the ground ~1.5 m a frame under an unstepped flyer).
        const DVec3 cl = S.cl.history.empty() ? b->toLocalD(S.p->position) : S.cl.history.back().local;
        const double off = (sl - cl).length() * b->scale, moved = (sl - S.start).length() * b->scale;
        c.note(TextFormat("server vs client %.3f m; the server's flyer moved %.1f m from the flight's start; %d corrections in the flight, "
                          "worst prediction error %.3f m; %d own states", off, moved, S.cl.corrections - S.correctionsAtStart, S.worstError, S.cl.states));
        c.lt("the server's flyer is at the client's (m; before #1922: at the host's parked camera, ~60 m off)", off, 0.05);
        c.gt("the server's flyer flew (m from the flight's start; 75 m flown)", moved, 70.0);
        c.lt("the client's prediction held (m worst error while flying)", S.worstError, 0.25);
        c.lt("corrections while flying (before #1922: one per state)", S.cl.corrections - S.correctionsAtStart, 2.5);
        c.check("the server says the flyer flies (state flag 2)", (S.cl.own.flags & 2) != 0, TextFormat("flags %u", (unsigned)S.cl.own.flags));
        NetClient_Disconnect(S.cl);
        return true;
    });
}

SCENARIO("net/creative-fly") { BuildFly(s); }
// The same on a slow link (#1922 was found on one): 100 ms each way, 15 ms jitter, 2 % loss (net.cpp's simulated
// conditions, #347). Lost input packets are covered by the next ones (each carries the newest 4 steps).
SCENARIO("net/creative-fly-lag") {
    BuildFly(s);
    s.env = { { "RAYLIBGAME_NET_LAG_MS", "100" }, { "RAYLIBGAME_NET_JITTER_MS", "15" }, { "RAYLIBGAME_NET_LOSS", "2" }, { "RAYLIBGAME_NET_SEED", "1922" } };
}
