#include "netguard.h"

#include <algorithm>
#include <cstdio>

// See netguard.h (#349).
namespace {
struct Limit { double rate, burst; }; // messages a second, and how many may come at once
constexpr double kKickAt = 40.0;      // violations (decaying one a second)
// What a client may send. Inputs: 60 packets a second at most (one a sim step), double for catch-up after a hitch.
Limit LimitOf(NetMsg k) {
    switch (k) {
        case NetMsg::Hello: return { 0.2, 2 };
        case NetMsg::Ping: return { 4, 8 };
        case NetMsg::Input: return { 130, 240 };
        case NetMsg::Chat: return { 1.5, 6 };
        case NetMsg::Loadout: return { 10, 20 };
        case NetMsg::EditRequest: return { 30, 60 };  // (#329 has its own stroke budget behind this)
        case NetMsg::PartRequest: return { 30, 60 };  // (#330 likewise)
        case NetMsg::SyncWant: return { 0.2, 3 };     // (#776) once a join (its bitmask: 4096 bytes is 32k entries)
        case NetMsg::AnimalHit: return { 10, 20 };    // (#334; a swing cooldown is checked behind this)
        case NetMsg::UseItem: return { 4, 8 };        // (#337: eating, a canister, a cell)
        case NetMsg::BedUse: return { 2, 4 };         // (#800: a right-click on a bed)
        case NetMsg::CraftRequest: return { 4, 10 };  // (#338: a click, or a few)
        case NetMsg::MachineAction: return { 8, 20 };  // (#332: clicks in a machine screen, I / O, watch on open and close)
        case NetMsg::PlantAction: return { 10, 20 };   // (#335: a fell (bushes and stones fall in one swing), I / O on farmland; the fell rate is checked behind this)
        case NetMsg::DropItems: return { 6, 12 };     // (#865: Q drops one per press, a click outside the window a stack)
        case NetMsg::CreativeGive: return { 8, 20 };  // (#866: Catalog clicks and middle-clicks, a drag across the Catalog)
        case NetMsg::Bye: return { 1, 2 };
        default: return { 0, 0 };                    // a server-to-client kind, or unknown: never from a client
    }
}
}

NetGuard::NetGuard() { for (int i = 0; i < 64; ++i) { tokens[i] = -1.0; tokenTime[i] = 0.0; } }

GuardVerdict NetGuard_Violation(NetGuard& g, double now) {
    if (g.last >= 0.0) g.violations = std::max(0.0, g.violations - (now - g.last)); // (decay, also when no message)
    g.last = std::max(g.last, now);
    g.violations += 1.0; ++g.dropped;
    return g.violations > kKickAt ? GuardVerdict::Kick : GuardVerdict::Drop;
}

GuardVerdict NetGuard_Admit(NetGuard& g, NetMsg kind, size_t bytes, bool wellFormed, double now, std::string& why) {
    const double dt = g.last < 0.0 ? 0.0 : std::max(0.0, now - g.last);
    g.violations = std::max(0.0, g.violations - dt);
    g.last = now;
    const int k = (int)kind;
    const Limit lim = (k >= 0 && k < 64) ? LimitOf(kind) : Limit{ 0, 0 };
    auto bad = [&](const char* w) { why = w; g.violations += 1.0; ++g.dropped; return g.violations > kKickAt ? GuardVerdict::Kick : GuardVerdict::Drop; };
    if (bytes > kGuardMaxBytes) return bad("oversized message");
    if (lim.burst <= 0.0) return bad("a message kind clients don't send");
    if (!wellFormed) return bad("malformed message");
    double& t = g.tokens[k];
    t = t < 0.0 ? lim.burst : std::min(lim.burst, t + std::max(0.0, now - g.tokenTime[k]) * lim.rate); // (since this kind's last top-up)
    g.tokenTime[k] = now;
    if (t < 1.0) return bad("too many messages");
    t -= 1.0;
    return GuardVerdict::Accept;
}

// --test-net-guard: a normal client never trips it; floods, oversized, unexpected and garbage messages are dropped and,
// kept up, end in a kick; a short burst (a hitch) passes; violations wear off.
bool NetGuard_Test() {
    bool ok = true;
    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: net guard: %s FAILED\n", what); ok = false; } };
    std::string why;
    { // 60 s of a normal client: 60 input packets, 1 ping a second, a chat line every 5 s, a loadout change every 10 s
        NetGuard g; int dropped = 0;
        for (int f = 0; f < 3600; ++f) {
            const double t = f / 60.0;
            dropped += NetGuard_Admit(g, NetMsg::Input, 180, true, t, why) != GuardVerdict::Accept;
            if (f % 60 == 0) dropped += NetGuard_Admit(g, NetMsg::Ping, 20, true, t, why) != GuardVerdict::Accept;
            if (f % 300 == 0) dropped += NetGuard_Admit(g, NetMsg::Chat, 60, true, t, why) != GuardVerdict::Accept;
            if (f % 600 == 0) dropped += NetGuard_Admit(g, NetMsg::Loadout, 30, true, t, why) != GuardVerdict::Accept;
        }
        check(dropped == 0, "a normal client is never dropped");
    }
    { // a 0.5 s hitch: the client's 30 queued input packets arrive at once
        NetGuard g; int dropped = 0;
        for (int f = 0; f < 60; ++f) NetGuard_Admit(g, NetMsg::Input, 180, true, f / 60.0, why);
        for (int i = 0; i < 30; ++i) dropped += NetGuard_Admit(g, NetMsg::Input, 180, true, 1.5, why) != GuardVerdict::Accept;
        check(dropped == 0, "a hitch's burst of inputs passes");
    }
    { // chat flood: 20 lines a second
        NetGuard g; int dropped = 0; bool kicked = false; double kickAt = -1.0;
        for (int i = 0; i < 400 && !kicked; ++i) {
            const GuardVerdict v = NetGuard_Admit(g, NetMsg::Chat, 40, true, i / 20.0, why);
            dropped += v != GuardVerdict::Accept; if (v == GuardVerdict::Kick) { kicked = true; kickAt = i / 20.0; }
        }
        check(dropped > 0 && kicked && kickAt > 1.0, "a chat flood is dropped, then kicked (not at once)");
        printf("bench: net guard: a 20/s chat flood: kicked after %.1f s\n", kickAt);
    }
    { // garbage: oversized, a server-only kind, malformed
        NetGuard g;
        check(NetGuard_Admit(g, NetMsg::Input, 10000, true, 0.0, why) == GuardVerdict::Drop && why == "oversized message", "oversized dropped");
        check(NetGuard_Admit(g, NetMsg::PlayerState, 40, true, 0.1, why) == GuardVerdict::Drop, "a server-only kind from a client dropped");
        check(NetGuard_Admit(g, NetMsg::Input, 40, false, 0.2, why) == GuardVerdict::Drop && why == "malformed message", "malformed dropped");
        bool kicked = false;
        for (int i = 0; i < 200 && !kicked; ++i) kicked = NetGuard_Admit(g, NetMsg::Count, 8, false, 0.3 + i * 0.01, why) == GuardVerdict::Kick;
        check(kicked, "a stream of garbage is kicked");
    }
    { // violations wear off: 30 bad messages, a minute of quiet, 30 more: no kick
        NetGuard g; bool kicked = false;
        for (int i = 0; i < 30; ++i) kicked |= NetGuard_Admit(g, NetMsg::PlayerState, 40, true, i * 0.01, why) == GuardVerdict::Kick;
        for (int i = 0; i < 30; ++i) kicked |= NetGuard_Admit(g, NetMsg::PlayerState, 40, true, 60.0 + i * 0.01, why) == GuardVerdict::Kick;
        check(!kicked, "violations decay");
    }
    printf("bench: net guard (#349): normal play, a hitch, a chat flood, garbage, decay -> %s\n", ok ? "PASS" : "FAIL");
    return ok;
}
