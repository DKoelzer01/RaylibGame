#ifndef NETGUARD_H
#define NETGUARD_H

// Per-connection message limits on the server (#349): what a client may send, how big, how often. Each message kind a
// client may send has a token bucket (a steady rate and a burst); others (the server's own kinds) aren't accepted from a
// client at all. Oversized, unexpected, malformed or over-rate messages are dropped and count as violations; violations
// decay (one a second), and past kKickAt the client is kicked ("flooding"). A real client never comes near the limits:
// inputs go out once a sim step (60 a second, in packets carrying the last 8), pings once a second.
// Standalone and testable: --test-net-guard runs a normal client and abusive ones through it without a network.

#include <cstddef>
#include <cstdint>
#include <string>
#include "netproto.h"

struct NetGuard {
    double tokens[64]; double tokenTime[64]; // per kind: tokens, and when they were last topped up
    double last = -1.0;
    double violations = 0.0;
    uint64_t dropped = 0;
    NetGuard();
};

enum class GuardVerdict { Accept, Drop, Kick };
constexpr size_t kGuardMaxBytes = 4096; // no client message is near this (inputs ~200 bytes, chat 256 + header)

// One message from a client: `kind` its id, `bytes` its size, `wellFormed` whether its reader accepted it. `why` says what
// was wrong (for the log) when it isn't Accept.
GuardVerdict NetGuard_Admit(NetGuard& g, NetMsg kind, size_t bytes, bool wellFormed, double now, std::string& why);
// Something wrong found later (an impossible input, ...): a violation, Kick when past the limit.
GuardVerdict NetGuard_Violation(NetGuard& g, double now);
bool NetGuard_Test();

#endif
