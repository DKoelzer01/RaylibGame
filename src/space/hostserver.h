#ifndef HOSTSERVER_H
#define HOSTSERVER_H

// Dedicated-server start (#311, second slice): `--host <save> [--port N] [--seed N] [--autosave-minutes N]` loads the
// save, or when the file does not exist starts a new world (from --seed N, else a random seed that is printed) and saves
// it there; an existing save keeps its own seed. It runs headless (--headless: a hidden window for the GL context),
// autosaves every 30 minutes by default (--autosave-minutes, 0 = off; user 2026-11-15) and on shutdown, and listens.
// With ENABLE_GNS it accepts connections over GameNetworkingSockets; without it the world still runs, listening is
// reported as unavailable.
//
// The connection lifecycle (#321) speaks the binary protocol (netproto.h):
//   client -> Hello (protocol, build commit, name, password)
//   host   -> Welcome (session, the world description: seed, bodies, fingerprint; settings; the server tick) + Shower,
//             or Reject (protocol, build, password, full, name already playing, garbled) and closes
//   client -> Ping every second; host -> Pong (round trip, the server tick). Nothing heard for 15 s: Bye(timed out).
//   either -> Bye (quit; the server closed). A player who comes back under the same name gets the same session.
// Players (#327): each session has a server-side Player stepped from its client's inputs (one a sim step) and sent to
// everyone as PlayerState at 30 Hz (unreliable); the owner reconciles its prediction, the others interpolate.

#include "dvec3.h"
#include <cstdint>
#include <string>
#include <vector>

class Scene;
struct NetMessage;

constexpr uint16_t kDefaultHostPort = 27015;
constexpr double kDefaultAutosaveMinutes = 30.0;

bool Host_Start(uint16_t port, std::string& status); // gNetRole = Host; false (with why) when not listening
void Host_SetAutosaveMinutes(double minutes);         // 0 = off
void Host_SetPassword(const std::string& password);   // "" none (--host-password)
void Host_SetMaxPlayers(int n);                       // default 8 (--max-players)
void Host_Update(Scene& world);                       // per frame: autosave, poll, the handshake, pings, timeouts
bool Host_SaveOnShutdown(Scene& world);               // before the scene goes away
void Host_Stop();                                     // tells every client the server closed
bool Host_Listening();
int Host_Joins();                                     // joins accepted so far
int Host_Autosaves();                                 // autosaves written so far
int Host_OnlinePlayers();                             // sessions connected now
// Tests (#327): a session's player as the server simulates it (body NetId, body-local feet); false: no such player.
bool Host_PlayerLocal(Scene& world, uint32_t session, uint64_t& body, DVec3& local);
// Tests: client connections opened in the same process; their messages are handed to the client, not to the host.
void Host_SetLocalClient(uint32_t conn);
std::vector<NetMessage> Host_TakeLocalClientMessages(uint32_t conn);

#endif
