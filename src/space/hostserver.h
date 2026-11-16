#ifndef HOSTSERVER_H
#define HOSTSERVER_H

// Dedicated-server start (#311, second slice): `--host <save> [--port N] [--seed N] [--autosave-minutes N]` loads the
// save, or when the file does not exist starts a new world (from --seed N, else a random seed that is printed) and saves
// it there; an existing save keeps its own seed. It runs headless (--headless: a hidden window for the GL context),
// autosaves every 30 minutes by default (--autosave-minutes, 0 = off; user 2026-11-15) and on shutdown, and listens.
// With ENABLE_GNS it accepts connections over GameNetworkingSockets; without it the world still runs, listening is
// reported as unavailable. This slice answers a join request only:
//   client -> host  "join"                     host -> client  "J" + JoinInfo_Write(...) (seed, bodies, fingerprint)
//                                                              "S <start> <end> <id>"     (the meteor shower, #316)
// Replication of edits, players and entities comes later (#304).

#include <cstdint>
#include <string>
#include <vector>

class Scene;
struct NetMessage;

constexpr uint16_t kDefaultHostPort = 27015;
constexpr double kDefaultAutosaveMinutes = 30.0;

bool Host_Start(uint16_t port, std::string& status); // gNetRole = Host; false (with why) when not listening
void Host_SetAutosaveMinutes(double minutes);         // 0 = off
void Host_Update(Scene& world);                       // per frame: autosave, poll, answer joins
bool Host_SaveOnShutdown(Scene& world);               // before the scene goes away
void Host_Stop();
bool Host_Listening();
int Host_Joins();                                     // join requests answered so far
int Host_Autosaves();                                 // autosaves written so far
// Tests: a client connection opened in the same process; its messages are handed to the test, not handled as the host.
void Host_SetLocalClient(uint32_t conn);
std::vector<NetMessage> Host_TakeLocalClientMessages();

#endif
