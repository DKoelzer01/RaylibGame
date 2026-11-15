#ifndef HOSTSERVER_H
#define HOSTSERVER_H

// Dedicated-server start (#311, second slice): `--host <save> [--port N]` loads the save (a missing file starts a new
// world from the seed and saves it there), runs it headless (--headless: a hidden window for the GL context) and listens.
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

bool Host_Start(uint16_t port, std::string& status); // gNetRole = Host; false (with why) when not listening
void Host_Update(Scene& world);                       // per frame: poll, answer joins
void Host_Stop();
bool Host_Listening();
int Host_Joins();                                     // join requests answered so far
// Tests: a client connection opened in the same process; its messages are handed to the test, not handled as the host.
void Host_SetLocalClient(uint32_t conn);
std::vector<NetMessage> Host_TakeLocalClientMessages();

#endif
