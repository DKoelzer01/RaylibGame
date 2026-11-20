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
//   client -> Ping every second; host -> Pong (round trip, the server tick). Nothing heard for 15 s: Bye(timed out)
//             (#1033: 120 s while a joined client is still loading, before its first input; Host_SilenceLimit).
//   either -> Bye (quit; the server closed). A player who comes back under the same name gets the same session.
// Players (#327): each session has a server-side Player stepped from its client's inputs (one a sim step) and sent to
// everyone as PlayerState at 30 Hz (unreliable); the owner reconciles its prediction, the others interpolate.

#include "dvec3.h"
#include <cstdint>
#include <functional>
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
enum class DeathDrops : uint8_t;
void Host_SetDeathDrops(DeathDrops d);                // (#804) --death-drops difficulty|keep|grave (default difficulty)
int Host_GravesMade();
int Host_DropsMade();                                 // (#865) drops made for joined players
int Host_CreativeGives();                             // (#866) creative gives done for joined players
int Host_CreativeRefused();                           // (#866) refused: not creative
void Host_Update(Scene& world);                       // per frame: autosave, poll, the handshake, pings, timeouts
bool Host_SaveOnShutdown(Scene& world);               // before the scene goes away
void Host_Stop();                                     // tells every client the server closed
bool Host_Listening();
int Host_Joins();                                     // joins accepted so far
int Host_Autosaves();                                 // autosaves written so far
int Host_OnlinePlayers();                             // sessions connected now
// Tests (#327): a session's player as the server simulates it (body NetId, body-local feet); false: no such player.
bool Host_PlayerLocal(Scene& world, uint32_t session, uint64_t& body, DVec3& local);
// Tests (#329): terrain edits applied for players so far; a player's inventory as the server keeps it.
int Host_EditsApplied();
double Host_InventoryAmount(uint32_t session, const std::string& item);
// Building (#330): parts placed or removed for players so far (tests); give a player items (the admin "give", tests):
// into the server's copy of their inventory, and the delta to them. False: no such player or item.
int Host_PartsApplied();
bool Host_Give(uint32_t session, const std::string& item, int count);
bool Host_SetCreative(const std::string& name, bool on); // (#823) /gamemode: an online player's game mode (false: not online)
// Tests (#346): how many of an item the server holds for a session (-1: no such player in the world); forget every
// player as a restart would (only with nobody online; their records stay on disk).
int Host_ItemCount(uint32_t session, const std::string& item);
double Host_SilenceLimit(bool joined, bool inGame); // (#1033) seconds of silence before a connection counts as gone
bool Host_TestSkillXp(uint32_t session, const std::string& skill, double xp); // (#920 tests)
bool Host_TestSetFood(uint32_t session, float food);                          // (#920 tests)
int Host_SkillSpends();                                                        // (#920 tests) allowed
int Host_SkillSpendsRefused();                                                 // (#920 tests) refused
bool Host_TestGive(uint32_t session, const std::string& item, int count); // (tests) into a remote player's inventory; count < 0 takes
double Host_PlayersApart(class Scene& world, uint32_t a, uint32_t b);          // (#865 tests) metres between two remote players
bool Host_TestPutPlayerAt(class Scene& world, uint32_t who, uint32_t at); // (#865 tests) one remote player onto another's spot
bool Host_ForgetPlayers();
// Join sync (#776, tests): syncs finished so far; bytes sent by them.
int Host_SyncsDone();
uint64_t Host_SyncBytesSent();
// The dedicated server's console and admin (#326, #345): who is (or was) here; kick by name (Bye kicked, with the reason as a chat
// line first); a chat line from the server (fromSession 0) or a player, to everyone (toSession 0) or one session; the message of
// the day (in the Welcome); a join check (bans, whitelist: false refuses with NetReject::Banned and why); a handler for
// players' chat (else it goes to everyone); joins, leaves and refusals for the audit log.
struct HostPlayer { uint32_t session = 0; std::string name; bool online = false; std::string address; };
std::vector<HostPlayer> Host_ListPlayers();
bool Host_Kick(const std::string& name, const std::string& reason);
void Host_Say(const std::string& text, uint32_t fromSession = 0, uint32_t toSession = 0);
void Host_SetMotd(const std::string& text);
void Host_SetJoinCheck(std::function<bool(const std::string& name, const std::string& address, std::string& why)> check);
void Host_SetChatHandler(std::function<void(uint32_t session, const std::string& name, const std::string& text)> handler);
void Host_SetEventHook(std::function<void(const char* what, uint32_t session, const std::string& detail)> hook);
// Tests (#349): messages the server's guard dropped (oversized, unexpected, malformed, over-rate, impossible), connections kicked.
int Host_GuardDropped();
int Host_GuardKicks();
// Tests: client connections opened in the same process; their messages are handed to the client, not to the host.
void Host_SetLocalClient(uint32_t conn);
std::vector<NetMessage> Host_TakeLocalClientMessages(uint32_t conn);

#endif
