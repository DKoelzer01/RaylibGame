#ifndef SERVERADMIN_H
#define SERVERADMIN_H

// Server administration (#345): roles (owner, op, player), bans (by name, and by IP for direct joins), a whitelist mode,
// and an audit log of admin actions and joins/leaves.
//
// (#1868) Owners and ops are a player's id (the UUID in their Hello and their player record), not a name. `op NAME` looks
// the name up among the players online and the player records and stores the id it finds; a name no one has played
// under is refused (the console may name the id: `op NAME UUID`). An entry without an id (an admin.txt from before)
// grants nothing until the name is op'd again; `ops` lists those as "unbound". Bans keep the name and, when known, the id.
// The whitelist lets a listed name in while it has no id (a friend's first join); with an id, only that player.
//
// State: <dir>/admin.txt, one entry a line, rewritten on every change:
//   owner NAME [id:UUID] | op NAME [id:UUID] | ban NAME [id:UUID] [reason] | banip IP [reason] | allow NAME [id:UUID] |
//   whitelist on
// Audit: <dir>/audit.log, time-stamped lines (who did what; joins, leaves and refusals).
//
// Console commands (the server console, and ops in chat): op NAME [UUID], deop NAME, ops, ban NAME [reason], unban NAME,
// banip IP [reason], unbanip IP, bans, whitelist on|off|add NAME|remove NAME|list, owner NAME [UUID] (console only).
// Admin_Command answers them.

#include <functional>
#include <string>

enum class AdminRole { Player, Op, Owner };

void Admin_Load(const std::string& dir);              // reads <dir>/admin.txt (none: empty lists)
// (#1086) A menu-hosted game's hosting player, by their id (UUID, sent only by their own client in Hello): always the
// owner, whatever name they or anyone else uses. "" (a dedicated server): none, admin.txt alone decides.
void Admin_SetHostOwner(const std::string& uuid);
// (#1868) How a name becomes a player's id: the server's players online, then its player records ("" no one by that name).
void Admin_SetResolver(std::function<std::string(const std::string& name)> resolve);
// A player's role, by their id: the host owner, else admin.txt's owners and ops. No id: a player. (The name is not used.)
AdminRole Admin_RoleOf(const std::string& name, const std::string& uuid = std::string());
const char* Admin_RoleName(AdminRole r);
// May this player join? False with the reason shown to them (banned by name, id or address; not on the whitelist).
bool Admin_CheckJoin(const std::string& name, const std::string& ip, std::string& why, const std::string& uuid = std::string());
// One admin command typed by `by` ("console" or a player's name, with their id `byUuid`). True when it was an admin
// command; `reply` says what happened (or why not: an op cannot change owners).
bool Admin_Command(const std::string& line, const std::string& by, std::string& reply, const std::string& byUuid = std::string());
void Admin_Audit(const std::string& text);            // a line in audit.log
bool Admin_Test();                                    // --test-server-admin

#endif
