#ifndef SERVERADMIN_H
#define SERVERADMIN_H

// Server administration (#345): roles (owner, op, player), bans (by name, and by IP for direct joins), a whitelist mode,
// and an audit log of admin actions and joins/leaves. Identity is the player name (there are no Steam IDs yet), so a
// password (server.cfg) is what keeps a name from being borrowed.
//
// State: <dir>/admin.txt, one entry a line, rewritten on every change:
//   owner NAME | op NAME | ban NAME [reason] | banip IP [reason] | allow NAME | whitelist on
// Audit: <dir>/audit.log, time-stamped lines (who did what; joins, leaves and refusals).
//
// Console commands (the server console, and later ops in chat): op NAME, deop NAME, ops, ban NAME [reason], unban NAME,
// banip IP [reason], unbanip IP, bans, whitelist on|off|add NAME|remove NAME|list. Admin_Command answers them.

#include <string>

enum class AdminRole { Player, Op, Owner };

void Admin_Load(const std::string& dir);              // reads <dir>/admin.txt (none: empty lists)
// (#1086) A menu-hosted game's hosting player, by their id (UUID, sent only by their own client in Hello): always the
// owner, whatever name they or anyone else uses. "" (a dedicated server): none, admin.txt alone decides.
void Admin_SetHostOwner(const std::string& uuid);
// A player's role: the host owner by `uuid` (when given and it matches), else by name from admin.txt.
AdminRole Admin_RoleOf(const std::string& name, const std::string& uuid = std::string());
const char* Admin_RoleName(AdminRole r);
// May this player join? False with the reason shown to them (banned, not on the whitelist).
bool Admin_CheckJoin(const std::string& name, const std::string& ip, std::string& why);
// One admin command typed by `by` ("console" or an op's name). True when it was an admin command; `reply` says what
// happened (or why not: an op cannot change owners).
bool Admin_Command(const std::string& line, const std::string& by, std::string& reply, const std::string& byUuid = std::string());
void Admin_Audit(const std::string& text);            // a line in audit.log
bool Admin_Test();                                    // --test-server-admin

#endif
