#include "serveradmin.h"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <map>
#include <set>
#include <sstream>

// See serveradmin.h (#345, #1868).
namespace {
std::string dir;
std::map<std::string, std::string> owners, ops, allowed; // lower-cased name -> the player's id ("": none yet)
struct Ban { std::string uuid, reason; };
std::map<std::string, Ban> bans;                         // lower-cased name -> its id (when known) and the reason
std::map<std::string, std::string> ipBans;               // ip -> reason
bool whitelist = false;
std::function<std::string(const std::string&)> resolver; // (#1868) name -> id

std::string Lower(std::string s) { for (char& c : s) c = (char)std::tolower((unsigned char)c); return s; }
std::string Path(const char* file) { return dir.empty() ? std::string(file) : dir + "/" + file; }
std::string Resolve(const std::string& name) { return resolver ? resolver(name) : std::string(); }
std::string IdField(const std::string& uuid) { return uuid.empty() ? std::string() : " id:" + uuid; }
bool HasId(const std::map<std::string, std::string>& m, const std::string& uuid) {
    if (uuid.empty()) return false;
    for (const auto& [n, u] : m) if (u == uuid) return true;
    return false;
}

void Save() {
    std::error_code ec;
    if (!dir.empty()) std::filesystem::create_directories(dir, ec);
    FILE* f = fopen(Path("admin.txt").c_str(), "w");
    if (!f) return;
    fputs("# RaylibGame server admin state (#345); rewritten by the server on every change. Owners and ops are by id (#1868).\n", f);
    for (const auto& [n, u] : owners) fprintf(f, "owner %s%s\n", n.c_str(), IdField(u).c_str());
    for (const auto& [n, u] : ops) fprintf(f, "op %s%s\n", n.c_str(), IdField(u).c_str());
    for (const auto& [n, b] : bans) fprintf(f, "ban %s%s%s%s\n", n.c_str(), IdField(b.uuid).c_str(), b.reason.empty() ? "" : " ", b.reason.c_str());
    for (const auto& [ip, r] : ipBans) fprintf(f, "banip %s%s%s\n", ip.c_str(), r.empty() ? "" : " ", r.c_str());
    for (const auto& [n, u] : allowed) fprintf(f, "allow %s%s\n", n.c_str(), IdField(u).c_str());
    if (whitelist) fputs("whitelist on\n", f);
    fclose(f);
}

// A role or whitelist entry for a name: the id given (the console), else the one the name resolves to. False: no player
// has played under that name and none was given.
bool Bind(const std::string& name, const std::string& given, std::string& uuid) {
    uuid = given.empty() ? Resolve(name) : given;
    return !uuid.empty();
}
}

void Admin_Load(const std::string& d) {
    dir = d;
    owners.clear(); ops.clear(); allowed.clear(); bans.clear(); ipBans.clear(); whitelist = false;
    FILE* f = fopen(Path("admin.txt").c_str(), "r");
    if (!f) return;
    char buf[512];
    int unbound = 0;
    while (fgets(buf, sizeof buf, f)) {
        std::istringstream in(buf);
        std::string kind, who, reason, uuid;
        if (!(in >> kind) || kind[0] == '#') continue;
        in >> who;
        std::getline(in, reason);
        reason.erase(0, reason.find_first_not_of(" \t")); while (!reason.empty() && (reason.back() == '\n' || reason.back() == '\r')) reason.pop_back();
        if (reason.rfind("id:", 0) == 0) { // (#1868) the player's id
            const size_t sp = reason.find(' ');
            uuid = reason.substr(3, sp == std::string::npos ? std::string::npos : sp - 3);
            reason = sp == std::string::npos ? std::string() : reason.substr(sp + 1);
        }
        const std::string n = Lower(who);
        if (kind == "owner") { owners[n] = uuid; unbound += uuid.empty(); }
        else if (kind == "op") { ops[n] = uuid; unbound += uuid.empty(); }
        else if (kind == "ban") bans[n] = { uuid, reason };
        else if (kind == "banip") ipBans[who] = reason;
        else if (kind == "allow") allowed[n] = uuid;
        else if (kind == "whitelist") whitelist = who == "on";
    }
    fclose(f);
    if (unbound) printf("admin: %d owner/op entr%s in admin.txt name%s no player id: %s grant nothing until op'd again (#1868; `ops` lists them)\n",
                        unbound, unbound == 1 ? "y" : "ies", unbound == 1 ? "s" : "", unbound == 1 ? "it" : "they");
}

static std::string hostOwnerId; // (#1086)
void Admin_SetHostOwner(const std::string& uuid) { hostOwnerId = uuid; }
void Admin_SetResolver(std::function<std::string(const std::string&)> resolve) { resolver = std::move(resolve); }
AdminRole Admin_RoleOf(const std::string&, const std::string& uuid) {
    if (uuid.empty()) return AdminRole::Player;
    if (uuid == hostOwnerId) return AdminRole::Owner; // the player hosting this game from their menu
    return HasId(owners, uuid) ? AdminRole::Owner : HasId(ops, uuid) ? AdminRole::Op : AdminRole::Player;
}
const char* Admin_RoleName(AdminRole r) { return r == AdminRole::Owner ? "owner" : r == AdminRole::Op ? "op" : "player"; }

bool Admin_CheckJoin(const std::string& name, const std::string& ip, std::string& why, const std::string& uuid) {
    const std::string n = Lower(name);
    auto banned = [&](const std::string& reason) { why = "banned" + (reason.empty() ? std::string() : ": " + reason); return false; };
    if (auto it = bans.find(n); it != bans.end()) return banned(it->second.reason);
    if (!uuid.empty()) for (const auto& [bn, b] : bans) if (b.uuid == uuid) return banned(b.reason); // (#1868) a new name doesn't lift it
    const std::string host = ip.substr(0, ip.rfind(':') == std::string::npos ? ip.size() : ip.rfind(':')); // "1.2.3.4:5678" -> the address
    if (auto it = ipBans.find(host); !host.empty() && it != ipBans.end()) return banned(it->second);
    if (whitelist) {
        const auto it = allowed.find(n);
        const bool listed = HasId(allowed, uuid) || (it != allowed.end() && it->second.empty()); // (an entry with an id: only that player)
        if (!listed && Admin_RoleOf(name, uuid) == AdminRole::Player) { why = "not on this server's whitelist"; return false; }
    }
    return true;
}

void Admin_Audit(const std::string& text) {
    std::error_code ec;
    if (!dir.empty()) std::filesystem::create_directories(dir, ec);
    FILE* f = fopen(Path("audit.log").c_str(), "a");
    if (!f) return;
    const std::time_t t = std::time(nullptr); char ts[32];
    std::strftime(ts, sizeof ts, "%Y-%m-%d %H:%M:%S", std::localtime(&t));
    fprintf(f, "[%s] %s\n", ts, text.c_str());
    fclose(f);
}

bool Admin_Command(const std::string& line, const std::string& by, std::string& reply, const std::string& byUuid) {
    std::istringstream in(line);
    std::string cmd, arg, rest;
    in >> cmd >> arg; std::getline(in, rest); rest.erase(0, rest.find_first_not_of(" \t"));
    const std::string a = Lower(arg);
    const bool console = by == "console";
    const AdminRole byRole = console ? AdminRole::Owner : Admin_RoleOf(by, byUuid);
    auto done = [&](const std::string& what) { reply = what; Admin_Audit(by + ": " + line + " -> " + what); Save(); return true; };
    auto need = [&](const char* usage) { reply = std::string("usage: ") + usage; return true; };
    auto list = [](const std::map<std::string, std::string>& m) { std::string s; for (auto& [n, u] : m) s += " " + n + (u.empty() ? " (unbound: op again)" : ""); return s; };
    // (#1868) an id named by the console (`op NAME UUID`); an op typing one is ignored (they can only op players who have joined)
    const std::string givenId = console ? rest.substr(0, rest.find(' ')) : std::string();
    const std::string noOne = a + " hasn't played on this server: roles are tied to a player's id, so op them once they've joined"
                                  " (or from the console: op NAME UUID)";
    if (cmd == "ops") {
        reply = "owners:" + list(owners) + "; ops:" + list(ops);
        return true;
    }
    if (cmd == "bans") {
        std::string s = "banned names:"; for (auto& [n, b] : bans) s += " " + n;
        s += "; banned addresses:"; for (auto& [ip, r] : ipBans) s += " " + ip;
        reply = s; return true;
    }
    if (cmd == "op" || cmd == "deop" || cmd == "ban" || cmd == "unban") {
        if (arg.empty()) return need((cmd + " NAME" + (cmd == "op" && console ? " [UUID]" : cmd == "ban" ? " [reason]" : "")).c_str());
        if (!console && byRole == AdminRole::Player) { reply = "only ops can do that"; return true; }
        if (!console && owners.count(a)) { reply = "an owner can only be changed from the server console"; return true; }
        if (cmd == "op") {
            std::string id;
            if (!Bind(a, givenId, id)) { reply = noOne; return true; }
            if (!console && HasId(owners, id)) { reply = "an owner can only be changed from the server console"; return true; }
            ops[a] = id;
            return done(a + " is now an op");
        }
        if (cmd == "deop") { const bool was = ops.erase(a) > 0; return done(was ? a + " is no longer an op" : a + " was not an op"); }
        if (cmd == "ban") {
            const std::string id = Resolve(a);
            if (!console && HasId(owners, id)) { reply = "an owner can only be changed from the server console"; return true; }
            bans[a] = { id, rest };
            ops.erase(a);
            for (auto it = ops.begin(); it != ops.end();) it = !id.empty() && it->second == id ? ops.erase(it) : std::next(it); // (under any name)
            return done(a + " is banned" + (rest.empty() ? "" : " (" + rest + ")"));
        }
        const bool was = bans.erase(a) > 0; return done(was ? a + " is unbanned" : a + " was not banned");
    }
    if (cmd == "banip" || cmd == "unbanip") {
        if (arg.empty()) return need((cmd + " IP" + (cmd == "banip" ? " [reason]" : "")).c_str());
        if (!console && byRole == AdminRole::Player) { reply = "only ops can do that"; return true; }
        if (cmd == "banip") { ipBans[arg] = rest; return done(arg + " is banned"); }
        const bool was = ipBans.erase(arg) > 0; return done(was ? arg + " is unbanned" : arg + " was not banned");
    }
    if (cmd == "owner") { // console only: the server's owner(s)
        if (!console) { reply = "only the server console can name owners"; return true; }
        if (arg.empty()) return need("owner NAME [UUID]");
        std::string id;
        if (!Bind(a, givenId, id)) { reply = noOne; return true; }
        owners[a] = id;
        return done(a + " is an owner");
    }
    if (cmd == "whitelist") {
        if (arg.empty() || a == "list") {
            std::string s = std::string("whitelist ") + (whitelist ? "on" : "off") + ":"; for (auto& [n, u] : allowed) s += " " + n;
            reply = s; return true;
        }
        if (!console && byRole == AdminRole::Player) { reply = "only ops can do that"; return true; }
        std::string who = Lower(rest); who = who.substr(0, who.find(' '));
        if (a == "on" || a == "off") { whitelist = a == "on"; return done(std::string("whitelist ") + (whitelist ? "on: only listed players (and ops) can join" : "off: anyone can join")); }
        if ((a == "add" || a == "remove") && who.empty()) return need("whitelist add|remove NAME");
        if (a == "add") { allowed[who] = Resolve(who); return done(who + " is on the whitelist"); } // (no id yet: by name until they've played)
        if (a == "remove") { const bool was = allowed.erase(who) > 0; return done(was ? who + " is off the whitelist" : who + " was not on it"); }
        return need("whitelist on|off|add NAME|remove NAME|list");
    }
    return false;
}

// --test-server-admin: the lists, the rules (ops can't touch owners; players can't do anything), roles by id (#1868), the
// join checks, and that admin.txt brings everything back after a restart.
bool Admin_Test() {
    namespace fs = std::filesystem;
    const fs::path tmp = fs::temp_directory_path() / "raylibgame_admin_test";
    std::error_code ec; fs::remove_all(tmp, ec);
    Admin_Load(tmp.string());
    // The players this server has seen (a test's player records): name -> id.
    const std::map<std::string, std::string> seen = { { "dee", "id-dee" }, { "alice", "id-alice" }, { "mallory", "id-mallory" }, { "carol", "id-carol" } };
    Admin_SetResolver([&seen](const std::string& n) { auto it = seen.find(n); return it == seen.end() ? std::string() : it->second; });
    bool ok = true;
    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: server admin: %s FAILED\n", what); ok = false; } };
    std::string r, why;
    check(Admin_Command("owner Dee", "console", r) && Admin_RoleOf("dee", "id-dee") == AdminRole::Owner, "console names an owner (by id)");
    check(Admin_Command("op Alice", "console", r) && Admin_RoleOf("ALICE", "id-alice") == AdminRole::Op, "op (names ignore case)");
    // (#1868) the role is the id's: the name with another id, or no id, gets nothing; the id under another name keeps it
    check(Admin_RoleOf("alice", "id-impostor") == AdminRole::Player && Admin_RoleOf("alice") == AdminRole::Player, "an op's name with another id is a player");
    check(Admin_RoleOf("alice2", "id-alice") == AdminRole::Op, "an op renamed keeps the role");
    check(Admin_Command("op zed", "console", r) && r.find("hasn't played") != std::string::npos && Admin_RoleOf("zed", "id-zed") == AdminRole::Player,
          "a name no one has played under can't be op'd");
    check(Admin_Command("op zed id-zed", "console", r) && Admin_RoleOf("zed", "id-zed") == AdminRole::Op && Admin_RoleOf("zed", "id-x") == AdminRole::Player,
          "the console may op a name by id");
    check(Admin_Command("deop zed", "console", r), "deop");
    check(Admin_Command("op bob", "mallory", r, "id-mallory") && r == "only ops can do that", "a player can't op");
    check(Admin_Command("op dee", "alice", r, "id-alice") && r.find("owner") != std::string::npos, "an op can't op an owner");
    check(Admin_Command("ban dee", "alice", r, "id-alice") && Admin_CheckJoin("dee", "", why, "id-dee"), "an op can't ban an owner");
    check(Admin_Command("ban mallory griefing", "alice", r, "id-alice") && !Admin_CheckJoin("Mallory", "", why, "id-mallory") && why == "banned: griefing", "an op bans; the reason is shown");
    check(!Admin_CheckJoin("mallory2", "", why, "id-mallory"), "a ban follows the id to a new name");
    check(Admin_Command("banip 10.0.0.9", "console", r) && !Admin_CheckJoin("eve", "10.0.0.9:51234", why) && Admin_CheckJoin("eve", "10.0.0.8:51234", why), "IP ban (port ignored)");
    check(Admin_Command("whitelist on", "console", r) && !Admin_CheckJoin("erin", "", why, "id-erin") && why == "not on this server's whitelist", "whitelist refuses others");
    check(Admin_CheckJoin("alice", "", why, "id-alice") && Admin_CheckJoin("dee", "", why, "id-dee"), "ops and owners pass the whitelist (by id)");
    check(!Admin_CheckJoin("alice", "", why, "id-impostor"), "an op's name with another id doesn't pass the whitelist");
    // (#1086) a menu-hosted game: its hosting player is the owner by their id, whatever the name
    Admin_SetHostOwner("0b1e2c3d-host");
    check(Admin_RoleOf("Host Player", "0b1e2c3d-host") == AdminRole::Owner, "the host by id is the owner (a name not in admin.txt)");
    check(Admin_RoleOf("Host Player", "9f9f9f9f-other") == AdminRole::Player && Admin_RoleOf("Host Player") == AdminRole::Player, "the host's name with another id is a player");
    check(Admin_Command("op carol", "Host Player", r, "0b1e2c3d-host") && Admin_RoleOf("carol", "id-carol") == AdminRole::Op, "the host by id may op");
    check(Admin_Command("deop carol", "Host Player", r, "9f9f9f9f-other") && r == "only ops can do that" && Admin_RoleOf("carol", "id-carol") == AdminRole::Op, "a spoof of the host's name may not");
    Admin_SetHostOwner("");
    check(Admin_RoleOf("Host Player", "0b1e2c3d-host") == AdminRole::Player, "a dedicated server (no host id): admin.txt alone");
    // The whitelist: a name not yet played under is let in by name; once it has an id, only that id.
    check(Admin_Command("whitelist add Frank", "alice", r, "id-alice") && Admin_CheckJoin("frank", "", why, "id-frank"), "whitelist add, before their first join");
    check(Admin_Command("whitelist add Carol", "alice", r, "id-alice") && Admin_CheckJoin("carol", "", why, "id-carol") && !Admin_CheckJoin("carol", "", why, "id-not-carol"),
          "whitelist add of a known player: by id");
    check(!Admin_Command("say hi", "console", r), "not an admin command");
    // A restart: everything comes back from admin.txt, with the ids.
    Admin_Load(tmp.string());
    check(Admin_RoleOf("dee", "id-dee") == AdminRole::Owner && Admin_RoleOf("alice", "id-alice") == AdminRole::Op && !Admin_CheckJoin("mallory", "", why, "id-mallory") &&
          !Admin_CheckJoin("x", "10.0.0.9:1", why) && Admin_CheckJoin("carol", "", why, "id-carol") && !Admin_CheckJoin("zed", "", why, "id-zed"), "admin.txt persists it all");
    check(Admin_Command("unban mallory", "alice", r, "id-alice") && Admin_Command("whitelist off", "alice", r, "id-alice") && Admin_CheckJoin("mallory", "", why, "id-mallory"), "unban, whitelist off");
    // An old admin.txt (names only): those owners and ops grant nothing until op'd again.
    if (FILE* f = fopen((tmp / "admin.txt").string().c_str(), "w")) { fputs("owner olddee\nop oldalice\n", f); fclose(f); }
    Admin_Load(tmp.string());
    check(Admin_RoleOf("olddee", "any-id") == AdminRole::Player && Admin_RoleOf("oldalice", "any-id") == AdminRole::Player, "name-only entries grant nothing");
    check(Admin_Command("ops", "console", r) && r.find("unbound") != std::string::npos, "ops lists them as unbound");
    // The audit log has the changes, attributed.
    int lines = 0; bool attributed = false;
    if (FILE* f = fopen((tmp / "audit.log").string().c_str(), "r")) { char b[512]; while (fgets(b, sizeof b, f)) { ++lines; attributed = attributed || strstr(b, "alice: ban mallory griefing") != nullptr; } fclose(f); }
    check(lines >= 8 && attributed, "audit log");
    printf("bench: server admin (#345, #1868): roles by id, bans (name, id, IP), whitelist, persistence, audit (%d lines) -> %s\n", lines, ok ? "PASS" : "FAIL");
    Admin_SetResolver(nullptr);
    fs::remove_all(tmp, ec);
    Admin_Load("");
    return ok;
}
