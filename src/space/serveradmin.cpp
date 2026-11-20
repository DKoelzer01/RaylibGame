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

// See serveradmin.h (#345).
namespace {
std::string dir;
std::set<std::string> owners, ops, allowed;      // lower-cased names
std::map<std::string, std::string> bans, ipBans; // name / ip -> reason
bool whitelist = false;

std::string Lower(std::string s) { for (char& c : s) c = (char)std::tolower((unsigned char)c); return s; }
std::string Path(const char* file) { return dir.empty() ? std::string(file) : dir + "/" + file; }

void Save() {
    std::error_code ec;
    if (!dir.empty()) std::filesystem::create_directories(dir, ec);
    FILE* f = fopen(Path("admin.txt").c_str(), "w");
    if (!f) return;
    fputs("# RaylibGame server admin state (#345); rewritten by the server on every change.\n", f);
    for (const std::string& n : owners) fprintf(f, "owner %s\n", n.c_str());
    for (const std::string& n : ops) fprintf(f, "op %s\n", n.c_str());
    for (const auto& [n, r] : bans) fprintf(f, "ban %s%s%s\n", n.c_str(), r.empty() ? "" : " ", r.c_str());
    for (const auto& [ip, r] : ipBans) fprintf(f, "banip %s%s%s\n", ip.c_str(), r.empty() ? "" : " ", r.c_str());
    for (const std::string& n : allowed) fprintf(f, "allow %s\n", n.c_str());
    if (whitelist) fputs("whitelist on\n", f);
    fclose(f);
}
}

void Admin_Load(const std::string& d) {
    dir = d;
    owners.clear(); ops.clear(); allowed.clear(); bans.clear(); ipBans.clear(); whitelist = false;
    FILE* f = fopen(Path("admin.txt").c_str(), "r");
    if (!f) return;
    char buf[512];
    while (fgets(buf, sizeof buf, f)) {
        std::istringstream in(buf);
        std::string kind, who, reason;
        if (!(in >> kind) || kind[0] == '#') continue;
        in >> who; std::getline(in, reason);
        reason.erase(0, reason.find_first_not_of(" \t")); while (!reason.empty() && (reason.back() == '\n' || reason.back() == '\r')) reason.pop_back();
        if (kind == "owner") owners.insert(Lower(who));
        else if (kind == "op") ops.insert(Lower(who));
        else if (kind == "ban") bans[Lower(who)] = reason;
        else if (kind == "banip") ipBans[who] = reason;
        else if (kind == "allow") allowed.insert(Lower(who));
        else if (kind == "whitelist") whitelist = who == "on";
    }
    fclose(f);
}

static std::string hostOwnerId; // (#1086)
void Admin_SetHostOwner(const std::string& uuid) { hostOwnerId = uuid; }
AdminRole Admin_RoleOf(const std::string& name, const std::string& uuid) {
    if (!uuid.empty() && uuid == hostOwnerId) return AdminRole::Owner; // the player hosting this game from their menu
    const std::string n = Lower(name);
    return owners.count(n) ? AdminRole::Owner : ops.count(n) ? AdminRole::Op : AdminRole::Player;
}
const char* Admin_RoleName(AdminRole r) { return r == AdminRole::Owner ? "owner" : r == AdminRole::Op ? "op" : "player"; }

bool Admin_CheckJoin(const std::string& name, const std::string& ip, std::string& why) {
    const std::string n = Lower(name);
    if (auto it = bans.find(n); it != bans.end()) { why = "banned" + (it->second.empty() ? std::string() : ": " + it->second); return false; }
    const std::string host = ip.substr(0, ip.rfind(':') == std::string::npos ? ip.size() : ip.rfind(':')); // "1.2.3.4:5678" -> the address
    if (auto it = ipBans.find(host); !host.empty() && it != ipBans.end()) { why = "banned" + (it->second.empty() ? std::string() : ": " + it->second); return false; }
    if (whitelist && !allowed.count(n) && !ops.count(n) && !owners.count(n)) { why = "not on this server's whitelist"; return false; }
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
    if (cmd == "ops") {
        std::string s = "owners:"; for (auto& n : owners) s += " " + n;
        s += "; ops:"; for (auto& n : ops) s += " " + n;
        reply = s; return true;
    }
    if (cmd == "bans") {
        std::string s = "banned names:"; for (auto& [n, r] : bans) s += " " + n;
        s += "; banned addresses:"; for (auto& [ip, r] : ipBans) s += " " + ip;
        reply = s; return true;
    }
    if (cmd == "op" || cmd == "deop" || cmd == "ban" || cmd == "unban") {
        if (arg.empty()) return need((cmd + " NAME" + (cmd == "ban" ? " [reason]" : "")).c_str());
        if (!console && byRole == AdminRole::Player) { reply = "only ops can do that"; return true; }
        if (!console && owners.count(a)) { reply = "an owner can only be changed from the server console"; return true; }
        if (cmd == "op") { ops.insert(a); return done(a + " is now an op"); }
        if (cmd == "deop") { const bool was = ops.erase(a) > 0; return done(was ? a + " is no longer an op" : a + " was not an op"); }
        if (cmd == "ban") { bans[a] = rest; ops.erase(a); return done(a + " is banned" + (rest.empty() ? "" : " (" + rest + ")")); }
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
        if (arg.empty()) return need("owner NAME");
        owners.insert(a); return done(a + " is an owner");
    }
    if (cmd == "whitelist") {
        if (arg.empty() || a == "list") {
            std::string s = std::string("whitelist ") + (whitelist ? "on" : "off") + ":"; for (auto& n : allowed) s += " " + n;
            reply = s; return true;
        }
        if (!console && byRole == AdminRole::Player) { reply = "only ops can do that"; return true; }
        std::string who = Lower(rest); who = who.substr(0, who.find(' '));
        if (a == "on" || a == "off") { whitelist = a == "on"; return done(std::string("whitelist ") + (whitelist ? "on: only listed players (and ops) can join" : "off: anyone can join")); }
        if ((a == "add" || a == "remove") && who.empty()) return need("whitelist add|remove NAME");
        if (a == "add") { allowed.insert(who); return done(who + " is on the whitelist"); }
        if (a == "remove") { const bool was = allowed.erase(who) > 0; return done(was ? who + " is off the whitelist" : who + " was not on it"); }
        return need("whitelist on|off|add NAME|remove NAME|list");
    }
    return false;
}

// --test-server-admin: the lists, the rules (ops can't touch owners; players can't do anything), the join checks, and
// that admin.txt brings everything back after a restart.
bool Admin_Test() {
    namespace fs = std::filesystem;
    const fs::path tmp = fs::temp_directory_path() / "raylibgame_admin_test";
    std::error_code ec; fs::remove_all(tmp, ec);
    Admin_Load(tmp.string());
    bool ok = true;
    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: server admin: %s FAILED\n", what); ok = false; } };
    std::string r, why;
    check(Admin_Command("owner Dee", "console", r) && Admin_RoleOf("dee") == AdminRole::Owner, "console names an owner");
    check(Admin_Command("op Alice", "console", r) && Admin_RoleOf("ALICE") == AdminRole::Op, "op (names ignore case)");
    check(Admin_Command("op bob", "mallory", r) && r == "only ops can do that" && Admin_RoleOf("bob") == AdminRole::Player, "a player can't op");
    check(Admin_Command("ban dee", "alice", r) && Admin_CheckJoin("dee", "", why), "an op can't ban an owner");
    check(Admin_Command("ban mallory griefing", "alice", r) && !Admin_CheckJoin("Mallory", "", why) && why == "banned: griefing", "an op bans; the reason is shown");
    check(Admin_Command("banip 10.0.0.9", "console", r) && !Admin_CheckJoin("eve", "10.0.0.9:51234", why) && Admin_CheckJoin("eve", "10.0.0.8:51234", why), "IP ban (port ignored)");
    check(Admin_Command("whitelist on", "console", r) && !Admin_CheckJoin("carol", "", why) && why == "not on this server's whitelist", "whitelist refuses others");
    check(Admin_CheckJoin("alice", "", why) && Admin_CheckJoin("dee", "", why), "ops and owners pass the whitelist");
    // (#1086) a menu-hosted game: its hosting player is the owner by their id, whatever the name; a joiner using the
    // host's name (or an owner's) without the id gets nothing from the name it borrowed beyond what admin.txt says
    Admin_SetHostOwner("0b1e2c3d-host");
    check(Admin_RoleOf("Host Player", "0b1e2c3d-host") == AdminRole::Owner, "the host by id is the owner (a name not in admin.txt)");
    check(Admin_RoleOf("Host Player", "9f9f9f9f-other") == AdminRole::Player && Admin_RoleOf("Host Player") == AdminRole::Player, "the host's name with another id is a player");
    check(Admin_Command("op mallory", "Host Player", r, "0b1e2c3d-host") && Admin_RoleOf("mallory") == AdminRole::Op, "the host by id may op");
    check(Admin_Command("deop mallory", "Host Player", r, "9f9f9f9f-other") && r == "only ops can do that" && Admin_RoleOf("mallory") == AdminRole::Op, "a spoof of the host's name may not");
    Admin_SetHostOwner("");
    check(Admin_RoleOf("Host Player", "0b1e2c3d-host") == AdminRole::Player, "a dedicated server (no host id): admin.txt alone");
    check(Admin_Command("whitelist add Carol", "alice", r) && Admin_CheckJoin("carol", "", why), "whitelist add");
    check(!Admin_Command("say hi", "console", r), "not an admin command");
    // A restart: everything comes back from admin.txt.
    Admin_Load(tmp.string());
    check(Admin_RoleOf("dee") == AdminRole::Owner && Admin_RoleOf("alice") == AdminRole::Op && !Admin_CheckJoin("mallory", "", why) &&
          !Admin_CheckJoin("x", "10.0.0.9:1", why) && Admin_CheckJoin("carol", "", why) && !Admin_CheckJoin("zed", "", why), "admin.txt persists it all");
    check(Admin_Command("unban mallory", "alice", r) && Admin_Command("whitelist off", "alice", r) && Admin_CheckJoin("mallory", "", why), "unban, whitelist off");
    // The audit log has the changes, attributed.
    int lines = 0; bool attributed = false;
    if (FILE* f = fopen((tmp / "audit.log").string().c_str(), "r")) { char b[512]; while (fgets(b, sizeof b, f)) { ++lines; attributed = attributed || strstr(b, "alice: ban mallory griefing") != nullptr; } fclose(f); }
    check(lines >= 8 && attributed, "audit log");
    printf("bench: server admin (#345): roles, bans (name, IP), whitelist, persistence, audit (%d lines) -> %s\n", lines, ok ? "PASS" : "FAIL");
    fs::remove_all(tmp, ec);
    Admin_Load("");
    return ok;
}
