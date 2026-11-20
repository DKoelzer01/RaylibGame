#ifndef SERVERCONSOLE_H
#define SERVERCONSOLE_H

// Dedicated server (#326): `--server [server.cfg]` runs --host from a config file, with a console on stdin, a log file
// and a clean shutdown. Pterodactyl (tools/pterodactyl/egg-raylibgame.json) starts it with that line, edits server.cfg
// through its config parser, and stops it by typing "stop".
//
// server.cfg: "key = value" lines, '#' comments. Keys: save, port, max_players, password, autosave_minutes, motd, seed
// (a new world only), log. A missing file is written with the defaults first. Command-line flags given after --server
// override it.
//
// Console (stdin, one command a line): help, status, list, save, say <text>, kick <name> [reason], stop.
// Shutdown: "stop", Ctrl+C, closing the console window, a logoff or shutdown (Windows), SIGTERM/SIGINT: the world is
// saved, every client hears "server closed", then the process exits (Windows gives a closing console ~5 s: enough).
// Log: everything the process prints also goes to the log file, each line stamped with the time.

#include <cstdint>
#include <string>

struct ServerConfig {
    std::string save = "server/world.sav";
    int port = 27015;
    int maxPlayers = 8;
    std::string password;
    double autosaveMinutes = 30.0;
    std::string motd;
    unsigned int seed = 0;          // 0: the save's own, or a random one for a new world
    std::string log = "server/server.log";
    std::string ownerId;            // (#1086) a menu-hosted game: the hosting player's id (UUID), always the owner whatever
                                    // name anyone uses; empty on a dedicated server (its admin.txt decides)
};

// Reads `path` (writing the defaults there first when it does not exist). False, with why, on a malformed line.
bool ServerConfig_Load(const std::string& path, ServerConfig& out, std::string& error);
std::string ServerConfig_Text(const ServerConfig& c); // the file's text (the defaults file, the "status" command)

void ServerLog_Start(const std::string& path);         // tee stdout (and stderr) into the file from here on
void ServerConsole_Start();                            // a thread reading stdin lines
bool ServerConsole_Next(std::string& line);            // the next command typed (non-blocking)

void ServerShutdown_Install();                         // Ctrl+C / close / logoff / SIGTERM -> a requested shutdown
void ServerShutdown_Request(const char* why);
bool ServerShutdown_Requested(const char** why = nullptr);
void ServerShutdown_Done();                            // saved and closed: a waiting console-close handler may return
// The last call of a --server process (after CloseWindow). Linux (#792): flushes the log tee and ends the process with
// _exit, because glibc's exit() flushes every stdio stream and would wait forever for stdin's lock, held by the
// console thread blocked reading a stdin that stays open (the panel's console). Windows: returns (a normal exit).
void ServerExit(int code);

#endif
