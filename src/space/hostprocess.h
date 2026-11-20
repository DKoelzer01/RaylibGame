#ifndef HOSTPROCESS_H
#define HOSTPROCESS_H

// The menu's "Host a world" (#341): the game starts its own dedicated server (this exe, `--server <cfg>`, headless) as
// a child process and joins it over loopback like any client, so hosting runs the same tested server as a dedicated
// host. The child's stdin is a pipe: "stop" (its console command, #326) saves the world and ends it; its output goes to
// a log file beside the save. Windows only (elsewhere Start fails with why).

#include <string>

// The child runs in a kill-on-close job object: if the game ends or crashes, the server ends with it (never an orphan
// holding the port).
bool HostProcess_Start(const std::string& configPath, const std::string& logPath, std::string& error);
bool HostProcess_PortFree(int port); // nothing listening on that UDP port here (a bind test)
bool HostProcess_Running();
// Sends "stop" and waits up to `seconds` for it to save and exit (then ends it). True if it exited by itself.
bool HostProcess_Stop(double seconds);

#endif
