#include "serverconsole.h"
#include "utils/logger.h"

#include <atomic>
#include <chrono>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <sstream>
#include <thread>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
// Only what the console handler needs (windows.h clashes with raylib's names).
extern "C" {
typedef int BOOL_;
typedef unsigned long DWORD_;
typedef BOOL_(__stdcall* HandlerRoutine_)(DWORD_);
__declspec(dllimport) BOOL_ __stdcall SetConsoleCtrlHandler(HandlerRoutine_, BOOL_);
}
#else
#include <unistd.h>
#endif

// See serverconsole.h (#326).
namespace {
std::mutex inputMutex;
std::deque<std::string> input;
std::atomic<bool> shutdownRequested{ false }, shutdownDone{ false };
std::atomic<const char*> shutdownWhy{ "" };
// The log tee (ServerLog_Start): the console's fd and the pipe's write end, so ServerExit can close the pipe and wait
// until the tee thread has written everything (teeDone).
int teeConsole = -1, teePipeWrite = -1;
std::atomic<bool> teeDone{ false };

std::string Trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
}

#ifdef _WIN32
BOOL_ __stdcall OnConsoleEvent(DWORD_ type) {
    // 0 Ctrl+C, 1 Ctrl+Break, 2 close, 5 logoff, 6 shutdown. For close/logoff/shutdown Windows ends the process when the
    // handler returns (or after ~5 s): wait for the save (up to 4.5 s).
    ServerShutdown_Request(type == 2 ? "the console closed" : type >= 5 ? "Windows is shutting down" : "Ctrl+C");
    if (type == 2 || type >= 5)
        for (int i = 0; i < 450 && !shutdownDone; ++i) std::this_thread::sleep_for(std::chrono::milliseconds(10));
    return 1;
}
#endif
void OnSignal(int sig) { ServerShutdown_Request(sig == SIGINT ? "SIGINT" : "SIGTERM"); }
}

std::string ServerConfig_Text(const ServerConfig& c) {
    std::ostringstream s;
    s << "# RaylibGame dedicated server (#326). key = value; '#' starts a comment. Restart the server to apply changes.\n"
      << "save = " << c.save << "\n"
      << "port = " << c.port << "\n"
      << "max_players = " << c.maxPlayers << "\n"
      << "password = " << c.password << "\n"
      << "autosave_minutes = " << c.autosaveMinutes << "\n"
      << "motd = " << c.motd << "\n"
      << "seed = " << c.seed << "\n"
      << "log = " << c.log << "\n";
    return s.str();
}

bool ServerConfig_Load(const std::string& path, ServerConfig& out, std::string& error) {
    namespace fs = std::filesystem;
    if (!fs::exists(path)) {
        std::error_code ec;
        if (fs::path(path).has_parent_path()) fs::create_directories(fs::path(path).parent_path(), ec);
        if (FILE* f = fopen(path.c_str(), "w")) { fputs(ServerConfig_Text(out).c_str(), f); fclose(f); printf("server: wrote the default config to %s\n", path.c_str()); }
        else { error = "cannot write " + path; return false; }
    }
    FILE* f = fopen(path.c_str(), "r");
    if (!f) { error = "cannot read " + path; return false; }
    char buf[1024]; int lineNo = 0;
    ServerConfig c;
    while (fgets(buf, sizeof buf, f)) {
        ++lineNo;
        std::string line = buf;
        if (const size_t h = line.find('#'); h != std::string::npos && line.find_first_not_of(" \t") == h) continue; // a comment line (a '#' inside a value stays)
        line = Trim(line);
        if (line.empty()) continue;
        const size_t eq = line.find('=');
        if (eq == std::string::npos) { error = path + ":" + std::to_string(lineNo) + ": expected key = value"; fclose(f); return false; }
        const std::string key = Trim(line.substr(0, eq)), val = Trim(line.substr(eq + 1));
        auto num = [&](double& v) { char* end = nullptr; v = strtod(val.c_str(), &end); return end && end != val.c_str() && *end == 0; };
        double v = 0.0;
        if (key == "save") c.save = val;
        else if (key == "password") c.password = val;
        else if (key == "motd") c.motd = val;
        else if (key == "log") c.log = val;
        else if (key == "port" && num(v) && v >= 1 && v <= 65535) c.port = (int)v;
        else if (key == "max_players" && num(v) && v >= 1 && v <= 64) c.maxPlayers = (int)v;
        else if (key == "autosave_minutes" && num(v) && v >= 0) c.autosaveMinutes = v;
        else if (key == "seed" && num(v) && v >= 0) c.seed = (unsigned int)v;
        else if (key == "port" || key == "max_players" || key == "autosave_minutes" || key == "seed") {
            error = path + ":" + std::to_string(lineNo) + ": bad value for " + key + ": " + val; fclose(f); return false;
        } else printf("server: %s:%d: unknown key '%s' (ignored)\n", path.c_str(), lineNo, key.c_str());
    }
    fclose(f);
    out = c;
    return true;
}

void ServerLog_Start(const std::string& path) {
    namespace fs = std::filesystem;
    std::error_code ec;
    if (fs::path(path).has_parent_path()) fs::create_directories(fs::path(path).parent_path(), ec);
    FILE* log = fopen(path.c_str(), "a");
    if (!log) { printf("server: cannot open the log %s\n", path.c_str()); return; }
    fflush(stdout); fflush(stderr);
#ifdef _WIN32
    int fds[2];
    if (_pipe(fds, 1 << 16, _O_BINARY) != 0) { fclose(log); return; }
    const int console = _dup(1);
    _dup2(fds[1], 1); _dup2(fds[1], 2);
    const int readFd = fds[0];
    auto readChunk = [readFd](char* b, int n) { return _read(readFd, b, (unsigned)n); };
    auto writeConsole = [console](const char* b, int n) { _write(console, b, (unsigned)n); };
#else
    int fds[2];
    if (pipe(fds) != 0) { fclose(log); return; }
    const int console = dup(1);
    dup2(fds[1], 1); dup2(fds[1], 2);
    teeConsole = console; teePipeWrite = fds[1];
    const int readFd = fds[0];
    auto readChunk = [readFd](char* b, int n) { return (int)read(readFd, b, (size_t)n); };
    auto writeConsole = [console](const char* b, int n) { (void)!write(console, b, (size_t)n); };
#endif
    setvbuf(stdout, nullptr, _IONBF, 0);
    std::thread([log, readChunk, writeConsole]() {
        char buf[4096]; bool lineStart = true;
        for (;;) {
            const int n = readChunk(buf, sizeof buf);
            if (n <= 0) break;
            writeConsole(buf, n);
            for (int i = 0; i < n; ++i) {
                if (lineStart) {
                    const std::time_t t = std::time(nullptr); char ts[32];
                    std::strftime(ts, sizeof ts, "[%Y-%m-%d %H:%M:%S] ", std::localtime(&t));
                    fputs(ts, log); lineStart = false;
                }
                fputc(buf[i], log);
                if (buf[i] == '\n') { lineStart = true; fflush(log); }
            }
        }
        fflush(log);
        teeDone = true;
    }).detach();
    printf("server: logging to %s\n", path.c_str());
}

void ServerConsole_Start() {
    std::thread([]() {
        std::string line;
        while (std::getline(std::cin, line)) {
            line = Trim(line);
            if (line.empty()) continue;
            std::lock_guard<std::mutex> lock(inputMutex);
            input.push_back(line);
        }
        // stdin closed (a detached server with no console): keep running; "stop" can't come any more, signals still can.
    }).detach();
}

bool ServerConsole_Next(std::string& line) {
    std::lock_guard<std::mutex> lock(inputMutex);
    if (input.empty()) return false;
    line = input.front(); input.pop_front();
    return true;
}

void ServerShutdown_Install() {
    std::signal(SIGINT, OnSignal);
    std::signal(SIGTERM, OnSignal);
#ifdef _WIN32
    SetConsoleCtrlHandler(OnConsoleEvent, 1);
#endif
}
void ServerShutdown_Request(const char* why) { if (!shutdownRequested.exchange(true)) shutdownWhy = why; }
bool ServerShutdown_Requested(const char** why) { if (why) *why = shutdownWhy; return shutdownRequested; }
void ServerShutdown_Done() { shutdownDone = true; }

void ServerExit(int code) {
#ifdef _WIN32
    (void)code;
#else
    Log_Flush();
    fflush(stdout); fflush(stderr);
    if (teeConsole >= 0) {
        // The pipe's write ends are fds 1, 2 and teePipeWrite: point 1 and 2 back at the console and close the last one,
        // so the tee thread reads the rest, then end of file (at most 2 s).
        dup2(teeConsole, 1); dup2(teeConsole, 2); close(teePipeWrite);
        for (int i = 0; i < 200 && !teeDone; ++i) std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    _exit(code);
#endif
}
