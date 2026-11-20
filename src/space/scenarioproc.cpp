// Child processes for `--scenario tag:TAG` (#870): each scenario runs in a game process of its own, with a hard timeout.
// Kept apart from scenario.cpp: windows.h and raylib.h clash (CloseWindow, Rectangle, ...).
#include <string>
#include <cstdlib>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

// Runs `cmdline` (stdout/stderr inherited), waits up to timeoutSec; past it the process is ended. Returns its exit code
// (-1 when it could not start, -2 when it timed out).
int ScenarioProc_Run(const std::string& cmdline, double timeoutSec) {
    STARTUPINFOA si{}; si.cb = sizeof si;
    PROCESS_INFORMATION pi{};
    std::string cmd = cmdline; // CreateProcessA may write to it
    if (!CreateProcessA(nullptr, cmd.data(), nullptr, nullptr, TRUE, BELOW_NORMAL_PRIORITY_CLASS, nullptr, nullptr, &si, &pi)) return -1;
    const DWORD ms = timeoutSec > 0.0 ? (DWORD)(timeoutSec * 1000.0) : INFINITE;
    int rc;
    if (WaitForSingleObject(pi.hProcess, ms) == WAIT_TIMEOUT) { TerminateProcess(pi.hProcess, 124); WaitForSingleObject(pi.hProcess, 5000); rc = -2; }
    else { DWORD code = 0; GetExitCodeProcess(pi.hProcess, &code); rc = (int)code; }
    CloseHandle(pi.hThread); CloseHandle(pi.hProcess);
    return rc;
}
std::string ScenarioProc_Exe() {
    char buf[MAX_PATH * 2] = {};
    GetModuleFileNameA(nullptr, buf, sizeof buf);
    return buf;
}
int ScenarioProc_Pid() { return (int)GetCurrentProcessId(); }
#else
#include <unistd.h>
#include <climits>
int ScenarioProc_Run(const std::string& cmdline, double timeoutSec) {
    const std::string c = timeoutSec > 0.0 ? "timeout " + std::to_string((int)timeoutSec) + " " + cmdline : cmdline;
    const int r = std::system(c.c_str());
    if (r == -1) return -1;
    const int code = (r >> 8) & 255;
    return code == 124 ? -2 : code;
}
std::string ScenarioProc_Exe() {
    char buf[PATH_MAX] = {};
    const ssize_t n = readlink("/proc/self/exe", buf, sizeof buf - 1);
    return n > 0 ? std::string(buf, (size_t)n) : std::string("./game");
}
int ScenarioProc_Pid() { return (int)getpid(); }
#endif
