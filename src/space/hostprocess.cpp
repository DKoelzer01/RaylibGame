#include "hostprocess.h"

// See hostprocess.h (#341). <windows.h> stays in this file (its names clash with raylib's).
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <winsock2.h>
#include <ws2tcpip.h>
#include <windows.h>
#include <iphlpapi.h>
#include <vector>

namespace {
PROCESS_INFORMATION child{};
HANDLE stdinWrite = nullptr;
HANDLE job = nullptr;   // the child is in it: when the game ends (or crashes) the job's last handle closes and it ends too
bool running = false;

// Is any UDP endpoint of this family on the port? (The system's own table: a test bind succeeded on Windows while a
// GameNetworkingSockets server held both 0.0.0.0 and [::] on that port, whatever the socket options.)
bool UdpPortListed(ULONG family, int port) {
    DWORD size = 0;
    GetExtendedUdpTable(nullptr, &size, FALSE, family, UDP_TABLE_BASIC, 0);
    std::vector<unsigned char> buf(size + 256);
    size = (DWORD)buf.size();
    if (GetExtendedUdpTable(buf.data(), &size, FALSE, family, UDP_TABLE_BASIC, 0) != NO_ERROR) return false;
    const u_short want = htons((u_short)port);
    if (family == AF_INET) {
        const auto* t = (const MIB_UDPTABLE*)buf.data();
        for (DWORD i = 0; i < t->dwNumEntries; ++i) if ((u_short)t->table[i].dwLocalPort == want) return true;
    } else {
        const auto* t = (const MIB_UDP6TABLE*)buf.data();
        for (DWORD i = 0; i < t->dwNumEntries; ++i) if ((u_short)t->table[i].dwLocalPort == want) return true;
    }
    return false;
}
}

bool HostProcess_PortFree(int port) { return !UdpPortListed(AF_INET, port) && !UdpPortListed(AF_INET6, port); }

bool HostProcess_Start(const std::string& configPath, const std::string& logPath, std::string& error) {
    if (HostProcess_Running()) { error = "a server is already running"; return false; }
    char exe[MAX_PATH] = {};
    if (!GetModuleFileNameA(nullptr, exe, MAX_PATH)) { error = "could not find the game's exe"; return false; }
    SECURITY_ATTRIBUTES sa{ sizeof(sa), nullptr, TRUE };
    HANDLE inRead = nullptr, inWrite = nullptr;
    if (!CreatePipe(&inRead, &inWrite, &sa, 0)) { error = "could not make the server's console pipe"; return false; }
    SetHandleInformation(inWrite, HANDLE_FLAG_INHERIT, 0); // (our end stays ours)
    HANDLE log = CreateFileA(logPath.c_str(), GENERIC_WRITE, FILE_SHARE_READ, &sa, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
    STARTUPINFOA si{};
    si.cb = sizeof(si);
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdInput = inRead;
    si.hStdOutput = log != INVALID_HANDLE_VALUE ? log : nullptr;
    si.hStdError = si.hStdOutput;
    std::string cmd = std::string("\"") + exe + "\" --server \"" + configPath + "\"";
    std::string mutableCmd = cmd;
    // Suspended until it is in the kill-on-close job: no moment where a crash of ours leaves it running alone.
    const BOOL ok = CreateProcessA(nullptr, mutableCmd.data(), nullptr, nullptr, TRUE, CREATE_NO_WINDOW | CREATE_SUSPENDED, nullptr, nullptr, &si, &child);
    CloseHandle(inRead);
    if (log != INVALID_HANDLE_VALUE) CloseHandle(log);
    if (!ok) { CloseHandle(inWrite); error = "could not start the server (" + std::to_string(GetLastError()) + ")"; return false; }
    if (!job) {
        job = CreateJobObjectA(nullptr, nullptr);
        JOBOBJECT_EXTENDED_LIMIT_INFORMATION li{};
        li.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
        if (job) SetInformationJobObject(job, JobObjectExtendedLimitInformation, &li, sizeof(li));
    }
    if (job) AssignProcessToJobObject(job, child.hProcess);
    ResumeThread(child.hThread);
    stdinWrite = inWrite;
    running = true;
    return true;
}

bool HostProcess_Running() {
    if (!running) return false;
    if (WaitForSingleObject(child.hProcess, 0) == WAIT_OBJECT_0) { // it ended (by itself, or a crash)
        CloseHandle(child.hProcess); CloseHandle(child.hThread);
        if (stdinWrite) CloseHandle(stdinWrite);
        stdinWrite = nullptr; child = {}; running = false;
    }
    return running;
}

bool HostProcess_Stop(double seconds) {
    if (!HostProcess_Running()) return true;
    DWORD n = 0;
    WriteFile(stdinWrite, "stop\n", 5, &n, nullptr);
    const bool exited = WaitForSingleObject(child.hProcess, (DWORD)(seconds * 1000.0)) == WAIT_OBJECT_0;
    if (!exited) TerminateProcess(child.hProcess, 1); // (it didn't save in time: the last autosave stands)
    HostProcess_Running(); // (clean up)
    if (running) { CloseHandle(child.hProcess); CloseHandle(child.hThread); if (stdinWrite) CloseHandle(stdinWrite); stdinWrite = nullptr; child = {}; running = false; }
    return exited;
}
#else
bool HostProcess_PortFree(int) { return true; }
bool HostProcess_Start(const std::string&, const std::string&, std::string& error) { error = "hosting from the menu is Windows-only for now"; return false; }
bool HostProcess_Running() { return false; }
bool HostProcess_Stop(double) { return true; }
#endif
