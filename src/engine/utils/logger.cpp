#include "logger.h"
#include <algorithm>
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <mutex>
#include <exception>
#include <new>
#include <csignal>
#include <string>
#include <thread>
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <psapi.h>
#include <dbghelp.h>
#else
#include <csignal>
#include <execinfo.h>
#include <fcntl.h>
#include <unistd.h>
#endif

// (#1526, #1527) The log is written through the process's own append-only handle and its own buffer, not a FILE*:
// - a crash handler can write without taking the CRT's FILE lock (another thread may hold it mid-line: the old handler
//   then deadlocked, or the log stopped mid-line with no crash entry);
// - every write lands at the end of the file, and the file is opened so the next run can rename it away (log.prev.txt):
//   a hung, killed process that still holds it writes into the old run's log, never into the middle of the new one.
namespace {
    std::mutex gMutex;
#ifdef _WIN32
    HANDLE gOut = INVALID_HANDLE_VALUE;
#else
    int gOut = -1;
#endif
    char gBuf[1 << 16]; // lines not yet written (gMutex)
    size_t gLen = 0;
    char gBuild[128] = "";
    char gCrashDir[260] = "saves/logs";
    const char* volatile gAbortWhy = nullptr; // what led to abort(), for its line: an assertion, LOG_FATAL, ...
    const auto gStart = std::chrono::steady_clock::now();
    std::atomic<int> gThreadCount{ 0 };
    thread_local char tName[16] = "";
    const char* const kNames[] = { "TRACE", "DEBUG", "INFO", "WARN", "ERROR", "FATAL", "OFF" };
    double gLastFlush = 0.0; // run time of the last flush (gMutex): lines reach the file within a second

    bool LogOpen() {
#ifdef _WIN32
        return gOut != INVALID_HANDLE_VALUE;
#else
        return gOut >= 0;
#endif
    }
    void RawWrite(const char* p, size_t n) {
        if (!LogOpen()) return;
#ifdef _WIN32
        while (n > 0) { DWORD w = 0; if (!WriteFile(gOut, p, (DWORD)std::min<size_t>(n, 1u << 20), &w, nullptr) || w == 0) return; p += w; n -= w; }
#else
        while (n > 0) { const ssize_t w = write(gOut, p, n); if (w <= 0) return; p += w; n -= (size_t)w; }
#endif
    }
    void FlushLocked() { RawWrite(gBuf, gLen); gLen = 0; }
    void AppendLocked(const char* p, size_t n) {
        if (gLen + n > sizeof(gBuf)) FlushLocked();
        if (n > sizeof(gBuf)) { RawWrite(p, n); return; }
        memcpy(gBuf + gLen, p, n); gLen += n;
    }
#if RAYLIBGAME_LOGGING
    // A crash: the buffered lines, then the crash's own. The lock only if it is free (the crashing thread may hold it).
    void CrashFlush() {
        if (gMutex.try_lock()) { FlushLocked(); gMutex.unlock(); }
        else { RawWrite(gBuf, gLen); gLen = 0; }
    }
    void CrashPrintf(const char* fmt, ...) {
        char b[768];
        va_list a; va_start(a, fmt); const int n = vsnprintf(b, sizeof b, fmt, a); va_end(a);
        if (n > 0) RawWrite(b, std::min<size_t>((size_t)n, sizeof b - 1));
    }
#endif
#ifdef _WIN32
    void FatalBox(const char* text) { if (Log_Interactive()) MessageBoxA(nullptr, text, "RaylibGame", MB_OK | MB_ICONERROR); }
#endif
#if defined(_WIN32) && RAYLIBGAME_LOGGING
    double CommitMB() { PROCESS_MEMORY_COUNTERS pmc{}; return GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof pmc) ? pmc.PagefileUsage / 1048576.0 : 0.0; }
    // "game.exe+0x1234" or "nvoglv64.dll+0x5678": the module an address is in, and the offset into it
    void ModuleOf(const void* at, char* out, size_t n) {
        HMODULE m = nullptr;
        char path[MAX_PATH] = "?";
        if (GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT, (LPCSTR)at, &m) && m)
            GetModuleFileNameA(m, path, sizeof path);
        const char* base = path;
        for (const char* c = path; *c; ++c) if (*c == '\\' || *c == '/') base = c + 1;
        snprintf(out, n, "%s+0x%llx", base, (unsigned long long)((uintptr_t)at - (uintptr_t)m));
    }
    // A minidump next to the other logs (open it with the exe and its symbols in WinDbg or Visual Studio).
    void WriteMinidump(EXCEPTION_POINTERS* ep) {
        char dir[260]; snprintf(dir, sizeof dir, "%s", gCrashDir);
        for (char* c = dir; *c; ++c) if (*c == '/' || *c == '\\') { const char k = *c; *c = 0; CreateDirectoryA(dir, nullptr); *c = k; }
        CreateDirectoryA(dir, nullptr);
        SYSTEMTIME t; GetLocalTime(&t);
        char path[400];
        snprintf(path, sizeof path, "%s/crash-%04d%02d%02d-%02d%02d%02d-%lu.dmp", dir, t.wYear, t.wMonth, t.wDay, t.wHour, t.wMinute, t.wSecond, (unsigned long)GetCurrentProcessId());
        HANDLE f = CreateFileA(path, GENERIC_WRITE, 0, nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
        if (f == INVALID_HANDLE_VALUE) { CrashPrintf("FATAL minidump: could not create %s\n", path); return; }
        MINIDUMP_EXCEPTION_INFORMATION mei{ GetCurrentThreadId(), ep, FALSE };
        const BOOL ok = MiniDumpWriteDump(GetCurrentProcess(), GetCurrentProcessId(), f, (MINIDUMP_TYPE)(MiniDumpWithThreadInfo | MiniDumpWithIndirectlyReferencedMemory),
                                          ep ? &mei : nullptr, nullptr, nullptr);
        CloseHandle(f);
        CrashPrintf(ok ? "FATAL minidump: %s\n" : "FATAL minidump failed: %s\n", path);
    }
    LONG WINAPI FlushOnCrash(EXCEPTION_POINTERS* ep) {
        // The crash itself in the log (#902, #1526): what, where (module + offset: resolve game.exe's with its symbols,
        // nm -C, image base 0x140000000), the stack (this filter runs on the crashing thread), the build, a minidump.
        static std::atomic<bool> once{ false };
        if (once.exchange(true)) return EXCEPTION_CONTINUE_SEARCH; // a second thread crashing meanwhile: the first one reports
        CrashFlush();
        const EXCEPTION_RECORD* er = ep ? ep->ExceptionRecord : nullptr;
        const void* at = er ? er->ExceptionAddress : nullptr;
        char where[MAX_PATH + 32];
        ModuleOf(at, where, sizeof where);
        const unsigned long code = er ? (unsigned long)er->ExceptionCode : 0ul;
        const char* name = code == EXCEPTION_ACCESS_VIOLATION ? " (access violation)" : code == EXCEPTION_STACK_OVERFLOW ? " (stack overflow)"
                         : code == EXCEPTION_INT_DIVIDE_BY_ZERO ? " (integer divide by zero)" : code == EXCEPTION_ILLEGAL_INSTRUCTION ? " (illegal instruction)" : "";
        CrashPrintf("\nFATAL crash: exception 0x%08lx%s at %p = %s (commit %.0f MB, thread %lu '%s', build %s)\n", code, name, at, where, CommitMB(),
                    (unsigned long)GetCurrentThreadId(), tName, gBuild[0] ? gBuild : "?");
        if (code == EXCEPTION_ACCESS_VIOLATION && er->NumberParameters >= 2)
            CrashPrintf("FATAL crash: %s address %p\n", er->ExceptionInformation[0] == 0 ? "reading" : er->ExceptionInformation[0] == 8 ? "executing" : "writing", (void*)er->ExceptionInformation[1]);
        void* frames[32];
        const USHORT n = RtlCaptureStackBackTrace(0, 32, frames, nullptr);
        CrashPrintf("FATAL stack:");
        for (USHORT i = 0; i < n; ++i) { ModuleOf(frames[i], where, sizeof where); CrashPrintf(" %s", where); }
        CrashPrintf("\n");
        if (code != EXCEPTION_STACK_OVERFLOW) WriteMinidump(ep); // (no stack left for dbghelp after an overflow)
        FatalBox("The game crashed and has to close. The details are in log.txt.");
        return EXCEPTION_CONTINUE_SEARCH;
    }
    // An exception nobody caught (std::bad_alloc at the commit limit, #837 / #902): std::terminate aborted with no line and
    // no Windows error event. Now: what it was, the commit charge, a message box.
    [[noreturn]] void OnTerminate() {
        const char* what = "an unknown exception";
        bool oom = false;
        if (std::exception_ptr e = std::current_exception()) {
            try { std::rethrow_exception(e); }
            catch (const std::bad_alloc&) { what = "std::bad_alloc (out of memory)"; oom = true; }
            catch (const std::exception& x) { what = x.what(); }
            catch (...) {}
        } else what = "std::terminate with no exception (a std::thread destroyed while still joinable?)";
        CrashFlush();
        CrashPrintf("\nFATAL uncaught exception: %s (commit %.0f MB, thread '%s', build %s)\n", what, CommitMB(), tName, gBuild[0] ? gBuild : "?");
        fprintf(stderr, "FATAL uncaught exception: %s\n", what);
        FatalBox(oom ? "The game ran out of memory and has to close. The details are in log.txt."
                     : "The game hit a fatal error and has to close. The details are in log.txt.");
        gAbortWhy = "the uncaught exception above";
        std::abort();
    }
    // (#1153) abort() (a failed assert or IM_ASSERT, LOG_FATAL, terminate) ended the game with no line here and no Windows
    // error event: it looked like a clean exit. Now it leaves a line (with what led to it) and the error box.
    void OnAbort(int) {
        CrashFlush();
        CrashPrintf("\nFATAL abort() after %s (commit %.0f MB, thread '%s', build %s)\n", gAbortWhy ? gAbortWhy : "no message (a library's own abort?)",
                    CommitMB(), tName, gBuild[0] ? gBuild : "?");
        fprintf(stderr, "FATAL abort()\n");
        if (!gAbortWhy) FatalBox("The game hit a fatal error and has to close. The details are in log.txt.");
        signal(SIGABRT, SIG_DFL);
        std::_Exit(3);
    }
    // (#1527) Ctrl+C / Ctrl+Break / closing the console (how scripts and timeouts stop a run): before, Windows ran exit()
    // under the game's threads, a still-joinable std::thread (the factory's) called std::terminate, the fatal box came
    // up and the process could hang on it. Now: one line, and the process ends at once (no save: it was killed).
    BOOL WINAPI OnConsoleStop(DWORD ev) {
        const char* what = ev == CTRL_C_EVENT ? "Ctrl+C" : ev == CTRL_BREAK_EVENT ? "Ctrl+Break" : ev == CTRL_CLOSE_EVENT ? "the console closed"
                         : ev == CTRL_LOGOFF_EVENT ? "a logoff" : "a Windows shutdown";
        CrashFlush();
        CrashPrintf("%9.3f stopped from outside (%s): ending at once, nothing saved\n",
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - gStart).count(), what);
        TerminateProcess(GetCurrentProcess(), 0xC000013A); // STATUS_CONTROL_C_EXIT, as Windows' own handler exits
        return TRUE;
    }
#elif RAYLIBGAME_LOGGING
    // Linux (#792): the same on a crash signal, then the default action (core dump / exit status) runs.
    void FlushOnCrash(int sig) {
        CrashFlush();
        void* frames[48]; // the stack, as module(+offset): addr2line -e <the unstripped game> resolves it
        const int n = backtrace(frames, 48);
        fprintf(stderr, "crash: signal %d, stack:\n", sig);
        backtrace_symbols_fd(frames, n, 2);
        CrashPrintf("\nFATAL crash: signal %d%s (thread '%s', build %s), stack:\n", sig, sig == SIGABRT && gAbortWhy ? " (abort)" : "", tName, gBuild[0] ? gBuild : "?");
        if (sig == SIGABRT && gAbortWhy) CrashPrintf("FATAL abort() after %s\n", gAbortWhy);
        if (gOut >= 0) backtrace_symbols_fd(frames, n, gOut);
        signal(sig, SIG_DFL);
        raise(sig);
    }
    void OnConsoleStop(int sig) { // (#1527) SIGINT / SIGTERM outside --server: one line, then end at once
        CrashFlush();
        CrashPrintf("%9.3f stopped from outside (%s): ending at once, nothing saved\n",
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - gStart).count(), sig == SIGINT ? "SIGINT" : "SIGTERM");
        _exit(128 + sig);
    }
#endif
#if RAYLIBGAME_LOGGING
    // (#1526) A failed assert (or IM_ASSERT): its message (expression, file, line) went to stderr only, which runs started
    // by scripts rarely keep. Now it goes into the log too, then abort() as before.
    [[noreturn]] void LogAssert(const char* expr, const char* file, unsigned line) {
        CrashFlush();
        CrashPrintf("\nFATAL assertion failed: %s, file %s, line %u (thread '%s')\n", expr, file, line, tName);
        fprintf(stderr, "Assertion failed: %s, file %s, line %u\n", expr, file, line);
        fflush(stderr);
        gAbortWhy = "the failed assertion above";
#ifdef _WIN32
        FatalBox("The game hit a fatal error (a failed assertion) and has to close. The details are in log.txt.");
#endif
        std::abort();
    }
#endif
}

#if RAYLIBGAME_LOGGING
#ifdef _WIN32
// MinGW's assert() calls _assert (or _wassert) through the import pointer __imp__assert. Defining that pointer here makes
// the linker take ours instead of msvcrt's import, for the game and every library linked into it (raylib, ImGui).
extern "C" {
static void __cdecl LogAssertA(const char* expr, const char* file, unsigned line) { LogAssert(expr, file, line); }
static void __cdecl LogAssertW(const wchar_t* expr, const wchar_t* file, unsigned line) {
    char e[512], f[260];
    snprintf(e, sizeof e, "%ls", expr); snprintf(f, sizeof f, "%ls", file);
    LogAssert(e, f, line);
}
void (__cdecl* __imp__assert)(const char*, const char*, unsigned) = LogAssertA;
void (__cdecl* __imp__wassert)(const wchar_t*, const wchar_t*, unsigned) = LogAssertW;
// The functions too (#2028): libmingw32's tls_atexit.o (pulled in by any thread_local with a destructor) calls _assert
// directly. Without a definition here the linker takes msvcrt's _assert member, which also defines __imp__assert and
// collides with the pointer above (#1802, #2015). With one, msvcrt's member is never pulled.
void __cdecl _assert(const char* expr, const char* file, unsigned line) { LogAssertA(expr, file, line); }
void __cdecl _wassert(const wchar_t* expr, const wchar_t* file, unsigned line) { LogAssertW(expr, file, line); }
}
#else
// glibc's assert() calls __assert_fail: the game's own definition comes first.
extern "C" [[noreturn]] void __assert_fail(const char* expr, const char* file, unsigned int line, const char* /*function*/) noexcept { LogAssert(expr, file, line); }
#endif
#endif

void Log_SetLevel(LogLevel l) { logdetail::gLevel.store((int)l, std::memory_order_relaxed); }
LogLevel Log_Level() { return (LogLevel)logdetail::gLevel.load(std::memory_order_relaxed); }
const char* Log_LevelName(LogLevel l) { return kNames[std::clamp((int)l, 0, 6)]; }

bool Log_ParseLevel(const char* name, LogLevel& out) {
    if (!name) return false;
    char s[16] = {};
    for (int i = 0; i < 15 && name[i]; ++i) s[i] = (char)((name[i] >= 'A' && name[i] <= 'Z') ? name[i] - 'A' + 'a' : name[i]);
    static const struct { const char* n; LogLevel l; } kParse[] = {
        { "trace", LogLevel::Trace }, { "debug", LogLevel::Debug }, { "info", LogLevel::Info }, { "warn", LogLevel::Warn },
        { "warning", LogLevel::Warn }, { "error", LogLevel::Error }, { "fatal", LogLevel::Fatal }, { "off", LogLevel::Off } };
    for (const auto& p : kParse) if (!strcmp(s, p.n)) { out = p.l; return true; }
    return false;
}

void Log_SetBuild(const char* build) {
    snprintf(gBuild, sizeof gBuild, "%s", build ? build : "");
#ifdef _WIN32
    const unsigned long pid = (unsigned long)GetCurrentProcessId();
#else
    const unsigned long pid = (unsigned long)getpid();
#endif
    LOG_INFO("Log", "build %s, pid %lu", gBuild, pid); // (which process wrote a log: #1527)
}
void Log_SetCrashDir(const char* dir) { snprintf(gCrashDir, sizeof gCrashDir, "%s", dir ? dir : "."); }

bool Log_Interactive() {
    if (const char* e = getenv("RAYLIBGAME_FATAL_BOX")) return e[0] == '1'; // 0: never a box, 1: always
#ifdef _WIN32
    // A console (a person ran it) or none at all (a double-click start): yes. A pipe or a file (scripts, agents, the
    // runner, timeouts): no one is there to press OK, and the box kept killed processes alive (#1527).
    const HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
    if (h == nullptr || h == INVALID_HANDLE_VALUE) return true;
    const DWORD t = GetFileType(h);
    return t == FILE_TYPE_CHAR || t == FILE_TYPE_UNKNOWN;
#else
    return isatty(1) != 0;
#endif
}

void Log_Open(const char* path) {
#if RAYLIBGAME_LOGGING
    std::lock_guard<std::mutex> lock(gMutex);
    if (LogOpen()) {
        FlushLocked();
#ifdef _WIN32
        CloseHandle(gOut); gOut = INVALID_HANDLE_VALUE;
#else
        close(gOut); gOut = -1;
#endif
    }
    { // (#1153) the last session's log is kept as <name>.prev<ext> (log.prev.txt): a failure, then a new start, kept nothing
        const std::string p = path;
        const size_t dot = p.find_last_of('.'), slash = p.find_last_of("/\\");
        const std::string prev = dot != std::string::npos && (slash == std::string::npos || dot > slash) ? p.substr(0, dot) + ".prev" + p.substr(dot) : p + ".prev";
        std::remove(prev.c_str());
        std::rename(path, prev.c_str()); // (works while an older process still has it open: it was opened with delete sharing)
    }
#ifdef _WIN32
    // Append-only (every write at the end), and shared for reading, writing and renaming (#1527: see the top of the file).
    gOut = CreateFileA(path, FILE_APPEND_DATA | SYNCHRONIZE, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, nullptr, CREATE_ALWAYS,
                       FILE_ATTRIBUTE_NORMAL, nullptr);
    ULONG stack = 64 * 1024; SetThreadStackGuarantee(&stack); // room for the crash handler after a stack overflow on this thread
    SetUnhandledExceptionFilter(FlushOnCrash);
    std::set_terminate(OnTerminate);
    signal(SIGABRT, OnAbort);
#else
    gOut = open(path, O_WRONLY | O_CREAT | O_TRUNC | O_APPEND | O_CLOEXEC, 0644);
    for (int sig : { SIGSEGV, SIGBUS, SIGFPE, SIGILL, SIGABRT }) signal(sig, FlushOnCrash);
#endif
#else
    (void)path;
#endif
}

void Log_InstallConsoleStop() {
    if (const char* e = getenv("RAYLIBGAME_CONSOLE_STOP")) if (e[0] == '0') return; // the old behaviour, for A/B
#if RAYLIBGAME_LOGGING
#ifdef _WIN32
    SetConsoleCtrlHandler(OnConsoleStop, TRUE);
#else
    signal(SIGINT, OnConsoleStop);
    signal(SIGTERM, OnConsoleStop);
#endif
#endif
}

void Log_SetThreadName(const char* name) { snprintf(tName, sizeof(tName), "%s", name); }

void Log_TestCrash(const char* kind) {
    const std::string k = kind ? kind : "";
    LOG_INFO("Log", "--test-crash %s", k.c_str());
    Log_Flush();
    if (k == "segv") { volatile int* p = nullptr; *p = 1; }
    else if (k == "thread-segv") { std::thread t([] { Log_SetThreadName("crasher"); volatile int* p = nullptr; *p = 1; }); t.join(); }
    else if (k == "assert") { // (assert() itself is compiled out of release builds: the CRT's entry, as assert() calls it)
#if RAYLIBGAME_LOGGING && defined(_WIN32)
        __imp__assert("--test-crash assert", __FILE__, __LINE__);
#elif RAYLIBGAME_LOGGING
        __assert_fail("--test-crash assert", __FILE__, __LINE__, __func__);
#endif
    }
    else if (k == "terminate") { std::thread t([] {}); } // destroyed while joinable: std::terminate
    else if (k == "fatal") LOG_FATAL("Log", "--test-crash fatal");
    else if (k == "wait") for (;;) std::this_thread::sleep_for(std::chrono::milliseconds(100));
    else if (k == "ctrl-c") { // (#1527) a thread still running (as the factory's), then Ctrl+C on this process after RAYLIBGAME_TEST_CRASH_DELAY s
        static std::thread busy([] { for (;;) std::this_thread::sleep_for(std::chrono::milliseconds(100)); });
        const char* d = getenv("RAYLIBGAME_TEST_CRASH_DELAY");
        std::this_thread::sleep_for(std::chrono::milliseconds((int)(1000.0 * (d ? atof(d) : 0.0))));
#ifdef _WIN32
        FreeConsole(); AllocConsole(); // a console of its own: the event reaches no other process (the shell that started it)
        SetConsoleCtrlHandler(nullptr, FALSE);
        Log_InstallConsoleStop(); // (a new console starts with no handlers: as a run started in a console has them)
        GenerateConsoleCtrlEvent(CTRL_C_EVENT, 0);
#else
        raise(SIGINT);
#endif
        for (;;) std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    LOG_WARN("Log", "--test-crash %s: unknown kind (segv, thread-segv, assert, terminate, fatal, wait)", k.c_str());
}

void Log_Flush() {
    std::lock_guard<std::mutex> lock(gMutex);
    FlushLocked();
}

void Log_Write(LogLevel l, const char* category, const char* fmt, ...) {
    char msg[2048];
    va_list args;
    va_start(args, fmt);
    vsnprintf(msg, sizeof(msg), fmt, args);
    va_end(args);
    size_t n = strlen(msg);
    while (n > 0 && (msg[n - 1] == '\n' || msg[n - 1] == '\r')) msg[--n] = '\0'; // the logger ends the line

    if (!tName[0]) snprintf(tName, sizeof(tName), "t%d", ++gThreadCount);
    const double run = std::chrono::duration<double>(std::chrono::steady_clock::now() - gStart).count();
    const auto wall = std::chrono::system_clock::now();
    const time_t secs = std::chrono::system_clock::to_time_t(wall);
    const int ms = (int)(std::chrono::duration_cast<std::chrono::milliseconds>(wall.time_since_epoch()).count() % 1000);
    struct tm tmv;
#ifdef _WIN32
    localtime_s(&tmv, &secs);
#else
    localtime_r(&secs, &tmv);
#endif
    char head[96];
    const int hn = snprintf(head, sizeof(head), "%9.3f %02d:%02d:%02d.%03d %-5s %-7s [%s] ", run, tmv.tm_hour, tmv.tm_min, tmv.tm_sec, ms,
                            Log_LevelName(l), tName, category ? category : "");
    const bool urgent = l >= LogLevel::Warn;
    {
        std::lock_guard<std::mutex> lock(gMutex);
        if (LogOpen()) {
            AppendLocked(head, (size_t)std::clamp(hn, 0, (int)sizeof(head) - 1));
            AppendLocked(msg, n);
            AppendLocked("\n", 1);
            if (urgent || run - gLastFlush > 1.0) { FlushLocked(); gLastFlush = run; }
        }
    }
    if (urgent) { fputs(head, stderr); fputs(msg, stderr); fputc('\n', stderr); }
}

void Log_FatalExit() {
    Log_Flush();
#ifdef _WIN32
    FatalBox("The game hit a fatal error and has to close. The details are in log.txt.");
#endif
#if RAYLIBGAME_LOGGING
    gAbortWhy = "the FATAL line above";
#endif
    std::abort();
}

void Log_StartupChecks() {
#ifdef _WIN32
    // A process started by a debugger gets the NT debug heap: every malloc/free is much slower, and chunk builds ran
    // ~26x slower in the playtests (#387, #509). _NO_DEBUG_HEAP=1 in the debugger's environment turns it off.
    const char* noDebugHeap = getenv("_NO_DEBUG_HEAP");
    if (IsDebuggerPresent() && !(noDebugHeap && noDebugHeap[0] == '1'))
        LOG_WARN("Log", "a debugger started this process with the Windows debug heap: allocation-heavy work (chunk builds) "
                        "runs many times slower; set _NO_DEBUG_HEAP=1 in the debugger's environment");
#endif
    LOG_INFO("Log", "log level %s (compiled in from %s)", Log_LevelName(Log_Level()), Log_LevelName((LogLevel)RAYLIBGAME_LOG_MIN));
}
