#include "logger.h"
#include <algorithm>
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <mutex>
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace {
    std::mutex gMutex;
    FILE* gFile = nullptr;
    const auto gStart = std::chrono::steady_clock::now();
    std::atomic<int> gThreadCount{ 0 };
    thread_local char tName[16] = "";
    const char* const kNames[] = { "TRACE", "DEBUG", "INFO", "WARN", "ERROR", "FATAL", "OFF" };
    double gLastFlush = 0.0; // run time of the last flush (gMutex): lines reach the file within a second
#if defined(_WIN32) && RAYLIBGAME_LOGGING
    // A crash: what is buffered goes out before Windows ends the process (no lock: the crashing thread may hold it).
    LONG WINAPI FlushOnCrash(EXCEPTION_POINTERS*) {
        if (gFile) fflush(gFile);
        return EXCEPTION_CONTINUE_SEARCH;
    }
#endif
}

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

void Log_Open(const char* path) {
#if RAYLIBGAME_LOGGING
    std::lock_guard<std::mutex> lock(gMutex);
    if (gFile) fclose(gFile);
    gFile = fopen(path, "w");
    if (gFile) setvbuf(gFile, nullptr, _IOFBF, 1 << 16);
#ifdef _WIN32
    SetUnhandledExceptionFilter(FlushOnCrash);
#endif
#else
    (void)path;
#endif
}

void Log_SetThreadName(const char* name) { snprintf(tName, sizeof(tName), "%s", name); }

void Log_Flush() {
    std::lock_guard<std::mutex> lock(gMutex);
    if (gFile) fflush(gFile);
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
    snprintf(head, sizeof(head), "%9.3f %02d:%02d:%02d.%03d %-5s %-7s [%s] ", run, tmv.tm_hour, tmv.tm_min, tmv.tm_sec, ms,
             Log_LevelName(l), tName, category ? category : "");
    const bool urgent = l >= LogLevel::Warn;
    std::lock_guard<std::mutex> lock(gMutex);
    if (gFile) {
        fputs(head, gFile); fputs(msg, gFile); fputc('\n', gFile);
        if (urgent || run - gLastFlush > 1.0) { fflush(gFile); gLastFlush = run; }
    }
    if (urgent) { fputs(head, stderr); fputs(msg, stderr); fputc('\n', stderr); }
}

void Log_FatalExit() {
    Log_Flush();
#ifdef _WIN32
    MessageBoxA(nullptr, "The game hit a fatal error and has to close. The details are in log.txt.", "RaylibGame", MB_OK | MB_ICONERROR);
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
