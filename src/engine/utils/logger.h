#pragma once
// Logging with levels (#514). Lines go to log.txt as
//   <run time s> <wall clock> <LEVEL> <thread> [Category] message
// e.g. "   12.480 18:02:03.117 INFO  main    [World] loaded saves/world_18.sav (seed 1633406125)".
//
//   LOG_TRACE / LOG_DEBUG / LOG_INFO / LOG_WARN / LOG_ERROR ("Category", "printf format", args...)
//   LOG_FATAL(...)  logs, flushes, shows a message box (Windows) and aborts.
//
// Runtime level: Log_SetLevel (--log-level, settings.txt "log_level"); default INFO in release builds (NDEBUG), DEBUG
// otherwise. Compile time (CMake): RAYLIBGAME_LOGGING=OFF turns every macro into nothing - the arguments are not
// evaluated and no format string reaches the binary - and no log.txt is written; RAYLIBGAME_LOG_MIN strips the levels
// below it the same way. WARN and above are flushed at once (and echoed to stderr) so they survive a crash.
// Thread-safe: worker threads log too (Log_SetThreadName names a thread in its lines).
// Benches and tests keep printing their "bench:" lines to stdout: scripts parse those, not log.txt.

#ifndef RAYLIBGAME_LOGGING
#define RAYLIBGAME_LOGGING 1
#endif
#ifndef RAYLIBGAME_LOG_MIN
#define RAYLIBGAME_LOG_MIN 0 // the lowest level compiled in: 0 TRACE ... 5 FATAL
#endif

#include <atomic>

enum class LogLevel : int { Trace = 0, Debug, Info, Warn, Error, Fatal, Off };

namespace logdetail {
    inline std::atomic<int> gLevel{
#ifdef NDEBUG
        (int)LogLevel::Info
#else
        (int)LogLevel::Debug
#endif
    };
    template <class... A> int Unused(const A&...); // never defined: only named inside sizeof (nothing is evaluated)
}

inline bool Log_Enabled(LogLevel l) { return (int)l >= logdetail::gLevel.load(std::memory_order_relaxed); }
void Log_SetLevel(LogLevel l);
LogLevel Log_Level();
bool Log_ParseLevel(const char* name, LogLevel& out); // trace, debug, info, warn(ing), error, fatal, off
const char* Log_LevelName(LogLevel l);
void Log_Open(const char* path);  // truncates; nothing is written before (or without) it
// (#1526) Crash reports: the build in every crash entry, and where a minidump goes (Windows; default saves/logs).
void Log_SetBuild(const char* build);
void Log_SetCrashDir(const char* dir);
// (#1527) A console Ctrl+C / Ctrl+Break / close (a timeout or a kill from a script) ends the process at once with a
// log line: no fatal box, no std::terminate from threads still running. A --server installs its own handler on top.
void Log_InstallConsoleStop();
bool Log_Interactive();           // a person may be watching: stdout is a console or there is none (a double-click start)
// --test-crash KIND (#1526): segv | thread-segv | assert | terminate | fatal | wait (sleeps until killed) | ctrl-c (#1527: a
// running thread, then Ctrl+C on itself after RAYLIBGAME_TEST_CRASH_DELAY s)
void Log_TestCrash(const char* kind);
void Log_SetThreadName(const char* name); // this thread's name in its lines (default t1, t2, ...)
void Log_Flush();                 // a bug report reads log.txt mid-run (#169)
void Log_StartupChecks();         // warnings every log should carry (e.g. the debugger's debug heap)
void Log_Write(LogLevel l, const char* category, const char* fmt, ...)
#if defined(__GNUC__)
    __attribute__((format(gnu_printf, 3, 4))) // gnu_printf: MinGW's plain printf check rejects %zu, which the UCRT handles
#endif
    ;
[[noreturn]] void Log_FatalExit(); // flush, message box, abort (also when logging is compiled out)

#define LOG_DISCARD_(...) ((void)sizeof(logdetail::Unused(__VA_ARGS__)))
#if RAYLIBGAME_LOGGING
#define LOG_AT_(lvl, ...) do { if (Log_Enabled(lvl)) Log_Write(lvl, __VA_ARGS__); } while (0)
#else
#define LOG_AT_(lvl, ...) LOG_DISCARD_(__VA_ARGS__)
#endif

#if RAYLIBGAME_LOG_MIN <= 0
#define LOG_TRACE(...) LOG_AT_(LogLevel::Trace, __VA_ARGS__)
#else
#define LOG_TRACE(...) LOG_DISCARD_(__VA_ARGS__)
#endif
#if RAYLIBGAME_LOG_MIN <= 1
#define LOG_DEBUG(...) LOG_AT_(LogLevel::Debug, __VA_ARGS__)
#else
#define LOG_DEBUG(...) LOG_DISCARD_(__VA_ARGS__)
#endif
#if RAYLIBGAME_LOG_MIN <= 2
#define LOG_INFO(...) LOG_AT_(LogLevel::Info, __VA_ARGS__)
#else
#define LOG_INFO(...) LOG_DISCARD_(__VA_ARGS__)
#endif
#if RAYLIBGAME_LOG_MIN <= 3
#define LOG_WARN(...) LOG_AT_(LogLevel::Warn, __VA_ARGS__)
#else
#define LOG_WARN(...) LOG_DISCARD_(__VA_ARGS__)
#endif
#if RAYLIBGAME_LOG_MIN <= 4
#define LOG_ERROR(...) LOG_AT_(LogLevel::Error, __VA_ARGS__)
#else
#define LOG_ERROR(...) LOG_DISCARD_(__VA_ARGS__)
#endif
#if RAYLIBGAME_LOGGING
#define LOG_FATAL(...) do { Log_Write(LogLevel::Fatal, __VA_ARGS__); Log_FatalExit(); } while (0)
#else
#define LOG_FATAL(...) do { LOG_DISCARD_(__VA_ARGS__); Log_FatalExit(); } while (0)
#endif
