#include "logger.h"
#include <iostream>
#include <cstdarg>
#include <cstring>

Logger::Logger(const std::string& filename) {
    // Open file in trunc mode to clear it on initialization
    logFile.open(filename, std::ios::out | std::ios::trunc);
}

// Flushing every line costs a syscall per message; only force it for errors/warnings so they
// survive a crash. Everything else is flushed by the stream buffer / on exit.
static bool isUrgent(const char* s) {
    return strstr(s, "ERROR") != nullptr || strstr(s, "WARN") != nullptr;
}

void Logger::log(const std::string& message) {
    std::lock_guard<std::mutex> lock(logMutex);
    if (!logFile.is_open()) {
        std::cerr << "Logger: Failed to write to log file. File is not open." << std::endl;
        return;
    }
    logFile << message << '\n';
    if (isUrgent(message.c_str())) logFile.flush();
}

void Logger::logf(const char* fmt, ...) {
    std::lock_guard<std::mutex> lock(logMutex);
    if (!logFile.is_open()) {
        std::cerr << "Logger: Failed to write to log file. File is not open." << std::endl;
        return;
    }
    constexpr size_t BUF_SIZE = 1024;
    char buf[BUF_SIZE];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buf, BUF_SIZE, fmt, args);
    va_end(args);
    logFile << buf;
    if (isUrgent(buf)) logFile.flush();
}
