#include "logger.h"
#include <iostream>
#include <cstdarg>

Logger::Logger(const std::string& filename) {
    // Open file in trunc mode to clear it on initialization
    logFile.open(filename, std::ios::out | std::ios::trunc);
}

void Logger::log(const std::string& message) {
    std::lock_guard<std::mutex> lock(logMutex);
    if (logFile.is_open()) {
        logFile << message << std::endl;
        logFile.flush();
    } else {
        // Handle error if file is not open
        std::cerr << "Logger: Failed to write to log file. File is not open." << std::endl;
    }
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
    logFile.flush();
}