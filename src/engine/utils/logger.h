#pragma once
#include <string>
#include <fstream>
#include <mutex>
#include <cstdarg>

class Logger {
public:
    // Initialize logger with a file path. Clears the log file on init.
    Logger(const std::string& filename);

    // Write a line to the log file (thread-safe).
    void log(const std::string& message);

    // printf-style logging
    void logf(const char* fmt, ...);

private:
    std::ofstream logFile;
    std::mutex logMutex;
};

extern Logger logger;