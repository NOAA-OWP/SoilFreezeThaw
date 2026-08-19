#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdlib.h>
#include <string>

#define LOG Logger::Log

enum class LogLevel {
    NONE    = 0,
    DEBUG   = 1,
    INFO    = 2,
    WARNING = 3,
    SEVERE  = 4,
    FATAL   = 5,
};

/**
 * Logger Class Used to Output Details of Current Application Flow
   All methods and variables are static so instantiating an object is unnecessary.
 */
class Logger {
  public:
    // Methods
    static void Log(LogLevel messageLevel, const char* message, ...);
    static void Log(LogLevel messageLevel, std::string message);
    static void Log(std::string message, LogLevel messageLevel = LogLevel::INFO);
    static bool IsLoggingEnabled(void);
    static LogLevel GetLogLevel(void);

  private:
    // Methods
    static std::string CreateDateString(void);
    static std::string CreateTimestamp(bool appendMS = true, bool iso = true);
    static bool CreateDirectory(const std::string& path);
    static std::string ConvertLogLevelToString(LogLevel level);
    static LogLevel ConvertStringToLogLevel(const std::string& logLevel);
    static bool DirectoryExists(const std::string& path);
    static std::string GetLogFilePath(void);
    static bool LogFileReady(bool appendMode=true);
    static void SetLogFilePath(void);
    static void SetLoggingFlag(void);
    static void SetLogLevel(void);
    static void SetLogModuleName(void);
    static void SetLogPreferences(void);
    static std::string ToUpper(const std::string& input);
    static std::string TrimString(const std::string& str);

    // Variables
    // Declaring these static so they exist in the class without needing to instantiate it. 
    static bool         loggerInitialized;
    static std::string  logFilePath;
    static bool         loggingEnabled;
    static std::fstream logFile;
    static LogLevel     logLevel;
    static std::string  moduleName;

    static std::shared_ptr<Logger> loggerInstance;
};

#endif
