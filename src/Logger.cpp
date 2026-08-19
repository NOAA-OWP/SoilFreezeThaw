#include "Logger.hpp"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib> // For getenv()
#include <cstring>
#include <fstream> // For file handling
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string> // For std::string
#include <sys/stat.h>
#include <sys/wait.h>
#include <thread>
#include <unordered_map>
#include <cstdarg>
#include <cstdio>
#include <vector>

const std::string  MODULE_NAME           = "SFT";
const std::string  LOG_DIR_NGENCERF      = "/ngencerf/data";      // ngenCERF log directory string if environement var empty.
const std::string  LOG_DIR_DEFAULT       = "run-logs";            // Default parent log directory string if env var empty  & ngencerf dosn't exist
const std::string  LOG_FILE_EXT          = "log";                 // Log file name extension
const std::string  DS                    = "/";                   // Directory separator
const unsigned int LOG_MODULE_NAME_LEN   = 8;                     // Width of module name for log entries

const std::string  EV_EWTS_LOGGING       = "NGEN_EWTS_LOGGING";   // Enable/disable of Error Warning and Trapping System  
const std::string  EV_NGEN_LOGFILEPATH   = "NGEN_LOG_FILE_PATH";  // ngen log file 
const std::string  EV_MODULE_LOGLEVEL    = "SFT_LOGLEVEL";        // This modules log level
const std::string  EV_MODULE_LOGFILEPATH = "SFT_LOGFILEPATH";     // This modules log full log filename

bool         Logger::loggerInitialized = false;
std::string  Logger::logFilePath = "";
bool         Logger::loggingEnabled = true;
std::fstream Logger::logFile;
LogLevel     Logger::logLevel = LogLevel::INFO; // Default Log Level
std::string  Logger::moduleName = "";

std::shared_ptr<Logger> Logger::loggerInstance;

using namespace std;

// String to LogLevel map
static const std::unordered_map<std::string, LogLevel> logLevelMap = {
    {"NONE",    LogLevel::NONE },
    {"0",       LogLevel::NONE },
    {"DEBUG",   LogLevel::DEBUG},
    {"1",       LogLevel::DEBUG},
    {"INFO",    LogLevel::INFO },
    {"2",       LogLevel::INFO },
    {"WARNING", LogLevel::WARNING},
    {"3",       LogLevel::WARNING},
    {"SEVERE",  LogLevel::SEVERE},
    {"4",       LogLevel::SEVERE},
    {"FATAL",   LogLevel::FATAL},
    {"5",       LogLevel::FATAL},
};

// Reverse map: LogLevel to String
static const std::unordered_map<LogLevel, std::string> logLevelToStringMap = {
    {LogLevel::NONE,    "NONE   "},
    {LogLevel::DEBUG,   "DEBUG  "},
    {LogLevel::INFO,    "INFO   "},
    {LogLevel::WARNING, "WARNING"},
    {LogLevel::SEVERE,  "SEVERE "},
    {LogLevel::FATAL,   "FATAL  "},
};

std::string Logger::ToUpper(const std::string& input) {
    std::string result = input;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c){ return std::toupper(c); });
    return result;
}

// Function to trim leading and trailing spaces
std::string Logger::TrimString(const std::string& str) {
    // Trim leading spaces
    size_t first = str.find_first_not_of(" \t\n\r\f\v");
    if (first == std::string::npos) {
        return ""; // No non-whitespace characters
    }

    // Trim trailing spaces
    size_t last = str.find_last_not_of(" \t\n\r\f\v");

    // Return the trimmed string
    return str.substr(first, last - first + 1);
}

LogLevel Logger::ConvertStringToLogLevel(const std::string& levelStr) {
    std::string level = TrimString(levelStr);
    if (!level.empty()) {
        // Convert string to LogLevel (supports both names and numbers)
        auto it = logLevelMap.find(level);
        if (it != logLevelMap.end()) {
            return it->second; // Found valid named or numeric log level
        }

        // Try parsing as an integer (for cases where an invalid numeric value is given)
        try {
            int levelNum = std::stoi(level);
            if (levelNum >= 0 && levelNum <= 5) {
                return static_cast<LogLevel>(levelNum);
            }
        } catch (...) {
            // Ignore errors (e.g., if std::stoi fails for non-numeric input)
        }
    }
    return LogLevel::NONE;
}

/**
 * Convert LogLevel to String Representation of Log Level
 * @param logLevel : LogLevel
 * @return String log level
 */
std::string Logger::ConvertLogLevelToString(LogLevel level) {
    auto it = logLevelToStringMap.find(level);
    if (it != logLevelToStringMap.end()) {
        return it->second; // Found valid named or numeric log level
    }
    return "NONE";
}

bool Logger::DirectoryExists(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        return false; // Cannot access path
    }
    return (info.st_mode & S_IFDIR) != 0;
}

/**
 * Create the directory checking both the call
 * to execute the command and the result of the command
 */
bool Logger::CreateDirectory(const std::string& path) {

    if (!DirectoryExists(path)) {
        std::string mkdir_cmd = "mkdir -p " + path;
        int status            = system(mkdir_cmd.c_str());

        if (status == -1) {
            std::cerr << "system() failed to run mkdir.\n";
            return false;
        } else if (WIFEXITED(status)) {
            int exitCode = WEXITSTATUS(status);
            if (exitCode != 0) {
                std::cerr << "mkdir command failed with exit code: " << exitCode << "\n";
                return false;
            }
        } else {
            std::cerr << "mkdir terminated abnormally.\n";
            return false;
        }
    }
    return true;
}

/**
 * Open log file and return open status. If already open,
 * ensure the write pointer is at the end of the file.
 *
 * return bool true if open and good, false otherwise
 */
bool Logger::LogFileReady(bool appendMode) {

    if (logFile.is_open() && logFile.good()) {
        logFile.seekp(0, std::ios::end); // Ensure write pointer is at the actual file end
        return true;
    } else {
        // Attempt to open
        if (!logFilePath.empty()) {
            logFile.open(
                logFilePath,
                ios::out | ((appendMode) ? ios::app : ios::trunc)
            ); // This will silently fail if already open.
            if (logFile.good()) {
                return true;
            }
        }
    }
    return false;
}

/**
 * Set the log file path name using the following pattern
 *  - Use the module log file if available (unset when first run by ngen), otherwise 
 *  - Use ngen log file if available, otherwise
 *  - Use /ngencerf/data/run-logs/<username>/<module>_<YYMMDDTHHMMSS> if available, otherwise
 *  - Use ~/run-logs/<YYYYMMDD>/<module>_<YYMMDDTHHMMSS>
 *  - Onced opened, save the full log path to the modules log environment variable so
 *    it is only opened once for each ngen run (vs for each catchment)
 */
void Logger::SetLogFilePath(void) {
    // Use module logfile path environment variable if it exists
    logFilePath = "";
    bool appendEntries = true;
    bool mdduleLogEnvExists = false;
    const char* envVar = std::getenv(EV_MODULE_LOGFILEPATH.c_str()); // Set once module has successfully opened a log file
    if (envVar != nullptr && envVar[0] != '\0') {
        logFilePath = envVar;
        mdduleLogEnvExists = true;
    } else {
        // log path not set yet.
        envVar = std::getenv(EV_NGEN_LOGFILEPATH.c_str()); // Currently set by ngen-cal but envision set for WCOSS at some point
        if (envVar != nullptr && envVar[0] != '\0') {
            logFilePath = envVar;
        } else {
            // ngen log file not defined. Create alternate log file
            appendEntries = false;
            std::string logFileDir;
            if (DirectoryExists(LOG_DIR_NGENCERF)) {
                logFileDir = LOG_DIR_NGENCERF + DS + LOG_DIR_DEFAULT;
            } else {
                const char* envHomeDir = std::getenv("HOME");
                std::string homeDir = (envHomeDir != nullptr && envHomeDir[0] != '\0') ? envHomeDir : "~";
                logFileDir = homeDir + DS + LOG_DIR_DEFAULT;
            }

            // Ensure parent log directory exists
            if (CreateDirectory(logFileDir)) {
                // Set full log directory path
                const char* envUsername = std::getenv("USER");
                if (envUsername != nullptr && envUsername[0] != '\0') {
                    std::string username = envUsername;
                    logFileDir = logFileDir + DS + username;
                } else {
                    logFileDir = logFileDir + DS + CreateDateString();
                }
                // Set the full path if log directory exists/created
                if (CreateDirectory(logFileDir))
                    logFilePath = logFileDir + DS + MODULE_NAME + "_" + CreateTimestamp(false, false) +
                                "." + LOG_FILE_EXT;
            }
        }
    }

    // Open log file. If path not found, it will try an alternate file. If neither
    // works false returned and logs will be written to stdout.
    if (LogFileReady(appendEntries)) {
        if (!mdduleLogEnvExists) { 
            setenv(EV_MODULE_LOGFILEPATH.c_str(), logFilePath.c_str(), 1);
            std::cout << "Module " << MODULE_NAME << " Log File: " << logFilePath << std::endl;
            LogLevel saveLevel = logLevel;
            logLevel = LogLevel::INFO; // Ensure this INFO message is always logged
            Log(logLevel, "Logging started. Log File Path: %s\n", logFilePath.c_str());
            logLevel = saveLevel;
        }
    } else {
        std::cout << "Unable to open log file ";
        if (!logFilePath.empty()) {
            std::cout << logFilePath;
            std::cout << " (Perhaps check permissions)" << std::endl;
        }
        std::cout << "Log entries will be writen to stdout" << std::endl;
    }
}

void Logger::SetLogLevel(void) {
    // Set the logger log level if environment var not found
    const char* envLogLevel = std::getenv(EV_MODULE_LOGLEVEL.c_str());
    if (envLogLevel != nullptr && envLogLevel[0] != '\0') {
        logLevel = ConvertStringToLogLevel(envLogLevel);
    }
    std::string llMsg = "Log level set to " + ConvertLogLevelToString(logLevel) + "\n";
    cout << MODULE_NAME << " " << llMsg;
    LogLevel saveLevel = logLevel;
    logLevel = LogLevel::INFO; // Ensure this INFO message is always logged
    Log(logLevel, llMsg);
    logLevel = saveLevel;

}

void Logger::SetLoggingFlag(void) {
    const char* envVar = std::getenv(EV_EWTS_LOGGING.c_str()); // Set once module has successfully opened a log file
    if (envVar != nullptr && envVar[0] != '\0') {
        std::string logState = ToUpper(TrimString(envVar));
        loggingEnabled = (logState == "ENABLED")?true:false;
    }
    std::cout << MODULE_NAME << " Logging " << ((loggingEnabled)?"ENABLED":"DISABLED") << std::endl;
}

bool Logger::IsLoggingEnabled(void) {
    return loggingEnabled;
}

void Logger::SetLogModuleName(void) {
    // Make sure the module name used for logging entries is all uppercase and 8 characters wide.
    moduleName            = MODULE_NAME;
    std::string upperName = moduleName.substr(0, LOG_MODULE_NAME_LEN); // Truncate to LOG_MODULE_NAME_LEN chars max
    std::transform(upperName.begin(), upperName.end(), upperName.begin(), ::toupper);

    std::ostringstream oss;
    oss << std::left << std::setw(8) << std::setfill(' ') << upperName;
    moduleName = oss.str();
}

/**
 * Configure Logger Preferences
 * @param logFile
 * @param level: LogLevel::INFO by Default
 * @return void
 */
void Logger::SetLogPreferences(void) {

    if (!loggerInitialized) {
        loggerInitialized = true; // Only call this once

        SetLoggingFlag();
        if (loggingEnabled) {
            SetLogModuleName();
            SetLogFilePath();
            SetLogLevel();
        }
    }
}

void Logger::Log(LogLevel messageLevel, const char* message, ...) {
    va_list args;
    va_start(args, message);

    // Make a copy to calculate required size
    va_list args_copy;
    va_copy(args_copy, args);
    int requiredLen = vsnprintf(nullptr, 0, message, args_copy);
    va_end(args_copy);

    if (requiredLen > 0) {
        std::vector<char> buffer(requiredLen + 1);  // +1 for null terminator
        vsnprintf(buffer.data(), buffer.size(), message, args);

        va_end(args);

        Log(std::string(buffer.data()), messageLevel);
    } else {
        va_end(args);  // still need to clean up
    }
}

/**
 * Log given message with defined parameters and generate message to pass on Console or File
 * @param message: Log Message
 * @param messageLevel: Log Level, LogLevel::INFO by default
 */
void Logger::Log(LogLevel messageLevel, std::string message) {
    Log(message, messageLevel);
}
void Logger::Log(std::string message, LogLevel messageLevel) {

    if (!loggerInitialized) SetLogPreferences(); // Cover case where Log is called before setup done

    // don't log if messageLevel < logLevel
    if (loggingEnabled && (messageLevel >= logLevel)) {
        std::string logType   = ConvertLogLevelToString(messageLevel);
        std::string logPrefix = CreateTimestamp() + " " + moduleName + " " + logType;

        // Log message, creating individual entries for a multi-line message
        std::istringstream logMsg(message);
        std::string line;
        if (LogFileReady()) {
            while (std::getline(logMsg, line)) {
                logFile << logPrefix + " " + line << std::endl;
            }
            logFile.flush();
        } else {
            // Log file not found. Write to stdout.
            while (std::getline(logMsg, line)) {
                cout << logPrefix + " " + line << std::endl;
            }
            cout << std::flush;
        }
    }
}

std::string Logger::CreateDateString(void) {
    std::time_t tt    = std::time(0);
    std::tm* timeinfo = std::gmtime(&tt); // Use std::localtime(&tt) if you want local time

    char buffer[11]; // Enough for "YYYY-MM-DD" + null terminator
    std::strftime(buffer, sizeof(buffer), "%F", timeinfo); // %F == %Y-%m-%d

    std::stringstream ss;
    ss << buffer;

    return ss.str();
}

std::string Logger::CreateTimestamp(bool appendMS, bool iso) {
    using namespace std::chrono;

    // Get current time point
    auto now        = system_clock::now();
    auto now_time_t = system_clock::to_time_t(now);

    // Get milliseconds
    auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;

    // Convert to UTC time
    std::tm utc_tm;
    gmtime_r(&now_time_t, &utc_tm);

    // Format date/time with strftime
    char buffer[32];
    if (iso) {
        std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S", &utc_tm);
    } else {
        std::strftime(buffer, sizeof(buffer), "%Y%m%dT%H%M%S", &utc_tm);
    }

    if (appendMS) {
        // Combine with milliseconds
        std::ostringstream oss;
        oss << buffer << '.' << std::setw(3) << std::setfill('0') << ms.count();
        return oss.str();
    }
    return std::string(buffer);
}

std::string Logger::GetLogFilePath() {
    return logFilePath;
}
