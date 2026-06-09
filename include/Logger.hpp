#ifndef SFT_LOGGER_HPP
#define SFT_LOGGER_HPP

#include "ewts/module_constants.hpp"
#include "ewts/logger.hpp"
#include "ewts/log_levels.hpp"

#define LOG(...) ::ewts::GetLogger(::ewts::modules::EWTS_ID_SFT).Log(__VA_ARGS__)
#define GetLogLevel() ::ewts::GetLogger(::ewts::modules::EWTS_ID_SFT).GetLogLevel()
#define IsLoggingEnabled() ::ewts::GetLogger(::ewts::modules::EWTS_ID_SFT).IsLoggingEnabled()

using ewts::EwtsInit;
using ewts::LogLevel;

inline constexpr const char* SFT_MODULE_ID = ewts::modules::EWTS_ID_SFT;

#endif /* SFT_LOGGER_HPP */
