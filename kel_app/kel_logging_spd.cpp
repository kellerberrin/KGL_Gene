
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Logging implementation uses the SPDLog library
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "kel_logging_spd.h"


namespace kel = kellerberrin;


kel::SPDLoggerImpl::SPDLoggerImpl(const std::string& module, const std::string& log_file) {

  spdlog::set_pattern(SPDLOG_DEFAULT_FORMAT);
  spdlog::set_level(spdlog::level::trace);
  std::vector<spdlog::sink_ptr> sinks;
  sinks.push_back(std::make_shared<spdlog::sinks::ansicolor_stderr_sink_mt>());
  sinks.push_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file));
  log_impl_ptr_ = std::make_unique<spdlog::logger>(module, sinks.begin(), sinks.end());

}


spdlog::source_loc kel::SPDLoggerImpl::spdlogLocation(const std::source_location &location) {

  return spdlog::source_loc{location.file_name(),
                            static_cast<std::int32_t>(location.line()),
                            location.function_name()};

}


spdlog::level::level_enum kel::SPDLoggerImpl::spdlogSeverity(ExecEnvLogger::LoggerSeverity severity) {

  switch(severity) {

    case ExecEnvLogger::LoggerSeverity::INFO:
      return spdlog::level::info;

    case ExecEnvLogger::LoggerSeverity::WARN:
      return spdlog::level::warn;

    case ExecEnvLogger::LoggerSeverity::ERROR:
      return spdlog::level::err;

    default:
    case ExecEnvLogger::LoggerSeverity::CRITICAL:
      return spdlog::level::critical;


  }

}


void kel::SPDLoggerImpl::formatImpl(std::string&& formatted_string, ExecEnvLogger::LoggerSeverity severity) noexcept {

  log_impl_ptr_->log(spdlogSeverity(severity), formatted_string);
  log_impl_ptr_->flush();

}

void kel::SPDLoggerImpl::locationImpl(const LogFormatLocation& format_location,
                                      const std::string& formatted_string,
                                      ExecEnvLogger::LoggerSeverity severity) noexcept {

  log_impl_ptr_->log(spdlogLocation(format_location.location()), spdlogSeverity(severity), fmt::runtime(formatted_string));
  log_impl_ptr_->flush();

}

