//
// Created by kellerberrin on 18/06/24.
//

#include "kel_logging_stream.h"
#include "kel_utility.h"
#include <chrono>
#include <iostream>


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Logging implementation uses the SPDLog library
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kel = kellerberrin;


kel::StreamLoggerImpl::StreamLoggerImpl(const std::string& module, const std::string& log_file) {

  padded_module_name_ = "[" + module + "] ";
  output_logfile_.open(log_file);
  if (not output_logfile_.good()) {

    std::cout << ANSI_RED << "Unable to open application log file: " << log_file << ANSI_RESET
              << ". Logging to console only." << std::endl;

  }

}


std::string kel::StreamLoggerImpl::severityANSI(ExecEnvLogger::LoggerSeverity severity, bool is_ansi) const {

  switch(severity) {

    case ExecEnvLogger::LoggerSeverity::INFO:
      return is_ansi ? ANSI_INFO : INFO;

    case ExecEnvLogger::LoggerSeverity::WARN:
      return is_ansi ? ANSI_WARN : WARN;

    case ExecEnvLogger::LoggerSeverity::ERROR:
      return is_ansi ? ANSI_ERROR : ERROR;

    default:
    case ExecEnvLogger::LoggerSeverity::CRITICAL:
      return is_ansi ? ANSI_CRITICAL : CRITICAL;

  }

}


void kel::StreamLoggerImpl::formatImpl(std::string&& formatted_string, ExecEnvLogger::LoggerSeverity severity) noexcept {

  std::string date_time = dateTime();
  std::cout << date_time <<  padded_module_name_ << severityANSI(severity, true) << formatted_string << '\n';
  if (output_logfile_.good()) {

    output_logfile_ << date_time <<  padded_module_name_ << severityANSI(severity, false) << formatted_string << '\n';

  }

}

void kel::StreamLoggerImpl::locationImpl(const LogFormatLocation& format_location,
                                      const std::string& formatted_string,
                                      ExecEnvLogger::LoggerSeverity severity) noexcept {

  std::string date_time = dateTime();
  std::string location_text = location(format_location);
  std::cout << date_time <<  padded_module_name_ << severityANSI(severity, true) << location_text << formatted_string << '\n';
  if (output_logfile_.good()) {

    output_logfile_ << date_time <<  padded_module_name_ << severityANSI(severity, false) << location_text << formatted_string << '\n';

  }

}


std::string kel::StreamLoggerImpl::dateTime() const
{

  const std::chrono::zoned_time cur_time{ std::chrono::current_zone(), std::chrono::system_clock::now() };
  return std::format("[{:%Y-%m-%d %X}] ", cur_time);

}

std::string kel::StreamLoggerImpl::location(const LogFormatLocation &format_location) const {

  std::string file_name = Utility::fileName(format_location.location().file_name());
  std::string location_text = std::string("[") + file_name + ":"
                            + std::to_string(format_location.location().line()) + "] ";
  return location_text;

}

