//
// Created by kellerberrin on 18/06/24.
//

#include "kel_logging_stream.h"
#include "kel_utility.h"

#include <chrono>
#include <format>
#include <iostream>
#include <syncstream>


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Logging implementation.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kel = kellerberrin;


kel::StreamLoggerImpl::StreamLoggerImpl(const std::string& module, const std::string& log_file) {

  padded_module_name_ = "[" + module + "] ";
  output_logfile_.open(log_file);
  if (not output_logfile_.good()) {

    std::osyncstream(std::cout) << ANSI_RED << "Unable to open application log file: " << log_file << ANSI_RESET
                  << ". Logging to console only." << std::endl;

  }

}


std::string kel::StreamLoggerImpl::severityANSI(ExecEnvLogger::LoggerSeverity severity, bool is_ansi) {

  auto plain = [](std::string_view sv) { return std::string(sv); };
  auto color = [](std::string_view color_sv, std::string_view text_sv) {
    return std::string(color_sv) + std::string(text_sv) + std::string(ANSI_RESET);
  };

  switch(severity) {

    case ExecEnvLogger::LoggerSeverity::INFO:
      return plain(INFO);

    case ExecEnvLogger::LoggerSeverity::WARN:
      return is_ansi ? color(ANSI_YELLOW, WARN) : plain(WARN);

    case ExecEnvLogger::LoggerSeverity::ERROR:
      return is_ansi ? color(ANSI_RED, ERROR) : plain(ERROR);

    default:
    case ExecEnvLogger::LoggerSeverity::CRITICAL:
      return is_ansi ? color(ANSI_BACKGROUND_RED, CRITICAL) : plain(CRITICAL);

  }

}


void kel::StreamLoggerImpl::formatImpl(const std::string& formatted_string, ExecEnvLogger::LoggerSeverity severity) noexcept {

  writeMessage(dateTime(), severity, formatted_string);

}

void kel::StreamLoggerImpl::locationImpl(const LogFormatLocation& format_location,
                                      const std::string& formatted_string,
                                      ExecEnvLogger::LoggerSeverity severity) noexcept {

  writeMessage(dateTime(), severity, location(format_location) + formatted_string);

}


void kel::StreamLoggerImpl::writeMessage(const std::string& date_time,
                                         ExecEnvLogger::LoggerSeverity severity,
                                         const std::string& body) noexcept {

  try {

    std::osyncstream synced_cout(std::cout);
    synced_cout << date_time << padded_module_name_ << severityANSI(severity, true) << body << '\n';

  } catch (...) {
    // Console output failure is not fatal; the log file is the authoritative record.
  }

  if (output_logfile_.good()) {

    output_logfile_ << date_time << padded_module_name_ << severityANSI(severity, false) << body << '\n';

  }

}


std::string kel::StreamLoggerImpl::dateTime()
{

  try {

    const std::chrono::zoned_time cur_time{ std::chrono::current_zone(), std::chrono::system_clock::now() };
    return std::format("[{:%Y-%m-%d %X}] ", cur_time);

  }
  catch (const std::exception& e) {

    // Fallback to UTC if the local timezone database is unavailable.
    auto now = std::chrono::system_clock::now();
    return std::format("[{:%Y-%m-%d %X}] ", std::chrono::zoned_time{"UTC", now});

  }

}

std::string kel::StreamLoggerImpl::location(const LogFormatLocation &format_location) {

  std::string file_name = Utility::fileName(format_location.location().file_name());
  std::string location_text = std::string("[") + file_name + ":"
                            + std::to_string(format_location.location().line()) + "] ";
  return location_text;

}

// Finalizes logging. flushes and closes any open files. Called by logging a critical msg
void kel::StreamLoggerImpl::flushStreamsImpl() noexcept {

  output_logfile_.close();

}
