//
// Created by kellerberrin on 20/10/23.
//

// Copyright 2023 Kellerberrin
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
// documentation files (the "Software"), to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
// and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
// OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
//


#include <spdlog/spdlog.h>  // Implement the logger using the spdlog library
#include "spdlog/sinks/basic_file_sink.h"

#include "kel_logging.h"

#include <iostream>


namespace kel = kellerberrin;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



namespace kellerberrin {   //  organization level namespace


class ExecEnvLoggerImpl {

public:

  ExecEnvLoggerImpl(const std::string &module, const std::string &log_file);
  ~ExecEnvLoggerImpl() = default;

  ExecEnvLoggerImpl(const ExecEnvLoggerImpl &) = delete;
  ExecEnvLoggerImpl(ExecEnvLoggerImpl &&) = delete;
  ExecEnvLoggerImpl &operator=(const ExecEnvLoggerImpl &) = delete;

  // These functions simply re-direct to the PIMPL implementation object.
  // These functions simply re-direct to the PIMPL implementation object.
  void formatImpl(std::string&& formatted_string, ExecEnvLogger::LoggerSeverity severity) noexcept;
  void locationImpl(const LogFormatLocation& format_location,
                    const std::string& formatted_string,
                    ExecEnvLogger::LoggerSeverity severity) noexcept;

private:

  std::unique_ptr<spdlog::logger> log_impl_ptr_;

  static constexpr const char *SPDLOG_DEFAULT_FORMAT{"%+"};  // Default format from spdlog

  [[nodiscard]] static spdlog::source_loc spdlogLocation(const std::source_location &location);
  [[nodiscard]] static spdlog::level::level_enum spdlogSeverity(ExecEnvLogger::LoggerSeverity severity);


};


} // Namespace.


kel::ExecEnvLoggerImpl::ExecEnvLoggerImpl(const std::string& module, const std::string& log_file) {

  spdlog::set_pattern(SPDLOG_DEFAULT_FORMAT);
  spdlog::set_level(spdlog::level::trace);
  std::vector<spdlog::sink_ptr> sinks;
  sinks.push_back(std::make_shared<spdlog::sinks::ansicolor_stderr_sink_mt>());
  sinks.push_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file));
  log_impl_ptr_ = std::make_unique<spdlog::logger>(module, sinks.begin(), sinks.end());

}


spdlog::source_loc kel::ExecEnvLoggerImpl::spdlogLocation(const std::source_location &location) {

  return spdlog::source_loc{location.file_name(),
                            static_cast<std::int32_t>(location.line()),
                            location.function_name()};

}


spdlog::level::level_enum kel::ExecEnvLoggerImpl::spdlogSeverity(ExecEnvLogger::LoggerSeverity severity) {

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


void kel::ExecEnvLoggerImpl::formatImpl(std::string&& formatted_string, ExecEnvLogger::LoggerSeverity severity) noexcept {

  log_impl_ptr_->log(spdlogSeverity(severity), formatted_string);
  log_impl_ptr_->flush();

}

void kel::ExecEnvLoggerImpl::locationImpl(const LogFormatLocation& format_location,
                                          const std::string& formatted_string,
                                          ExecEnvLogger::LoggerSeverity severity) noexcept {

  log_impl_ptr_->log(spdlogLocation(format_location.location()), spdlogSeverity(severity), fmt::runtime(formatted_string));
  log_impl_ptr_->flush();

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kel::ExecEnvLogger::ExecEnvLogger(const std::string& module, const std::string& log_file) {

  log_impl_ptr_ = std::make_unique<ExecEnvLoggerImpl>(module, log_file);

}

kel::ExecEnvLogger::~ExecEnvLogger() {

  log_impl_ptr_ = nullptr;

}

void kel::ExecEnvLogger::formatImpl(std::string&& formatted_string, LoggerSeverity severity) noexcept {

  log_impl_ptr_->formatImpl(std::forward<std::string&&>(formatted_string), severity);

}

void kel::ExecEnvLogger::locationImpl(const LogFormatLocation& format_location,
                                      const std::string& formatted_string,
                                      LoggerSeverity severity) noexcept {

  log_impl_ptr_->locationImpl(format_location, formatted_string, severity);

}

bool kel::ExecEnvLogger::warnMessageLimits() {

  std::lock_guard<std::mutex> lock(limit_mutex_);

  ++warn_message_count_;
  if (warn_message_count_ == max_warn_messages_) {

    formatImpl(std::format("Maximum warning messages: {} issued.", max_warn_messages_), LoggerSeverity::WARN);
    formatImpl(std::format("Further warning messages will be suppressed."),  LoggerSeverity::WARN);

  }

  if (max_warn_messages_ == 0 or warn_message_count_ <= max_warn_messages_) {

    return true;

  }

  return false;

}

bool kel::ExecEnvLogger::errorMessageLimits() {

  std::lock_guard<std::mutex> lock(limit_mutex_);

  ++error_message_count_;
  if (max_error_messages_ > 0 and error_message_count_ > max_error_messages_) {

    formatImpl(std::format("Maximum error messages: {} issued.", max_error_messages_), LoggerSeverity::ERROR);
    formatImpl(std::format("Forced Program exit. May terminate abnormally."), LoggerSeverity::ERROR);
    std::exit(EXIT_FAILURE);

  }

  return true;

}

