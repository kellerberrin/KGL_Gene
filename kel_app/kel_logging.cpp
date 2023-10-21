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
  void info(const std::string& formatted_string) noexcept;
  void warn(const LogFormatLocation &format_location, const std::string& formatted_string) noexcept;
  void error(const LogFormatLocation &format_location, const std::string& formatted_string) noexcept;
  void critical(const LogFormatLocation &format_location, const std::string& formatted_string) noexcept;

  spdlog::logger& getLoggerImpl() { return *log_impl_ptr_; }

private:

  std::unique_ptr<spdlog::logger> log_impl_ptr_;

  static constexpr const char *SPDLOG_DEFAULT_FORMAT{"%+"};  // Default format from spdlog

  [[nodiscard]] static spdlog::source_loc spdLocation(const std::source_location &location);

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


spdlog::source_loc kel::ExecEnvLoggerImpl::spdLocation(const std::source_location &location) {

  return spdlog::source_loc{location.file_name(),
                            static_cast<std::int32_t>(location.line()),
                            location.function_name()};

}


void kel::ExecEnvLoggerImpl::info(const std::string& formatted_message) noexcept {

  log_impl_ptr_->log(spdlog::level::info, fmt::runtime(formatted_message));
  log_impl_ptr_->flush();

}


void kel::ExecEnvLoggerImpl::warn(const LogFormatLocation& format_location, const std::string& formatted_string) noexcept {

  log_impl_ptr_->log(spdLocation(format_location.location()), spdlog::level::warn, fmt::runtime(formatted_string));
  log_impl_ptr_->flush();

}

void kel::ExecEnvLoggerImpl::error(const LogFormatLocation& format_location, const std::string& formatted_string) noexcept {


  log_impl_ptr_->log(spdLocation(format_location.location()), spdlog::level::err, fmt::runtime(formatted_string));
  log_impl_ptr_->flush();


}

// Critical always displays the calling function.
void kel::ExecEnvLoggerImpl::critical(const LogFormatLocation& format_location, const std::string& formatted_string) noexcept {

  log_impl_ptr_->log(spdLocation(format_location.location()), spdlog::level::critical, fmt::runtime(formatted_string));
  log_impl_ptr_->log(spdLocation(format_location.location()), spdlog::level::critical, "Forced Program exit. May terminate abnormally.");
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


// These functions simply re-direct to the PIMPL implementation object.
void kel::ExecEnvLogger::infoImpl(const std::string& formatted_string) noexcept {

  log_impl_ptr_->info(formatted_string);

}

void kel::ExecEnvLogger::warnImpl(const LogFormatLocation& format_location, const std::string& formatted_string) noexcept {

  log_impl_ptr_->warn(format_location, formatted_string);

}

void kel::ExecEnvLogger::errorImpl(const LogFormatLocation& format_location, const std::string& formatted_string) noexcept {

  log_impl_ptr_->error(format_location, formatted_string);

}

void kel::ExecEnvLogger::criticalImpl(const LogFormatLocation& format_location, const std::string& formatted_string) noexcept {

  log_impl_ptr_->critical(format_location, formatted_string);

}


bool kel::ExecEnvLogger::warnMessageLimits() {

  std::lock_guard<std::mutex> lock(limit_mutex_);

  ++warn_message_count_;
  if (warn_message_count_ == max_warn_messages_) {

    log_impl_ptr_->getLoggerImpl().log(spdlog::level::warn, "Maximum warning messages: {} issued.", max_warn_messages_);
    log_impl_ptr_->getLoggerImpl().log(spdlog::level::warn, "Further warning messages will be suppressed.");
    log_impl_ptr_->getLoggerImpl().flush();

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

    log_impl_ptr_->getLoggerImpl().log(spdlog::level::err, "Maximum error messages: {} issued.", max_error_messages_);
    log_impl_ptr_->getLoggerImpl().log(spdlog::level::err, "Forced Program exit. May terminate abnormally.");
    log_impl_ptr_->getLoggerImpl().flush();
    std::exit(EXIT_FAILURE);

  }

  return true;

}

