//
// Created by kellerberrin on 20/10/23.
//

#include "kel_logging.h"
#include "kel_logging_stream.h"

#include <iostream>


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace kellerberrin {   //  organization level namespace

ExecEnvLogger::ExecEnvLogger(const std::string &module, const std::string &log_file) {


  log_impl_ptr_ = std::make_unique<StreamLoggerImpl>(module, log_file);

}

ExecEnvLogger::~ExecEnvLogger() {

  formatImpl(std::format("Message summary; INFO: {}, WARN: {}, ERROR: {}",
                         info_message_count_.load(std::memory_order_relaxed),
                         warn_message_count_.load(std::memory_order_relaxed),
                         error_message_count_.load(std::memory_order_relaxed)), LoggerSeverity::INFO);
  log_impl_ptr_ = nullptr;

}

void ExecEnvLogger::formatImpl(const std::string &formatted_string, LoggerSeverity severity) noexcept {

  std::lock_guard lock(message_mutex_);
  log_impl_ptr_->formatImpl(formatted_string, severity);

}

void ExecEnvLogger::locationImpl(const LogFormatLocation &format_location,
                                      const std::string &formatted_string,
                                      LoggerSeverity severity) noexcept {

  std::lock_guard lock(message_mutex_);
  log_impl_ptr_->locationImpl(format_location, formatted_string, severity);

}

bool ExecEnvLogger::warnMessageLimits() noexcept {

  const size_t count = ++warn_message_count_;
  const size_t max_messages = max_warn_messages_.load(std::memory_order_relaxed);
  if (max_messages != UNBOUNDED_MESSAGES and count == max_messages) {

    formatImpl(std::format("Maximum warning messages: {} issued.", max_messages), LoggerSeverity::WARN);
    formatImpl(std::format("Further warning messages will be suppressed."), LoggerSeverity::WARN);

  }

  return max_messages == UNBOUNDED_MESSAGES or count <= max_messages;

}

bool ExecEnvLogger::errorMessageLimits() noexcept {

  const size_t count = ++error_message_count_;
  const size_t max_messages = max_error_messages_.load(std::memory_order_relaxed);
  return max_messages == UNBOUNDED_MESSAGES or count <= max_messages;

}

// Finalizes logging. flushes and closes any open files. Called by logging a critical msg
void ExecEnvLogger::flushStreamsImpl() noexcept {

  log_impl_ptr_->flushStreamsImpl();

}


} // End namespace
