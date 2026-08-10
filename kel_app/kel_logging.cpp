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

  log_impl_ptr_->formatImpl(formatted_string, severity);

}

void ExecEnvLogger::locationImpl(const LogFormatLocation &format_location,
                                      const std::string &formatted_string,
                                      LoggerSeverity severity) noexcept {

  log_impl_ptr_->locationImpl(format_location, formatted_string, severity);

}

bool ExecEnvLogger::warnMessageLimits() noexcept {

  ++warn_message_count_;
  if (warn_message_count_ == max_warn_messages_ and max_warn_messages_ != UNBOUNDED_MESSAGES) {

    formatImpl(std::format("Maximum warning messages: {} issued.",
                                         max_warn_messages_.load(std::memory_order_relaxed)),
                                         LoggerSeverity::WARN);
    formatImpl(std::format("Further warning messages will be suppressed."), LoggerSeverity::WARN);

  }

  return max_warn_messages_ == UNBOUNDED_MESSAGES or warn_message_count_ <= max_warn_messages_;

}

bool ExecEnvLogger::errorMessageLimits() noexcept {

  ++error_message_count_;
  return max_error_messages_ == UNBOUNDED_MESSAGES or error_message_count_ <= max_error_messages_;

}


} // End namespace

