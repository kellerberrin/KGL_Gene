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
                         info_message_count_,
                         warn_message_count_,
                         error_message_count_), LoggerSeverity::INFO);
  log_impl_ptr_ = nullptr;

}

void ExecEnvLogger::formatImpl(std::string &&formatted_string, LoggerSeverity severity) noexcept {

  log_impl_ptr_->formatImpl(std::forward<std::string &&>(formatted_string), severity);

}

void ExecEnvLogger::locationImpl(const LogFormatLocation &format_location,
                                      const std::string &formatted_string,
                                      LoggerSeverity severity) noexcept {

  log_impl_ptr_->locationImpl(format_location, formatted_string, severity);

}

bool ExecEnvLogger::warnMessageLimits() {

  ++warn_message_count_;
  if (warn_message_count_ == max_warn_messages_) {

    formatImpl(std::format("Maximum warning messages: {} issued.", max_warn_messages_), LoggerSeverity::WARN);
    formatImpl(std::format("Further warning messages will be suppressed."), LoggerSeverity::WARN);

  }

  if (max_warn_messages_ == 0 or warn_message_count_ <= max_warn_messages_) {

    return true;

  }

  return false;

}

bool ExecEnvLogger::errorMessageLimits() {

  ++error_message_count_;
  if (max_error_messages_ > 0 and error_message_count_ > max_error_messages_) {

    formatImpl(std::format("Maximum error messages: {} issued.", max_error_messages_), LoggerSeverity::ERROR);
    formatImpl(std::format("Forced Program exit. May terminate abnormally."), LoggerSeverity::ERROR);
    std::exit(EXIT_FAILURE);

  }

  return true;

}

} // End namespace

