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


#include "kel_logging_new.h"
#include "spdlog/sinks/basic_file_sink.h"

#include <iostream>


namespace kel = kellerberrin;


kel::NewLogger::NewLogger(const std::string& module, const std::string& log_file) {

  spdlog::set_pattern(SPDLOG_DEFAULT_FORMAT);
  spdlog::set_level(spdlog::level::trace);
  std::vector<spdlog::sink_ptr> sinks;
  sinks.push_back(std::make_shared<spdlog::sinks::ansicolor_stderr_sink_mt>());
  sinks.push_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file));
  log_impl_ptr_ = std::make_unique<spdlog::logger>(module, sinks.begin(), sinks.end());

}


bool kel::NewLogger::warnMessageLimits() {

  std::lock_guard<std::mutex> lock(limit_mutex_);

  ++warn_message_count_;
  if (warn_message_count_ == max_warn_messages_) {

    log_impl_ptr_->log(spdlog::level::warn, "Maximum warning messages: {} issued.", max_warn_messages_);
    log_impl_ptr_->log(spdlog::level::warn, "Further warning messages will be suppressed.");
    log_impl_ptr_->flush();

  }

  if (max_warn_messages_ == 0 or warn_message_count_ <= max_warn_messages_) {

    return true;

  }

  return false;

}

bool kel::NewLogger::errorMessageLimits() {

  std::lock_guard<std::mutex> lock(limit_mutex_);

  ++error_message_count_;
  if (max_error_messages_ > 0 and error_message_count_ > max_error_messages_) {

    log_impl_ptr_->log(spdlog::level::err, "Maximum error messages: {} issued.", max_error_messages_);
    log_impl_ptr_->log(spdlog::level::err, "Forced Program exit. May terminate abnormally.");
    log_impl_ptr_->flush();
    std::exit(EXIT_FAILURE);

  }

  return true;

}

