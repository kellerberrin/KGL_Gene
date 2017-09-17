// MIT License
//
// Copyright (c) 2017
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
//
//
// Created by kellerberrin on 6/09/17.
//

#ifndef KGL_LOGGING_H
#define KGL_LOGGING_H


#include <memory>
#include "spdlog/spdlog.h"  // Implement the logger using the spdlog library


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class Logger {

public:

  Logger(const std::string& module, const std::string& log_file);
  ~Logger() = default;
  Logger(const Logger&) = delete;
  Logger(Logger&&) = delete;
  Logger& operator=(const Logger&) = delete;

  enum class Severity { Trace, Info, Warn, Error, Critical };

  void SetLevel(Severity level) noexcept;
  void SetFormat(const std::string& message) noexcept;

  static constexpr const int NO_EXIT_ON_ERROR = -1;  // No exit on error messages
  void SetMaxErrorMessages(int max_messages) { max_error_messages_ = max_messages; }

  template<typename M, typename... Args> void trace(M& message, Args... args) noexcept;
  template<typename M, typename... Args> void info(M& message, Args... args) noexcept;
  template<typename M, typename... Args> void warn(M& message, Args... args) noexcept;
  template<typename M, typename... Args> void error(M& message, Args... args) noexcept;
  template<typename M, typename... Args> void critical(M& message, Args... args) noexcept;

private:

  static constexpr const char* SPDLOG_DEFAULT_FORMAT{"%+"};  // Default format from spdlog

  std::unique_ptr<spdlog::logger> plog_impl_;
  std::atomic<int> max_error_messages_{100};     // Defaults to 100 error messages
  std::atomic<int> error_message_count_{0};     // number of error messages issued.

};

template<typename M, typename... Args> void Logger::trace(M& message, Args... args) noexcept {

  plog_impl_->trace(message, args...);
  plog_impl_->flush();

}

template<typename M, typename... Args> void Logger::info(M& message, Args... args) noexcept {

  plog_impl_->info(message, args...);
  plog_impl_->flush();

}

template<typename M, typename... Args> void Logger::warn(M& message, Args... args) noexcept {

  plog_impl_->warn(message, args...);
  plog_impl_->flush();

}

template<typename M, typename... Args> void Logger::error(M& message, Args... args) noexcept {

  plog_impl_->error(message, args...);
  plog_impl_->flush();

  ++error_message_count_;

  if (max_error_messages_ >= 0 and error_message_count_ >= max_error_messages_) {

    plog_impl_->error("Maximum error messages: {} issued.", max_error_messages_);
    plog_impl_->error("Program exits.");
    std::exit(EXIT_FAILURE);

  }

}

template<typename M, typename... Args> void Logger::critical(M& message, Args... args) noexcept {

  plog_impl_->critical(message, args...);
  plog_impl_->flush();
  plog_impl_->error("Critical error; program exits.");
  std::exit(EXIT_FAILURE);

}

}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_LOGGING_H


