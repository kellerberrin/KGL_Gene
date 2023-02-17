//
// Created by kellerberrin on 6/09/17.
//

#ifndef KGL_LOGGING_H
#define KGL_LOGGING_H


#include <memory>
#include <spdlog/spdlog.h>  // Implement the logger using the spdlog library


namespace kellerberrin {   //  organization level namespace


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
  void SetMaxwarningMessages(int max_messages) { max_warn_messages_ = max_messages; }


  template<typename... Args> void trace(const std::string& message, Args... args) noexcept;
  template<typename... Args> void info(const std::string& message, Args... args) noexcept;
  template<typename... Args> void warn(const std::string& message, Args... args) noexcept;
  template<typename... Args> void error(const std::string& message, Args... args) noexcept;
  template<typename... Args> void critical(const std::string& message, Args... args) noexcept;


private:

  static constexpr const char* SPDLOG_DEFAULT_FORMAT{"%+"};  // Default format from spdlog
  static constexpr const char* DEFAULT_FLOAT_FORMAT = "0:.2f";

  std::unique_ptr<spdlog::logger> plog_impl_;
  std::atomic<int> max_error_messages_{100};     // Defaults to 100 error messages
  std::atomic<int> error_message_count_{0};     // number of error messages issued.
  std::atomic<int> max_warn_messages_{100};     // Defaults to 100 warning messages
  std::atomic<int> warn_message_count_{0};     // number of warning messages issued.


};


template<typename... Args> void Logger::trace(const std::string& message, Args... args) noexcept {

  plog_impl_->trace(message, args...);
  plog_impl_->flush();

}

template<typename... Args> void Logger::info(const std::string& message, Args... args) noexcept {

  plog_impl_->info(message, args...);
  plog_impl_->flush();

}

template<typename... Args> void Logger::warn(const std::string& message, Args... args) noexcept {

  if (max_warn_messages_ < 0 or warn_message_count_ <= max_warn_messages_) {

    plog_impl_->warn(message, args...);
    plog_impl_->flush();

  }

  if (max_warn_messages_ >= 0 and warn_message_count_ == max_warn_messages_) {

    plog_impl_->warn("Maximum warning messages: {} issued.", max_warn_messages_);
    plog_impl_->warn("Further warning messages will be suppressed.");

  }

  ++warn_message_count_;

}

template<typename... Args> void Logger::error(const std::string& message, Args... args) noexcept {

  plog_impl_->error(message , args...);
  plog_impl_->flush();

  ++error_message_count_;

  if (max_error_messages_ >= 0 and error_message_count_ >= max_error_messages_) {

    plog_impl_->error("Maximum error messages: {} issued.", max_error_messages_);
    plog_impl_->error("Program exits.");
    std::exit(EXIT_FAILURE);

  }

}

template<typename... Args> void Logger::critical(const std::string& message, Args... args) noexcept {


  plog_impl_->critical(message, args...);
  plog_impl_->flush();
  plog_impl_->error("Critical error; program exits.");
  std::exit(EXIT_FAILURE);

}

} // end namespace.

#endif //KGL_LOGGING_H


