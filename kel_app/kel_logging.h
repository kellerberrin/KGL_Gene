
#ifndef KEL_LOGGING_H
#define KEL_LOGGING_H


#include <memory>
#include <iostream>
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

  void setLevel(Severity level) noexcept;
  void setFormat(const std::string& message) noexcept;

  void setMaxErrorMessages(int max_messages) { max_error_messages_ = max_messages; }
  void setMaxWarningMessages(int max_messages) { max_warn_messages_ = max_messages; }


  template<typename... Args> void trace(const std::string& message, Args... args) noexcept;
  template<typename... Args> void info(const std::string& message, Args... args) noexcept;
  template<typename... Args> void warn(const std::string& message, Args... args) noexcept;
  template<typename... Args> void error(const std::string& message, Args... args) noexcept;
  template<typename... Args> void critical(const std::string& message, Args... args) noexcept;


private:

  static constexpr const char* SPDLOG_DEFAULT_FORMAT{"%+"};  // Default format from spdlog
  static constexpr const char* DEFAULT_FLOAT_FORMAT = "0:.2f";

  std::unique_ptr<spdlog::logger> plog_impl_;
  // Negative values such as -1 allow unlimited warning and error messages.
  int max_error_messages_{100};     // Defaults to 100 error messages
  std::atomic<int> error_message_count_{0};     // Number of error messages issued.
  int max_warn_messages_{100};     // Defaults to 100 warning messages
  std::atomic<int> warn_message_count_{0};     // Number of warning messages issued.


};


template<typename... Args> void Logger::trace(const std::string& message, Args... args) noexcept {

  plog_impl_->trace(fmt::runtime(message), args...);
  plog_impl_->flush();

}

template<typename... Args> void Logger::info(const std::string& message, Args... args) noexcept {

  plog_impl_->info(fmt::runtime(message), args...);
  plog_impl_->flush();

}

template<typename... Args> void Logger::warn(const std::string& message, Args... args) noexcept {

  if (max_warn_messages_ < 0 or warn_message_count_ <= max_warn_messages_) {

    plog_impl_->warn(fmt::runtime(message), args...);
    plog_impl_->flush();

  }

  if (max_warn_messages_ >= 0 and warn_message_count_ == max_warn_messages_) {

    plog_impl_->warn("Maximum warning messages: {} issued.", max_warn_messages_);
    plog_impl_->warn("Further warning messages will be suppressed.");

  }

  ++warn_message_count_;

}

template<typename... Args> void Logger::error(const std::string& message, Args... args) noexcept {

  plog_impl_->error(fmt::runtime(message) , args...);
  plog_impl_->flush();

  ++error_message_count_;

  if (max_error_messages_ >= 0 and error_message_count_ == max_error_messages_) {

    plog_impl_->error("Maximum error messages: {} issued.", max_error_messages_);
    plog_impl_->error("Program exits.");
    plog_impl_->flush();
    std::exit(EXIT_FAILURE);

  }

}

template<typename... Args> void Logger::critical(const std::string& message, Args... args) noexcept {

  plog_impl_->critical(fmt::runtime(message), args...);
  plog_impl_->critical("Program exits.");
  plog_impl_->flush();
  std::exit(EXIT_FAILURE);

}

} // end namespace.

#endif //KEl_LOGGING_H


