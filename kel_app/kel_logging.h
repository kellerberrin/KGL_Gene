
#ifndef KEL_LOGGING_H
#define KEL_LOGGING_H


#include <spdlog/spdlog.h>  // Implement the logger using the spdlog library


#include <memory>
#include <iostream>
#include <source_location>
#include <string_view>



namespace kellerberrin {   //  organization level namespace



[[nodiscard]] constexpr spdlog::source_loc getSourceLocation(const std::source_location &location) {

  return spdlog::source_loc{location.file_name(),
                            static_cast<std::int32_t>(location.line()),
                            location.function_name()};

}

class FormatLocation {

public:

  template<typename String>
  FormatLocation(const String &format,
                 const std::source_location &location = std::source_location::current())
      : format_{format}, spdlog_locaton_{getSourceLocation(location)} {}
  ~FormatLocation() = default;

  [[nodiscard]] const std::string& format() const {return format_; }
  [[nodiscard]] const spdlog::source_loc& sourceLocation() const { return spdlog_locaton_; }

private:

  std::string format_;
  spdlog::source_loc spdlog_locaton_;

};


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


  template<typename... Args> void trace(const std::string& message, Args&&... args) noexcept;
  template<typename... Args> void info(const std::string& message, Args&&... args) noexcept;
  template<typename... Args> void warn(const std::string& message, Args&&... args) noexcept;
  template<typename... Args> void warn_location(FormatLocation format_location, Args &&...args) noexcept;
  template<typename... Args> void error(const std::string& message, Args&&... args) noexcept;
  template<typename... Args> void critical(const std::string& message, Args&&... args) noexcept;


private:

  static constexpr const char* SPDLOG_DEFAULT_FORMAT{"%+"};  // Default format from spdlog
  static constexpr const char* DEFAULT_FLOAT_FORMAT = "0:.2f";

  std::unique_ptr<spdlog::logger> log_impl_ptr_;
  // Negative values such as -1 allow unlimited warning and error messages.
  int max_error_messages_{100};     // Defaults to 100 error messages
  std::atomic<int> error_message_count_{0};     // Number of error messages issued.
  int max_warn_messages_{100};     // Defaults to 100 warning messages
  std::atomic<int> warn_message_count_{0};     // Number of warning messages issued.


};


template<typename... Args> void Logger::trace(const std::string& message, Args&&... args) noexcept {

  log_impl_ptr_->trace(fmt::runtime(message), std::forward<Args>(args)...);
  log_impl_ptr_->flush();

}

template<typename... Args> void Logger::info(const std::string& message, Args&&... args) noexcept {

  log_impl_ptr_->info(fmt::runtime(message), std::forward<Args>(args)...);
  log_impl_ptr_->flush();

}

template<typename... Args> void Logger::warn(const std::string& message, Args&&... args) noexcept {

  if (max_warn_messages_ < 0 or warn_message_count_ <= max_warn_messages_) {

    log_impl_ptr_->warn(fmt::runtime(message), std::forward<Args>(args)...);
    log_impl_ptr_->flush();

  }

  if (max_warn_messages_ >= 0 and warn_message_count_ == max_warn_messages_) {

    log_impl_ptr_->warn("Maximum warning messages: {} issued.", max_warn_messages_);
    log_impl_ptr_->warn("Further warning messages will be suppressed.");

  }

  ++warn_message_count_;

}

template<typename... Args>
void Logger::warn_location(FormatLocation format_location, Args&&... args) noexcept {

  if (max_warn_messages_ < 0 or warn_message_count_ <= max_warn_messages_) {

    log_impl_ptr_->log(format_location.sourceLocation(),
                       spdlog::level::warn,
                       fmt::runtime(format_location.format()),
                       std::forward<Args>(args)...);
    log_impl_ptr_->flush();

  }

  if (max_warn_messages_ >= 0 and warn_message_count_ == max_warn_messages_) {

    log_impl_ptr_->warn("Maximum warning messages: {} issued.", max_warn_messages_);
    log_impl_ptr_->warn("Further warning messages will be suppressed.");

  }

  ++warn_message_count_;

}


template<typename... Args> void Logger::error(const std::string& message, Args&&... args) noexcept {

  log_impl_ptr_->error(fmt::runtime(message) , std::forward<Args>(args)...);
  log_impl_ptr_->flush();

  ++error_message_count_;

  if (max_error_messages_ >= 0 and error_message_count_ == max_error_messages_) {

    log_impl_ptr_->error("Maximum error messages: {} issued.", max_error_messages_);
    log_impl_ptr_->error("Program exits.");
    log_impl_ptr_->flush();
    std::exit(EXIT_FAILURE);

  }

}

template<typename... Args> void Logger::critical(const std::string& message, Args&&... args) noexcept {

  log_impl_ptr_->critical(fmt::runtime(message), std::forward<Args>(args)...);
  log_impl_ptr_->critical("Program exits.");
  log_impl_ptr_->flush();
  std::exit(EXIT_FAILURE);

}

} // end namespace.

#endif //KEl_LOGGING_H


