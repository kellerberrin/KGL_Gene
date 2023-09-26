
#ifndef KEL_LOGGING_H
#define KEL_LOGGING_H


#include <spdlog/spdlog.h>  // Implement the logger using the spdlog library


#include <memory>
#include <iostream>
#include <source_location>
#include <string_view>
#include <mutex>



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
      : format_{format}, spdlog_locaton_{getSourceLocation(location)} {

        function_format_ = "[" + std::string(location.function_name()) + "] " + format;

  }
  ~FormatLocation() = default;

  [[nodiscard]] const std::string& format() const {return format_; }
  [[nodiscard]] const std::string& functionFormat() const {return function_format_; }
  [[nodiscard]] const spdlog::source_loc& sourceLocation() const { return spdlog_locaton_; }

private:

  std::string format_;
  std::string function_format_;
  spdlog::source_loc spdlog_locaton_;

};


class Logger {

public:

  Logger(const std::string& module, const std::string& log_file);
  ~Logger() = default;
  Logger(const Logger&) = delete;
  Logger(Logger&&) = delete;
  Logger& operator=(const Logger&) = delete;

  void setMaxErrorMessages(size_t max_messages) { max_error_messages_ = max_messages; } // Zero (0) is unlimited.
  void setMaxWarningMessages(size_t max_messages) { max_warn_messages_ = max_messages; } // Zero (0) is unlimited.

// LOGGER_SOURCE_AUGMENTATION displays all source information for all messages types.
//#define LOGGER_SOURCE_AUGMENTATION 1
#ifdef LOGGER_SOURCE_AUGMENTATION

  template<typename... Args> void info(FormatLocation format_location, Args&&... args) noexcept;

#else

  template<typename... Args> void info(const std::string& message, Args&&... args) noexcept;

#endif


  template<typename... Args> void warn(FormatLocation format_location, Args &&...args) noexcept;
  template<typename... Args> void error(FormatLocation format_location, Args&&... args) noexcept;
  template<typename... Args> void critical(FormatLocation format_location, Args&&... args) noexcept;


private:

  std::unique_ptr<spdlog::logger> log_impl_ptr_;
  size_t max_error_messages_{100};     // Defaults to 100 error messages, zero (0) is unlimited.
  size_t error_message_count_{0};     // Number of error messages issued.
  size_t max_warn_messages_{100};     // Defaults to 100 warning messages, zero (0) is unlimited.
  size_t warn_message_count_{0};     // Number of warning messages issued.
  std::mutex limit_mutex_;

  static constexpr const char* SPDLOG_DEFAULT_FORMAT{"%+"};  // Default format from spdlog

  bool warnMessageLimits(); // Stops issuing messages after max warn messages reached.
  bool errorMessageLimits(); // Forces program termination after max error messages reached.

};

#ifdef LOGGER_SOURCE_AUGMENTATION

template<typename... Args>
void Logger::info(FormatLocation format_location, Args&&... args) noexcept {

  log_impl_ptr_->log(format_location.sourceLocation(),
                     spdlog::level::info,
                     fmt::runtime(format_location.functionFormat()),
                     std::forward<Args>(args)...);
  log_impl_ptr_->flush();

}


template<typename... Args>
void Logger::warn(FormatLocation format_location, Args&&... args) noexcept {

  if (warnMessageLimits()) {

    log_impl_ptr_->log(format_location.sourceLocation(),
                       spdlog::level::warn,
                       fmt::runtime(format_location.functionFormat()),
                       std::forward<Args>(args)...);
    log_impl_ptr_->flush();

  }

}

template<typename... Args>
void Logger::error(FormatLocation format_location, Args&&... args) noexcept {

  if (errorMessageLimits()) {

    log_impl_ptr_->log(format_location.sourceLocation(),
                       spdlog::level::err,
                       fmt::runtime(format_location.functionFormat()),
                       std::forward<Args>(args)...);
    log_impl_ptr_->flush();

  }

}


#else // LOGGER_SOURCE_AUGMENTATION


template<typename... Args>
void Logger::info(const std::string& message, Args&&... args) noexcept {

  log_impl_ptr_->log(spdlog::level::info, fmt::runtime(message), std::forward<Args>(args)...);
  log_impl_ptr_->flush();

}


template<typename... Args>
void Logger::warn(FormatLocation format_location, Args&&... args) noexcept {

  if (warnMessageLimits()) {

    log_impl_ptr_->log(format_location.sourceLocation(),
                       spdlog::level::warn,
                       fmt::runtime(format_location.format()),
                       std::forward<Args>(args)...);
    log_impl_ptr_->flush();

  }

}

template<typename... Args>
void Logger::error(FormatLocation format_location, Args&&... args) noexcept {

  if (errorMessageLimits()) {

    log_impl_ptr_->log(format_location.sourceLocation(),
                       spdlog::level::err,
                       fmt::runtime(format_location.format()),
                       std::forward<Args>(args)...);
    log_impl_ptr_->flush();

  }

}

#endif  // LOGGER_SOURCE_AUGMENTATION

// Critical always displays the calling function.
template<typename... Args>
void Logger::critical(FormatLocation format_location, Args&&... args) noexcept {

  log_impl_ptr_->log(format_location.sourceLocation(),
                     spdlog::level::critical,
                     fmt::runtime(format_location.functionFormat()),
                     std::forward<Args>(args)...);
  log_impl_ptr_->log(format_location.sourceLocation(),
                     spdlog::level::critical,
                     "Forced Program exit. May terminate abnormally.");
  log_impl_ptr_->flush();
  std::exit(EXIT_FAILURE);

}


} // end namespace.

#endif //KEl_LOGGING_H


