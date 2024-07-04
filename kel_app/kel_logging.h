//
// Created by kellerberrin on 20/10/23.
//

#ifndef KEL_EXECENV_LOGGING_H
#define KEL_EXECENV_LOGGING_H

#include <source_location>
#include <format>
#include <memory>
#include <mutex>


namespace kellerberrin {   //  organization level namespace



/////////////////////////////////////////////////////////////////////
//
// An auxiliary object to retrieve std::source_location information.
//
/////////////////////////////////////////////////////////////////////


class LogFormatLocation {

public:

  // Important note. This constructor must be a template in order to modify the order of argument substitution.
  template<typename String>
  LogFormatLocation(const String &format,
                    const std::source_location &location = std::source_location::current())
  : format_(format), location_(location) {}

  ~LogFormatLocation() = default;

  [[nodiscard]] const std::string &format() const { return format_; }

  [[nodiscard]] const std::source_location &location() const { return location_; }

private:

  const std::string format_;
  const std::source_location location_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The logger class used by the ExecEnv application execution environment which
// provides logging functionality to all application source files.
// The logging syntax uses the standard idiom; std::format("Arg1: {}, Arg2: {} ... Argn: {}", arg1, arg1, ..., argn).
// Four levels of message are provided, info, warn, error and critical.
// If a preset number of warn() messages are issued, further warnings are suppressed.
// If a preset number of error() messages are issued, the application terminates.
// If critical() is called the message is output and the application terminates immediately.
// Message logging is thread-safe, however an abrupt logger initiated application termination may
// cause a seg-fault in multi-threaded code.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Forward Declaration of the implementation PIMPL object (implemented using iostreams).
class StreamLoggerImpl;

// The logger syntax uses the standard idiom; std::format("Arg1: {}, Arg2: {} ... Argn: {}", arg1, arg1, ..., argn).
// The logger is available to all application files that have "#include exec_env.h" in the include chain.
class ExecEnvLogger {

public:

  ExecEnvLogger(const std::string& module, const std::string& log_file);
  ~ExecEnvLogger();
  ExecEnvLogger(const ExecEnvLogger&) = delete;
  ExecEnvLogger(ExecEnvLogger&&) = delete;
  ExecEnvLogger& operator=(const ExecEnvLogger&) = delete;

  void setMaxErrorMessages(size_t max_messages) { max_error_messages_ = max_messages; } // Zero (0) is unlimited.
  void setMaxWarningMessages(size_t max_messages) { max_warn_messages_ = max_messages; } // Zero (0) is unlimited.
  enum class LoggerSeverity { INFO, WARN, ERROR, CRITICAL };

// Select between message location information or compile-time argument checking.
//#define EXECENV_LOGGER_INFO_LOCATION 1
#ifdef EXECENV_LOGGER_INFO_LOCATION

  template<typename... Args> void info(LogFormatLocation format_location, Args &&...args) noexcept {

    std::lock_guard<std::mutex> lock(message_mutex_);

    info_message_count_++;
    std::string formatted_message = std::vformat(format_location.format(), std::make_format_args(args...));
    locationImpl(format_location, formatted_message, LoggerSeverity::INFO);

  }

#else

  template<typename... Args> void info(std::format_string<Args...> format, Args&&... args) noexcept {

    std::lock_guard<std::mutex> lock(message_mutex_);

    info_message_count_++;
    formatImpl(std::format(format, std::forward<Args>(args)...), LoggerSeverity::INFO);

  }

#endif

// Select between message location information or compile-time argument checking.
#define EXECENV_LOGGER_WARN_LOCATION 1
#ifdef EXECENV_LOGGER_WARN_LOCATION

  template<typename... Args> void warn(LogFormatLocation format_location, Args &&...args) noexcept {

    std::lock_guard<std::mutex> lock(message_mutex_);

    if (warnMessageLimits()) {

      std::string formatted_message = std::vformat(format_location.format(), std::make_format_args(args...));
      locationImpl(format_location, formatted_message, LoggerSeverity::WARN);

    }

  }

#else

  template<typename... Args> void warn(std::format_string<Args...> format, Args&&... args) noexcept {

    std::lock_guard<std::mutex> lock(message_mutex_);

    if (warnMessageLimits()) {

      formatImpl(std::format(format, std::forward<Args>(args)...), LoggerSeverity::WARN);

    }

  }

#endif

// Select between message location information or compile-time argument checking.
#define EXECENV_LOGGER_ERROR_LOCATION 1
#ifdef EXECENV_LOGGER_ERROR_LOCATION


  template<typename... Args> void error(LogFormatLocation format_location, Args&&... args) noexcept {

    std::lock_guard<std::mutex> lock(message_mutex_);

    if (errorMessageLimits()) {

      std::string formatted_message = std::vformat(format_location.format(), std::make_format_args(args...));
      locationImpl(format_location, formatted_message, LoggerSeverity::ERROR);

    }

  }


#else

  template<typename... Args> void error(std::format_string<Args...> format, Args&&... args) noexcept {

    std::lock_guard<std::mutex> lock(message_mutex_);

    if (errorMessageLimits()) {

      formatImpl(std::format(format, std::forward<Args>(args)...), LoggerSeverity::ERROR);

    }

  }

#endif

// Select between message location information or compile-time argument checking.
#define EXECENV_LOGGER_CRITICAL_LOCATION 1
#ifdef EXECENV_LOGGER_CRITICAL_LOCATION

  template<typename... Args> void critical(LogFormatLocation format_location, Args&&... args) noexcept {

    std::lock_guard<std::mutex> lock(message_mutex_);

    std::string formatted_message = std::vformat(format_location.format(), std::make_format_args(args...));
    locationImpl(format_location, formatted_message, LoggerSeverity::CRITICAL);
    formatImpl(std::format("Forced Program exit. May terminate abnormally."), LoggerSeverity::CRITICAL);
    std::exit(EXIT_FAILURE);

  }

#else

  template<typename... Args> void critical(std::format_string<Args...> format, Args&&... args) noexcept {

    std::lock_guard<std::mutex> lock(message_mutex_);

    formatImpl(std::format(format, std::forward<Args>(args)...), LoggerSeverity::CRITICAL);
    formatImpl(std::format("Forced Program exit. May terminate abnormally."), LoggerSeverity::CRITICAL);
    std::exit(EXIT_FAILURE);

  }

#endif


private:

  std::unique_ptr<StreamLoggerImpl> log_impl_ptr_; // The PIMPL logging implementation object.
  size_t info_message_count_{0};     // Number of info messages issued.
  size_t max_error_messages_{100};     // Defaults to 100 error messages, zero (0) is unlimited.
  size_t error_message_count_{0};     // Number of error messages issued.
  size_t max_warn_messages_{100};     // Defaults to 100 warning messages, zero (0) is unlimited.
  size_t warn_message_count_{0};     // Number of warning messages issued.
  std::mutex message_mutex_;          // Messaging is thread safe.

  bool warnMessageLimits(); // Stops issuing messages after max warn messages reached.
  bool errorMessageLimits(); // Forces program termination after max error messages reached.

  // These functions simply re-direct to the PIMPL implementation object.
  void formatImpl(std::string&& formatted_string, LoggerSeverity severity) noexcept;
  void locationImpl(const LogFormatLocation& format_location,
                    const std::string& formatted_string,
                    LoggerSeverity severity) noexcept;

};



} // Namespace.


#endif  // KEL_EXECENV_LOGGING_H
