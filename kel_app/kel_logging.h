//
// Created by kellerberrin on 20/10/23.
//

#ifndef KEL_EXECENV_LOGGING_H
#define KEL_EXECENV_LOGGING_H

#include <atomic>
#include <cstddef>
#include <cstdlib>
#include <format>
#include <memory>
#include <mutex>
#include <source_location>
#include <string>
#include <utility>


namespace kellerberrin {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////
//
// An auxiliary object to retrieve std::source_location information.
//
////////////////////////////////////////////////////////////////////////////////////////

class LogFormatLocation {

public:

  // Important note. This constructor must be a template in order to modify the order of argument substitution.
  template<typename String>
  LogFormatLocation(const String &format,
                    const std::source_location &location = std::source_location::current())
  : format_(format), location_(location) {}

  ~LogFormatLocation() = default;

  [[nodiscard]] const std::string &format() const noexcept { return format_; }
  [[nodiscard]] const std::source_location &location() const noexcept { return location_; }

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
// If a preset number of warn() messages are issued, further warnings are suppressed (disabled with UNBOUNDED_MESSAGES).
// If a preset number of error() messages are issued, the application terminates (disabled with UNBOUNDED_MESSAGES).
// If critical() is called the message is output and the application terminates immediately.
// Message logging is thread-safe, however an abrupt logger initiated application termination may
// cause a seg-fault in multi-threaded code.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Forward Declaration of the implementation PIMPL object (implemented using iostreams).
class StreamLoggerImpl;

// The logger syntax uses the standard idiom; std::format("Arg1: {}, Arg2: {} ... Argn: {}", arg1, arg1, ..., argn).
// The logger is available to all application files that have "#include exec_env.h" in the include chain.
// Syntax is:
// ExecEnv::log().info("Arg1: {}, Arg2: {} ... Argn: {}", arg1, arg1, ..., argn);
// ExecEnv::log().warn("Arg1: {}, Arg2: {} ... Argn: {}", arg1, arg1, ..., argn);
// ExecEnv::log().error("Arg1: {}, Arg2: {} ... Argn: {}", arg1, arg1, ..., argn);
// ExecEnv::log().critical("Arg1: {}, Arg2: {} ... Argn: {}", arg1, arg1, ..., argn);

class ExecEnvLogger {

public:

  ExecEnvLogger(const std::string& module, const std::string& log_file);
  ~ExecEnvLogger();

  ExecEnvLogger(const ExecEnvLogger&) = delete;
  ExecEnvLogger(ExecEnvLogger&&) = delete;
  ExecEnvLogger& operator=(const ExecEnvLogger&) = delete;
  ExecEnvLogger& operator=(ExecEnvLogger&&) = delete;

  constexpr static size_t UNBOUNDED_MESSAGES = 0; // Zero (0) is unlimited messages for WARN abd ERROR
  void setMaxErrorMessages(size_t max_messages) { max_error_messages_ = max_messages; } // Zero (0) is unlimited.
  void setMaxWarningMessages(size_t max_messages) { max_warn_messages_ = max_messages; } // Zero (0) is unlimited.

  enum class LoggerSeverity { INFO, WARN, ERROR, CRITICAL };

// Select between message location information or compile-time argument checking.
// In general; WARN, ERROR and CRITICAL will report the code location "[module_file:line_no]" when emitting
// a message. However, disabling the location #define will re-compile with rigorous compile-time argument checking.
// Very useful for debugging any possible runtime errors.

//#define EXECENV_LOGGER_INFO_LOCATION 1    // Defaults to no location information (compile time argument checking).
#define EXECENV_LOGGER_WARN_LOCATION 1      // Defaults to source file and line location information (no argument checking).
#define EXECENV_LOGGER_ERROR_LOCATION 1     // Defaults to source file and line location information (no argument checking).
#define EXECENV_LOGGER_CRITICAL_LOCATION 1  // Defaults to source file and line location information (no argument checking).

#ifdef EXECENV_LOGGER_INFO_LOCATION

  template<typename... Args> void info(const LogFormatLocation& format_location, Args &&...args) noexcept {

    info_message_count_++;
    locationFormat(format_location, LoggerSeverity::INFO, std::forward<Args>(args)...);

  }

#else

  template<typename... Args> void info(std::format_string<Args...> format, Args&&... args) noexcept {

    info_message_count_++;
    logFormat(format, LoggerSeverity::INFO, std::forward<Args>(args)...);

  }


#endif

// Select between message location information or compile-time argument checking.
#ifdef EXECENV_LOGGER_WARN_LOCATION

  template<typename... Args> void warn(const LogFormatLocation& format_location, Args &&...args) noexcept {

    if (warnMessageLimits()) {

      locationFormat(format_location, LoggerSeverity::WARN, std::forward<Args>(args)...);

    }

  }

#else

  template<typename... Args> void warn(std::format_string<Args...> format, Args&&... args) noexcept {

    if (warnMessageLimits()) {

      logFormat(format, LoggerSeverity::WARN, std::forward<Args>(args)...);

    }

  }

#endif
  
// Select between message location information or compile-time argument checking.
#ifdef EXECENV_LOGGER_ERROR_LOCATION

  template<typename... Args> void error(const LogFormatLocation& format_location, Args&&... args) noexcept {

    locationFormat(format_location, LoggerSeverity::ERROR, std::forward<Args>(args)...);
    if (not errorMessageLimits()) {

      logFormat("Maximum error messages: {} issued.",
                 LoggerSeverity::ERROR,
                 max_error_messages_.load(std::memory_order_relaxed));
      logFormat("Forced Program exit. May terminate abnormally.", LoggerSeverity::ERROR);
      std::exit(EXIT_FAILURE);

    }

  }

#else

  template<typename... Args> void error(std::format_string<Args...> format, Args&&... args) noexcept {

    logFormat(format, LoggerSeverity::ERROR, std::forward<Args>(args)...);
    if (not errorMessageLimits()) {

      logFormat("Maximum error messages: {} issued.", LoggerSeverity::ERROR, max_error_messages_);
      logFormat("Forced Program exit. May terminate abnormally.", LoggerSeverity::ERROR);
      std::exit(EXIT_FAILURE);

    }

  }

#endif

// Select between message location information or compile-time argument checking.
#ifdef EXECENV_LOGGER_CRITICAL_LOCATION

  template<typename... Args> [[noreturn]] void critical(const LogFormatLocation& format_location, Args&&... args) noexcept {

    locationFormat(format_location, LoggerSeverity::CRITICAL, std::forward<Args>(args)...);
    logFormat("Forced Program exit. May terminate abnormally.", LoggerSeverity::CRITICAL);
    std::exit(EXIT_FAILURE);

  }

#else

  template<typename... Args> void critical(std::format_string<Args...> format, Args&&... args) noexcept {

    logFormat(format, LoggerSeverity::CRITICAL, std::forward<Args>(args)...);
    logFormat("Forced Program exit. May terminate abnormally.", LoggerSeverity::CRITICAL);
    std::exit(EXIT_FAILURE);

  }

#endif


private:

  std::unique_ptr<StreamLoggerImpl> log_impl_ptr_; // The PIMPL logging implementation object.
  // Message counters are thread safe.
  std::atomic<std::size_t> info_message_count_{0};
  std::atomic<std::size_t> warn_message_count_{0};
  std::atomic<std::size_t> error_message_count_{0};

  std::atomic<std::size_t> max_error_messages_{100};
  std::atomic<std::size_t> max_warn_messages_{100};
  // Messaging output is thread safe.
  std::mutex message_mutex_;

  bool warnMessageLimits() noexcept; // Stops issuing messages after max warn messages reached.
  bool errorMessageLimits() noexcept; // Forces program termination after max error messages reached.

  // This function captures any formatting errors and argument mis-matches at compile time.
  template<typename... Args> void logFormat( std::format_string<Args...> format,
                                             LoggerSeverity severity,
                                             Args &&...args) noexcept;
  // This function automatically appends source file and line number ("[kel_utility.cpp:46]") to the message.
  template<typename... Args> void locationFormat( const LogFormatLocation& format_location,
                                                  LoggerSeverity severity,
                                                  Args&&... args) noexcept;

  // These functions simply re-direct to the PIMPL implementation object.
  void formatImpl(const std::string& formatted_string, LoggerSeverity severity) noexcept;
  void locationImpl(const LogFormatLocation& format_location,
                    const std::string& formatted_string,
                    LoggerSeverity severity) noexcept;

};

template<typename... Args> void ExecEnvLogger::locationFormat( const LogFormatLocation& format_location,
                                                               LoggerSeverity severity,
                                                               Args&&... args) noexcept {

  try {

    std::string formatted_message = std::vformat(format_location.format(), std::make_format_args(args...));
    std::lock_guard lock(message_mutex_);
    locationImpl(format_location, formatted_message, severity);

  } catch (const std::exception& e) {

    formatImpl(std::format("Unexpected exception logging - error: {}", e.what()), LoggerSeverity::ERROR);

  } catch (...) {

    formatImpl(std::format("Unexpected exception logging"), LoggerSeverity::ERROR);

  }

}

template<typename... Args> void ExecEnvLogger::logFormat( std::format_string<Args...> format,
                                                          LoggerSeverity severity,
                                                          Args &&...args) noexcept {

  try {

    std::string formatted_message = std::format(format, std::forward<Args>(args)...);
    std::lock_guard<std::mutex> lock(message_mutex_);
    formatImpl(formatted_message, severity);

  } catch (const std::exception& e) {

    formatImpl(std::format("Unexpected exception logging - error: {}", e.what()), LoggerSeverity::ERROR);

  } catch (...) {

    formatImpl(std::format("Unexpected exception logging"), LoggerSeverity::ERROR);

  };

}

} // Namespace.


#endif  // KEL_EXECENV_LOGGING_H
