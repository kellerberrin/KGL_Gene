//
// Created by kellerberrin on 20/10/23.
//

#ifndef KEL_EXECENV_LOGGING_H
#define KEL_EXECENV_LOGGING_H


#include <format>
#include <memory>
#include <source_location>
#include <mutex>


namespace kellerberrin {   //  organization level namespace

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
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// An object to automatically retrieve std::source_location information for warn, error and critical messages.
class LogFormatLocation {

public:

  template<typename String>
  LogFormatLocation(const String &format,
                    const std::source_location &location = std::source_location::current())
  : format_{format}, location_(location) {}
  ~LogFormatLocation() = default;

  [[nodiscard]] const std::string& format() const {return format_; }
  [[nodiscard]] const std::source_location& location() const { return location_; }

private:

  const std::string format_;
  const std::source_location location_;

};

// Forward Declaration of the implementation PIMPL object (currently implemented using the spdlog library).
class ExecEnvLoggerImpl;

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


  template<typename... Args> void info(std::string message, Args&&... args) noexcept;
  template<typename... Args> void warn(LogFormatLocation format_location, Args &&...args) noexcept;
  template<typename... Args> void error(LogFormatLocation format_location, Args&&... args) noexcept;
  template<typename... Args> void critical(LogFormatLocation format_location, Args&&... args) noexcept;


private:

  std::unique_ptr<ExecEnvLoggerImpl> log_impl_ptr_;
  size_t max_error_messages_{100};     // Defaults to 100 error messages, zero (0) is unlimited.
  size_t error_message_count_{0};     // Number of error messages issued.
  size_t max_warn_messages_{100};     // Defaults to 100 warning messages, zero (0) is unlimited.
  size_t warn_message_count_{0};     // Number of warning messages issued.
  std::mutex limit_mutex_;          // Messaging is thread safe.

  bool warnMessageLimits(); // Stops issuing messages after max warn messages reached.
  bool errorMessageLimits(); // Forces program termination after max error messages reached.

  // These functions simply re-direct to the PIMPL implementation object.
  void infoImpl(const std::string& formatted_string) noexcept;
  void warnImpl(const LogFormatLocation& format_location, const std::string& formatted_string) noexcept;
  void errorImpl(const LogFormatLocation& format_location, const std::string& formatted_string) noexcept;
  void criticalImpl(const LogFormatLocation& format_location, const std::string& formatted_string) noexcept;

};



template<typename... Args>
void ExecEnvLogger::info(std::string message, Args&&... args) noexcept {

  std::string formatted_message = std::vformat(message, std::make_format_args(args...));
  infoImpl(formatted_message);

}


template<typename... Args>
void ExecEnvLogger::warn(LogFormatLocation format_location, Args&&... args) noexcept {

  if (warnMessageLimits()) {

    std::string formatted_message = std::vformat(format_location.format(), std::make_format_args(args...));
    warnImpl(format_location, formatted_message);

  }

}

template<typename... Args>
void ExecEnvLogger::error(LogFormatLocation format_location, Args&&... args) noexcept {

  if (errorMessageLimits()) {

    std::string formatted_message = std::vformat(format_location.format(), std::make_format_args(args...));
    errorImpl(format_location, formatted_message);

  }

}

template<typename... Args>
void ExecEnvLogger::critical(LogFormatLocation format_location, Args&&... args) noexcept {

  std::string formatted_message = std::vformat(format_location.format(), std::make_format_args(args...));
  criticalImpl(format_location, formatted_message);
  std::exit(EXIT_FAILURE);

}



} // Namespace.


#endif  // KEL_EXECENV_LOGGING_H
