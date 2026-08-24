//
// Created by kellerberrin on 18/06/24.
//

#ifndef KEL_LOGGING_STREAM_H
#define KEL_LOGGING_STREAM_H


#include "kel_logging.h"

#include <fstream>
#include <string_view>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Messages format <DD-MM-YY HH SS.000> <Module> <ANSI Status \ANSI> <message>
// Note that ANSI colouration is only used on std::cout stream messages only, not file messages.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin {   //  organization level namespace



/// StreamLoggerImpl is the PIMPL logging implementation. It writes each message to std::cout
/// (with ANSI colouration) and, when a valid log file was opened, to that file (no ANSI codes).
///
/// Instances are owned and invoked exclusively by ExecEnvLogger, which serialises all calls
/// with its message mutex, so this class need not be internally thread-safe. std::osyncstream
/// is used for std::cout to ensure complete log lines are emitted atomically when other code
/// also writes to stdout.
class StreamLoggerImpl {

public:

  /// Opens the log file (if possible) and prepares the padded module name prefix.
  StreamLoggerImpl(const std::string &module, const std::string &log_file);
  ~StreamLoggerImpl() = default;

  StreamLoggerImpl(const StreamLoggerImpl &) = delete;
  StreamLoggerImpl(StreamLoggerImpl &&) = delete;
  StreamLoggerImpl &operator=(const StreamLoggerImpl &) = delete;

  /// Writes a fully formatted message (with severity prefix) to cout and the log file.
  void formatImpl(const std::string &formatted_string, ExecEnvLogger::LoggerSeverity severity) noexcept;

  /// Writes a formatted message, prepending the source location, to cout and the log file.
  void locationImpl(const LogFormatLocation &format_location,
                    const std::string &formatted_string,
                    ExecEnvLogger::LoggerSeverity severity) noexcept;

private:

  std::ofstream output_logfile_; // Output file is non ANSI (no color); std::cout is assumed to be an ANSI terminal.
  std::string padded_module_name_; // Leading and trailing spaces are pre-added for efficiency.

  static constexpr std::string_view ANSI_YELLOW = "\033[1;33m";
  static constexpr std::string_view ANSI_RED = "\033[1;31m";
  static constexpr std::string_view ANSI_BACKGROUND_RED = "\033[1;41m";
  static constexpr std::string_view ANSI_RESET = "\033[0m";

  static constexpr std::string_view INFO = "[INFO] ";
  static constexpr std::string_view WARN = "[WARN] ";
  static constexpr std::string_view ERROR = "[ERROR] ";
  static constexpr std::string_view CRITICAL = "[CRITICAL] ";

  [[nodiscard]] static std::string severityANSI(ExecEnvLogger::LoggerSeverity severity, bool is_ansi = false);
  /// Returns the current date and time formatted as "[YYYY-MM-DD HH:MM:SS] ".
  [[nodiscard]] static std::string dateTime();
  [[nodiscard]] static std::string location(const LogFormatLocation &format_location);

  void writeMessage(const std::string &date_time,
                    ExecEnvLogger::LoggerSeverity severity,
                    const std::string &body) noexcept;


};

} // Namespace.




#endif //KEL_LOGGING_STREAM_H
