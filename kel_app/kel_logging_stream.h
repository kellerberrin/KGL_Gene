//
// Created by kellerberrin on 18/06/24.
//

#ifndef KEL_LOGGING_STREAM_H
#define KEL_LOGGING_STREAM_H



//
// Created by kellerberrin on 20/10/23.
//


#include "kel_logging.h"

#include <fstream>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Messages format <DD-MM-YY HH SS.000> <Module> <ANSI Status \ANSI> <message>
// Note that ANSI colouration is only used on std::cout stream messages only, not file messages.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin {   //  organization level namespace



class StreamLoggerImpl {

public:

  StreamLoggerImpl(const std::string &module, const std::string &log_file);
  ~StreamLoggerImpl() = default;

  StreamLoggerImpl(const StreamLoggerImpl &) = delete;
  StreamLoggerImpl(StreamLoggerImpl &&) = delete;
  StreamLoggerImpl &operator=(const StreamLoggerImpl &) = delete;

  void formatImpl(std::string &&formatted_string, ExecEnvLogger::LoggerSeverity severity) noexcept;

  void locationImpl(const LogFormatLocation &format_location,
                    const std::string &formatted_string,
                    ExecEnvLogger::LoggerSeverity severity) noexcept;

private:

  std::ofstream output_logfile_; // Output file is non ANSI (no color); std::cout is assumed to be an ANSI terminal.
  std::string padded_module_name_; // Leading and trailing spaces are pre-added for efficiency.

  inline static constexpr std::string ANSI_YELLOW = "\033[1;33m";
  inline static constexpr std::string ANSI_RED = "\033[1;31m";
  inline static constexpr std::string ANSI_BACKGROUND_RED = "\033[1;41m";
  inline static constexpr std::string ANSI_RESET = "\033[0m";

  inline static constexpr std::string INFO = "[INFO] ";
  inline static constexpr std::string WARN = "[WARN] ";
  inline static constexpr std::string ERROR = "[ERROR] ";
  inline static constexpr std::string CRITICAL = "[CRITICAL] ";

  inline static const std::string ANSI_INFO = INFO;
  inline static const std::string ANSI_WARN = ANSI_YELLOW + WARN + ANSI_RESET;
  inline static const std::string ANSI_ERROR = ANSI_RED + ERROR + ANSI_RESET;
  inline static const std::string ANSI_CRITICAL = ANSI_BACKGROUND_RED + CRITICAL + ANSI_RESET;

  [[nodiscard]] std::string severityANSI(ExecEnvLogger::LoggerSeverity severity, bool is_ansi) const;
  [[nodiscard]] std::string dateTime() const;
  [[nodiscard]] std::string location(const LogFormatLocation &format_location) const;


};

} // Namespace.





#endif //KEL_LOGGING_STREAM_H
