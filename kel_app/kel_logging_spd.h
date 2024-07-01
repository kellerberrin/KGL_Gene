//
// Created by kellerberrin on 18/06/24.
//

#ifndef KEL_LOGGING_SPD_H
#define KEL_LOGGING_SPD_H


//
// Created by kellerberrin on 20/10/23.
//


#include <spdlog/spdlog.h>  // Implement the logger using the spdlog library
#include "spdlog/sinks/basic_file_sink.h"

#include <iostream>

#include "kel_logging.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin {   //  organization level namespace


class SPDLoggerImpl {

public:

  SPDLoggerImpl(const std::string &module, const std::string &log_file);
  ~SPDLoggerImpl() = default;

  SPDLoggerImpl(const SPDLoggerImpl &) = delete;
  SPDLoggerImpl(SPDLoggerImpl &&) = delete;
  SPDLoggerImpl &operator=(const SPDLoggerImpl &) = delete;

  // These functions simply re-direct to the PIMPL implementation object.
  // These functions simply re-direct to the PIMPL implementation object.
  void formatImpl(std::string&& formatted_string, ExecEnvLogger::LoggerSeverity severity) noexcept;
  void locationImpl(const LogFormatLocation& format_location,
                    const std::string& formatted_string,
                    ExecEnvLogger::LoggerSeverity severity) noexcept;

private:

  std::unique_ptr<spdlog::logger> log_impl_ptr_;

  static constexpr const char *SPDLOG_DEFAULT_FORMAT{"%+"};  // Default format from spdlog

  [[nodiscard]] static spdlog::source_loc spdlogLocation(const std::source_location &location);
  [[nodiscard]] static spdlog::level::level_enum spdlogSeverity(ExecEnvLogger::LoggerSeverity severity);


};


} // Namespace.




#endif //KEL_LOGGING_SPD_H
