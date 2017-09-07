//
// Created by kellerberrin on 6/09/17.
//
#include <memory>

#ifndef KGL_LOGGING_H
#define KGL_LOGGING_H



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

Logger& log(void) noexcept;  // Global function returns a static instance of the logger.

class Logger {

public:

  enum class Severity { Trace, Info, Warn, Error, Critical };

  void SetLevel(Severity level) noexcept;
  void SetFormat(const std::string& message) noexcept;

  template<typename M, typename... Args> void trace(M& message, Args... args) noexcept;
  template<typename M, typename... Args> void info(M& message, Args... args) noexcept;
  template<typename M, typename... Args> void warn(M& message, Args... args) noexcept;
  template<typename M, typename... Args> void error(M& message, Args... args) noexcept;
  template<typename M, typename... Args> void critical(M& message, Args... args) noexcept;

  friend Logger& log(void) noexcept;

private:

  Logger(const std::string& log_file);
  ~Logger() = default;
  Logger(const Logger&) = delete;
  Logger(Logger&&) = delete;
  Logger& operator=(const Logger&) = delete;

  class LogImpl;
  std::unique_ptr<LogImpl> plog_impl_;

};

int test_logging();


void async_example();
void syslog_example();
void android_example();
void user_defined_example();
void err_handler_example();

namespace spd = spdlog;


}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_LOGGING_H


