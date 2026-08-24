// Copyright 2023 Kellerberrin
//

#ifndef KGL_EXEC_ENV_H
#define KGL_EXEC_ENV_H

#include <string>
#include <memory>
#include <vector>
#include <mutex>

#include "kel_logging.h"

/// The ExecEnv class sets up the application runtime environment as a series of static
/// variables and member functions. The class is the first and only statement in main().
///
/// It provides the application logger and the parsed command line arguments. Access to the
/// logger and command-line state is internally synchronized, so reads are safe from multiple
/// threads once the environment has been initialized.
///
/// The class is not instantiable and only exposes static members.

namespace kellerberrin {   //  organization level namespace

class ExecEnv {

public:

  ExecEnv() = delete;
  ~ExecEnv() = delete;


  /// Definition of this template application function is in "kel_exec_env_app.h"
  template<class Environment> static int runApplication(int argc, char const ** argv);

  /// Returns a reference to the application logger. Terminates if the logger is uninitialized.
  /// All access is synchronized; safe to call from multiple threads once initialized.
  static ExecEnvLogger& log();
  /// Signal handler for ctrl-C. Logs a warning and terminates the program.
  /// Note: logging from a POSIX signal handler is not strictly async-signal-safe.
  static void ctrlC(int);
  /// Returns the full command line as a single string. Synchronized with setCommandTokens.
  static std::string commandLine() noexcept;
  /// Stores the raw command line tokens. Synchronized; intended to be called before threads are launched.
  static void setCommandTokens(int argc, char const ** argv) noexcept;
  /// Returns the raw command line tokens. The returned reference is valid only while no thread
  /// calls setCommandTokens().
  static const std::vector<std::string>& getCommandTokens() noexcept { return command_tokens_; }

  /// Creates and configures a new application logger.
  static std::unique_ptr<ExecEnvLogger> createLogger(const std::string& module,
                                                     const std::string& log_file,
                                                     size_t max_error_message,
                                                     size_t max_warning_messages);

private:

  inline static std::mutex env_mutex_;
  inline static std::vector<std::string> command_tokens_;
  inline static std::unique_ptr<ExecEnvLogger> log_ptr_;

};




}   // end namespace

#endif //KGL_EXEC_ENV_H
