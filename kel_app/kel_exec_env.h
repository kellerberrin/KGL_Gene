// Copyright 2023 Kellerberrin
//

#ifndef KGL_EXEC_ENV_H
#define KGL_EXEC_ENV_H

#include <string>
#include <memory>
#include <vector>

#include "kel_logging.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This class sets up the application runtime environment as a series of static variables
// and member functions. The class is the first and only statement in main().
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin {   //  organization level namespace

class ExecEnv {

public:

  ExecEnv() = delete;
  ~ExecEnv() = delete;


  // Definition of this template application function is in "kel_exec_env_app.h"
  template<class Environment> static int runApplication(int argc, char const ** argv);

  static ExecEnvLogger& log();
  static void ctrlC(int);
  static std::string commandLine();
  static void setCommandTokens(int argc, char const ** argv);
  static const std::vector<std::string>& getCommandTokens() { return command_tokens_; }

  static std::unique_ptr<ExecEnvLogger> createLogger(const std::string& module,
                                                     const std::string& log_file,
                                                     size_t max_error_message,
                                                     size_t max_warning_messages);

private:

  inline static std::vector<std::string> command_tokens_;
  inline static std::unique_ptr<ExecEnvLogger> log_ptr_;

};




}   // end namespace

#endif //KGL_EXEC_ENV_H
