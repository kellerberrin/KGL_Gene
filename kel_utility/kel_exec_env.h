//
// Created by kellerberrin on 30/09/17.
//

#ifndef KGL_EXEC_ENV_H
#define KGL_EXEC_ENV_H

#include <string>
#include <memory>
#include "kel_logging.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Singleton. This class sets up the application runtime environment as a series of static variables
// and member functions. The class is never instantiated and is the first and only statement in main() (see kgl_main.cc).
// Only uses C++ 14 features so that application logging can be provided to Cuda *.cu source files.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ImplLogger;

namespace kellerberrin {   //  organization level namespace

class ExecEnv {

public:

  ExecEnv()=delete;
  ~ExecEnv()=delete;


  // Definition of this template application function is in "kel_exec_env_app.h"
  template<class Environment> static int runApplication(int argc, char const ** argv);

  static const std::string& commandLine() { return command_line_; }
  static Logger& log() { return *log_ptr_; }
  static ImplLogger& impllog() { return *impl_log_ptr_; }

  static void ctrlC(int);
  static void getCommandLine(int argc, char const ** argv);
  static void createLogger(const std::string& module,
                           const std::string& log_file,
                           int max_error_message,
                           int max_warning_messages);

private:

  static std::string command_line_;
  static std::unique_ptr<Logger> log_ptr_;
  static std::unique_ptr<ImplLogger> impl_log_ptr_;


};




}   // end namespace

#endif //KGL_EXEC_ENV_H
