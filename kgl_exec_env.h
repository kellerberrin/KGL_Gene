//
// Created by kellerberrin on 30/09/17.
//

#ifndef KGL_EXEC_ENV_H
#define KGL_EXEC_ENV_H

#include <string>
#include <memory>
#include "kgl_logging.h"



namespace kellerberrin {   //  organization level namespace
  namespace genome {   // project level namespace

// Singleton. This class sets up the application runtime environment as a series of static variables
// and member functions. The class is never instantiated and is the first statement in main() (see kgl_main.cc).

class ExecEnv {

public:

  ExecEnv()=delete;
  ~ExecEnv()=delete;


  static Logger& log();
  static void getElpasedTime(double& Clock, double& System, double& User);
  static std::string filePath(const std::string& file_name, const std::string& path); // utility for "path/file"
  static const std::string& commandLine() { return command_line_; }
  static void createLogger(const std::string& module,
                           const std::string& log_file,
                           int max_error_message,
                           int max_warning_messages);

protected:

  static void commandLine(const std::string command_line) { command_line_ = command_line; }

private:

  static std::string command_line_;
  static std::unique_ptr<Logger> log_ptr_;

};



  }   // namespace genome
}   // namespace kellerberrin

#endif //KGL_EXEC_ENV_H
