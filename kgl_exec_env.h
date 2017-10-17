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

// This class sets up the application runtime environment as a series of static variables and member functions.
// The class is never instantiated and is the first statement in main() (see kgl_main.cc).

class ExecEnv {

public:

  ExecEnv()=delete;
  ~ExecEnv()=delete;


  static Logger& log();
  static void createLogger(const std::string& module, const std::string& log_file);
  static void getElpasedTime(double& Clock, double& System, double& User);

public:

  static std::unique_ptr<Logger> log_ptr_;

};



  }   // namespace genome
}   // namespace kellerberrin

#endif //KGL_EXEC_ENV_H
