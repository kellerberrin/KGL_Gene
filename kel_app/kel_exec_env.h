// Copyright 2023 Kellerberrin
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
// documentation files (the "Software"), to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
// and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
// OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
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
