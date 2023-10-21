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

#include "kel_exec_env.h"
#include "kel_utility.h"

#include <iostream>


// Define namespace alias
namespace kel = kellerberrin;

kel::ExecEnvLogger& kel::ExecEnv::log() {

  if (not log_ptr_) {
    std::cerr << "Critical - attempt to log to uninitialized logger.\n";
    std::cerr << "Program exits." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return *log_ptr_;

}


std::unique_ptr<kel::ExecEnvLogger> kel::ExecEnv::createLogger(const std::string& module,
                                                               const std::string& log_file,
                                                               size_t max_error_messages,
                                                               size_t max_warning_messages) {

  std::unique_ptr<kel::ExecEnvLogger> log_ptr;

  log_ptr = std::make_unique<ExecEnvLogger>(module, log_file);

  log_ptr->setMaxErrorMessages(max_error_messages);

  log_ptr->setMaxWarningMessages(max_warning_messages);

  return log_ptr;

}


void kel::ExecEnv::ctrlC(int) {

  ExecEnv::log().warn("Control-C. Program terminates. Output files may be corrupt. Multi-threaded code may hang.");
  std::exit(EXIT_FAILURE);

}

void kel::ExecEnv::setCommandTokens(int argc, char const ** argv) {

  for (int idx = 0; idx < argc; ++idx) {

    command_tokens_.emplace_back(argv[idx]);;

  }

}

std::string kel::ExecEnv::commandLine() {

  std::string command_line;

  for (auto const& token : command_tokens_) {

    command_line += token;
    command_line += ' ';

  }

  return command_line;

}

