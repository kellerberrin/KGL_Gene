///
// Created by kellerberrin on 30/09/17.
//

#include "kgl_exec_env.h"


// Define namespace alias
namespace kgl = kellerberrin::genome;


// Static private member declarations.
std::unique_ptr<kgl::Logger> kgl::ExecEnv::log_ptr_;

std::string kgl::ExecEnv::command_line_;

// Public static member functions.
kgl::Logger& kgl::ExecEnv::log() { return *log_ptr_; }

void kgl::ExecEnv::createLogger(const std::string& module,
                                const std::string& log_file,
                                int max_error_messages,
                                int max_warning_messages) {

  kgl::ExecEnv::log_ptr_ = std::make_unique<kgl::Logger>(module, log_file);

  kgl::ExecEnv::log_ptr_->SetMaxErrorMessages(max_error_messages);

  kgl::ExecEnv::log_ptr_->SetMaxwarningMessages(max_warning_messages);

}


void kgl::ExecEnv::ctrlC(int) {

  ExecEnv::log().warn("Control-C. Program terminates. Output files may be corrupt. Multi-threaded code may hang.");
  std::exit(EXIT_FAILURE);

}


void kgl::ExecEnv::getCommandLine(int argc, char const ** argv) {

  std::string command_line;

  for (int idx = 0; idx < argc; ++idx) {

    command_line += argv[idx];
    command_line += ' ';

  }

  command_line_ = command_line;

}


