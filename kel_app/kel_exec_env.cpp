///
// Created by kellerberrin on 30/09/17.
//

#include "kel_exec_env.h"
#include "kel_utility.h"


// Define namespace alias
namespace kel = kellerberrin;



std::unique_ptr<kel::Logger> kel::ExecEnv::createLogger( const std::string& module,
                                                         const std::string& log_file,
                                                         int max_error_messages,
                                                         int max_warning_messages) {

  std::unique_ptr<kel::Logger> log_ptr;

  log_ptr = std::make_unique<Logger>(module, log_file);

  log_ptr->SetMaxErrorMessages(max_error_messages);

  log_ptr->SetMaxwarningMessages(max_warning_messages);

  return log_ptr;

}


void kel::ExecEnv::ctrlC(int) {

  ExecEnv::log().warn("Control-C. Program terminates. Output files may be corrupt. Multi-threaded code may hang.");
  std::exit(EXIT_FAILURE);

}

void kel::ExecEnv::getCommandLine(int argc, char const ** argv) {

  std::string command_line;

  for (int idx = 0; idx < argc; ++idx) {

    command_line += argv[idx];
    command_line += ' ';

  }

  command_line_ = command_line;

}


