///
// Created by kellerberrin on 30/09/17.
//

#include "kel_exec_env.h"
#include <boost/timer/timer.hpp>
#define BOOST_FILESYSTEM_NO_DEPRECATED // Recommended by boost filesystem documentation.

// Define namespace alias
namespace bt = boost::timer;


// Define namespace alias
namespace kel = kellerberrin;



void kel::ExecEnv::createLogger(const std::string& module,
                                const std::string& log_file,
                                int max_error_messages,
                                int max_warning_messages) {

  log_ptr_ = std::make_unique<Logger>(module, log_file);

  log_ptr_->SetMaxErrorMessages(max_error_messages);

  log_ptr_->SetMaxwarningMessages(max_warning_messages);

}

// Hide the boost cpu timer in an anonymous namespace.
namespace {  bt::cpu_timer cpu_timer; }

void kel::ExecEnv::getElapsedTime(double &Clock, double &User, double &System) {

  Clock = 0; User = 0; System = 0;
  bt::cpu_times elapsedtime = cpu_timer.elapsed();
  Clock = elapsedtime.wall / 1e09; // Convert from nanoseconds to seconds
  User = elapsedtime.user / 1e09;
  System = elapsedtime.system / 1e09;

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


