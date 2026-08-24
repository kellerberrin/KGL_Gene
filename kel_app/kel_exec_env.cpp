// Copyright 2023 Kellerberrin


#include "kel_exec_env.h"
#include "kel_utility.h"

#include <iostream>
#include <mutex>


// Define namespace alias
namespace kel = kellerberrin;

kel::ExecEnvLogger& kel::ExecEnv::log() {

  std::lock_guard lock(env_mutex_);
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

  auto log_ptr(std::make_unique<ExecEnvLogger>(module, log_file));

  log_ptr->setMaxErrorMessages(max_error_messages);

  log_ptr->setMaxWarningMessages(max_warning_messages);

  return log_ptr;

}


void kel::ExecEnv::ctrlC(int) {

  log().warn("Control-C. Program terminates. Output files may be corrupt. Multi-threaded code may hang.");
  std::exit(EXIT_FAILURE);

}

void kel::ExecEnv::setCommandTokens(int argc, char const ** argv) noexcept {

  std::lock_guard lock(env_mutex_);
  command_tokens_.clear();
  command_tokens_.reserve(static_cast<size_t>(argc));

  for (int idx = 0; idx < argc; ++idx) {

    command_tokens_.emplace_back(argv[idx]);

  }

}

std::string kel::ExecEnv::commandLine() noexcept {

  std::lock_guard lock(env_mutex_);
  std::string command_line;

  for (auto const& token : command_tokens_) {

    command_line += token;
    command_line += ' ';

  }

  return command_line;

}
