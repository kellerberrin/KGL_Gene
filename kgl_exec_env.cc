///
// Created by kellerberrin on 30/09/17.
//


#include <iostream>
#include <seqan/arg_parse.h>
#include "kgl_exec_env.h"


// Define namespace alias
namespace kgl = kellerberrin::genome;


// Static private member declarations.
std::unique_ptr<kgl::Logger> kgl::ExecEnv::log_ptr_;

// Public static member functions.
kgl::Logger& kgl::ExecEnv::log() { return *log_ptr_; }

void kgl::ExecEnv::createLogger(const std::string& module, const std::string& log_file) {

  kgl::ExecEnv::log_ptr_ = std::make_unique<kgl::Logger>(module, log_file);

}




