//
// Created by kellerberrin on 30/09/17.
//


#include <iostream>
#include <seqan/arg_parse.h>
#include "kgl_exec_env.h"


// Define namespace alias
namespace kgl = kellerberrin::genome;




// Static member declarations.
kgl::ExecEnv::Args kgl::ExecEnv::args_;
std::unique_ptr<kgl::Logger> kgl::ExecEnv::log_ptr_;
kgl::Logger& kgl::ExecEnv::log() { return *log_ptr_; }
const kgl::ExecEnv::Args& kgl::ExecEnv::args() { return args_; }
constexpr const char* kgl::ExecEnv::MODULE_NAME;
constexpr const char* kgl::ExecEnv::VERSION;



