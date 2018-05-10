//
// Created by kellerberrin on 4/05/18.
//

#include "kgd_deploid_app.h"
#include "kgl_utility.h"


namespace kgd = kellerberrin::deploid;


// Static private member declarations.
kgd::DeploidArgs kgd::DeploidExecEnv::args_;

// Public static member functions.
const kgd::DeploidArgs& kgd::DeploidExecEnv::getArgs() { return args_; }

// Constants for the executable.
constexpr const char* kgd::DeploidExecEnv::MODULE_NAME;
constexpr const char* kgd::DeploidExecEnv::VERSION;


void kgd::DeploidExecEnv::executeApp() {

  kgd::deploidMain();

}

