//
// Created by kellerberrin on 4/05/18.
//

#include "kgl_deploid_app.h"
#include "kgl_utility.h"


namespace kgl = kellerberrin::genome;


// Static private member declarations.
kgl::DeploidArgs kgl::DeploidExecEnv::args_;

// Public static member functions.
const kgl::DeploidArgs& kgl::DeploidExecEnv::getArgs() { return args_; }

// Constants for the executable.
constexpr const char* kgl::DeploidExecEnv::MODULE_NAME;
constexpr const char* kgl::DeploidExecEnv::VERSION;


void kgl::DeploidExecEnv::executeApp() {

  const DeploidArgs& args = getArgs();

  kgl::deploidMain(args.argc, args.argv);

}

