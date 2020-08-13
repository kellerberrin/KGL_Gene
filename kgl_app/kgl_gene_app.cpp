//
// Created by kellerberrin on 10/11/17.
//

#include "kgl_gene_app.h"
#include "kgl_properties.h"
#include "kgl_package.h"


namespace kgl = kellerberrin::genome;



void kgl::GeneExecEnv::executeApp() {

  // Command line arguments
  const CmdLineArgs &args = getArgs();
  // XML program runtime options, this defines the program runtime.
  const RuntimeProperties &runtime_options = getRuntimeOptions();
  // Disassemble the XML runtime into a series of data and analysis operations.
  const ExecutePackage execute_package(runtime_options, args.workDirectory);
  // Individually executes the specified XML components (the package).
  // Executes the application logic and performs requested analysis.
  execute_package.executeActive();

}
