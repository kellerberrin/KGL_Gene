//
// Created by kellerberrin on 10/11/17.
//

#include "kgl_gene_app.h"
#include "kgl_package.h"


namespace kgl = kellerberrin::genome;


/// Executes the gene application by reading XML options and running the analysis packages.
void kgl::GeneExecEnv::executeApp() {

  // Command line arguments
  const CmdLineArgs &args = getArgs();
  // Read the XML program options.
  const RuntimeProperties runtime_options(args.reference_directory,args.options_file,args.option_file_out, args.parsed_option_out);
  // Disassemble the XML runtime into a series of data and analysis operations.
  const ExecutePackage execute_package(runtime_options, args.reference_directory);
  // Individually executes the specified XML components (the package).
  // Executes the application logic and performs requested analysis.
  execute_package.executeActive();

}