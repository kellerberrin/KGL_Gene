//
// Created by kellerberrin on 10/11/17.
//

#include "kgl_gene_app.h"
#include "kgl_package.h"


namespace kgl = kellerberrin::genome;


/// Executes the gene application by reading XML options and running the analysis packages.
void kgl::GeneExecEnv::executeApp() {

  // Get command line arguments
  const CmdLineArgs &args = getArgs();

  // Read the XML program options.
  // The XML program options are located relative to the reference directory.
  // Then Disassemble the XML runtime into a series of data and analysis operations called a package.
  const ExecutePackage execute_package( args.options_file,
                                        args.option_file_out,
                                        args.parsed_option_out ,
                                        args.reference_directory);

  // Then for each active package, execute the application logic and perform requested analysis.
  execute_package.executeActive();

}