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

 // Read the XML program options.
  runtime_options_.setWorkDirectory(args_.workDirectory);
  if (not runtime_options_.readProperties(args_.options_file)) {

      std::string options_file_path = Utility::filePath(args_.options_file, args_.workDirectory);
      ExecEnv::log().critical("parseCommandLine; could not read specified runtime properties file: {}", options_file_path);

  }

  // Disassemble the XML runtime into a series of data and analysis operations.
  const ExecutePackage execute_package(runtime_options_, args.workDirectory);
  // Individually executes the specified XML components (the package).
  // Executes the application logic and performs requested analysis.
  execute_package.executeActive();

}
