//
// Created by kellerberrin on 4/05/18.
//

#ifndef KGL_PHYLOGENETIC_APP_H
#define KGL_PHYLOGENETIC_APP_H


#include "kgl_properties.h"
#include "kgl_genome_types.h"

#include "kel_exec_env.h"

namespace kellerberrin::genome {   //  organization::project level namespace


// Holds the Commandline Arguments.
struct CmdLineArgs {

  std::string workDirectory{"./"};
  std::string logFile{"kgl_phylo.log"};
  std::string options_file{"runtime_options.xml"};
  size_t max_error_count{1000};
  size_t max_warn_count{1000};

};

// The Runtime environment.
class GeneExecEnv {

public:

  GeneExecEnv()=delete;
  ~GeneExecEnv()=delete;


// The following 5 static members are required for all applications.
  inline static constexpr const char* VERSION = "0.9";
  inline static constexpr const char* MODULE_NAME = "kglGene";
  static void executeApp(); // Application mainline.
  [[nodiscard]] static bool parseCommandLine(int argc, char const ** argv);  // Parse command line arguments.
  [[nodiscard]] static std::unique_ptr<ExecEnvLogger> createLogger(); // Create application logger.


  [[nodiscard]] inline static const CmdLineArgs& getArgs() { return args_; }
  [[nodiscard]] inline static const RuntimeProperties& getRuntimeOptions() { return runtime_options_; }

private:

  inline static CmdLineArgs args_;
  inline static RuntimeProperties runtime_options_;


};



} //  end namespace



#endif //KGL_PHYLOGENETIC_APP_H
