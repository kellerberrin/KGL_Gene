//
// Created by kellerberrin on 4/05/18.
//

#ifndef KGL_PHYLOGENETIC_APP_H
#define KGL_PHYLOGENETIC_APP_H


#include "kgl_genome_types.h"
#include "kel_exec_env.h"
#include "kgl_properties.h"


namespace kellerberrin::genome {   //  organization::project level namespace


// Holds the Commandline Arguments.
struct Phylogenetic {

  std::string workDirectory{"./"};
  std::string logFile{"kgl_phylo.log"};
  std::string options_file{"runtime_options.xml"};
  int max_error_count{1000};
  int max_warn_count{1000};
  bool verbose{false};

};

// The Phylogenetic Runtime environment.
class PhylogeneticExecEnv {

public:

  PhylogeneticExecEnv()=default;
  ~PhylogeneticExecEnv()=default;

  [[nodiscard]] static const Phylogenetic& getArgs();
  [[nodiscard]] static const RuntimeProperties& getRuntimeOptions();

// The following 4 static members are required for all applications.
  inline static constexpr const char* VERSION = "0.1";
  inline static constexpr const char* MODULE_NAME = "kgl_phylo";
  static void executeApp(); // Application mainline.
  [[nodiscard]] static bool parseCommandLine(int argc, char const ** argv);  // Parse command line arguments.

private:

  inline static Phylogenetic args_;
  inline static RuntimeProperties runtime_options_;


};



} //  end namespace



#endif //KGL_PHYLOGENETIC_APP_H
