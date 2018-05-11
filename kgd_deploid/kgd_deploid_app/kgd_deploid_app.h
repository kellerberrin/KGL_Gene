//
// Created by kellerberrin on 4/05/18.
//

#ifndef KGL_DEPLOID_APP_H
#define KGL_DEPLOID_APP_H


#include "kgl_genome_types.h"
#include "kgl_exec_env.h"


namespace kellerberrin {    // organization level namespace
namespace deploid {          // project level namespace


// Holds the kgd_deploid Arguments.
struct DeploidArgs {

  std::string vcfFile{"kgd_deploid.vcf"};
  std::string plafFile{"kgl_plaf.txt"};
  std::string outputTemplate{"kgl_deploid_out"};
  std::string panelFile{NOT_SPECIFIED};
  std::string excludeFile{NOT_SPECIFIED};
  std::string paintFile{NOT_SPECIFIED};
  std::string vcfOutFile{NOT_SPECIFIED};
  size_t maxStrains{5};
  size_t MCMCSamples{800};
  size_t MCMCSampleRate{5};
  long MCMCRandomSeed{1};
  double MCMCBurnRate{0.5};
  bool noPanelFlag{true};
  bool identityByDescentFlag{false};
  bool identityByDescentPainting{false};
  bool inbreedingProbabilitiesFlag{false};
  std::vector<double> initialStrainProportions;

// Logger args.
  std::string workDirectory{"./"};
  std::string logFile{"kgd_deploid.log"};
  int max_error_count{1000};
  int max_warn_count{1000};
  bool verbose{false};

  static constexpr const char* NOT_SPECIFIED = "NOT_SPECIFIED";

};

// Singleton. The kgd_deploid Runtime environment.
class DeploidExecEnv {

public:

  DeploidExecEnv()=delete;
  ~DeploidExecEnv()=delete;

// The following static members are required for all applications.
  static constexpr const char* VERSION = "0.1";
  static constexpr const char* MODULE_NAME = "kgd_deploid";
  static void executeApp(); // Application mainline.
  static bool parseCommandLine(int argc, char const ** argv);  // Parse command line arguments.

  static const DeploidArgs& getArgs();

private:

  static DeploidArgs args_;

};



}   // organization level namespace
}   // project level namespace



#endif //KGL_DEPLOID_APP_H
