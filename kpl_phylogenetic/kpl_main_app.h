//
// Created by kellerberrin on 8/12/22.
//

#ifndef KPL_MAIN_APP_H
#define KPL_MAIN_APP_H


#include "kel_exec_env.h"
#include "kpl_strom.h"


namespace kellerberrin::phylogenetic {   //  organization level namespace


// Holds the Commandline Arguments.
struct CmdLineArgs {

  std::string workDirectory{"./"};
  std::string logFile{"kpl_phylo.log"};
  int max_error_count{1000};
  int max_warn_count{1000};

};

// The Runtime environment.
class PhyloExecEnv {

public:

  PhyloExecEnv()=delete;
  ~PhyloExecEnv()=delete;

  [[nodiscard]] inline static const CmdLineArgs& getArgs() { return args_; }

  // The following 5 static members are required for all applications.
  inline static constexpr const char* VERSION = "0.1";
  inline static constexpr const char* MODULE_NAME = "kpl_phyloTree";
  static void executeApp(); // Application mainline.
  [[nodiscard]] static bool parseCommandLine(int argc, char const ** argv);  // Parse command line arguments.
  [[nodiscard]] static std::unique_ptr<ExecEnvLogger> createLogger(); // Create application logger.

private:

  inline static CmdLineArgs args_;
  inline static Strom strom_app_;

};



} //  end namespace


#endif //KPL_MAIN_APP_H
