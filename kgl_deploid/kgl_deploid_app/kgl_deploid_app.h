//
// Created by kellerberrin on 4/05/18.
//

#ifndef KGL_DEPLOID_APP_H
#define KGL_DEPLOID_APP_H


#include "kgl_genome_types.h"
#include "kgl_exec_env.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


// Holds the kgl_deploid Arguments.
struct DeploidArgs {

  std::string workDirectory{"./"};
  std::string logFile{"kgl_deploid.log"};
  int max_error_count{1000};
  int max_warn_count{1000};
  bool verbose{false};
  int argc;
  char const ** argv;

};

// Singleton. The kgl_deploid Runtime environment.
class DeploidExecEnv {

public:

  DeploidExecEnv()=delete;
  ~DeploidExecEnv()=delete;

// The following 4 static members are required for all applications.
  static constexpr const char* VERSION = "0.1";
  static constexpr const char* MODULE_NAME = "kgl_deploid";
  static void executeApp(); // Application mainline.
  static bool parseCommandLine(int argc, char const ** argv);  // Parse command line arguments.

private:

  static const DeploidArgs& getArgs();
  static DeploidArgs args_;

};



} //  organization level namespace
}  // project level namespace



#endif //KGL_DEPLOID_APP_H
