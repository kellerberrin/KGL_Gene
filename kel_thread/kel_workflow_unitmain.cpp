//
// Created by kellerberrin on 31/12/22.
//


#include "kel_workflow_unittest.h"


namespace kel = kellerberrin;


void unit_test() {

  kel::OrderedWorkFunctions::unit_test_moveable();

}


// Holds the Commandline Arguments.
struct CmdLineArgs {

  std::string workDirectory{"./"};
  std::string logFile{"thread_unit_test.log"};
  int max_error_count{1000};
  int max_warn_count{1000};

};

// The Runtime environment.
class ThreadExecEnv {

public:

  ThreadExecEnv() = delete;

  ~ThreadExecEnv() = delete;

  // The following 5 static members are required for all applications.
  inline static constexpr const char *VERSION = "0.1";
  inline static constexpr const char *MODULE_NAME = "Thread_Unit_Test";

  // Application mainline.
  static void executeApp() { unit_test(); }

  // No command line args, so dummy function.
  [[nodiscard]] static bool parseCommandLine(int /* argc */, char const ** /* argv */) { return true; }

  // Create application logger.
  [[nodiscard]] static std::unique_ptr<kel::Logger> createLogger() {
    return kel::ExecEnv::createLogger(MODULE_NAME, args_.logFile, args_.max_error_count, args_.max_warn_count);
  }

private:

  inline static CmdLineArgs args_;

};



/// The mainline.
int main(int argc, char const ** argv)
{

  return kel::ExecEnv::runApplication<ThreadExecEnv>(argc, argv);

}

