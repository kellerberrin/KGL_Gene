// Copyright 2023 Kellerberrin
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
// documentation files (the "Software"), to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
// and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
// OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//


#include "kel_workflow_unittest.h"


namespace kel = kellerberrin;


void unit_test() {

  kel::SynchQueueUnitTest::synchMoveable();
  kel::AsynchQueueUnitTest::asynchMoveable();

}


// Holds the Commandline Arguments.
struct CmdLineArgs {

  std::string workDirectory{"./"};
  std::string logFile{"thread_unit_test.log"};
  int max_error_count{1000};
  int max_warn_count{1000};

};

// The Static Runtime environment.
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

