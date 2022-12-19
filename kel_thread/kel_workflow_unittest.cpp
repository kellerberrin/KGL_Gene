//
// Created by kellerberrin on 18/12/22.
//

#include "kel_workflow_queue.h"
#include "kel_exec_env_app.h"
#include "kel_logging.h"


namespace kel = kellerberrin;


void unit_test() {

  const size_t test_thread_count = 5;
  struct InputObject {

    size_t count_;
    std::string count_string_;

  };

  struct OutputObject {

    size_t count_;
    std::string count_string_;

  };

  using InputType = std::shared_ptr<InputObject>;
  using OutputType = std::shared_ptr<OutputObject>;


  auto callback = [](InputType /*p*/) -> OutputType { return std::make_shared<OutputObject>(); };
  kel::WorkflowThreads<InputType, OutputType> work_flow;
  work_flow.registerCallback(nullptr, nullptr, test_thread_count, callback);

  for (size_t i = 0; i < 1000; ++i) {

    work_flow.push(std::make_unique<InputObject>());

  }

  // Stop the processing threads.
  work_flow.joinThreads();

  size_t out_size{0};
  OutputType out_obj = work_flow.waitAndPop();
  while(out_obj) {

    ++out_size;
    kel::ExecEnv::log().info("Out objects processed: {}, Out Queue Size: {}", out_size, work_flow.outputQueue().size());
    out_obj = work_flow.waitAndPop();

  }

  kel::ExecEnv::log().info("Final - Out objects processed: {}, Out Queue Size: {}", out_size, work_flow.outputQueue().size());

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

  ThreadExecEnv()=delete;
  ~ThreadExecEnv()=delete;

  // The following 5 static members are required for all applications.
  inline static constexpr const char* VERSION = "0.1";
  inline static constexpr const char* MODULE_NAME = "Thread_Unit_Test";
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

