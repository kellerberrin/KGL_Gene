//
// Created by kellerberrin on 18/12/22.
//

#include "kel_workflow_queue.h"
#include "kel_exec_env_app.h"
#include "kel_logging.h"


namespace kel = kellerberrin;


const size_t test_thread_count = 5;
struct InputObject {

  size_t count_;
  std::string count_string_;

};

struct IntermediateObject {

  size_t count_;
  size_t in_count_;
  std::string count_string_;

};

struct OutputObject {

  size_t count_;
  size_t in_count_;
  std::string count_string_;

};

using InputType = std::shared_ptr<InputObject>;
using OutputType = std::shared_ptr<OutputObject>;

void unit_test_copyable() {


  auto callback = [](size_t /*i*/, InputType /*p*/) -> OutputType { return std::make_shared<OutputObject>(); };
  kel::WorkflowQueues<InputType, OutputType> work_flow;
  work_flow.registerProcessingFn(nullptr, nullptr, test_thread_count, callback, 7);

  for (size_t i = 0; i < 1000; ++i) {

    work_flow.push(std::make_unique<InputObject>());

  }

  // Stop the processing threads by pushing a stop token
  work_flow.push(nullptr);

  size_t out_size{0};
  // Dequeue the output objects until the output stop token is encountered.
  OutputType out_obj = work_flow.waitAndPop();
  while(out_obj) {

    ++out_size;
    kel::ExecEnv::log().info("Copy - Out objects processed: {}, Out Queue Size: {}", out_size, work_flow.outputQueue().size());
    out_obj = work_flow.waitAndPop();

  }

  kel::ExecEnv::log().info("Final - Copy - Out objects processed: {}, Out Queue Size: {}", out_size, work_flow.outputQueue().size());
  kel::ExecEnv::log().info("Final - Copy - In Queue Size: {}", work_flow.inputQueue().size());

}



void unit_test_moveable() {


  using InputType = std::unique_ptr<InputObject>;
  using IntermediateType = std::unique_ptr<IntermediateObject>;
  using OutputType = std::unique_ptr<OutputObject>;

  std::unique_ptr<kel::WorkflowSingle<InputType>> work_flow(std::make_unique<kel::WorkflowSingle<InputType>>(nullptr));
  std::shared_ptr<kel::WorkflowSingle<IntermediateType>> intermediate_queue(std::make_shared<kel::WorkflowSingle<IntermediateType>>(nullptr));
  std::shared_ptr<kel::WorkflowSingle<OutputType>> output_queue(std::make_shared<kel::WorkflowSingle<OutputType>>(nullptr));

  auto interqueue = [](std::shared_ptr<kel::WorkflowSingle<IntermediateType>> intermediate_queue
      , IntermediateType intermediate_stop
      , InputType input_stop
      ,InputType input_item) -> void {

    if (input_item == input_stop) {

      intermediate_queue->push(std::move(intermediate_stop));

    } else {

      intermediate_queue->push(std::make_unique<IntermediateObject>());

    }


  };

  auto outqueue = [](std::shared_ptr<kel::WorkflowSingle<OutputType>> output_queue
      , OutputType output_stop
      , IntermediateType intermediate_stop
      , IntermediateType intermediate_item) -> void {

    if (intermediate_item == intermediate_stop) {

      output_queue->push(std::move(output_stop));

    } else {

      output_queue->push(std::make_unique<OutputObject>());

    }


  };

  intermediate_queue->registerProcessingFn(test_thread_count, outqueue, output_queue, nullptr, nullptr);
  work_flow->registerProcessingFn(test_thread_count, interqueue, std::move(intermediate_queue), nullptr, nullptr);

  for (size_t i = 0; i < 1000; ++i) {

    work_flow->push(std::make_unique<InputObject>());

  }

  // Stop the processing threads.
  work_flow->push(nullptr);

  size_t out_size{0};
  OutputType out_obj = output_queue->waitAndPop();
  while(out_obj) {

    ++out_size;
    kel::ExecEnv::log().info("Move - Out objects processed: {}, Out Queue Size: {}", out_size, output_queue->inputQueue().size());
    out_obj = output_queue->waitAndPop();

  }

  kel::ExecEnv::log().info("Final - Move - Out objects processed: {}, Out Queue Size: {}", out_size, output_queue->inputQueue().size());
  kel::ExecEnv::log().info("Final - Move - In Queue Size: {}", work_flow->inputQueue().size());

}


void unit_test() {

  unit_test_copyable();
  unit_test_moveable();

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

