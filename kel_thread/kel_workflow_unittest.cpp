//
// Created by kellerberrin on 18/12/22.
//

#include "kel_workflow_queue.h"
#include "kel_exec_env_app.h"
#include "kel_logging.h"


namespace kel = kellerberrin;


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

const size_t high_tide = 5000;
const size_t low_tide = 1000;
const size_t mon_freq_ms = 100;

const size_t input_thread_count = 50;
const size_t intermediate_thread_count = 50;
const size_t output_thread_count = 50;

const size_t iterations = 1000000;
const size_t work_iterations = 1000000; // work for each thread

using InputType = std::unique_ptr<InputObject>;
using OutputType = std::unique_ptr<OutputObject>;
using IntermediateType = std::unique_ptr<IntermediateObject>;

// using QueueType = kel::MtQueue;

using InQueue = kel::WorkflowQueue<InputType, kel::BoundedMtQueue>;
using MedQueue = kel::WorkflowQueue<IntermediateType, kel::BoundedMtQueue>;
using OutQueue = kel::WorkflowQueue<OutputType, kel::BoundedMtQueue>;


class WorkFunctions {

public:

  WorkFunctions() = default;
  ~WorkFunctions() = default;

  void interqueue_work_fn(std::shared_ptr<MedQueue> intermediate_queue, IntermediateType intermediate_stop, InputType input_stop, InputType input_item) {

    if (input_item == input_stop) {

      intermediate_queue->push(std::move(intermediate_stop));

    } else {

      auto med_ptr = std::make_unique<IntermediateObject>();
      for (size_t i = 0; i < work_iterations; ++i) {

        // Random Work
        med_ptr->count_ += input_item->count_ ^ i;

      }

      intermediate_queue->push(std::move(med_ptr));

    }


  }


  void outqueue_work_fn(std::shared_ptr<OutQueue> output_queue, OutputType output_stop, IntermediateType intermediate_stop, IntermediateType intermediate_item) {

    if (intermediate_item == intermediate_stop) {

      output_queue->push(std::move(output_stop));

    } else {

      auto out_ptr = std::make_unique<OutputObject>();
      for (size_t i = 0; i < work_iterations; ++i) {

        // Random Work
        out_ptr->count_ += intermediate_item->count_ ^ i;

      }

      output_queue->push(std::move(out_ptr));

    }


  }

  void push_input_thread(std::shared_ptr<InQueue> input_queue) {

    // Push input objects onto the input queue
    for (size_t i = 0; i < iterations; ++i) {

      auto input_ptr = std::make_unique<InputObject>();
      input_ptr->count_ = i;
      input_queue->push(std::move(input_ptr));

    }
    // Stop the input processing by pushing a stop token.
    input_queue->push(nullptr);


  }


  void retrieve_output_thread(std::shared_ptr<OutQueue> output_queue) {

    size_t out_size{0};
    size_t check_sum{0};

    OutputType out_obj = output_queue->waitAndPop();
    while (out_obj) {

      ++out_size;
      check_sum += out_obj->count_;

      out_obj = output_queue->waitAndPop();

    }

    kel::ExecEnv::log().info("Final - Move - Out Queue Items: {}, CheckSum: {}", out_size, check_sum);

  }

};


void unit_test_moveable() {

  WorkFunctions work_functions;
  // Create the 3 work queues
  auto input_queue_impl_ptr = std::make_unique<kel::BoundedMtQueue<InputType>>(high_tide, low_tide, "Input_Queue", mon_freq_ms);
  auto input_queue = std::make_shared<InQueue>(nullptr, std::move(input_queue_impl_ptr));

  auto intermediate_queue_impl_ptr = std::make_unique<kel::BoundedMtQueue<IntermediateType>>(high_tide, low_tide, "Intermediate_Queue", mon_freq_ms);
  auto intermediate_queue = std::make_shared<MedQueue>(nullptr, std::move(intermediate_queue_impl_ptr));

  auto output_queue_impl_ptr = std::make_unique<kel::BoundedMtQueue<OutputType>>(high_tide, low_tide, "Output_Queue", mon_freq_ms);
  auto output_queue = std::make_shared<OutQueue>(nullptr, std::move(output_queue_impl_ptr));

  // Link the queues together
  input_queue->registerProcessingFn(input_thread_count, &WorkFunctions::interqueue_work_fn, &work_functions, intermediate_queue, nullptr, nullptr);
  intermediate_queue->registerProcessingFn(intermediate_thread_count, &WorkFunctions::outqueue_work_fn, &work_functions, output_queue, nullptr, nullptr);

  // Asynchronously add objects to the beginning of the linked queues.
  std::thread input_thread(&WorkFunctions::push_input_thread, &work_functions, input_queue);

  // Retrieve objects from the out queue until we encounter a stop token.
  work_functions.retrieve_output_thread(output_queue);

  kel::ExecEnv::log().info("Final - Move - Out Queue Size: {}", output_queue->inputQueue().size());
  kel::ExecEnv::log().info("Final - Move - In Queue Size: {}", input_queue->inputQueue().size());

  input_thread.join();

}



void unit_test() {

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

