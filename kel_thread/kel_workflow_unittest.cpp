//
// Created by kellerberrin on 18/12/22.
//

#include "kel_workflow_queue.h"
#include "kel_ordered_workflow_queue.h"
#include "kel_exec_env_app.h"
#include "kel_logging.h"

#include <map>

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

const size_t input_thread_count = 100;
const size_t intermediate_thread_count = 100;
const size_t output_thread_count = 20;

const size_t iterations = 10000000;
const size_t report_iterations = 100000;
const size_t work_iterations = 1000000; // work for each thread

using InputType = std::unique_ptr<InputObject>;
using OutputType = std::unique_ptr<OutputObject>;
using IntermediateType = std::unique_ptr<IntermediateObject>;

template <typename T> using QueueType = kel::MtQueue<T>;
//template <typename T> using QueueType = kel::BoundedMtQueue<T>;


//using InQueue = kel::WorkflowQueue<InputType, QueueType>;
//using MedQueue = kel::WorkflowQueue<IntermediateType, QueueType>;
//using OutQueue = kel::WorkflowQueue<OutputType, QueueType>;

template <typename T> using OrderedQueueType = kel::BoundedMtQueue<std::pair<kel::WorkFlowObjectCounter, T>>;
//using OrderedQueueType = kel::BoundedMtQueue;
using InQueue = kel::WorkflowOrderedQueue<InputType, IntermediateType, OrderedQueueType>;
using MedQueue = kel::WorkflowOrderedQueue<IntermediateType, OutputType, OrderedQueueType>;
//using InQueue = kel::WorkflowOrderedQueue<InputType, IntermediateType>;
//using MedQueue = kel::WorkflowOrderedQueue<IntermediateType, OutputType>;



class WorkFunctions {

public:

  WorkFunctions() = default;
  ~WorkFunctions() = default;

  IntermediateType interqueue_work_fn(InputType input_item) {

    auto med_ptr = std::make_unique<IntermediateObject>();
    volatile u_int64_t work_count{0};
    for (size_t i = 0; i < work_iterations; ++i) {

      // Random Work
      work_count = input_item->count_ + 1;

    }

    med_ptr->count_ = work_count;
    med_ptr->count_ = input_item->count_;

    return med_ptr;

  }


  OutputType outqueue_work_fn(IntermediateType intermediate_item) {


    auto out_ptr = std::make_unique<OutputObject>();
    volatile u_int64_t work_count{0};
    for (size_t i = 0; i < work_iterations; ++i) {
      // Random Work
      work_count = intermediate_item->count_ + 1;

    }

    out_ptr->count_ = work_count;
    out_ptr->count_ = intermediate_item->count_;

    return out_ptr;


  }


  void input_to_ouput_thread(std::shared_ptr<InQueue> input_queue, std::shared_ptr<MedQueue> med_queue) {

    auto item = input_queue->waitAndPop();
    while (true) {

      if (item) {

        med_queue->push(std::move(item));

      } else {

        med_queue->push(std::move(item));
        break;

      }

      item = input_queue->waitAndPop();

    }

  }


  void push_input_thread(std::shared_ptr<InQueue> input_queue) {

    std::atomic<size_t> object_counter{0};
    // Push input objects onto the input queue
    for (size_t i = 1; i <= iterations; ++i) {

      auto input_ptr = std::make_unique<InputObject>();
      input_ptr->count_ = i;
      input_queue->push(std::move(input_ptr));
      if (i % report_iterations == 0) {

        kel::ExecEnv::log().info("Input Objects processed: {}", i);

      }

    }
    // Stop the input processing by pushing a stop token.
    input_queue->push(nullptr);
    kel::ExecEnv::log().info("Input Objects All Queued");


  }


  void retrieve_output_thread(std::shared_ptr<MedQueue> output_queue) {

    size_t out_size{0};
    size_t check_sum{0};

    auto out_obj = output_queue->waitAndPop();
    while (out_obj) {

      ++out_size;
      check_sum += out_obj->count_;
      if (out_size != out_obj->count_) {

        kel::ExecEnv::log().info("Expected Object Id: {} Actual Object Id: {}", out_size, out_obj->count_);
        break;

      }
      if (out_size % report_iterations == 0) {

        kel::ExecEnv::log().info("Output Objects processed: {}", out_size);

      }

      out_obj = output_queue->waitAndPop();

    }

    kel::ExecEnv::log().info("Final - Move - Out Queue Items: {}, CheckSum: {}", out_size, check_sum);

  }

};


void unit_test_moveable() {

  WorkFunctions work_functions;
  // Create the 3 work queues
//  auto input_queue_impl_ptr = std::make_unique<QueueType<InputType>>(high_tide, low_tide, "Input_Queue", mon_freq_ms);
//  auto input_queue = std::make_shared<InQueue>(nullptr, std::move(input_queue_impl_ptr));

//  auto intermediate_queue_impl_ptr = std::make_unique<QueueType<IntermediateType>>(high_tide, low_tide, "Intermediate_Queue", mon_freq_ms);
//  auto intermediate_queue = std::make_shared<MedQueue>(nullptr, std::move(intermediate_queue_impl_ptr));

//  auto output_queue_impl_ptr = std::make_unique<QueueType<OutputType>>(high_tide, low_tide, "Output_Queue", mon_freq_ms);
//  auto output_queue = std::make_shared<OutQueue>(nullptr, std::move(output_queue_impl_ptr));

    auto input_queue_impl_ptr = std::make_unique<OrderedQueueType<InputType>>(high_tide, low_tide, "Input_Queue", mon_freq_ms);
    auto input_queue = std::make_shared<InQueue>(nullptr, nullptr, std::move(input_queue_impl_ptr));
    auto intermediate_queue_impl_ptr = std::make_unique<OrderedQueueType<IntermediateType>>(high_tide, low_tide, "Intermediate_Queue", mon_freq_ms);
    auto intermediate_queue = std::make_shared<MedQueue>(nullptr, nullptr, std::move(intermediate_queue_impl_ptr));
//    auto output_queue = std::make_shared<OutQueue>(nullptr);


  // Work functions.
  input_queue->registerProcessingFn(input_thread_count, &WorkFunctions::interqueue_work_fn, &work_functions);
  intermediate_queue->registerProcessingFn(intermediate_thread_count, &WorkFunctions::outqueue_work_fn, &work_functions);

  // Asynchronously add objects to the beginning of the linked queues.
  std::thread input_thread(&WorkFunctions::push_input_thread, &work_functions, input_queue);
  // Asynchronously transfer objects between queues (including stop token).
  std::thread spool_thread(&WorkFunctions::input_to_ouput_thread, &work_functions, input_queue, intermediate_queue);
  // Retrieve objects from the out queue until we encounter a stop token.
  work_functions.retrieve_output_thread(intermediate_queue);

//  kel::ExecEnv::log().info("Final - Move - Out Queue Size: {}", intermediate_queue->ObjectQueue().size());
//  kel::ExecEnv::log().info("Final - Move - In Queue Size: {}", input_queue->ObjectQueue().size());

  input_thread.join();
  spool_thread.join();

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

