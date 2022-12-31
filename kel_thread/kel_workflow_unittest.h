//
// Created by kellerberrin on 31/12/22.
//

#ifndef KEL_WORKFLOW_UNITTEST_H
#define KEL_WORKFLOW_UNITTEST_H




#include "kel_workflow_queue.h"
#include "kel_ordered_workflow_queue.h"
#include "kel_exec_env_app.h"
#include "kel_logging.h"

#include <map>


namespace kellerberrin {  //  organization level namespace


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


using InputType = std::unique_ptr<InputObject>;
using OutputType = std::unique_ptr<OutputObject>;
using IntermediateType = std::unique_ptr<IntermediateObject>;

template<typename T> using QueueType = MtQueue<T>;
//template <typename T> using QueueType = kel::BoundedMtQueue<T>;


//using InQueue = kel::WorkflowQueue<InputType, QueueType>;
//using MedQueue = kel::WorkflowQueue<IntermediateType, QueueType>;
//using OutQueue = kel::WorkflowQueue<OutputType, QueueType>;

template<typename T> using OrderedQueueType = BoundedMtQueue<std::pair<WorkFlowObjectCounter, T>>;
//using OrderedQueueType = kel::BoundedMtQueue;
using InQueue = WorkflowOrderedQueue<InputType, IntermediateType, OrderedQueueType>;
using MedQueue = WorkflowOrderedQueue<IntermediateType, OutputType, OrderedQueueType>;
//using InQueue = kel::WorkflowOrderedQueue<InputType, IntermediateType>;
//using MedQueue = kel::WorkflowOrderedQueue<IntermediateType, OutputType>;


const size_t iterations = 10000000;
const size_t report_iterations = 100000;
const size_t work_iterations = 1000000; // work for each thread


const size_t high_tide = 5000;
const size_t low_tide = 1000;
const size_t mon_freq_ms = 100;

const size_t input_thread_count = 100;
const size_t intermediate_thread_count = 100;
const size_t output_thread_count = 20;


class OrderedWorkFunctions {

public:

  OrderedWorkFunctions() = default;
  ~OrderedWorkFunctions() = default;

  static void unit_test_moveable();

  [[nodiscard]] IntermediateType interqueue_work_fn(InputType input_item);

  [[nodiscard]] OutputType outqueue_work_fn(IntermediateType intermediate_item);

  void input_to_ouput_thread(std::shared_ptr<InQueue> input_queue, std::shared_ptr<MedQueue> med_queue);

  void push_input_thread(std::shared_ptr<InQueue> input_queue);

  void retrieve_output_thread(std::shared_ptr<MedQueue> output_queue);


};



} // end namespace.

#endif //KEL_WORKFLOW_UNITTEST_H
