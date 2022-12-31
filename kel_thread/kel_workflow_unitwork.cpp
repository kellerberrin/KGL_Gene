//
// Created by kellerberrin on 31/12/22.
//

#include "kel_workflow_unittest.h"

namespace kel = kellerberrin;



kel::IntermediateType kel::OrderedWorkFunctions::interqueue_work_fn(InputType input_item) {

  auto med_ptr = std::make_unique<kel::IntermediateObject>();
  volatile u_int64_t work_count{0};
  for (size_t i = 0; i < work_iterations; ++i) {

    // Random Work
    work_count = input_item->count_ + 1;

  }

  med_ptr->count_ = work_count;
  med_ptr->count_ = input_item->count_;

  return med_ptr;

}


kel::OutputType kel::OrderedWorkFunctions::outqueue_work_fn(IntermediateType intermediate_item) {


  auto out_ptr = std::make_unique<kel::OutputObject>();
  volatile u_int64_t work_count{0};
  for (size_t i = 0; i < work_iterations; ++i) {
    // Random Work
    work_count = intermediate_item->count_ + 1;

  }

  out_ptr->count_ = work_count;
  out_ptr->count_ = intermediate_item->count_;

  return out_ptr;


}


void kel::OrderedWorkFunctions::input_to_ouput_thread(std::shared_ptr<InQueue> input_queue, std::shared_ptr<MedQueue> med_queue) {

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


void kel::OrderedWorkFunctions::push_input_thread(std::shared_ptr<InQueue> input_queue) {

  std::atomic<size_t> object_counter{0};
  // Push input objects onto the input queue
  for (size_t i = 1; i <= iterations; ++i) {

    auto input_ptr = std::make_unique<InputObject>();
    input_ptr->count_ = i;
    input_queue->push(std::move(input_ptr));
    if (i % report_iterations == 0) {

      ExecEnv::log().info("Input Objects processed: {}", i);

    }

  }
  // Stop the input processing by pushing a stop token.
  input_queue->push(nullptr);
  ExecEnv::log().info("Input Objects All Queued");


}


void kel::OrderedWorkFunctions::retrieve_output_thread(std::shared_ptr<MedQueue> output_queue) {

  size_t out_size{0};
  size_t check_sum{0};

  auto out_obj = output_queue->waitAndPop();
  while (out_obj) {

    ++out_size;
    check_sum += out_obj->count_;
    if (out_size != out_obj->count_) {

      ExecEnv::log().info("Expected Object Id: {} Actual Object Id: {}", out_size, out_obj->count_);
      break;

    }
    if (out_size % report_iterations == 0) {

      ExecEnv::log().info("Output Objects processed: {}", out_size);

    }

    out_obj = output_queue->waitAndPop();

  }

  ExecEnv::log().info("Final - Move - Out Queue Items: {}, CheckSum: {}", out_size, check_sum);

}


void kel::OrderedWorkFunctions::unit_test_moveable() {

  kel::OrderedWorkFunctions work_functions;
  // Create the 3 work queues
  //  auto input_queue_impl_ptr = std::make_unique<QueueType<InputType>>(high_tide, low_tide, "Input_Queue", mon_freq_ms);
  //  auto input_queue = std::make_shared<InQueue>(nullptr, std::move(input_queue_impl_ptr));

  //  auto intermediate_queue_impl_ptr = std::make_unique<QueueType<IntermediateType>>(high_tide, low_tide, "Intermediate_Queue", mon_freq_ms);
  //  auto intermediate_queue = std::make_shared<MedQueue>(nullptr, std::move(intermediate_queue_impl_ptr));

  //  auto output_queue_impl_ptr = std::make_unique<QueueType<OutputType>>(high_tide, low_tide, "Output_Queue", mon_freq_ms);
  //  auto output_queue = std::make_shared<OutQueue>(nullptr, std::move(output_queue_impl_ptr));

  auto input_queue_impl_ptr = std::make_unique<kel::OrderedQueueType<kel::InputType>>(high_tide, low_tide, "Input_Queue", mon_freq_ms);
  auto input_queue = std::make_shared<kel::InQueue>(nullptr, nullptr, std::move(input_queue_impl_ptr));
  auto intermediate_queue_impl_ptr = std::make_unique<kel::OrderedQueueType<kel::IntermediateType>>(high_tide, low_tide, "Intermediate_Queue", mon_freq_ms);
  auto intermediate_queue = std::make_shared<kel::MedQueue>(nullptr, nullptr, std::move(intermediate_queue_impl_ptr));
  //    auto output_queue = std::make_shared<OutQueue>(nullptr);


  // Work functions.
  input_queue->registerProcessingFn(input_thread_count, &kel::OrderedWorkFunctions::interqueue_work_fn, &work_functions);
  intermediate_queue->registerProcessingFn(intermediate_thread_count, &kel::OrderedWorkFunctions::outqueue_work_fn, &work_functions);

  // Asynchronously add objects to the beginning of the linked queues.
  std::thread input_thread(&kel::OrderedWorkFunctions::push_input_thread, &work_functions, input_queue);
  // Asynchronously transfer objects between queues (including stop token).
  std::thread spool_thread(&kel::OrderedWorkFunctions::input_to_ouput_thread, &work_functions, input_queue, intermediate_queue);
  // Retrieve objects from the out queue until we encounter a stop token.
  work_functions.retrieve_output_thread(intermediate_queue);

  //  kel::ExecEnv::log().info("Final - Move - Out Queue Size: {}", intermediate_queue->ObjectQueue().size());
  //  kel::ExecEnv::log().info("Final - Move - In Queue Size: {}", input_queue->ObjectQueue().size());

  input_thread.join();
  spool_thread.join();

}

