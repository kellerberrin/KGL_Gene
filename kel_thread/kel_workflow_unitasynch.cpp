//
// Created by kellerberrin on 1/01/23.
//

#include "kel_workflow_unittest.h"

namespace kel = kellerberrin;



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kel::AsynchQueueUnitTest::asynchInputWork( std::shared_ptr<AsynchMedQueue> med_queue
                                              , InputType input_item) {

  if (not input_item) {

    ExecEnv::log().info("Asynch input objects processed before stop token: {}", input_count_.load());
    med_queue->push(nullptr);
    return;

  }

  ++input_count_;

  if (input_count_ % report_iterations == 0) {

    ExecEnv::log().info("Asynch Input Objects processed: {}", input_count_.load());

  }

  auto med_ptr = std::make_unique<kel::IntermediateObject>();
  volatile u_int64_t work_count{0};
  for (size_t i = 0; i < work_iterations; ++i) {

    // Random Work
    work_count = input_item->count_ + 1;

  }

  med_ptr->count_ = work_count;
  med_ptr->count_ = input_item->count_;

  // Queue to the intermediate queue.
  med_queue->push(std::move(med_ptr));

}


void kel::AsynchQueueUnitTest::asynchIntermediateWork(  std::shared_ptr<const AsynchMedQueue> med_queue
                                                      , std::shared_ptr<AsynchOutQueue> output_queue
                                                      , IntermediateType intermediate_item) {



  if (not intermediate_item) {

    ExecEnv::log().info("Asynch Intermediate Objects processed before stop token: {}, Queue size: {}", intermediate_count_.load(), med_queue->ObjectQueue().size());
    output_queue->push(nullptr);
    return;

  }

  ++intermediate_count_;

  if (intermediate_count_ % report_iterations == 0) {

    ExecEnv::log().info("Asynch Intermediate Objects processed: {}", intermediate_count_.load());

  }

  auto out_ptr = std::make_unique<kel::OutputObject>();
  volatile u_int64_t work_count{0};
  for (size_t i = 0; i < work_iterations; ++i) {
    // Random Work
    work_count = intermediate_item->count_ + 1;

  }

  out_ptr->count_ = work_count;
  out_ptr->count_ = intermediate_item->count_;

  // Queue to the output queue.
  output_queue->push(std::move(out_ptr));

}



void kel::AsynchQueueUnitTest::pushInput(std::shared_ptr<AsynchInQueue> input_queue) {

  std::vector<InputType> input_vector;
  input_vector.reserve(iterations + 1);
  // Push input objects onto the input queue
  for (size_t i = 1; i <= iterations; ++i) {

    auto input_ptr = std::make_unique<InputObject>();
    input_ptr->count_ = i;
    input_vector.push_back(std::move(input_ptr));
    if (i % report_iterations == 0) {

      ExecEnv::log().info("Input Objects created: {}", i);

    }

  }
  ExecEnv::log().info("Finished creating input objects: {}", iterations);
  // Stop the input processing by pushing a stop token.
  size_t objects{0};
  for (auto& input_ptr : input_vector) {

    input_queue->push(std::move(input_ptr));
    ++objects;

  }
  // Stop the queue.
  input_queue->push(nullptr);

  ExecEnv::log().info("Pushed all objects: {} onto input work queue.", objects);

}


void kel::AsynchQueueUnitTest::asynchOutputWork(std::shared_ptr<const AsynchOutQueue> output_queue, OutputType output_item) {


  if (output_item) {

    ++output_count_;
    check_sum_ += output_item->count_;

    if (output_count_ % report_iterations == 0) {

      ExecEnv::log().info("Asynch Output Objects processed: {}", output_count_.load());

    }

    volatile u_int64_t work_count{0};
    for (size_t i = 0; i < work_iterations; ++i) {
      // Random Work
      work_count = output_item->count_ + 1;

    }

    output_item->count_ = work_count;

  } else {

    ExecEnv::log().info("Asynch Final - Move - Out Queue Items: {}, CheckSum: {}, Queue Size: {}", output_count_.load(), check_sum_, output_queue->ObjectQueue().size());

  }

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kel::AsynchQueueUnitTest::asynchMoveable() {

  kel::AsynchQueueUnitTest work_functions;

  // Create the 3 work queues
  // The input queue.
  auto input_queue_impl_ptr = std::make_unique<AsynchQueueType<InputType>>(high_tide, low_tide, "Input_Queue", mon_freq_ms);
  auto input_queue = std::make_shared<AsynchInQueue>(nullptr, std::move(input_queue_impl_ptr));

  // The middle queue.
  auto intermediate_queue_impl_ptr = std::make_unique<AsynchQueueType<IntermediateType>>(high_tide, low_tide, "Intermediate_Queue", mon_freq_ms);
  auto intermediate_queue = std::make_shared<AsynchMedQueue>(nullptr, std::move(intermediate_queue_impl_ptr));

  // The output queue
  auto output_queue_impl_ptr = std::make_unique<AsynchQueueType<OutputType>>(high_tide, low_tide, "Output_Queue", mon_freq_ms);
  auto output_queue = std::make_shared<AsynchOutQueue>(nullptr, std::move(output_queue_impl_ptr));

  // The input queue.
//  auto input_queue = std::make_shared<AsynchInQueue>(nullptr);

  // The middle queue.
//  auto intermediate_queue = std::make_shared<AsynchMedQueue>(nullptr);

  // The output queue
//  auto output_queue = std::make_shared<AsynchOutQueue>(nullptr);


  // Assign work functions to the queues.
  input_queue->registerProcessingFn(input_thread_count, &kel::AsynchQueueUnitTest::asynchInputWork, &work_functions, intermediate_queue);
  intermediate_queue->registerProcessingFn(intermediate_thread_count, &kel::AsynchQueueUnitTest::asynchIntermediateWork, &work_functions, intermediate_queue, output_queue);
  output_queue->registerProcessingFn(output_thread_count, &kel::AsynchQueueUnitTest::asynchOutputWork, &work_functions, output_queue);

  // Asynchronously add objects to the beginning of the linked queues.
  std::thread input_thread(&kel::AsynchQueueUnitTest::pushInput, &work_functions, input_queue);
  input_thread.join();

  // Wait until the queues have stopped processing.
  output_queue->waitUntilStopped();
  input_queue->waitUntilStopped();
  intermediate_queue->waitUntilStopped();

  ExecEnv::log().info("Asynch Workflow Unit Test exits. Input Shared Count: {}, Intermediate Count: {}, Output Count: {}"
                      , input_queue.use_count(), intermediate_queue.use_count(), output_queue.use_count());

}

