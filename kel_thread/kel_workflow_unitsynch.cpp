//
// Created by kellerberrin on 31/12/22.
//

#include "kel_workflow_unittest.h"

namespace kel = kellerberrin;



kel::InputType kel::SynchQueueUnitTest::synchRequeueWork(InputType input_item) {

  // Check for stop token.
  if (not input_item) {

    return nullptr;

  }

  volatile u_int64_t work_count{0};
  for (size_t i = 0; i < work_iterations; ++i) {

    // Random Work
    work_count = input_item->count_ + 1;

  }

  input_item->count_ = work_count;

  return input_item;

}


kel::IntermediateType kel::SynchQueueUnitTest::synchInputWork(InputType input_item) {

  // Check for stop token.
  if (not input_item) {

    return nullptr;

  }

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


kel::OutputType kel::SynchQueueUnitTest::synchIntermediateWork(IntermediateType intermediate_item) {


  // Check for stop token.
  if (not intermediate_item) {

    return nullptr;

  }

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


void kel::SynchQueueUnitTest::connectInputIntermediate(std::shared_ptr<InQueue> input_queue, std::shared_ptr<MedQueue> med_queue) {

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


void kel::SynchQueueUnitTest::queueInputObjects(std::shared_ptr<InQueue> input_queue) {

  std::atomic<size_t> object_counter{0};
  // Push input objects onto the input queue
  for (size_t i = 1; i <= iterations; ++i) {

    auto input_ptr = std::make_unique<InputObject>();
    input_ptr->count_ = i;
    input_queue->push(std::move(input_ptr));
    if (i % report_iterations == 0) {

      ExecEnv::log().info("Synch Input Objects processed: {}", i);

    }

  }
  // Stop the input processing by pushing a stop token.
  input_queue->push(nullptr);
  ExecEnv::log().info("Synch Input Objects All Queued");


}


void kel::SynchQueueUnitTest::requeueObjects(std::shared_ptr<ReQueue> input_queue) {

  std::atomic<size_t> object_counter{0};
  // Push input objects onto the input queue
  for (size_t i = 1; i <= iterations; ++i) {

    auto input_ptr = std::make_unique<InputObject>();
    input_ptr->count_ = i;
    input_queue->push(std::move(input_ptr));
    if (i % report_iterations == 0) {

      ExecEnv::log().info("Synch Requeue Objects processed: {}", i);

    }

  }
  // Stop the input processing by pushing a stop token.
  input_queue->push(nullptr);
  ExecEnv::log().info("Synch Requeue Objects All Queued");


}



void kel::SynchQueueUnitTest::retrieveOutputObjects(std::shared_ptr<MedQueue> output_queue) {

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

      ExecEnv::log().info("Synch Output Objects processed: {}", out_size);

    }

    out_obj = output_queue->waitAndPop();

  }

  ExecEnv::log().info("Synch Final - Move - Out Queue Items: {}, CheckSum: {}", out_size, check_sum);

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kel::SynchQueueUnitTest::synchMoveable() {

  kel::SynchQueueUnitTest work_functions;


  auto input_requeue = std::make_shared<kel::ReQueue>(nullptr);
  input_requeue->registerProcessingFn(input_thread_count, &kel::SynchQueueUnitTest::synchRequeueWork, &work_functions);

  // Asynchronously add objects to the beginning of the linked queues.
  std::thread requeue_thread(&kel::SynchQueueUnitTest::requeueObjects, &work_functions, input_requeue);
  requeue_thread.join();
  // Wait on queue.
  input_requeue->waitUntilStopped();
  // And check the contents.
  ExecEnv::log().info("Synchronous requeue, input queue size: {}, output queue size: {}", input_requeue->inputQueue().size(), input_requeue->outputQueue().size());

  // Reactivate the queue.
  input_requeue->registerProcessingFn(input_thread_count, &kel::SynchQueueUnitTest::synchRequeueWork, &work_functions);

  // Asynchronously add objects to the beginning of the linked queues.
  std::thread rerequeue_thread(&kel::SynchQueueUnitTest::requeueObjects, &work_functions, input_requeue);
  rerequeue_thread.join();
  // Wait on queue.
  input_requeue->waitUntilStopped();
  // And check the contents.
  ExecEnv::log().info("Synchronous requeue, input queue size: {}, output queue size: {}", input_requeue->inputQueue().size(), input_requeue->outputQueue().size());


  auto input_queue_impl_ptr = std::make_unique<kel::OrderedQueueType<kel::InputType>>(high_tide, low_tide, "Input_Queue", mon_freq_ms);
  auto input_queue = std::make_shared<kel::InQueue>(nullptr, std::move(input_queue_impl_ptr));
  auto intermediate_queue_impl_ptr = std::make_unique<kel::OrderedQueueType<kel::IntermediateType>>(high_tide, low_tide, "Intermediate_Queue", mon_freq_ms);
  auto intermediate_queue = std::make_shared<kel::MedQueue>(nullptr, std::move(intermediate_queue_impl_ptr));

  // Work functions.
  input_queue->registerProcessingFn(input_thread_count, &kel::SynchQueueUnitTest::synchInputWork, &work_functions);
  intermediate_queue->registerProcessingFn(intermediate_thread_count, &kel::SynchQueueUnitTest::synchIntermediateWork, &work_functions);

  // Asynchronously add objects to the beginning of the linked queues.
  std::thread input_thread(&kel::SynchQueueUnitTest::queueInputObjects, &work_functions, input_queue);
  // Asynchronously transfer objects between queues (including stop token).
  std::thread spool_thread(&kel::SynchQueueUnitTest::connectInputIntermediate, &work_functions, input_queue, intermediate_queue);
  // Retrieve objects from the out queue until we encounter a stop token.
  work_functions.retrieveOutputObjects(intermediate_queue);

  input_thread.join();
  spool_thread.join();

}

