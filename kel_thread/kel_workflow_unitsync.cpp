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



kel::InputType kel::SynchQueueUnitTest::synchRequeueWork(InputType input_item) {

  // Check for stop token.
  if (not input_item) {

    return nullptr;

  }

  // Make work nonsense.
  volatile u_int64_t work_count{0};
  size_t save_count = input_item->count_;
  for (size_t i = 0; i < work_iterations; ++i) {

    // Random Work
    work_count = input_item->count_ + 1;

  }

  input_item->count_ = work_count;
  input_item->count_ = save_count;

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

  // Create the workflow with a bounded tidal input queue and nullptr stop token.
  auto bounded_queue_ptr = std::make_unique<BoundedSyncInput<InputType>>();
  WorkflowSyncBounded<InputType, OutputType> workflow_bounded(nullptr, std::move(bounded_queue_ptr));

  // Simple workflow lambda increments the object count and queues input objects.
  // Workflow functions must handle stop tokens as a special case.
  auto bounded_lambda = [](InputType t)->OutputType {

    // Work functions must always handle input stop tokens.
    if (not t) return nullptr; // Just re-queue the null pointer on the output queue.
    return std::make_unique<OutputObject>(t->count_);  // Transfer the input count into the output object.

  };

  // Activate the workflow with 20 threads. This is best done after creating a bounded workflow and before any input
  // objects are processed. Otherwise, the input queue will block after it reaches high tide.
  workflow_bounded.activateWorkflow(20, bounded_lambda);

  // Now place some objects in the input queue.
  for (size_t i = 1; i <= 1000000; ++i) {

    workflow_bounded.push(std::make_unique<InputObject>(i));

  }
  // Push a stop token onto the input. This shuts down all active threads.
  workflow_bounded.push(nullptr);

  // Simple lambda to check the output queue object ordering.
  auto check_lambda = [&workflow_bounded]() -> void {

    size_t out_order{0};
    auto out_obj = workflow_bounded.waitAndPop();

  // Check for the stop token.
    while(out_obj) {

      ++out_order;
      // Check the ordering of the output objects.
      if (out_order != out_obj->count_) break;

      out_obj = workflow_bounded.waitAndPop();

    }

  };

  // Do the check asynchronously.
  std::thread check_thread(check_lambda);

  // Wait until all output objects have been checked and dequeued.
  check_thread.join();

  // Input and output queue sizes should be empty here.
  ExecEnv::log().info( "Monitored Bounded Queue, input queue size: {}, output queue size: {}"
      , workflow_bounded.inputQueue().size(), workflow_bounded.outputQueue().size());



  // An unbounded input queue workflow. The nullptr has been specified as the stop token.
  WorkflowSync<InputType, InputType> workflow_queue(nullptr);

  // Asynchronously place some objects in the queue.
  // This can be done before the queue is activated because an unbounded queue will not block.
  auto fill_lambda = [&workflow_queue]()->void{

    for (size_t i = 0; i < 1000000; ++i) {

      auto in_obj = std::make_unique<InputObject>();
      // Tag the input object with an ordering.
      in_obj->count_ = i + 1;
      workflow_queue.push(std::move(in_obj));

    }
    // Queue the stop token.
    workflow_queue.push(nullptr);

  };
  std::thread fill_thread(fill_lambda);

  // Activate the workflow with 20 threads. Just move input objects to the output queue.
  workflow_queue.activateWorkflow(20, [](InputType t)->InputType{ return t; });

  // Dequeue from the output queue.
  // The stop token will be transferred to the output queue as the final object, and we can check for that.
  size_t out_order{0};
  auto out_obj = workflow_queue.waitAndPop();
  while(out_obj) {

    ++out_order;
    // Check the ordering of the output objects.
    if (out_order != out_obj->count_) {

      ExecEnv::log().error("Unbounded Sync Queue output object ordering incorrect; expected: {}, dequeued: {}", out_order, out_obj->count_);
      break;

    }

    out_obj = workflow_queue.waitAndPop();

  }
  // Check the queue sizes. Input and Output queues should be empty.
  ExecEnv::log().info("Unbounded Sync Queue, input queue size: {}, output queue size: {}",
                      workflow_queue.inputQueue().size(), workflow_queue.outputQueue().size());

  fill_thread.join();

  kel::SynchQueueUnitTest work_functions;

  auto mt_queue_ptr = std::make_unique<kel::ReQueueType<InputType>>("Monitored Mt queue", mon_freq_ms);
  auto input_requeue = std::make_shared<kel::ReQueue>(nullptr, std::move(mt_queue_ptr));
  input_requeue->activateWorkflow(input_thread_count, &kel::SynchQueueUnitTest::synchRequeueWork, &work_functions);

  // Asynchronously add objects to the beginning of the linked queues.
  std::thread requeue_thread(&kel::SynchQueueUnitTest::requeueObjects, &work_functions, input_requeue);
  requeue_thread.join();
  // Wait on queue.
  input_requeue->waitUntilStopped();
  // And check the contents.
  ExecEnv::log().info("Synchronous requeue, input queue size: {}, output queue size: {}", input_requeue->inputQueue().size(), input_requeue->outputQueue().size());

  // Dequeue all the objects.
  InputType obj = input_requeue->waitAndPop();
  size_t expected{0};
  while(obj) {

    if (++expected != obj->count_) {

      ExecEnv::log().info("Synchronous requeue, expected: {}, actual: {}", expected, obj->count_);
      break;

    }

    obj = input_requeue->waitAndPop();

  }

  ExecEnv::log().info("Synchronous requeue before reactivate, workflow state: {}", (input_requeue->workflowState() == SyncWorkflowState::STOPPED ? "Stopped" : "Active"));
  // Reactivate the queue.
  input_requeue->activateWorkflow(input_thread_count, &kel::SynchQueueUnitTest::synchRequeueWork, &work_functions);
  ExecEnv::log().info("Synchronous requeue after reactivate");

  // Asynchronously add objects to the beginning of the linked queues.
  std::thread rerequeue_thread(&kel::SynchQueueUnitTest::requeueObjects, &work_functions, input_requeue);
  rerequeue_thread.join();
  // Wait on queue.
  input_requeue->waitUntilStopped();
  // And check the contents.
  ExecEnv::log().info("Synchronous requeue, input queue size: {}, output queue size: {}", input_requeue->inputQueue().size(), input_requeue->outputQueue().size());
  // End queue lifetime.
  input_requeue = nullptr;

  auto input_queue_impl_ptr = std::make_unique<kel::OrderedQueueType<kel::InputType>>(high_tide, low_tide, "Input_Queue", mon_freq_ms);
  auto input_queue = std::make_shared<kel::InQueue>(nullptr, std::move(input_queue_impl_ptr));
  auto intermediate_queue_impl_ptr = std::make_unique<kel::OrderedQueueType<kel::IntermediateType>>(high_tide, low_tide, "Intermediate_Queue", mon_freq_ms);
  auto intermediate_queue = std::make_shared<kel::MedQueue>(nullptr, std::move(intermediate_queue_impl_ptr));

  // Work functions.
  input_queue->activateWorkflow(input_thread_count, &kel::SynchQueueUnitTest::synchInputWork, &work_functions);
  intermediate_queue->activateWorkflow(intermediate_thread_count, &kel::SynchQueueUnitTest::synchIntermediateWork, &work_functions);

  // Asynchronously add objects to the beginning of the linked queues.
  std::thread input_thread(&kel::SynchQueueUnitTest::queueInputObjects, &work_functions, input_queue);
  // Asynchronously transfer objects between queues (including stop token).
  std::thread spool_thread(&kel::SynchQueueUnitTest::connectInputIntermediate, &work_functions, input_queue, intermediate_queue);
  // Retrieve objects from the out queue until we encounter a stop token.
  work_functions.retrieveOutputObjects(intermediate_queue);

  input_thread.join();
  spool_thread.join();

}

