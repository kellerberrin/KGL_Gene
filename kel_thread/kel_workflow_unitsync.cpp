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

#include "kel_workflow_sync.h"

// All the workflow objects are in the kellerberrin (Noongar - 'Hill where the Kingfisher bird is found.') namespace.
namespace kel = kellerberrin;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Example of a synchronous (preserves object ordering) workflow with an unbounded input queue.
// Workflows with unbounded input queues will accept objects indefinitely (subject to memory constraints) without blocking.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Simple example object just contains a counter to check object ordering.
struct Example1Object{

  explicit Example1Object(size_t count) : count_(count) {}

  size_t count_{0};

};
// Moving objects around on the workflow via a pointer is generally optimal.
using Example1Type = std::unique_ptr<Example1Object>;

// The actual work performed by the workflow threads.
auto task_lambda = [](Example1Type t) ->Example1Type {

  // Always check for the stop token.
  if (not t) return t; // Just place on the output queue. This will be the last object on the queue.

  // Do some example work, volatile will not be optimized away by the compiler.
  volatile u_int64_t work_count{0};
  for (size_t i = 0; i < 1000000; ++i) {

    work_count = i + 1;

  }
  work_count = work_count + t->count_;

  return t; // Requeue on the output queue.

};


void syncUnboundedExample1() {

  // An unbounded synchronised workflow. The nullptr has been specified as the stop token.
  // The same object type is used for both input and output queues.
  kel::WorkflowSync<Example1Type, Example1Type> workflow_unbounded(nullptr);

  // Activate the workflow with 20 threads. Perform the example work and move input objects to the output queue.
  workflow_unbounded.activateWorkflow(20, task_lambda);

  // Asynchronously place some objects in the queue.
  // This can also be done before the queue is activated because an unbounded queue will not block.
  auto fill_lambda = [&workflow_unbounded]()->void {

    for (size_t i = 1; i <= 1000000; ++i) {

      workflow_unbounded.push(std::make_unique<Example1Object>(i));

    }
    // Queue the stop token as the last object.
    workflow_unbounded.push(nullptr);

  };
  std::thread fill_thread(fill_lambda);

  // Dequeue from the output queue.
  // The stop token will be transferred to the output queue as the final object, and we check for that.
  size_t object_order{0};
  auto out_obj = workflow_unbounded.waitAndPop();
  while(out_obj) {

    ++object_order;
    // Check the ordering of the output objects.
    if (object_order != out_obj->count_) {

      break;

    }

    out_obj = workflow_unbounded.waitAndPop();

  }
  // Check the queue sizes. Input and Output queues should now be empty.
  kel::ExecEnv::log().info( "Unbounded Sync Queue, input queue size: {}, output queue size: {}"
                          , workflow_unbounded.inputQueue().size(), workflow_unbounded.outputQueue().size());

  fill_thread.join();

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Example of a synchronous (preserves object ordering) workflow with a bounded input queue.
// Bounded queues are useful because they automatically balance thread CPU usage between producer and consumer
// threads by blocking producer threads once the queue high-tide level is reached until the consumer threads drain
// the queue back to low-tide level.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Simple example input object contains a counter to check object ordering.
struct Example2InputObject{

  explicit Example2InputObject(size_t count) : count_(count) {}

  size_t count_{0};

};
// Moving objects around on the workflow via a pointer is generally optimal.
using Example2InputType = std::unique_ptr<Example2InputObject>;

// Simple example output object contains a string counter to check object ordering.
struct Example2OutputObject{

  explicit Example2OutputObject(size_t count)  { str_count_ = std::to_string(count); }

  [[nodiscard]] size_t count() const { return std::stoull(str_count_); }

  std::string str_count_;

};
// Moving objects around on the workflow via a pointer is generally optimal.
using Example2OutputType = std::unique_ptr<Example2OutputObject>;

// The workflow function is implemented as a member of the Example2Task object.
struct Example2Task {

  [[nodiscard]] Example2OutputType task (Example2InputType t) {

    // Work functions must always handle input stop tokens.
    if (not t) return nullptr; // Just re-queue the null pointer as a stop token on the output queue.

    ++objects_processed_;

    return std::make_unique<Example2OutputObject>(t->count_);  // Transfer the input count to the output object.

  }

  // Will be simultaneously accessed by all workflow threads so must be std::atomic.
  std::atomic<size_t> objects_processed_{0};

};


void syncBoundedExample2() {

  // Create the workflow with a bounded tidal input queue and a nullptr stop token.
  kel::WorkflowSyncBounded<Example2InputType, Example2OutputType> workflow_bounded(nullptr);

  // An instance of the Example2Task object.
  Example2Task task_object;
  // Activate the workflow with 20 threads. This is best done before any input
  // objects are processed. Otherwise, the bounded input queue will block after it reaches high tide.
  workflow_bounded.activateWorkflow(20, &Example2Task::task, &task_object);

  // Now place some objects in the input queue tagged with their position in the queue.
  for (size_t i = 1; i <= 1000000; ++i) {

    workflow_bounded.push(std::make_unique<Example2InputObject>(i));

  }
  // Push a stop token onto the input.
  // This will shut down all active threads after all the preceding input objects been processed.
  workflow_bounded.push(nullptr);

  // Lambda to dequeue objects from the output queue and check the object ordering.
  auto check_lambda = [&workflow_bounded]() -> void {

    size_t out_order{0};
    auto out_obj = workflow_bounded.waitAndPop();
    // Check for the stop token on the output queue.
    while(out_obj) {

      ++out_order;
      // Check the ordering of the output objects.
      if (out_order != out_obj->count()) break;

      // Next output object.
      out_obj = workflow_bounded.waitAndPop();

    }

  };

  // Do the sequential output object check asynchronously.
  std::thread check_thread(check_lambda);

  // Wait until all output objects have been checked and dequeued.
  check_thread.join();

  // Input and output queues are empty here.
  kel::ExecEnv::log().info( "Bounded Sync Queue, input queue size: {}, output queue size: {}, objects processed: {}"
      , workflow_bounded.inputQueue().size(), workflow_bounded.outputQueue().size(), task_object.objects_processed_.load());

}

