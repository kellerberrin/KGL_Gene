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
#include "kel_workflow_async.h"

#include <iostream>

// All the workflow objects are in the kellerberrin namespace.
namespace kel = kellerberrin;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Asynchronous workflows do not preserve object input order.
// The exception is the stop token which is guaranteed to be the last object processed.
//
// An example of using asynchronous workflow queues. The example below uses 3 workflows connected together
// by a task object (ExampleAsyncTask). The workflows are implemented as bounded queues (WorkflowAsyncBounded).
// This is useful as it automatically load balances the CPU as the workflows perform their tasks.
// (Read the documentation on bounded tidal queues.) Another advantage is that memory usage is minimised
// as the bounded queues can only grow to a size of high-tide (default is 10000 objects).
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Toy object to be processed.
class ExampleAsync{

public:

  explicit ExampleAsync(size_t count) : count_(count) {}

private:

  size_t count_{0};

};
// Moving objects around on the workflows via a pointer is generally optimal.
using ExampleAsyncType = std::unique_ptr<ExampleAsync>;

// Forward definition of the async bounded workflows.
using ExampleWorkflowType = std::shared_ptr<kel::WorkflowAsyncBounded<ExampleAsyncType>>;

// An object that encapsulates the CPU tasks that the workflows will perform on each workflow object.
struct ExampleAsyncTask {

  // Do some work and place the object on the next workflow.
  // Note that the workflow object is always the last argument to the task function.
  void workflowTask (ExampleWorkflowType workflow_ptr, size_t work_iterations, ExampleAsyncType t) {

    // Always check for stop tokens.
    if (t) {

      workIterations(work_iterations); // Perform some work.

    }

    // Forward all objects including the final stop token to the next workflow.
     workflow_ptr->push(std::move(t));

  }

  // Do some work and place the object in an output queue.
  // Note that the workflow object is always the last argument to the task function.
  void queueTask (size_t work_iterations, ExampleAsyncType t) {

    // Do not push the stop token onto the output queue.
    if (t) {

      workIterations(work_iterations);
      output_queue_.push(std::move(t));

    }

  }

  // Simple make work function to soak up some CPU cycles.
  static void workIterations(size_t iterations) {

    // Do some example work, volatile will not be optimized away by the compiler.
    volatile u_int64_t work_count{0};
    for (size_t i = 0; i < iterations; ++i) {

      work_count = i + 1;

    }
    work_count = work_count + 1;

  }

  // Will be simultaneously accessed by all output workflow threads.
  // This is OK, MtQueue and BoundedMtQueue can be accessed by multiple consumer and producer threads.
  kel::MtQueue<ExampleAsyncType> output_queue_;

};


void asyncExample1() {

  // Workflows can be concatenated together to perform multi-stage processing. Stop tokens are defined as nullptr.
  auto input_workflow_ptr = std::make_shared<kel::WorkflowAsyncBounded<ExampleAsyncType>>(nullptr);
  auto intermediate_workflow_ptr = std::make_shared<kel::WorkflowAsyncBounded<ExampleAsyncType>>(nullptr);
  auto output_workflow_ptr = std::make_shared<kel::WorkflowAsyncBounded<ExampleAsyncType>>(nullptr);

  // An instance of the task object defined above.
  ExampleAsyncTask async_task;

  // Activate the workflows with 100 threads each and assign different CPU demand work functions to the queues.
  // More threads have been allocated than are available in hardware.
  // However, CPU usage will still balance to approximately 100% of the available hardware CPU.
  // This is the advantage of using bounded queues in a workflow.
  input_workflow_ptr->activateWorkflow(100, &ExampleAsyncTask::workflowTask, &async_task, intermediate_workflow_ptr, 100000);
  intermediate_workflow_ptr->activateWorkflow(100, &ExampleAsyncTask::workflowTask, &async_task, output_workflow_ptr, 1000000);
  output_workflow_ptr->activateWorkflow(100, &ExampleAsyncTask::queueTask, &async_task, 10000);

  // Asynchronously add objects to the input of the linked queues.
  auto fill_lambda = [input_workflow_ptr]() {

    for (size_t i = 1; i <= 1000000; ++i) {

      input_workflow_ptr->push(std::make_unique<ExampleAsync>(i));

    }
    // Stop token
    input_workflow_ptr->push(nullptr);

  };
  std::thread fill_thread(fill_lambda);

  output_workflow_ptr->waitUntilStopped();
  // The workflows are empty and stopped, so we can join on the fill thread.
  fill_thread.join();

  // Check that the workflows are empty and there are 1000000 objects in the final output queue.
  std::cout << "Input workflow size: " << input_workflow_ptr->objectQueue().size()
            << ", Intermediate workflow size: " <<  intermediate_workflow_ptr->objectQueue().size()
            << ", Output workflow size: " << output_workflow_ptr->objectQueue().size()
            << ", Final output queue size: " << async_task.output_queue_.size() << std::endl;


}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Example of a synchronous (preserves object ordering) workflow with an unbounded input queue.
// Workflows with unbounded input queues will accept objects indefinitely (subject to memory constraints) without blocking.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Simple example object just contains a counter to check object ordering.
struct SyncExample1{

  explicit SyncExample1(size_t count) : count_(count) {}

  size_t count_{0};

};
// Moving objects around on the workflow via a pointer is generally optimal.
using SyncExample1Type = std::unique_ptr<SyncExample1>;

// The actual work performed by the workflow threads.
auto task_lambda = [](SyncExample1Type t) ->std::optional<SyncExample1Type> {

  // Check for the last stop token.
  if (not t) return t;
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
  kel::WorkflowSync<SyncExample1Type, SyncExample1Type> workflow_unbounded(nullptr);
  workflow_unbounded.activateWorkflow(20, task_lambda);

  // Place some objects in the workflow.
  auto fill_lambda = [&workflow_unbounded]()->void {

    for (size_t i = 1; i <= 1000000; ++i) {

      workflow_unbounded.push(std::make_unique<SyncExample1>(i));

    }
    // Stop the workflow by pushing a stop token.
    workflow_unbounded.push(nullptr);

  };
  std::thread fill_thread(fill_lambda);

  // Dequeue from the output queue until the stop token is received.
  size_t count{0};
  auto out_obj =  workflow_unbounded.waitAndPop();
  while(out_obj) {
    // Next object.
    count++;
    out_obj =  workflow_unbounded.waitAndPop();

  }

  fill_thread.join();

  // Check the queue sizes. Input and Output queues should now be empty.
  std::cout << "Unbounded Sync Queue, input queue size: " <<  workflow_unbounded.inputQueue().size()
            << ", output queue size: " << workflow_unbounded.outputQueue().size()
            << ", objects processed: " << count << std::endl;


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
struct SyncExample2{

  explicit SyncExample2(size_t count) : count_(count) {}

  size_t count_{0};

};
// Moving objects around on the workflow via a pointer is generally optimal.
using SyncExample2Type = std::unique_ptr<SyncExample2>;

// Simple example output object contains a string counter to check object ordering.
struct SyncExample2Output{

  explicit SyncExample2Output(size_t count)  { str_count_ = std::to_string(count); }

  [[nodiscard]] size_t count() const { return std::stoull(str_count_); }

  std::string str_count_;

};
// Moving objects around on the workflow via a pointer is generally optimal.
using SyncExample2OutputType = std::unique_ptr<SyncExample2Output>;

// The workflow function is implemented as a member of the SyncExample2Task object.
struct SyncExample2Task {

  [[nodiscard]] SyncExample2OutputType task (SyncExample2Type t) {

    // Work functions must always handle input stop tokens.
    if (not t) return nullptr; // Just re-queue the null pointer as a stop token on the output queue.

    ++objects_processed_;

    return std::make_unique<SyncExample2Output>(t->count_);  // Transfer the input count to the output object.

  }

  // Will be simultaneously accessed by all workflow threads so must be std::atomic.
  std::atomic<size_t> objects_processed_{0};

};


void syncBoundedExample2() {

  // Create the workflow with a bounded tidal input queue and a nullptr stop token.
  auto workflow_ptr = std::make_shared<kel::WorkflowSync<SyncExample2Type, SyncExample2OutputType>>(nullptr);

  // An instance of the SyncExample2Task object.
  SyncExample2Task task_object;
  // Activate the workflow with 20 threads. This is best done before any input
  // objects are processed. Otherwise, the bounded input queue will block after it reaches high tide.
  workflow_ptr->activateWorkflow(20, &SyncExample2Task::task, &task_object);

  auto fill_lambda = [workflow_ptr]()->void { // Now place some objects in the input queue tagged with their position in the queue.

    for (size_t i = 1; i <= 1000000; ++i) {

      workflow_ptr->push(std::make_unique<SyncExample2>(i));

    }
    // Push a stop token onto the input.
    // This will shut down all active threads after all the preceding input objects been processed.
    workflow_ptr->push(nullptr);

  };

  // Do the sequential output object check asynchronously.
  std::thread fill_thread(fill_lambda);


  // Lambda to dequeue objects from the output queue and check the object ordering.
  auto check_lambda = [workflow_ptr]() -> void {

    size_t out_order{0};
    auto out_obj = workflow_ptr->waitAndPop();
    // Check for the stop token on the output queue.
    while(out_obj) {

      ++out_order;
      // Check the ordering of the output objects.
      if (out_order != out_obj->count()) break;

      // Next output object.
      out_obj = workflow_ptr->waitAndPop();

    }

  };

  // Do the sequential output object check asynchronously.
  std::thread check_thread(check_lambda);

  // Wait until all output objects have been checked and dequeued.
  workflow_ptr->waitUntilStopped();

  // Input and output queues are empty here.
  std::cout << "Bounded Sync Queue, input queue size: " << workflow_ptr->inputQueue().size()
            << ", output queue size: " << workflow_ptr->outputQueue().size()
            << ", objects processed: " << task_object.objects_processed_.load() << std::endl;

  check_thread.join();
  fill_thread.join();

}



int main(int /* argc */, char const ** /* argv */)
{

  for (size_t i = 0; i < 1000000; ++i) {

    std::cout << "*** syncUnboundedExample1()" << std::endl;
    syncUnboundedExample1();
    std::cout << "*** syncBoundedExample2()" << std::endl;
    syncBoundedExample2();
    std::cout << "*** asyncExample1()" << std::endl;
    asyncExample1();

  }

    return EXIT_SUCCESS;

}

