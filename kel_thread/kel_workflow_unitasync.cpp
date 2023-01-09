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

#include "kel_workflow_async.h"
#include <iostream>

namespace kel = kellerberrin;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Asynchronous workflows do not preserve object input order.
//
// An example of using asynchronous workflow queues. The example below uses 3 workflows connected together
// by a task object (ExampleAsyncTask). The workflows are implemented as bounded queues (WorkflowAsyncBounded).
// This is useful as it automatically load balances the CPU as the workflows perform their tasks.
// (Read the documentation on bounded tidal queues.) Another advantage is that memory usage is minimised
// as the bounded queues can only grow to a size of high-tide (default is 10000 objects).
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The toy object to be processed.
struct ExampleAsync{

  explicit ExampleAsync(size_t count) : count_(count) {}

  size_t count_{0};

};
// Moving objects around on the workflows via a pointer is generally optimal.
using ExampleAsyncType = std::unique_ptr<ExampleAsync>;

// Forward definition of the bounded workflows.
using ExampleWorkflowType = std::shared_ptr<kel::WorkflowAsyncBounded<ExampleAsyncType>>;

// The workflow function is implemented as a member of the ExampleAsyncTask object.
struct ExampleAsyncTask {

  // Do some work and place the object on the next workflow.
  // Note that the workflow object is always the last argument to the task function.
  void workflowTask (ExampleWorkflowType workflow_ptr, size_t work_iterations, ExampleAsyncType t) {

    // Always check for stop tokens.
     if (not t) {

       ++stop_tokens_;
       workflow_ptr->push(std::move(t)); // Just re-queue the stop token on the next workflow.

     }
    // Perform some work.
    workIterations(work_iterations);
    workflow_ptr->push(std::move(t));  // Transfer to the next workflow.

  }

  // Do some work and place the object in an output queue.
  // Note that the workflow object is always the last argument to the task function.
  void queueTask (size_t work_iterations, ExampleAsyncType t) {

    // Ignore the stop token.
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
  std::atomic<size_t> stop_tokens_{0};

};


void asyncExample1() {

  // We can concatenate workflows together to perform 3 stage processing.
  auto input_workflow_ptr = std::make_shared<kel::WorkflowAsyncBounded<ExampleAsyncType>>(nullptr);
  auto intermediate_workflow_ptr = std::make_shared<kel::WorkflowAsyncBounded<ExampleAsyncType>>(nullptr);
  auto output_workflow_ptr = std::make_shared<kel::WorkflowAsyncBounded<ExampleAsyncType>>(nullptr);

  ExampleAsyncTask async_task;
  // Activate the workflows with 100 threads each and assign different CPU demand work functions to the queues.
  // More threads have been allocated than are available in hardware.
  // However, CPU usage will still balance to approximately 100% of the available hardware CPU.
  // This is the advantage of using bounded queues in a workflow.
  input_workflow_ptr->activateWorkflow(100, &ExampleAsyncTask::workflowTask, &async_task, intermediate_workflow_ptr, 100000);
  intermediate_workflow_ptr->activateWorkflow(100, &ExampleAsyncTask::workflowTask, &async_task, output_workflow_ptr, 10000);
  output_workflow_ptr->activateWorkflow(100, &ExampleAsyncTask::queueTask, &async_task, 1000);

  // Asynchronously add objects to the input of the linked queues.
  auto fill_lambda = [input_workflow_ptr]() {

    for (size_t i = 1; i <= 1000000; ++i) {

      input_workflow_ptr->push(std::make_unique<ExampleAsync>(i));

    }
    // Stop token
    input_workflow_ptr->push(nullptr);

  };

  std::cout << "Asynchronous concatenated workflow example begins ....." << std::endl;

  std::thread fill_thread(fill_lambda);
  // Wait until processing is finished,
  input_workflow_ptr->waitUntilStopped();
  intermediate_workflow_ptr->waitUntilStopped();
  output_workflow_ptr->waitUntilStopped();

  std::cout << "Asynchronous concatenated workflow example ends." << std::endl;

  // Check that the workflows are empty and there are 1000000 objects in the final output queue.
  std::cout << "Input workflow size: " << input_workflow_ptr->objectQueue().size()
            << ", Intermediate workflow size: " <<  intermediate_workflow_ptr->objectQueue().size()
            << ", Output workflow size: " << output_workflow_ptr->objectQueue().size()
            << ", Final output queue size: " << async_task.output_queue_.size() << std::endl;

  // The workflows are empty and stopped, so we can join on the fill thread.
  fill_thread.join();

  if (output_workflow_ptr->objectQueue().size() != 0) {

    auto obj = output_workflow_ptr->waitAndPop();
    if (not obj) {

      std::cout << "Final object is a stop token, stop tokens moved: " << async_task.stop_tokens_.load() << std::endl;

    } else {

      std::cout << "Final object is a NOT stop token" << std::endl;

    }


  }

}

