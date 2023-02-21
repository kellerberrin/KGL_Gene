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

#ifndef KEL_WORKFLOW_ASYNC_H
#define KEL_WORKFLOW_ASYNC_H

#include "kel_queue_mt_safe.h"
#include "kel_queue_tidal.h"

#include <functional>
#include <vector>
#include <thread>
#include <optional>

namespace kellerberrin {  //  organization level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A threaded workflow for std::move constructable objects (such as std::unique_ptr).
// There is no guarantee that objects (except the stop token) are processed in the same order as they were pushed onto the queue.
// The workflow is stopped and no threads are active after a stop token is pushed onto the workflow.
// The stop token is guaranteed to be the last object processed.
// This workflow can be readily 'ganged' together with other WorkflowAsync to provide multi thread, multi stage processing.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The objects must be std::move constructable. In addition, the objects should be comparable
// to enable the detection of a stop token placed on the workflow queue.
// The object queue can be specified as QueueMtSafe (unbounded) or QueueTidal (bounded tidal).
template<typename QueuedObj, template <typename> typename Queue = QueueMtSafe>
requires (std::move_constructible<QueuedObj> && std::equality_comparable<QueuedObj>)
class WorkflowAsync
{

  using WorkProc = std::function<void(QueuedObj)>;

public:

  // The constructor requires that a stop token is specified.
  // If the Object is a pointer (a typical case is InputObject = std::unique_ptr<T>) then this will be nullptr.
  // The queue will be either a QueueMtSafe (unbounded) or BondedMtQueue (a bounded tidal queue).
  explicit WorkflowAsync(QueuedObj stop_token, std::unique_ptr<Queue<QueuedObj>> queue_ptr = std::make_unique<Queue<QueuedObj>>())
    : stop_token_(std::move(stop_token)) , queue_ptr_(std::move(queue_ptr)) {}
  ~WorkflowAsync() {

       // This stop token is not processed by the workflow function.
    queue_ptr_->push(std::move(stop_token_));
    joinAndDeleteThreads();

  }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe.
  // If the work function is a non-static class member function then the first ...args should be a pointer (MyClass* this) to the class instance.
  // Calling this function on an active workflow queue will return false.
  template<typename F, typename... Args>
  bool activateWorkflow(size_t threads, F&& f, Args&&... args)
  {

    if (active_threads_ > 0) {

      return false;

    }
    workflow_callback_ = std::bind_front(std::forward<F>(f), std::forward<Args>(args)...);
    queueThreads(threads);

    return true;

  }


  // This will block if the workflow is not active.
  void push(QueuedObj input_obj) {

    queue_ptr_->push(std::move(input_obj));

  }

  // This will block until a stop token is pushed onto the queue.
  void joinAndDeleteThreads() {

    // Join all the threads
    for(auto& thread : threads_) {

      thread.join();

    }

    // Delete the joined threads.
    threads_.clear();

  }


  // Underlying object queue access routine.
  [[nodiscard]] const Queue<QueuedObj>& objectQueue() const { return *queue_ptr_; }

private:

  QueuedObj stop_token_;
  std::unique_ptr<Queue<QueuedObj>> queue_ptr_;
  std::vector<std::jthread> threads_;
  std::atomic<uint32_t> active_threads_{0};
  WorkProc workflow_callback_;

  void queueThreads(size_t threads)
  {

    // Always have at least one worker thread.
    threads = threads < 1 ? 1 : threads;
    // Queue the worker threads,
    for (size_t i = 0; i < threads; ++i) {

      threads_.emplace_back(&WorkflowAsync::threadProlog, this);

    }
    active_threads_.store(threads);

  }

  // Shuts down the all the active threads if a stop token is found on the object queue.
  void threadProlog() {

    while(true) {

      // Get the next work item from the underlying object queue.
      QueuedObj work_item = queue_ptr_->waitAndPop();

      if (work_item == stop_token_) {

        // If not the last thread then re-queue the stop token and terminate.
        if (active_threads_.fetch_sub(1) != 1) {

         queue_ptr_->push(std::move(work_item));

        } else {

          // This is guaranteed to be the last active thread.
          // Call the workflow function with the stop token.
          // The stop token is guaranteed to be the last object processed before the workflow is STOPPED.
          workflow_callback_(std::move(work_item));

        }

        return; // Thread terminates and can be joined.

      }
      else {

        // The thread performs work with the dequeued work item.
        workflow_callback_(std::move(work_item));

      }

    }

  }


};


// Convenience Asynchronous workflow alias.

// Implemented with a bounded tidal queue.
template<typename WorkObj> using BoundedAsync = QueueTidal<WorkObj>;
template<typename WorkObj> using WorkflowAsyncBounded = WorkflowAsync<WorkObj, BoundedAsync>;

// Implemented with an unbounded queue.
template<typename WorkObj> using AsyncQueue = QueueMtSafe<WorkObj>;
template<typename WorkObj> using WorkflowAsyncUnbounded = WorkflowAsync<WorkObj, AsyncQueue>;


}   // end namespace

#endif //KEL_WORKFLOW_ASYNC_H
