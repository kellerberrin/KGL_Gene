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

#include "kel_mt_queue.h"
#include "kel_bound_queue.h"

#include <functional>
#include <vector>
#include <thread>
#include <optional>

namespace kellerberrin {  //  organization level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A threaded workflow for std::move constructable objects (such as std::unique_ptr).
// There is no guarantee that objects (except the stop token) are processed in the same order as they were pushed onto the queue.
// The workflow is STOPPED and no threads are active after a stop token is pushed onto the workflow.
// The supplied processing function must check for the stop token (a nullptr in the std::unique_ptr case).
// The stop token is guaranteed to be the last object processed.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The workflow is STOPPED if there are no active work threads.
// The workflow is STOPPED on object creation and before activating the workflow with a task function.
// The workflow is in a SHUTDOWN state if a stop token has been received.
// SHUTDOWN is a temporary state that will transition to STOPPED.
// The workflow is ACTIVE if a task function has been supplied and all threads are ready for processing.

enum class AsyncWorkflowState { ACTIVE, SHUTDOWN, STOPPED};

// The objects must be std::move constructable. In addition, the objects should be comparable
// to enable the detection of a stop token placed on the workflow queue.
// The object queue can be specified as MtQueue (unbounded) or BoundedMtQueue (bounded tidal).
template<typename QueuedObj, template <typename> typename Queue = MtQueue>
requires (std::move_constructible<QueuedObj> && std::equality_comparable<QueuedObj>)
class WorkflowBase
{

  using WorkProc = std::function<void(QueuedObj)>;

public:

  // The constructor requires that a top token is specified.
  // If the Object is a pointer (a typical case is InputObject = std::unique_ptr<T>) then this will be nullptr.
  // The queue will be either a MtQueue (unbounded) or BondedMtQueue (a bounded tidal queue).
  explicit WorkflowBase(QueuedObj stop_token, std::unique_ptr<Queue<QueuedObj>> queue_ptr = std::make_unique<Queue<QueuedObj>>())
    : stop_token_(std::move(stop_token))
    , queue_ptr_(std::move(queue_ptr)) {}
  ~WorkflowBase() {

    // If any active threads then push the stop token onto the workflow queue.
    if (workflow_state_ != AsyncWorkflowState::STOPPED) {

      // This stop token is not processed by the workflow function.
      queue_ptr_->push(std::move(stop_token_));

    }
    joinAndDeleteThreads();

  }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe.
  // The callback lambda is not mutable and great care (thread safe!) must be taken when modifying the arguments within the supplied function.
  // If the work function is a non-static class member function then the first ...args should be a pointer (MyClass* this) to the class instance.
  // This function is not multi-threaded and must be called after workflow queue creation for workflow processing to be ACTIVE.
  // If the queue has been STOPPED, this function can be called with different workflow functions and thread counts.
  // Calling this function on an active workflow queue will return false.
  template<typename F, typename... Args>
  bool activateWorkflow(size_t threads, F&& f, Args&&... args) noexcept
  {

    { // mutex scope
      std::scoped_lock<std::mutex> lock(state_mutex_);

      if (workflow_state_ != AsyncWorkflowState::STOPPED) {

        return false;

      }

    } // ~mutex scope
    joinAndDeleteThreads();
    workflow_callback_ = std::bind_front(f, args...);
    queueThreads(threads);
    {
      std::scoped_lock<std::mutex> lock(state_mutex_);
      workflow_state_ = AsyncWorkflowState::ACTIVE;
    }
    active_condition_.notify_one();

    return true;

  }


  // This will block if the workflow is not active.
  void push(QueuedObj input_obj) {

    auto in_obj_opt = pushNoBlock(std::move(input_obj));

    while (in_obj_opt.has_value()) {

      waitUntilActive();
      in_obj_opt = pushNoBlock(std::move(in_obj_opt.value()));

    }

  }

  // Does not block, the input object is returned to the caller if the workflow is not active.
  [[nodiscard]] std::optional<QueuedObj> pushNoBlock(QueuedObj input_obj) {

    // Input objects cannot be pushed onto a stopped/shutdown workflow.
    { // mutex scope
      std::scoped_lock<std::mutex> lock(state_mutex_);

      if ( workflow_state_ != AsyncWorkflowState::ACTIVE) {

        return std::make_optional<QueuedObj>(std::move(input_obj));

      }

      if (input_obj == stop_token_) {

        workflow_state_ = AsyncWorkflowState::SHUTDOWN;

      }

    } // ~mutex.

    queue_ptr_->push(std::move(input_obj));

    // The input object has been consumed.
    return std::nullopt;

  }

  // Underlying object queue access routine.
  // All const public members of MtQueue and BoundedMtQueue are thread safe.
  [[nodiscard]] const Queue<QueuedObj>& objectQueue() const { return *queue_ptr_; }

  // This function is not thread safe! A race condition potentially exists.
  // If a producer thread has pushed a stop token, another thread calling this function may still return ACTIVE.
  // Consider using waitOnStopped() or waitOnActive() instead.
  [[nodiscard]] AsyncWorkflowState workflowState() const { return workflow_state_.load(); }

  // Calling thread(s) wait on a std::condition_variable until the queue is STOPPED.
  void waitUntilStopped() const {

    {
      std::unique_lock<std::mutex> lock(state_mutex_);
      stopped_condition_.wait(lock, [this]{ return workflow_state_ == AsyncWorkflowState::STOPPED; });
    }
    stopped_condition_.notify_one(); // Multiple threads may be blocked.

  }

  // Calling thread(s) wait on a std::condition_variable until the workflow queue is ACTIVE.
  void waitUntilActive() const {

    {
      std::unique_lock<std::mutex> lock(state_mutex_);
      active_condition_.wait(lock, [this]{ return workflow_state_ == AsyncWorkflowState::ACTIVE; });
    }
    active_condition_.notify_one(); // In case multiple threads are blocked.

  }


private:

  QueuedObj stop_token_;
  std::unique_ptr<Queue<QueuedObj>> queue_ptr_;
  std::vector<std::jthread> threads_;
  std::atomic<uint32_t> active_threads_{0};
  WorkProc workflow_callback_;
  std::atomic<AsyncWorkflowState> workflow_state_{AsyncWorkflowState::STOPPED };
  mutable std::mutex state_mutex_;
  mutable std::condition_variable stopped_condition_;
  mutable std::condition_variable active_condition_;


  // Remove an object from the underlying workflow object queue.
  // Blocks the calling thread if the queue is empty.
  // Can be concurrently called by multiple consumer threads.
  [[nodiscard]] QueuedObj waitAndPop() {

    return queue_ptr_->waitAndPop();

  }

  void queueThreads(size_t threads)
  {

    // Always have at least one worker thread.
    threads = threads < 1 ? 1 : threads;
    // Queue the worker threads,
    for (size_t i = 0; i < threads; ++i) {

      threads_.emplace_back(&WorkflowBase::threadProlog, this);

    }
    active_threads_.store(threads);

  }

  // Shuts down the all the active threads if a stop token is found on the object queue.
  void threadProlog() {

    while(true) {

      // Get the next work item from the underlying object queue.
      QueuedObj work_item = waitAndPop();

      if (work_item == stop_token_) {

        // If not the last thread then re-queue the stop token and terminate.
        if (active_threads_.fetch_sub(1) != 1) {

         queue_ptr_->push(std::move(work_item));

        } else {

          // This is guaranteed to be the last active thread.
          // Call the workflow function with the stop token.
          // The stop token is guaranteed to be the last object processed before the workflow is STOPPED.
          workflow_callback_(std::move(work_item));

          {
            std::scoped_lock<std::mutex> lock(state_mutex_);
            workflow_state_ = AsyncWorkflowState::STOPPED;
          }
          stopped_condition_.notify_one();

        }

        return; // Thread terminates and can be joined.

      }
      else {

        // The thread performs work with the dequeued work item.
        workflow_callback_(std::move(work_item));

      }

    }

  }

  void joinAndDeleteThreads() {

    // Join all the threads
    for(auto& thread : threads_) {

      thread.join();

    }

    // Delete the joined threads.
    threads_.clear();

  }

};


// Convenience Asynchronous workflow typedefs.

// Implemented with a bounded tidal queue.
template<typename WorkObj> using BoundedAsync = BoundedMtQueue<WorkObj>;
template<typename WorkObj> using WorkflowAsyncBounded = WorkflowBase<WorkObj, BoundedAsync>;

// Implemented with an unbounded queue.
template<typename WorkObj> using AsyncQueue = MtQueue<WorkObj>;
template<typename WorkObj> using WorkflowAsync = WorkflowBase<WorkObj, AsyncQueue>;


}   // end namespace

#endif //KEL_WORKFLOW_ASYNC_H
