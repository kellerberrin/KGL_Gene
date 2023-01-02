//
// Created by kellerberrin on 14/12/22.
//

#ifndef KEL_WORKFLOW_ASYNCH_H
#define KEL_WORKFLOW_ASYNCH_H

#include <functional>
#include <vector>
#include <map>
#include <set>
#include <thread>

#include "kel_mt_queue.h"
#include "kel_bound_queue.h"


namespace kellerberrin {  //  organization level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A threaded workflow for std::move constructable objects (std::unique_ptr).
// There is no guarantee that objects are processed in the same order as they were pushed onto the queue.
// The supplied processing function must be able to handle the stop token (generally a null pointer).
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class AsynchWorkflowState { ACTIVE, STOPPED};

template<typename QueuedObj, template <typename> typename Queue = MtQueue>
requires (std::move_constructible<QueuedObj> && std::equality_comparable<QueuedObj>)
class WorkflowAsynchQueue
{

public:

  using WorkProc = std::function<void(QueuedObj)>;

  WorkflowAsynchQueue(QueuedObj stop_token, std::unique_ptr<Queue<QueuedObj>> queue_ptr = std::make_unique<Queue<QueuedObj>>())
    : stop_token_(std::move(stop_token))
    , queue_ptr_(std::move(queue_ptr)) {}
  ~WorkflowAsynchQueue() {
    // Deactivate and join all the workflow threads.
    stopProcessing();
    // Explicitly remove the workflow lambda to prevent circular references from any lambda pointer arguments.
    workflow_callback_ = nullptr;

  }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe.
  // The callback lambda is not mutable and great care (thread safe!) must be taken when modifying the arguments within the supplied function.
  // If the work function is a non-static class member function then the first ...args should be a pointer (MyClass* this) to the class instance.
  template<typename F, typename... Args>
  void activateWorkflow(size_t threads, F&& f, Args&&... args) noexcept
  {

    workflow_callback_ = [f, args...](QueuedObj t)->void { std::invoke(f, args..., std::move(t)); };
    queueThreads(threads);
    workflow_state_ = AsynchWorkflowState::ACTIVE;

  }

  // Queue objects onto the workflow queue.
  // Can be concurrently called by multiple producer threads.
  // If the underlying object queue is a bounded (load balancing) queue then this will block
  // until space becomes available on the queue. Pushing a stop token onto the workflow queue
  // shuts down all the active threads and the workflow queue is in a STOPPED condition.
  void push(QueuedObj input_obj) {

    queue_ptr_->push(std::move(input_obj));

  }

  // Remove an object from the underlying workflow object queue.
  // Blocks the calling thread if the queue is empty.
  // Can be concurrently called by multiple consumer threads.
  [[nodiscard]] QueuedObj waitAndPop() {

    return queue_ptr_->waitAndPop();

  }

  // Underlying object queue access routine.
  [[nodiscard]] const Queue<QueuedObj>& ObjectQueue() const { return *queue_ptr_; }

  // Workflow is STOPPED if there are no active work threads.
  // The Workflow queue is STOPPED on object creation and before registering the workflow function and creating active threads.
  // The Workflow queue is also STOPPED after a stop token is placed on the object queue and there are no longer any active threads.
  [[nodiscard]] AsynchWorkflowState workflowState() const { return workflow_state_; }

  // Calling thread(s) wait on a std::condition_variable until the queue is STOPPED.
  void waitUntilStopped() const {

    {
      std::unique_lock<std::mutex> lock(mutex_);
      stopped_condition_.wait(lock, [this]{ return workflow_state_ == AsynchWorkflowState::STOPPED; });
    }
    stopped_condition_.notify_one(); // In case multiple threads are blocked.

  }


private:

  QueuedObj stop_token_;
  std::unique_ptr<Queue<QueuedObj>> queue_ptr_;
  std::vector<std::thread> threads_;
  std::atomic<uint32_t> active_threads_{0};
  WorkProc workflow_callback_;
  AsynchWorkflowState workflow_state_{ AsynchWorkflowState::STOPPED };
  mutable std::mutex mutex_;
  mutable std::condition_variable stopped_condition_;


  void queueThreads(size_t threads)
  {

    // Remove any existing threads.
    stopProcessing();

    // Always have at least one worker thread queued.
    threads = threads < 1 ? 1 : threads;

    // Queue the worker threads,
    for(size_t i = 0; i < threads; ++i) {

      threads_.emplace_back(&WorkflowAsynchQueue::threadProlog, this);
      ++active_threads_;

    }

  }

  void threadProlog() {

    while(true) {

      // Get the next work item from the underlying object queue.
      QueuedObj work_item = waitAndPop();

      if (work_item == stop_token_) {

        // If not the last thread then re-queue the stop token and terminate.
        if (--active_threads_ != 0) {

          push(std::move(work_item));

        } else {

          // This is guaranteed to be the last active thread.

          // Call the workflow function with the stop token
          // to notify the processing function logic that the workflow queue will be STOPPED.
          workflow_callback_(std::move(work_item));

          // Explicitly remove the workflow function to prevent circular references from pointer arguments.
          workflow_callback_ = nullptr;

          // Notify any threads waiting on the workflow STOPPED condition.
          {
            std::lock_guard<std::mutex> lock(mutex_);
            workflow_state_ = AsynchWorkflowState::STOPPED;
          }
          stopped_condition_.notify_one();

        }
        break; // Thread terminates and can be joined.

      } else {

        // The thread performs work with the dequeued work item.
        workflow_callback_(std::move(work_item));

      }

    }

  }

  void stopProcessing() {

    // If any active threads then push the stop token onto the input queue.
    if (active_threads_ != 0) {

      push(std::move(stop_token_));

    }

    // Join all the threads
    for(auto& thread : threads_) {

      thread.join();

    }

    // Delete the deactivated threads.
    threads_.clear();

  }

};


}   // end namespace

#endif //KEL_WORKFLOW_ASYNCH_H
