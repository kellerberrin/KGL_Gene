// Copyright 2023 Kellerberrin
//

#ifndef KEL_WORKFLOW_ASYNC_H
#define KEL_WORKFLOW_ASYNC_H

#include "kel_exec_env.h"
#include "kel_queue_mt_safe.h"
#include "kel_queue_tidal.h"

#include <atomic>
#include <concepts>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace kellerberrin {  //  organization level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A threaded workflow for std::move constructable objects (such as std::unique_ptr).
// There is no guarantee that objects (except the stop token) are processed in the same order as they were pushed onto the queue.
// The workflow is stopped and no threads are active after a stop token is pushed onto the workflow.
// The stop token is guaranteed to be the last object processed.
// This workflow can be readily concatenated together with other WorkflowAsync to provide multi thread, multi stage processing.
//
// Shutdown protocol: a worker that pops the stop token re-queues it and terminates, so the single token
// circulates until the last worker, which delivers it to the user callback (each concatenated stage therefore
// 'processes' the stop token exactly once). Exceptions thrown by the callback are caught and logged; a throwing
// callback never kills a worker thread.
//
// push() before activateWorkflow() simply queues the object; it is processed once workers exist. For a bounded
// QueueTidal the push blocks at high tide, for the unbounded QueueMtSafe it never blocks - the queue's state, not
// the workflow's activity, governs whether push() waits.
//
// The destructor pushes the stop token and joins the workers - but only if the workflow was
// activated (active_threads_ > 0). A workflow that was constructed but never activated would
// otherwise hang the destructor pushing into a full QueueTidal with no consumers. If a bounded
// queue is stalled at high tide after activation (consumers all died) the push still blocks
// forever; keep at least one consumer alive until shutdown.
//
// If activateWorkflow() fails partway through thread creation (resource exhaustion) the partially started workers
// are shut down with stop tokens and joined, and the exception is rethrown. When the stop token is move-only this
// rollback consumes it, so a workflow whose activation failed must not be re-activated (construct a new one);
// copyable stop tokens are preserved.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The objects must be std::move constructable. In addition, the objects should be comparable
// to enable the detection of a stop token placed on the workflow queue.
// The object queue can be specified as QueueMtSafe (unbounded) or QueueTidal (bounded tidal).
template<typename QueuedObj, template <typename> typename Queue = QueueMtSafe>
requires (std::move_constructible<QueuedObj> && std::equality_comparable<QueuedObj>)
class WorkflowAsync
{

  // std::move_only_function (C++23) accepts move-only callables that std::function cannot.
  using WorkProc = std::move_only_function<void(QueuedObj)>;
  // The callback is held in a shared_ptr and accessed atomically. This eliminates the
  // re-activation data race at the primitive level: a worker that has loaded the shared_ptr
  // holds a stable reference even if activateWorkflow() installs a new callback concurrently.
  using WorkProcPtr = std::shared_ptr<WorkProc>;

public:

  // The constructor requires that a stop token is specified.
  // If the Object is a pointer (a typical case is InputObject = std::unique_ptr<T>) then this will be nullptr.
  // The queue will be either a QueueMtSafe (unbounded) or BoundedMtQueue (a bounded tidal queue).
  explicit WorkflowAsync(QueuedObj stop_token, std::unique_ptr<Queue<QueuedObj>> queue_ptr = std::make_unique<Queue<QueuedObj>>())
    : stop_token_(std::move(stop_token)) , queue_ptr_(std::move(queue_ptr)) {}
  ~WorkflowAsync() noexcept {

    // Only push a stop token if there are active threads to consume it. Pushing onto an
    // inactive, full QueueTidal workflow would otherwise block the destructor forever.
    // This guards the common never-activated case; for an activated-but-stalled bounded
    // queue the caller must still keep at least one consumer alive until shutdown.
    try {

      if (active_threads_ > 0) {

        // This stop token is not processed by the workflow function.
        queue_ptr_->push(std::move(stop_token_));
        joinAndDeleteThreads();

      }

    } catch (const std::exception& e) {

      ExecEnv::log().error("~WorkflowAsync() caught exception: {}", e.what());

    } catch (...) {

      ExecEnv::log().error("~WorkflowAsync() caught unknown exception");

    }

  }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe.
  // If the work function is a non-static class member function then the first ...args should be a pointer (MyClass* this) to the class instance.
  // Calling this function on an active workflow queue will return false.
  template<typename F, typename... Args>
  bool activateWorkflow(size_t threads, F&& f, Args&&... args)
  {

    std::scoped_lock lock{activation_mutex_};

    if (active_threads_ > 0) {

      return false;

    }

    // A previous run may have stopped (active_threads_ == 0) while its last worker was still finishing the
    // stop-token callback. Join the stale threads first so workflow_callback_ is never overwritten while a
    // worker is still executing it (otherwise a data race on the callback).
    for (auto& thread : threads_) {

      thread.join();

    }
    threads_.clear();

    auto callback = std::bind_front(std::forward<F>(f), std::forward<Args>(args)...);

    // Store the callback atomically in a shared_ptr so worker threads can load a
    // consistent snapshot without data races during re-activation.
    workflow_callback_.store(std::make_shared<WorkProc>(std::move(callback)),
                             std::memory_order_release);
    queueThreads(threads);

    return true;

  }


  // Enqueue an object onto the workflow queue.
  // Objects pushed before activateWorkflow() are processed once the workers start.
  // For a bounded QueueTidal this blocks at high tide; for the unbounded QueueMtSafe it never blocks.
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
  [[nodiscard]] const Queue<QueuedObj>& objectQueue() const noexcept { return *queue_ptr_; }

private:

  QueuedObj stop_token_;
  std::unique_ptr<Queue<QueuedObj>> queue_ptr_;
  std::vector<std::jthread> threads_;
  std::atomic<uint32_t> active_threads_{0};
  std::mutex activation_mutex_;
  // Shared by all worker threads. Accessed via atomic load/store to avoid data
  // races during re-activation.
  std::atomic<WorkProcPtr> workflow_callback_;

  void queueThreads(size_t threads)
  {

    // Always have at least one worker thread.
    threads = threads < 1 ? 1 : threads;
    // Publish the worker count before any worker can pop a stop token, so the token-circulation
    // protocol in threadProlog() always decrements from an accurate count.
    active_threads_.store(threads);

    size_t started = 0;
    try {

      // Queue the worker threads,
      for (; started < threads; ++started) {

        threads_.emplace_back(&WorkflowAsync::threadProlog, this);

      }

    } catch (...) {

      // Thread creation failed partway. The partially started workers are all blocked on the empty queue,
      // so wake each of them with a stop token, then join them: leaving joinable threads behind would make
      // the vector destructor hang (or std::terminate for std::thread).
      active_threads_.store(started);
      for (size_t i = 0; i < started; ++i) {

        if constexpr (std::copy_constructible<QueuedObj>) {

          queue_ptr_->push(stop_token_);

        } else {

          // A move-only token (e.g. a null unique_ptr) still matches threadProlog()'s comparison after being
          // moved from, because the moved-from value equals itself. Copyable tokens are left untouched.
          queue_ptr_->push(std::move(stop_token_));

        }

      }
      for (auto& thread : threads_) {

        thread.join();

      }
      threads_.clear();
      active_threads_.store(0);
      throw;

    }

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
          try {
            auto callback = workflow_callback_.load(std::memory_order_acquire);
            (*callback)(std::move(work_item));
          } catch (const std::exception& e) {
            // Swallow exceptions during stop-token processing to keep shutdown safe.
            ExecEnv::log().error("WorkflowAsync stop-token processing threw exception: {}", e.what());
          } catch (...) {
            ExecEnv::log().error("WorkflowAsync stop-token processing threw an unknown exception.");
          }

        }

        return; // Thread terminates and can be joined.

      }
      else {

        // The thread performs work with the dequeued work item.
        try {
          auto callback = workflow_callback_.load(std::memory_order_acquire);
          (*callback)(std::move(work_item));
        } catch (const std::exception& e) {
          // Swallow user-callback exceptions to keep worker threads alive.
          ExecEnv::log().error("WorkflowAsync task threw exception: {}", e.what());
        } catch (...) {
          ExecEnv::log().error("WorkflowAsync task threw an unknown exception.");
        }

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
