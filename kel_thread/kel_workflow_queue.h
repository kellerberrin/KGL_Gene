//
// Created by kellerberrin on 14/12/22.
//

#ifndef KEL_WORKFLOW_QUEUES_H
#define KEL_WORKFLOW_QUEUES_H

#include <functional>
#include <vector>
#include <thread>

#include "kel_mt_queue.h"
#include "kel_bound_queue.h"


namespace kellerberrin {  //  organization level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A threaded workflow for std::move constructable objects (std::unique_ptr).
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename T, template <typename> typename Queue>
requires (std::move_constructible<T> && std::equality_comparable<T>)
class WorkflowQueue
{

public:

  using WorkProc = std::function<void(T)>;

  WorkflowQueue(T stop_token, std::unique_ptr<Queue<T>> queue_ptr)
    : queue_ptr_(std::move(queue_ptr))
    , stop_token_(std::move(stop_token)) {}
  ~WorkflowQueue() { stopProcessing(); }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe.
  template<typename F, typename... Args>
  void registerProcessingFn(size_t threads, F&& f, Args&&... args) noexcept
  {

//    std::forward<_Callable>(__f), std::forward<_Args>(__args)...)
    workflow_callback_ = [f, args...](T t)->void { std::invoke(f, args..., std::move(t)); };
    queueThreads(threads);

  }

  // Input stop tokens are rejected.
  void push(T&& input_obj) {

    queue_ptr_->push(std::move(input_obj));

  }

  [[nodiscard]] T waitAndPop() {

    return queue_ptr_->waitAndPop();

  }

  // Queue state access routines.
  [[nodiscard]] const Queue<T>& ObjectQueue() const { return *queue_ptr_; }


private:

  std::unique_ptr<Queue<T>> queue_ptr_;
  T stop_token_;
  std::vector<std::thread> threads_;
  std::atomic<uint32_t> active_threads_{0};
  WorkProc workflow_callback_;

  void queueThreads(size_t threads)
  {

    // Remove any existing threads.
    stopProcessing();

    // Always have at least one worker thread queued.
    threads = threads < 1 ? 1 : threads;

    // Queue the worker threads,
    for(size_t i = 0; i < threads; ++i) {

      threads_.emplace_back(&WorkflowQueue::threadProlog, this);
      ++active_threads_;

    }

  }

  void threadProlog() {

    while(true) {

      T work_item = waitAndPop();

      if (work_item == stop_token_) {

        // If the last thread then do not queue a stop token.
        if (--active_threads_ != 0) {

          push(std::move(work_item));

        } else {

          workflow_callback_(std::move(work_item));

        }
        break; // Thread terminates,

      } else {

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

    threads_.clear();

  }

};


}   // end namespace

#endif //KEL_WORKFLOW_QUEUES_H
