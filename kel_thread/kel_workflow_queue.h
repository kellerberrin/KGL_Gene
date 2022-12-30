//
// Created by kellerberrin on 14/12/22.
//

#ifndef KEL_WORKFLOW_QUEUES_H
#define KEL_WORKFLOW_QUEUES_H

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


template<typename QueuedObj, template <typename> typename Queue = MtQueue>
requires (std::move_constructible<QueuedObj> && std::equality_comparable<QueuedObj>)
class WorkflowQueue
{

public:

  using WorkProc = std::function<void(QueuedObj)>;

  WorkflowQueue(QueuedObj stop_token, std::unique_ptr<Queue<QueuedObj>> queue_ptr = std::make_unique<Queue<QueuedObj>>())
    : stop_token_(std::move(stop_token)
    , queue_ptr_(std::move(queue_ptr))) {}
  ~WorkflowQueue() { stopProcessing(); }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe.
  // If the work function is a non-static class member function then the first ...args should be a pointer (Class* this) to the class instance.
  template<typename F, typename... Args>
  void registerProcessingFn(size_t threads, F&& f, Args&&... args) noexcept
  {

    workflow_callback_ = [f, args...](QueuedObj t)->void { std::invoke(f, args..., std::move(t)); };
    queueThreads(threads);

  }

  // Input stop tokens are rejected.
  void push(QueuedObj input_obj) {

    queue_ptr_->push(std::move(input_obj));

  }

  [[nodiscard]] QueuedObj waitAndPop() {

    return queue_ptr_->waitAndPop();

  }

  // Queue state access routines.
  [[nodiscard]] const Queue<QueuedObj>& ObjectQueue() const { return *queue_ptr_; }


private:

  std::unique_ptr<Queue<QueuedObj>> queue_ptr_;
  QueuedObj stop_token_;
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

      QueuedObj work_item = waitAndPop();

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
