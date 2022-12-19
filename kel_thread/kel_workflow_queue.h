//
// Created by kellerberrin on 14/12/22.
//

#ifndef KEL_WORKFLOW_QUEUE_H
#define KEL_WORKFLOW_QUEUE_H

#include <condition_variable>
#include <functional>
#include <iostream>
#include <future>
#include <vector>
#include <thread>
#include <queue>
#include <algorithm>


#include "kel_mt_queue.h"


namespace kellerberrin {   //  organization level namespace



template<typename T, typename U> class WorkflowThreads
{

  using Proc = std::function<U(T)>;

public:

  WorkflowThreads() = default;
  ~WorkflowThreads() noexcept {

    joinThreads();

  }

  // Convenience routines, default is available hardware threads minus 1, minimum 1 thread.
  [[nodiscard]] static size_t defaultThreads() { return std::max<size_t>(std::thread::hardware_concurrency() - 1, 1); }
  [[nodiscard]] static size_t defaultThreads(size_t job_size) { return (job_size > 0 ? std::min<size_t>(defaultThreads(), job_size) : 1); }

  // Assumes the function has a void return type, does not return a future.
  template<typename F, typename... Args>
  void registerCallback(T input_stop, U output_stop, size_t threads, F&& f, Args&&... args) noexcept
  {

    input_stop_ = input_stop;
    output_stop_ = output_stop;
    workflow_callback_ = [f, args...](T t_obj)->U { return f(t_obj, args...); };
    queueThreads(threads);

  }

  [[nodiscard]] size_t threadCount() const { return threads_.size(); }

  void push(T input_obj) {

    input_queue_.push(input_obj);

  }

  [[nodiscard]] U waitAndPop() {

    return output_queue_.waitAndPop();

  }

  void joinThreads() {

    push(input_stop_);

    for(auto& thread : threads_) {

      thread.join();

    }

    threads_.clear();

  // There are no active threads processing the input queue.

    output_queue_.push(output_stop_);

  }

  const MtQueue<T>& inputQueue() const { return input_queue_; }
  const MtQueue<U>& outputQueue() const { return output_queue_; }

private:

  MtQueue<T> input_queue_;
  T input_stop_;
  MtQueue<U> output_queue_;
  U output_stop_;
  Proc workflow_callback_;
  std::vector<std::thread> threads_;

  void queueThreads(size_t threads)
  {

    // Always have at least one worker thread queued.
    threads = std::max<size_t>(threads, 1);

    // Queue the worker threads,
    for(size_t i = 0; i < threads; ++i) {

      threads_.emplace_back(&WorkflowThreads::threadProlog, this);

    }

  }

  void threadProlog() {

    while(true) {

      T workItem = input_queue_.waitAndPop();

      if (workItem == input_stop_) {

        push(input_stop_);

        break;

      }

      output_queue_.push(workflow_callback_(workItem));

    }

  }

};




}   // end namespace

#endif //KEL_WORKFLOW_QUEUE_H
