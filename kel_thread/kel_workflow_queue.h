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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A threaded workflow for copyable objects, specifically NOT std::unique_ptr. Use std::shared_ptr with this object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin {   //  organization level namespace


template<typename T, typename U>
requires (std::copy_constructible<T> && std::equality_comparable<T>) && (std::copy_constructible<U> && std::equality_comparable<U>)
class WorkflowThreads
{

  using Proc = std::function<U(T)>;

public:

  WorkflowThreads() = default;
  ~WorkflowThreads() noexcept {

    stopProcessing();

  }

  // Convenience routines, default is available hardware threads minus 1, minimum 1 thread.
  [[nodiscard]] static size_t defaultThreads() { return std::max<size_t>(std::thread::hardware_concurrency() - 1, 1); }
  [[nodiscard]] static size_t defaultThreads(size_t job_size) { return (job_size > 0 ? std::min<size_t>(defaultThreads(), job_size) : 1); }

  // Assumes the function has a void return type, does not return a future.
  template<typename F, typename... Args>
  void registerProcessingFn(T input_stop, U output_stop, size_t threads, F f, Args... args) noexcept
  {

    input_stop_ = input_stop;
    output_stop_ = output_stop;
    workflow_callback_ = [f, args...](T t)->U { return f(t, args...); };
    queueThreads(threads);

  }


  void push(T input_obj) {

    input_queue_.push(input_obj);

  }

  [[nodiscard]] U waitAndPop() {

    return output_queue_.waitAndPop();

  }

  void stopProcessing() {

    push(input_stop_);

    for(auto& thread : threads_) {

      thread.join();

    }

    threads_.clear();

  // There are no active threads processing the input queue.

    output_queue_.push(output_stop_);

  }

  [[nodiscard]] size_t threadCount() const { return threads_.size(); }
  [[nodiscard]] const MtQueue<T>& inputQueue() const { return input_queue_; }
  [[nodiscard]] const MtQueue<U>& outputQueue() const { return output_queue_; }

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

      } else {

        output_queue_.push(workflow_callback_(workItem));

      }

    }

  }

};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A threaded workflow for std::move constructable objects. Use std::unique_ptr with this object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename T, typename U>
requires (std::move_constructible<T> && std::equality_comparable<T>) && (std::move_constructible<U> && std::equality_comparable<U>)
class WorkflowMove
{

  using Proc = std::function<U(T)>;

public:

  WorkflowMove() = default;
  ~WorkflowMove() noexcept {

    stopProcessing();

  }

  // Convenience routines, default is available hardware threads minus 1, minimum 1 thread.
  [[nodiscard]] static size_t defaultThreads() { return std::max<size_t>(std::thread::hardware_concurrency() - 1, 1); }
  [[nodiscard]] static size_t defaultThreads(size_t job_size) { return (job_size > 0 ? std::min<size_t>(defaultThreads(), job_size) : 1); }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe.
  template<typename F, typename... Args>
  void registerProcessingFn(T&& input_stop, U&& output_stop, size_t threads, F&& f, Args&&... args) noexcept
  {

    input_stop_ = std::move(input_stop);
    output_stop_ = std::move(output_stop);
    workflow_callback_ = [f, args...](T&& t)->U { return f(std::move(t), args...); };
    queueThreads(threads);

  }


  // Input stop tokens are rejected.
  void push(T&& input_obj) {

    if (input_obj == input_stop_) {

      return;

    }

    input_queue_.push(std::move(input_obj));

  }

  [[nodiscard]] U waitAndPop() {

    return output_queue_.waitAndPop();

  }

  void stopProcessing() {

    // Push the stop token onto the input queue.
    input_queue_.push(std::move(input_stop_));

    // Join the active threads
    for(auto& thread : threads_) {

      thread.join();

    }

    // Remove the inactive threads.
    threads_.clear();
    // Ensure the input queue is empty
    while(not input_queue_.empty()) {

      T t = input_queue_.waitAndPop();

    }

    // Push the stop token onto the output queue.
    output_queue_.push(std::move(output_stop_));

  }

  [[nodiscard]] size_t threadCount() const { return threads_.size(); }
  [[nodiscard]] const MtQueue<T>& inputQueue() const { return input_queue_; }
  [[nodiscard]] const MtQueue<U>& outputQueue() const { return output_queue_; }

private:

  MtQueue<T> input_queue_;
  T input_stop_;
  MtQueue<U> output_queue_;
  U output_stop_;
  Proc workflow_callback_;
  std::vector<std::jthread> threads_;

  void queueThreads(size_t threads)
  {

    // Always have at least one worker thread queued.
    threads = std::max<size_t>(threads, 1);

    // Queue the worker threads,
    for(size_t i = 0; i < threads; ++i) {

      threads_.emplace_back(&WorkflowMove::threadProlog, this);

    }

  }

  void threadProlog() {

    while(true) {

      T workItem = input_queue_.waitAndPop();

      if (workItem == input_stop_) {

        input_queue_.push(std::move(workItem));
        break;

      } else {

        U result = std::move(workflow_callback_(std::move(workItem)));
        output_queue_.push(std::move(result));

      }

    }

  }

};


}   // end namespace

#endif //KEL_WORKFLOW_QUEUE_H
