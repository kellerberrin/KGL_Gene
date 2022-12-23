//
// Created by kellerberrin on 14/12/22.
//

#ifndef KEL_WORKFLOW_QUEUES_H
#define KEL_WORKFLOW_QUEUES_H

#include <functional>
#include <vector>
#include <thread>

#include "kel_mt_queue.h"


namespace kellerberrin {  //  organization level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A threaded workflow for std::move constructable objects (std::unique_ptr).
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename T, typename U>
requires (std::move_constructible<T> && std::equality_comparable<T>) && (std::move_constructible<U> && std::equality_comparable<U>)
class WorkflowQueues
{

  using Proc = std::function<U(T)>;

public:

  WorkflowQueues() = default;
  ~WorkflowQueues() noexcept { stopProcessing(); }

  // Convenience routines, default is available hardware threads minus 1, minimum 1 thread.
  [[nodiscard]] static size_t defaultThreads() { return std::max<size_t>(std::thread::hardware_concurrency() - 1, 1); }
  [[nodiscard]] static size_t defaultThreads(size_t job_size) { return (job_size > 0 ? std::min<size_t>(defaultThreads(), job_size) : 1); }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe.
  template<typename F, typename... Args>
  void registerProcessingFn(T&& input_stop, U&& output_stop, size_t threads, F&& f, Args&&... args) noexcept
  {

    input_stop_ = std::move(input_stop);
    output_stop_ = std::move(output_stop);
    workflow_callback_ = [f, args...](T&& t)->U { return f(args..., std::move(t)); };
    queueThreads(threads);

  }

  // Input stop tokens are rejected.
  void push(T&& input_obj) {

    input_queue_.push(std::move(input_obj));

  }

  [[nodiscard]] U waitAndPop() {

    return output_queue_.waitAndPop();

  }


  // Queue state access routines.
  [[nodiscard]] const MtQueue<T>& inputQueue() const { return input_queue_; }
  [[nodiscard]] const MtQueue<U>& outputQueue() const { return output_queue_; }

private:

  MtQueue<T> input_queue_;
  T input_stop_;
  MtQueue<U> output_queue_;
  U output_stop_;
  Proc workflow_callback_;
  std::vector<std::thread> threads_;
  std::atomic<uint32_t> active_threads_{0};


  void queueThreads(size_t threads)
  {

    // Remove any existing threads.
    stopProcessing();

    // Always have at least one worker thread queued.
    threads = std::max<size_t>(threads, 1);

    // Queue the worker threads,
    for(size_t i = 0; i < threads; ++i) {

      threads_.emplace_back(&WorkflowQueues::threadProlog, this);
      ++active_threads_;

    }

  }

  void threadProlog() {

    while(true) {

      T work_item = input_queue_.waitAndPop();

      if (work_item == input_stop_) {

        // If the last thread then queue a stop token on the output queue.
        if (--active_threads_ == 0) {

          // Push the stop token onto the output queue.
          output_queue_.push(std::move(output_stop_));

        } else { // Push a stop token on the input queue to stop another active thread.

          input_queue_.push(std::move(work_item));

        }
        break;

      } else {

        U u = std::move(workflow_callback_(std::move(work_item)));
        output_queue_.push(std::move(u));

      }

    }

  }

  void stopProcessing() {

    // If any active threads then push the stop token onto the input queue.
    if (active_threads_ != 0) {

      input_queue_.push(std::move(input_stop_));

    }

    // Join all the threads
    for(auto& thread : threads_) {

      thread.join();

    }

    threads_.clear();

  }

};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A threaded workflow for std::move constructable objects (std::unique_ptr).
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename T>
requires (std::move_constructible<T> && std::equality_comparable<T>)
class WorkflowBase
{

public:

  explicit WorkflowBase(T stop_token) : stop_token_(std::move(stop_token)) {}
  ~WorkflowBase() { stopProcessing(); };

  // Convenience routines, default is available hardware threads minus 1, minimum 1 thread.
  [[nodiscard]] static size_t defaultThreads() { return std::max<size_t>(std::thread::hardware_concurrency() - 1, 1); }
  [[nodiscard]] static size_t defaultThreads(size_t job_size) { return (job_size > 0 ? std::min<size_t>(defaultThreads(), job_size) : 1); }


  // Input stop tokens are rejected.
  virtual void push(T&& input_obj) = 0;


//protected:

  MtQueue<T> queue_;
  T stop_token_;
  std::vector<std::thread> threads_;
  std::atomic<uint32_t> active_threads_{0};

  void stopProcessing() {

    // If any active threads then push the stop token onto the input queue.
    if (active_threads_ != 0) {

      queue_.push(std::move(stop_token_));

    }

    // Join all the threads
    for(auto& thread : threads_) {

      thread.join();

    }

    threads_.clear();

  }

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A threaded workflow for std::move constructable objects (std::unique_ptr).
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename T, typename U>
requires (std::move_constructible<T> && std::equality_comparable<T>) && (std::move_constructible<U> && std::equality_comparable<U>)
class WorkflowDouble : public WorkflowBase<T>
{

  using Proc = std::function<U(T)>;

public:

  WorkflowDouble(T stop_token, std::unique_ptr<WorkflowBase<U>>&& output_queue_ptr) : WorkflowBase<T>(std::move(stop_token)), output_queue_ptr_(std::move(output_queue_ptr))  {}
  ~WorkflowDouble() = default;


  // Note that the variadic args... are presented to ALL active threads and must be thread safe.
  template<typename F, typename... Args>
  void registerProcessingFn(size_t threads, F&& f, Args&&... args) noexcept
  {

    workflow_callback_ = [f, args...](T&& t)->U { return f(args..., std::move(t)); };
    queueThreads(threads);

  }

  void push(T&& input_obj) final {

    WorkflowBase<T>::queue_.push(std::move(input_obj));

  }

  [[nodiscard]] U waitAndPop() {

    return output_queue_ptr_->queue_.waitAndPop();

  }


  // Queue state access routines.
  [[nodiscard]] const MtQueue<T>& inputQueue() const { return WorkflowBase<T>::queue_; }
  [[nodiscard]] const MtQueue<U>& outputQueue() const { return output_queue_ptr_->queue_; }

private:

  std::unique_ptr<WorkflowBase<U>> output_queue_ptr_;
  Proc workflow_callback_;

  void queueThreads(size_t threads)
  {

    // Remove any existing threads.
    WorkflowBase<T>::stopProcessing();

    // Always have at least one worker thread queued.
    threads = threads < 1 ? 1 : threads;

    // Queue the worker threads,
    for(size_t i = 0; i < threads; ++i) {

      WorkflowBase<T>::threads_.emplace_back(&WorkflowDouble::threadProlog, this);
      ++WorkflowBase<T>::active_threads_;

    }

  }

  void threadProlog() {

    while(true) {

      T work_item = WorkflowBase<T>::queue_.waitAndPop();

      if (work_item == WorkflowBase<T>::stop_token_) {

        // If the last thread then queue a stop token on the output queue.
        if (--WorkflowBase<T>::active_threads_ == 0) {

          // Push the stop token onto the output queue.
          output_queue_ptr_->queue_.push(std::move(output_queue_ptr_->stop_token_));

        } else { // Push a stop token on the input queue to stop another active thread.

          WorkflowBase<T>::queue_.push(std::move(work_item));

        }
        break;

      } else {

        U u = std::move(workflow_callback_(std::move(work_item)));
        output_queue_ptr_->queue_.push(std::move(u));

      }

    }

  }

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A threaded workflow for std::move constructable objects (std::unique_ptr).
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename T>
requires (std::move_constructible<T> && std::equality_comparable<T>)
class WorkflowSingle : public WorkflowBase<T>
{

  using Proc = std::function<void(T)>;
  using Prolog = void (WorkflowSingle<T>::*)();

public:

  explicit WorkflowSingle(T stop_token) : WorkflowBase<T>(std::move(stop_token)), prolog_proc_(&WorkflowSingle<T>::threadProlog) {}
  ~WorkflowSingle() = default;

  // Note that the variadic args... are presented to ALL active threads and must be thread safe.
  template<typename F, typename... Args>
  void registerProcessingFn(size_t threads, F&& f, Args&&... args) noexcept
  {

    workflow_callback_ = [f, args...](T&& t)->void { return f(args..., std::move(t)); };
    queueThreads(threads);

  }

  // Input stop tokens are rejected.
  void push(T&& input_obj) {

    WorkflowBase<T>::queue_.push(std::move(input_obj));

  }

  [[nodiscard]] T waitAndPop() {

    return WorkflowBase<T>::queue_.waitAndPop();

  }

  // Queue state access routines.
  [[nodiscard]] const MtQueue<T>& inputQueue() const { return WorkflowBase<T>::queue_; }

private:

  Proc workflow_callback_;
  Prolog prolog_proc_;

  void queueThreads(size_t threads)
  {

    // Remove any existing threads.
    WorkflowBase<T>::stopProcessing();

    // Always have at least one worker thread queued.
    threads = std::max<size_t>(threads, 1);

    // Queue the worker threads,
    for(size_t i = 0; i < threads; ++i) {

      WorkflowBase<T>::threads_.emplace_back(prolog_proc_, this);
      ++WorkflowBase<T>::active_threads_;

    }

  }

  void threadProlog() {

    while(true) {

      T work_item = WorkflowBase<T>::queue_.waitAndPop();

      if (work_item == WorkflowBase<T>::stop_token_) {

        // If the last thread then do not queue a stop token.
        if (--WorkflowBase<T>::active_threads_ != 0) {

          WorkflowBase<T>::queue_.push(std::move(work_item));

        } else {

          workflow_callback_(std::move(work_item));

        }
        break;

      } else {

        workflow_callback_(std::move(work_item));

      }

    }

  }

};



}   // end namespace

#endif //KEL_WORKFLOW_QUEUES_H
