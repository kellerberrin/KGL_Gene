// Copyright 2023 Kellerberrin
//

#ifndef KEL_THREAD_POOL_H
#define KEL_THREAD_POOL_H

#include "kel_queue_mt_safe.h"

#include <functional>
#include <future>
#include <vector>
#include <thread>


namespace kellerberrin {  //  organization level namespace


///////////////////////////////////////////////////////////////////////////////////////////
//
// A general purpose thread pool class returns a std::future with user work function results.
//
///////////////////////////////////////////////////////////////////////////////////////////

template<typename... Args>
concept move_constructable_variadic = std::conjunction_v<std::is_move_constructible<Args>...>;

class WorkflowThreads
{

  using ThreadFunc = std::move_only_function<void(void)>; // Function args can be moved (std::unique_ptr).
  using ThreadFuncPtr = std::unique_ptr<ThreadFunc>; // Functions can be heterogeneous (type erasure).

public:

  WorkflowThreads() = default;
  explicit WorkflowThreads(size_t threads) { queueThreads(threads); }
  ~WorkflowThreads() { joinThreads(); }

  // Convenience routines, default is available hardware threads minus 1, minimum 1 thread.
  [[nodiscard]] static size_t defaultThreads() { return std::max<size_t>(std::thread::hardware_concurrency() - 1, 1); }
  [[nodiscard]] static size_t defaultThreads(size_t job_size) { return (job_size > 0 ? std::min<size_t>(defaultThreads(), job_size) : 1); }


  // A task is a work function and associated arguments. Tasks can be heterogeneous.
  // Returns a std::future holding the work function return value.
  // Note that any exceptions thrown by the work function will also be returned in the std::future and
  // can be captured enclosing the external future.get() in a try/catch block.
  template<typename F, typename... Args>
  requires std::invocable<F, Args...> && move_constructable_variadic<Args...>
  [[nodiscard]] auto enqueueFuture(F&& f,  Args&&... args) -> std::future<std::invoke_result_t<F, Args...>>
  {

    using return_type = std::invoke_result_t<F, Args...>;

    // Create a callable by binding the function task arguments. Then package the task.
    // Anonymize the task callable type (type erasure) using a unique pointer.
    // Push the pointer onto the thread-safe work queue, and return a future.
    auto callable = std::bind_front(std::forward<F>(f), std::forward<Args>(args)...);
    auto task = std::packaged_task<return_type()>(callable);
    std::future<return_type> future = task.get_future();
    ThreadFuncPtr work_callable_ptr = std::make_unique<ThreadFunc>(std::move(task));
    work_queue_.push(std::move(work_callable_ptr));

    return future;

  }

  // Assumes the work function has a void return type and therefore does not return a future.
  template<typename F, typename... Args>
  requires std::invocable<F, Args...> && move_constructable_variadic<Args...>
  void enqueueVoid(F&& f, Args&&... args)
  {

    auto callable = std::bind_front(std::forward<F>(f), std::forward<Args>(args)...);
    ThreadFuncPtr work_callable_ptr = std::make_unique<ThreadFunc>(std::move(callable));
    work_queue_.push(std::move(work_callable_ptr));

  }


  // Force all threads to join().
  void joinThreads() {

    work_queue_.push(nullptr);

    for(auto& thread : threads_) {

      thread.join();

    }

    threads_.clear();
    work_queue_.clear();

  }

  // This function is only valid if there are no active threads (threadCount() == 0).
  bool queueThreads(size_t threads)
  {

    if (not threads_.empty()) {

      return false;

    }

    // Always have at least one worker thread queued.
    threads = std::max<size_t>(threads, 1);

    // Queue the worker threads,
    for(size_t i = 0; i < threads; ++i) {

      threads_.emplace_back(&WorkflowThreads::threadProlog, this);

    }

    return true;

  }

  [[nodiscard]] size_t threadCount() const { return threads_.size(); }

private:

  std::vector<std::thread> threads_;
  QueueMtSafe<ThreadFuncPtr> work_queue_;

  void threadProlog() {

    while(true)
    {

      ThreadFuncPtr work_callable_ptr = work_queue_.waitAndPop();

      if (not work_callable_ptr) {

        work_queue_.push(nullptr);
        break;

      }

      std::invoke(*work_callable_ptr);

    }

  }


};


} // namespace


#endif //KEL_THREAD_POOL_H
