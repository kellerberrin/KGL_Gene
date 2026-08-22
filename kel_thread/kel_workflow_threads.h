// Copyright 2023 Kellerberrin
//

#ifndef KEL_WORKFLOW_THREADS_H
#define KEL_WORKFLOW_THREADS_H

#include "kel_queue_mt_safe.h"

#include <algorithm>
#include <exception>
#include <functional>
#include <future>
#include <optional>
#include <expected>
#include <thread>


namespace kellerberrin {  // organization level namespace


///////////////////////////////////////////////////////////////////////////////////////////
//
// A general purpose thread pool class that returns a std::future with user work function
// results. Work is submitted as std::packaged_task objects so that exceptions thrown by
// the user callback are captured in the future returned by enqueueFuture().
//
// The thread pool is stopped by pushing "stop tokens" (empty optionals) onto the work
// queue, one per worker. Each worker that sees a stop token terminates.
//
// Exceptions from void tasks (enqueueVoid) are caught by the worker threads and stored.
// They are NOT rethrown by joinThreads() so that shutdown remains safe for destructors and
// legacy callers. Call capturedException() after joinThreads() to observe the first void-task
// error, or call shutdownAndRethrow() if you prefer explicit exception propagation.
//
///////////////////////////////////////////////////////////////////////////////////////////

/// Concept: every argument type (after reference decay) must be move constructible so it
/// can be bound into a packaged_task / move_only_function.
template<typename... Args>
concept move_constructible_variadic = (std::move_constructible<std::decay_t<Args>> && ...);

class WorkflowThreads
{
  // Type-erased move-only callable. We use this rather than packaged_task<void()>
  // directly so that enqueueVoid() does not pay for a shared state it never exposes.
  using WorkTask = std::move_only_function<void()>;
  using WorkTaskOpt = std::optional<WorkTask>;

public:

  WorkflowThreads() = default;
  explicit WorkflowThreads(size_t threads) { queueThreads(threads); }
  ~WorkflowThreads() noexcept { joinThreads(); }

  // Non-copyable, non-movable.
  WorkflowThreads(const WorkflowThreads&) = delete;
  WorkflowThreads& operator=(const WorkflowThreads&) = delete;
  WorkflowThreads(WorkflowThreads&&) = delete;
  WorkflowThreads& operator=(WorkflowThreads&&) = delete;

  // Convenience routines. The default number of worker threads is the hardware concurrency
  // minus one, with a minimum of one. If hardware_concurrency() returns 0 we fall back to a
  // single thread rather than underflowing the unsigned arithmetic.
  [[nodiscard]] static size_t defaultThreads() noexcept {

    const size_t concurrency = std::thread::hardware_concurrency();
    return concurrency == 0 ? 1 : std::max<size_t>(concurrency - 1, 1);

  }

  [[nodiscard]] static size_t defaultThreads(size_t job_size) noexcept {

    return job_size > 0 ? std::min<size_t>(defaultThreads(), job_size) : 1;

  }

  /// Starts the requested number of worker threads.
  ///
  /// @param threads  Number of worker threads to start (clamped to at least 1).
  /// @return true if threads were started, false if the pool is already running.
  bool queueThreads(size_t threads)
  {

    if (not workers_.empty()) {

      ExecEnv::log().error("Attempt to queue threads: {} to active thread pool", threads);
      return false;

    }

    threads = std::max<size_t>(threads, 1);
    for (size_t i = 0; i < threads; ++i) {

      workers_.push(std::thread(&WorkflowThreads::workerLoop, this));

    }

    return true;

  }

  /// Signals all worker threads to stop and waits for them to finish.
  ///
  /// This function is noexcept: any exception that occurs while stopping threads or any
  /// exception captured from a void task is swallowed. Call capturedException() afterwards
  /// to inspect void-task errors, or call shutdownAndRethrow() to propagate them.
  void joinThreads() noexcept
  {

    try {
      // Push one stop token per worker. Unlike a single re-queued token, this is robust
      // against concurrent work submissions and guarantees every worker wakes up.

      const size_t token_count = workers_.size();
      for (size_t i = 0; i < token_count; ++i) {

        work_queue_.push(WorkTaskOpt{});

      }

      while (not workers_.empty()) {

        auto worker = workers_.waitAndPop();
        worker.join();

      }

      work_queue_.clear();

    } catch (...) {

      ExecEnv::log().error("Thread pool joinThreads threw exception");

    }

  }

  /// Enables asynchronous queue monitoring.
  void launchMonitor(std::string queue_name, size_t sample_frequency_ms, bool monitor_stalled = true)
  {

    work_queue_.monitor().launchStats(sample_frequency_ms, std::move(queue_name), monitor_stalled);

  }

  [[nodiscard]] bool running() const
  {

    return not workers_.empty();

  }

  [[nodiscard]] size_t threadCount() const
  {

    return workers_.size();

  }


  /// Enqueues a work function and its arguments, returning a future that will hold the result
  /// Any exceptions by the callable are returned in the std::future.
  template<typename F, typename... Args>
  requires std::invocable<F, Args...> && move_constructible_variadic<F, Args...>
  [[nodiscard]] auto enqueueFuture(F&& f, Args&&... args) -> std::future<std::invoke_result_t<F, Args...>>
  {

    using return_type = std::invoke_result_t<F, Args...>;

    auto callable = std::bind_front(std::forward<F>(f), std::forward<Args>(args)...);
    auto typed_task = std::packaged_task<return_type()>(std::move(callable));
    std::future<return_type> future = typed_task.get_future();
    // Wrap the typed packaged_task inside a void() move_only_function. This lets all tasks
    // share the same queue type while preserving the per-task future and exception state.
    work_queue_.push(WorkTaskOpt{[task = std::move(typed_task)] () mutable ->void { task(); }});

    return future;

  }

  /// Enqueues a work function with a void (or ignored) return type. No future is returned.
  /// Exceptions thrown by the callable are caught when executed in workerLoop().
  template<typename F, typename... Args>
  requires std::invocable<F, Args...> && move_constructible_variadic<F, Args...>
  void enqueueVoid(F&& f, Args&&... args)
  {

    auto callable = std::bind_front(std::forward<F>(f), std::forward<Args>(args)...);
    work_queue_.push(WorkTaskOpt{[task = std::move(callable)] () mutable ->void { task(); }});

  }

private:

  QueueMtSafe<std::thread> workers_;
  QueueMtSafe<WorkTaskOpt> work_queue_;


  /// Worker thread main loop. Each iteration waits for a task and executes it.
  /// Tasks queued with futures will return exceptions in the future.
  /// Void task exceptions are caught here.
  void workerLoop()
  {

    while (true) {

      WorkTaskOpt task_opt = work_queue_.waitAndPop();

      if (not task_opt.has_value()) { break; }

      try {

        (*task_opt)();

      } catch (const std::exception& e) {

        ExecEnv::log().error("Thread pool void task threw exception: {}", e.what());

      } catch (...) {

        ExecEnv::log().error("Thread pool void task threw unknown exception.");

      }

    }

  }

};


} // namespace


#endif // KEL_WORKFLOW_THREADS_H
