// Copyright 2023 Kellerberrin
//

#ifndef KEL_WORKFLOW_THREADS_H
#define KEL_WORKFLOW_THREADS_H

#include "kel_exec_env.h"
#include "kel_queue_mt_safe.h"

#include <algorithm>
#include <concepts>
#include <exception>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>


namespace kellerberrin {  // organization level namespace


///////////////////////////////////////////////////////////////////////////////////////////
//
// A general purpose thread pool class that returns a std::future with user work function
// results. Work is submitted as std::packaged_task objects so that exceptions thrown by
// the user callback are captured in the future returned by enqueueFuture().
//
// The thread pool is stopped by pushing nullptr stop tokens onto the work queue, one per
// worker. Each worker that sees a stop token terminates.
//
// The work item is held in a std::unique_ptr<move_only_function<void()>>. std::unique_ptr's
// move is always noexcept, so a throwing move can never accidentally turn a real task into
// a stop token (the moved-from-empty-optional risk that an std::optional payload would face
// if the contained type's move constructor could throw).
//
// Exceptions from void tasks (enqueueVoid) are caught by the worker threads and logged;
// they are not stored or rethrown. joinThreads() is noexcept and never propagates an
// exception, so shutdown remains safe for destructors and legacy callers.
//
// Enqueueing work on a stopped pool (after joinThreads()) never runs the task:
// enqueueFuture() returns a future already resolved with std::logic_error and
// enqueueVoid() logs an error and drops the task.
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
  // nullptr is used as a stop token for worker threads. std::unique_ptr is noexcept-movable,
  // so a throwing move cannot accidentally turn a real task into a stop token.
  using WorkTaskPtr = std::unique_ptr<WorkTask>;

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
  ///
  /// If std::thread construction fails partway through (resource exhaustion) the
  /// exception propagates and the pool is rolled back to the Stopped state.
  bool queueThreads(size_t threads)
  {

    // Serialize pool lifecycle operations (queueThreads / joinThreads / enqueue*).
    std::scoped_lock lock(pool_mutex_);

    if (state_ == PoolState::Running) {

      ExecEnv::log().error("Attempt to queue threads: {} to active thread pool", threads);
      return false;

    }

    threads = std::max<size_t>(threads, 1);
    state_ = PoolState::Running;
    size_t started = 0;
    try {

      for (; started < threads; ++started) {

        workers_.push(std::thread(&WorkflowThreads::workerLoop, this));

      }

    } catch (...) {

      // Roll back the partially started pool. The freshly started workers have not
      // processed any tasks, so they cannot re-enter the pool; joining while holding
      // the lock is safe here. If the rollback itself fails (resource exhaustion) the
      // pool is left Running and joinThreads() can still shut it down.
      for (size_t i = 0; i < started; ++i) {

        work_queue_.push(nullptr);

      }
      while (not workers_.empty()) {

        auto worker = workers_.waitAndPop();
        worker.join();

      }
      state_ = PoolState::Stopped;
      throw;

    }

    return true;

  }

  /// Signals all worker threads to stop and waits for them to finish.
  ///
  /// This function is noexcept: any exception that occurs while stopping threads is
  /// swallowed. Void-task exceptions are logged by the worker threads and are not
  /// observable here.
  void joinThreads() noexcept
  {

    try {

      // Push one stop token per worker. Unlike a single re-queued token, this is robust
      // against concurrent work submissions and guarantees every worker wakes up.
      // Holding the lock here excludes concurrent enqueue, so no task can be placed
      // after the stop tokens; state_ is set to Stopped before the lock is released so
      // that enqueue during the join phase is rejected rather than silently dropped.
      {
        std::scoped_lock lock(pool_mutex_);

        const size_t token_count = workers_.size();
        for (size_t i = 0; i < token_count; ++i) {

          // A failed push would strand a worker and hang the join loop below. The queue
          // is unbounded, so a push only fails on transient resource exhaustion; retry
          // until the token is accepted.
          for (;;) {

            try {

              work_queue_.push(nullptr);
              break;

            } catch (...) {

              ExecEnv::log().error("Thread pool joinThreads: stop token push failed, retrying");

            }

          }

        }

        state_ = PoolState::Stopped;

      }

      // Join without holding the lock: a worker task may itself enqueue more work
      // (e.g. a workflow stage), and holding the lock here would deadlock it.
      while (not workers_.empty()) {

        auto worker = workers_.waitAndPop();
        worker.join();

      }

      {
        std::scoped_lock lock(pool_mutex_);
        work_queue_.clear();
      }

    } catch (...) {

      ExecEnv::log().error("Thread pool joinThreads threw exception");

    }

  }

  /// Enables asynchronous queue monitoring.
  void launchMonitor(std::string queue_name, size_t sample_frequency_ms, bool monitor_stalled = true)
  {

    work_queue_.monitor().launchStats(sample_frequency_ms, std::move(queue_name), monitor_stalled);

  }

  [[nodiscard]] bool running() const noexcept
  {

    return not workers_.empty();

  }

  [[nodiscard]] size_t threadCount() const noexcept
  {

    return workers_.size();

  }


  /// Enqueues a work function and its arguments, returning a future that will hold the result
  /// Any exceptions by the callable are returned in the std::future.
  /// If the pool has been stopped (see joinThreads()) the task is not run; the returned
  /// future is already resolved with std::logic_error so get() throws rather than hanging.
  /// Note: arguments are bound by value (std::bind_front decays and copies/moves them);
  /// pass large objects by rvalue to move them into the task.
  template<typename F, typename... Args>
  requires std::invocable<F, Args...> && move_constructible_variadic<F, Args...>
  [[nodiscard]] auto enqueueFuture(F&& f, Args&&... args) -> std::future<std::invoke_result_t<F, Args...>>
  {

    using return_type = std::invoke_result_t<F, Args...>;
    static_assert(std::is_void_v<return_type> or std::is_object_v<return_type>,
                  "enqueueFuture requires a void or object return type; std::packaged_task "
                  "does not support reference return types");

    std::future<return_type> future;

    {
      // Serialize against joinThreads() so a task cannot be enqueued after the pool has
      // stopped (which would leave the future permanently unresolved).
      std::scoped_lock lock(pool_mutex_);

      if (state_ != PoolState::Running) {

        // Do not run the user's callable on a stopped pool. Return a future that is already
        // satisfied with std::logic_error so the caller's get() throws immediately instead
        // of blocking forever.
        ExecEnv::log().error("WorkflowThreads::enqueueFuture called on a stopped thread pool");
        std::packaged_task<return_type()> failed_task{[]() -> return_type {
          throw std::logic_error("WorkflowThreads::enqueueFuture called on a stopped thread pool");
        }};
        future = failed_task.get_future();
        failed_task();

      }
      else {

        auto callable = std::bind_front(std::forward<F>(f), std::forward<Args>(args)...);
        auto typed_task = std::packaged_task<return_type()>(std::move(callable));
        future = typed_task.get_future();

        // Wrap the typed packaged_task inside a void() move_only_function, then store it on
        // the heap. This lets all tasks share the same queue type while preserving the
        // per-task future and exception state, and the noexcept-movable unique_ptr avoids
        // the moved-from-empty-optional problem on a throwing move.
        work_queue_.push(std::make_unique<WorkTask>([task = std::move(typed_task)] () mutable ->void { task(); }));

      }

    }

    return future;

  }

  /// Enqueues a work function with a void (or ignored) return type. No future is returned.
  /// Exceptions thrown by the callable are caught and logged when executed in workerLoop().
  /// If the pool has been stopped (see joinThreads()) an error is logged and the task is dropped.
  /// Note: arguments are bound by value (std::bind_front decays and copies/moves them);
  /// pass large objects by rvalue to move them into the task.
  template<typename F, typename... Args>
  requires std::invocable<F, Args...> && move_constructible_variadic<F, Args...>
  void enqueueVoid(F&& f, Args&&... args)
  {

    {
      // Serialize against joinThreads() so a task cannot be enqueued after the pool has
      // stopped (which would silently drop the work).
      std::scoped_lock lock(pool_mutex_);

      if (state_ != PoolState::Running) {

        ExecEnv::log().error("WorkflowThreads::enqueueVoid called on a stopped thread pool");
        return;

      }

      auto callable = std::bind_front(std::forward<F>(f), std::forward<Args>(args)...);
      work_queue_.push(std::make_unique<WorkTask>([task = std::move(callable)] () mutable ->void { task(); }));

    }

  }

private:

  // Pool lifecycle. Always accessed under pool_mutex_.
  enum class PoolState { Stopped, Running };

  QueueMtSafe<std::thread> workers_;
  QueueMtSafe<WorkTaskPtr> work_queue_;
  // Serializes pool lifecycle operations: queueThreads(), joinThreads() and the enqueue
  // paths. Guards against starting threads twice, enqueueing after shutdown, and the
  // check-then-act race between those operations.
  std::mutex pool_mutex_;
  PoolState state_{PoolState::Stopped};


  /// Worker thread main loop. Each iteration waits for a task and executes it.
  /// Tasks queued with futures will return exceptions in the future.
  /// Void task exceptions are caught here.
  void workerLoop()
  {

    while (true) {

      WorkTaskPtr task_opt = work_queue_.waitAndPop();

      if (not task_opt) { break; }

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
