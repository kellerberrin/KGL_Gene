// Copyright 2023 Kellerberrin
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
// documentation files (the "Software"), to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
// and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
// OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
//

#ifndef KEL_THREAD_POOL_H
#define KEL_THREAD_POOL_H

#include "kel_queue_mt_safe.h"
#include "kel_movefunction.h"

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

// Thread pool state.
enum class WorkflowThreadState { ACTIVE, STOPPED};

class WorkflowThreads
{

// Convert to std::move_only_function when compiling to c++23 standard.
//  using ThreadFunc = std::move_only_function<void(void)>;
  using ThreadFunc = MoveFunction<void(void)>;
  using ThreadFuncPtr = std::unique_ptr<const ThreadFunc>;

public:

  WorkflowThreads() = default;
  explicit WorkflowThreads(size_t threads) { queueThreads(threads); }
  ~WorkflowThreads() { joinThreads(); }

  // Convenience routines, default is available hardware threads minus 1, minimum 1 thread.
  [[nodiscard]] static size_t defaultThreads() { return std::max<size_t>(std::thread::hardware_concurrency() - 1, 1); }
  [[nodiscard]] static size_t defaultThreads(size_t job_size) { return (job_size > 0 ? std::min<size_t>(defaultThreads(), job_size) : 1); }

  // Assumes the work function has a void return type and therefore does not return a future.
  template<typename F, typename... Args>
  requires std::invocable<F, Args...> && move_constructable_variadic<Args...>
  void enqueueVoid(F&& f, Args&&... args)
  {

    auto callable = std::bind_front(std::forward<F>(f), std::forward<Args>(args)...);
    work_queue_.push(std::make_unique<ThreadFunc>(std::move(callable)));

  }

  // Returns a std::future holding the work function return value.
  template<typename F, typename... Args>
  requires std::invocable<F, Args...> && move_constructable_variadic<Args...>
  [[nodiscard]] auto enqueueFuture(F&& f,  Args&&... args) -> std::future<std::invoke_result_t<F, Args...>>
  {

    using return_type = std::invoke_result_t<F, Args...>;

    auto callable = std::bind_front(std::forward<F>(f), std::forward<Args>(args)...);
    auto task = std::packaged_task<return_type()>(callable);
    std::future<return_type> future = task.get_future();
    work_queue_.push(std::make_unique<ThreadFunc>(std::move(task)));

    return future;

  }


  void joinThreads() {

    work_queue_.push(nullptr);

    for(auto& thread : threads_) {

      thread.join();

    }

    threads_.clear();
    work_queue_.clear();
    workflow_thread_state_ = WorkflowThreadState::STOPPED;

  }

  bool queueThreads(size_t threads)
  {

    if (workflow_thread_state_ != WorkflowThreadState::STOPPED) {

      return false;

    }

    // Always have at least one worker thread queued.
    threads = std::max<size_t>(threads, 1);

    // Queue the worker threads,
    for(size_t i = 0; i < threads; ++i) {

      threads_.emplace_back(&WorkflowThreads::threadProlog, this);

    }

    workflow_thread_state_ = WorkflowThreadState::ACTIVE;

    return true;

  }

  [[nodiscard]] WorkflowThreadState threadState() const { return workflow_thread_state_; }
  [[nodiscard]] const QueueMtSafe<ThreadFuncPtr>& workQueue() const { return work_queue_; }
  [[nodiscard]] size_t threadCount() const { return threads_.size(); }

private:

  std::vector<std::thread> threads_;
  QueueMtSafe<ThreadFuncPtr> work_queue_;
  WorkflowThreadState workflow_thread_state_{WorkflowThreadState::STOPPED};

  void threadProlog() {

    while(true)
    {

      ThreadFuncPtr workItem = work_queue_.waitAndPop();

      if (workItem == nullptr) {

        work_queue_.push(nullptr);
        break;

      }

      (*workItem)();

    }

  }


};


} // namespace


#endif //KEL_THREAD_POOL_H
