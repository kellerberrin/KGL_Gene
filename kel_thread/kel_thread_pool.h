//
// Created by kellerberrin on 25/7/20.
//

#ifndef KEL_THREAD_POOL_H
#define KEL_THREAD_POOL_H

#include "kel_mt_queue.h"

#include <condition_variable>
#include <functional>
#include <iostream>
#include <future>
#include <vector>
#include <thread>
#include <queue>
#include <algorithm>


namespace kellerberrin {  //  organization level namespace


///////////////////////////////////////////////////////////////////////////////////////////
//
// A general purpose thread pool class.
//
///////////////////////////////////////////////////////////////////////////////////////////

class ThreadPool
{

  using Proc = std::function<void(void)>;

public:

  explicit ThreadPool(size_t threads) { queueThreads(threads); }
  ~ThreadPool() noexcept { joinThreads(); }

  // Convenience routines, default is available hardware threads minus 1, minimum 1 thread.
  [[nodiscard]] static size_t defaultThreads() { return std::max<size_t>(std::thread::hardware_concurrency() - 1, 1); }
  [[nodiscard]] static size_t defaultThreads(size_t job_size) { return (job_size > 0 ? std::min<size_t>(defaultThreads(), job_size) : 1); }

  // Assumes the function has a void return type, does not return a future.
  template<typename F, typename... Args>
  void enqueueWork(F&& f, Args&&... args) noexcept(std::is_nothrow_invocable<decltype(&MtQueue<Proc>::push)>::value)
  {

    auto fb = std::bind(std::forward<F>(f), std::forward<Args>(args)...);
    work_queue_.push([fb]() { fb(); });

  }


  // Returns a std::future holding the function return value.
  template<typename F, typename... Args>
  [[nodiscard]] auto enqueueTask(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type>
  {

    using return_type = typename std::result_of<F(Args...)>::type;

    auto task = std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
    std::future<return_type> future = task->get_future();
    work_queue_.push([task]{ (*task)(); });
    return future;

  }

  [[nodiscard]] size_t threadCount() const { return threads_.size(); }

private:

  std::vector<std::thread> threads_;
  MtQueue<Proc> work_queue_;

  void queueThreads(size_t threads)
  {

    if (threads < 1) {

      ExecEnv::log().warn("ThreadPool::queueThreads; attempted to initialize with zero (0) threads, initialized with 1 thread.");
      threads = 1;

    }

    // Queue the worker threads,
    for(size_t i = 0; i < threads; ++i) {

      threads_.emplace_back(&ThreadPool::threadProlog, this);

    }

  }

  void threadProlog() {

    while(true)
    {

      auto workItem = work_queue_.waitAndPop();

      if (workItem == nullptr) {

        work_queue_.push(nullptr);
        break;

      }

      workItem();

    }

  }

  void joinThreads() {

    work_queue_.push(nullptr);

    for(auto& thread : threads_) {

      thread.join();

    }

  }

};


} // namespace


#endif //KEL_THREAD_POOL_H
