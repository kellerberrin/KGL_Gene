//
// Created by kellerberrin on 25/7/20.
//
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

#include "kel_mt_queue.h"

#include <condition_variable>
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

class WorkflowThreads
{

  using Proc = std::function<void(void)>;

public:

  explicit WorkflowThreads(size_t threads) { queueThreads(threads); }
  ~WorkflowThreads() noexcept { joinThreads(); }

  // Convenience routines, default is available hardware threads minus 1, minimum 1 thread.
  [[nodiscard]] static size_t defaultThreads() { return std::max<size_t>(std::thread::hardware_concurrency() - 1, 1); }
  [[nodiscard]] static size_t defaultThreads(size_t job_size) { return (job_size > 0 ? std::min<size_t>(defaultThreads(), job_size) : 1); }

  // Assumes the work function has a void return type and therefore does not return a future.
  template<typename F, typename... Args>
  void enqueueWork(F&& f, Args&&... args) noexcept(std::is_nothrow_invocable<decltype(&MtQueue<Proc>::push)>::value)
  {

    auto task = std::bind(std::forward<F>(f), std::forward<Args>(args)...);
    work_queue_.push([task]{ task(); });

  }


  // Returns a std::future holding the work function return value.
  template<typename F, typename... Args>
  [[nodiscard]] auto enqueueTask(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type>
  {

    using return_type = typename std::result_of<F(Args...)>::type;

    auto task_ptr = std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
    std::future<return_type> future = task_ptr->get_future();
    work_queue_.push([task_ptr]{ (*task_ptr)(); });
    return future;

  }

  [[nodiscard]] size_t threadCount() const { return threads_.size(); }

private:

  std::vector<std::thread> threads_;
  MtQueue<Proc> work_queue_;

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

    while(true)
    {

      Proc workItem = work_queue_.waitAndPop();

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
