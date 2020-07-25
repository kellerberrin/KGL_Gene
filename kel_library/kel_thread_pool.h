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
// A general thread pool class.
//
///////////////////////////////////////////////////////////////////////////////////////////

class ThreadPool
{

public:

  explicit ThreadPool(size_t threads = std::thread::hardware_concurrency() - 1) { startThreads(threads); }
  ~ThreadPool() noexcept
  {

    work_queue_.push(nullptr);

    for(auto& thread : threads_) {

      thread.join();

    }

  }

  using Proc = std::function<void(void)>;

  template<typename F, typename... Args>
  void enqueue_work(F&& f, Args&&... args) noexcept(std::is_nothrow_invocable<decltype(&MtQueue<Proc>::push)>::value)
  {

    auto fb = std::bind(std::forward<F>(f), std::forward<Args>(args)...);
    work_queue_.push([=]() { fb(); });

  }

  template<typename F, typename... Args>
  auto enqueue_task(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type>
  {

    using return_type = typename std::result_of<F(Args...)>::type;
    auto task = std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
    std::future<return_type> res = task->get_future();
    work_queue_.push([task](){ (*task)(); });
    return res;

  }

private:

  std::vector<std::thread> threads_;
  MtQueue<Proc> work_queue_;

  void startThreads(size_t threads)
  {

    for(size_t i = 0; i < threads; ++i)

      threads_.emplace_back(std::thread([this]() {

        while(true)
        {

          auto workItem = work_queue_.waitAndPop();

          if (workItem == nullptr) {

            work_queue_.push(nullptr);
            break;

          }

          workItem();

        }
      }));
  }

};

} // namespace


#endif //KEL_THREAD_POOL_H
