//
// Created by kellerberrin on 6/09/17.
//

#ifndef KEL_MT_QUEUE_H
#define KEL_MT_QUEUE_H

#include "kel_exec_env.h"

#include <queue>
#include <mutex>
#include <condition_variable>
#include <memory>

namespace kellerberrin {   //  organization level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Thread safe queue for multiple consumer and producer threads.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template<typename T> class MtQueue {

public:

  MtQueue() = default;
  ~MtQueue() = default;
  MtQueue(const MtQueue&) = delete;
  MtQueue(MtQueue&&) = delete;
  MtQueue& operator=(const MtQueue&) = delete;

  void push(T value) {

    // Scope for automatic locking/unlocking.
    {

      std::scoped_lock lock(mutex_);

      data_queue_.push(std::move(value));
      ++size_;
      ++activity_;

    }

    // The notification is sent after the queue is unlocked so that any threads on waitAndPop() can immediately execute.
    data_cond_.notify_one();


  }

  [[nodiscard]] T waitAndPop() {

    std::unique_lock<std::mutex> lock(mutex_);
    // wait on non-empty
    data_cond_.wait(lock, [this]{ return not empty(); });

    --size_;
    ++activity_;

    T value(std::move(data_queue_.front()));
    data_queue_.pop();

    // Unlock the mutex.
    lock.unlock();

    // Notify any other waiting threads after the queue is unlocked.
    data_cond_.notify_one();

    return value;

  }

  [[nodiscard]] bool empty() const { return size_ == 0; }

  [[nodiscard]] size_t size() const { return size_; }

  [[nodiscard]] size_t activity() const { return activity_; }

private:

  std::mutex mutex_;
  std::queue<T> data_queue_;
  std::condition_variable data_cond_;
  std::atomic<size_t> size_{0};
  std::atomic<size_t> activity_{0};


};


}   // end namespace


#endif // KGL_MtQueue_H
