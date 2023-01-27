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

#ifndef KEL_MT_QUEUE_H
#define KEL_MT_QUEUE_H

#include "kel_queue_monitor.h"

#include <queue>
#include <mutex>
#include <condition_variable>
#include <memory>

namespace kellerberrin {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Thread safe queue for multiple consumer and producer threads.
// This queue can potentially grow without bound if producer threads can enqueue faster than consumer threads can dequeue.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template<typename T> requires std::move_constructible<T> class MtQueue {

public:

  MtQueue() = default;
  // Create the MtQueue with an asynchronous queue monitor that checks for 'stalled' queues and collects queue statistics.
  MtQueue(std::string queue_name, size_t sample_frequency)  {

    monitor_ptr_ = std::make_unique<MtQueueMonitor<T>>();
    monitor_ptr_->launchStats(this, sample_frequency, queue_name);

  }

  ~MtQueue() = default;

  // Copy constructors are removed.
  MtQueue(const MtQueue&) = delete;
  MtQueue(MtQueue&&) = delete;
  MtQueue& operator=(const MtQueue&) = delete;

  // Enqueue function can be called by multiple threads, this queue can potentially grow without bound.
  // These threads will only block on other producer thread enqueue activity.
  void push(T value) {

    // Scope for RAII locking/unlocking.
    {

      std::scoped_lock lock(mutex_);

      data_queue_.push(std::move(value));
      ++size_;
      ++activity_;

    }

    // The notification is sent after the queue is unlocked so that any threads on waitAndPop() can immediately execute.
    data_cond_.notify_one();


  }

  // Dequeue function can be called by multiple threads.
  // These threads will only block if the queue is empty or on other consumer thread dequeue activity.
  [[nodiscard]] T waitAndPop() {

    std::unique_lock<std::mutex> lock(mutex_);
    // Wait on non-empty queue.
    data_cond_.wait(lock, [this]{ return not empty(); });

    --size_;
    ++activity_;

    T value(std::move(data_queue_.front()));
    data_queue_.pop();

    // Unlock the mutex.
    lock.unlock();

    // Notify waiting threads after the queue is unlocked.
    data_cond_.notify_one();

    return value;

  }

  // All of these functions are thread safe.
  [[nodiscard]] bool empty() const { return size_ == 0; }
  [[nodiscard]] size_t size() const { return size_; }
  [[nodiscard]] size_t activity() const { return activity_; }

private:

  std::mutex mutex_;
  std::queue<T> data_queue_; // Implemented as a std::deque.
  std::condition_variable data_cond_;
  std::atomic<size_t> size_{0};
  std::atomic<size_t> activity_{0};

  // Held in a pointer for explicit object lifetime.
  std::unique_ptr<MtQueueMonitor<T>> monitor_ptr_;

};


}   // end namespace


#endif // KGL_MtQueue_H