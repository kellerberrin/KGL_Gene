//

#ifndef KEL_QUEUE_MT_SAFE_H
#define KEL_QUEUE_MT_SAFE_H

#include "kel_queue_monitor.h"

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <memory>
#include <mutex>
#include <queue>
#include <string>

namespace kellerberrin {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Thread safe queue for multiple consumer and producer threads.
// This queue can potentially grow without bound if producer threads can enqueue faster than consumer threads can dequeue.
// Objects on the queue must be move_constructible<T> (T=std::unique_ptr<QueuedObject>).
//
// NOTE: clear() destroys every queued object. Do not enqueue objects that require explicit
// lifecycle handling (e.g. std::thread) and call clear() while they are in the queue; destroying
// such an object is undefined behaviour. waitAndPop() move-constructs its result; if T's move
// constructor can throw, the exception propagates and the queue state is left unchanged.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template<typename T>
requires std::move_constructible<T>
class QueueMtSafe {

public:

  QueueMtSafe() { monitor_ptr_ = std::make_unique<MonitorMtSafe<T>>(this); }
  // Create the QueueMtSafe with an asynchronous queue monitor that checks for 'stalled' queues and collects queue statistics.
  QueueMtSafe(std::string queue_name, size_t sample_frequency)  {

    monitor_ptr_ = std::make_unique<MonitorMtSafe<T>>(this);
    monitor_ptr_->launchStats(sample_frequency, queue_name);

  }

  // Destroy the queue monitor, which joins its statistics thread. The queue itself has no
  // shutdown mechanism; waiting consumers must be released via a stop-token sentinel value.
  ~QueueMtSafe() noexcept { monitor_ptr_ = nullptr; }

  // Non-copyable, non-movable.
  QueueMtSafe(const QueueMtSafe&) = delete;
  QueueMtSafe(QueueMtSafe&&) = delete;
  QueueMtSafe& operator=(const QueueMtSafe&) = delete;
  QueueMtSafe& operator=(QueueMtSafe&&) = delete;

  // Enqueue function can be called by multiple threads, this queue can potentially grow without bound.
  // These threads will only block on other producer thread enqueue activity.
  void push(T value) {

    // Scope for RAII locking/unlocking.
    {

      std::scoped_lock lock(data_mutex_);

      data_queue_.push(std::move(value));
      ++size_;
      ++activity_;

    }

    // The notification is sent after the queue is unlocked so that any threads on waitAndPop() can immediately execute.
    data_cond_.notify_one();


  }


  // Dequeue function can be called by multiple threads.
  // These threads will only block if the queue is empty or on other consumer thread dequeue activity.
  // The counters are updated only after the element has been successfully moved out, so a
  // throwing T move constructor leaves the queue state consistent.
  [[nodiscard]] T waitAndPop() {

    std::unique_lock<std::mutex> lock(data_mutex_);
    // Wait on non-empty queue.
    data_cond_.wait(lock, [this]{ return not empty(); });

    T value(std::move(data_queue_.front()));
    data_queue_.pop();

    --size_;
    ++activity_;

    // Unlock the mutex.
    lock.unlock();

    // Notify waiting threads after the queue is unlocked.
    data_cond_.notify_one();

    return value;

  }

  // Unconditionally empty the queue. Destroys every queued object.
  void clear() {

    {
      std::scoped_lock<std::mutex> lock(data_mutex_);
      // Swap with a local empty queue (std::queue::swap is noexcept for std::deque), then
      // destroy the old elements when the local goes out of scope. Safer than assignment.
      std::queue<T>().swap(data_queue_);
      size_ = 0;
    }

    data_cond_.notify_all();

  }

  [[nodiscard]] MonitorMtSafe<T>& monitor() const noexcept { return *monitor_ptr_; }
  // All of these functions are thread safe.
  [[nodiscard]] bool empty() const noexcept { return size_ == 0; }
  [[nodiscard]] size_t size() const noexcept { return size_; }
  [[nodiscard]] size_t activity() const noexcept { return activity_; }

private:

  std::queue<T> data_queue_; // Implemented as a std::deque.
  std::mutex data_mutex_;
  std::condition_variable data_cond_;
  std::atomic<size_t> size_{0};
  std::atomic<size_t> activity_{0};
  // Held in a pointer for explicit object lifetime.
  std::unique_ptr<MonitorMtSafe<T>> monitor_ptr_;

};


}   // end namespace


#endif // KEL_QUEUE_MT_SAFE_H
