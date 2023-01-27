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

#ifndef KEL_BOUND_QUEUE_H
#define KEL_BOUND_QUEUE_H


#include "kel_mt_queue.h"
#include "kel_queue_monitor.h"
#include "kel_exec_env.h"

#include <iostream>
#include <cstddef>
#include <array>
#include <optional>
#include <queue>
#include <utility>

namespace kellerberrin {   //  organization level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The bounded multi thread queue has a maximum of high_tide elements and a low tide
// when the producer(s) can again start pushing elements after a high tide event.
// Producer threads can push elements onto the queue until high tide is reached, the producer threads are then
// blocked by a condition variable. The queue is then drained by consumer threads (ebbing tide) until the low tide
// level is reached when the producer threads are unblocked and can once again push elements onto the queue (flood tide).
// This automatically balances CPU usage between producer threads and consumer threads.
// It also stops excessive memory usage if producers can queue records faster than consumers can remove them.
//
// Note, if there are no active consumers then this queue will block forever on high tide, a condition known as a 'stalled'
// queue. This stalled condition can be optionally monitored. The monitor also generates statistics on 'high-tide',
// 'low-tide', 'flood-tide', 'ebb-tide' and average queue size. These statistics can be used to efficiently allocate
// CPU resources (threads) between producers and consumers.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The bounded (tidal) queue implemented using a std::queue.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
class RingBuffer {

public:

  explicit RingBuffer(size_t buf_size) : buffer_size_(buf_size), data_(buf_size) {}
  ~RingBuffer() = default;

  [[nodiscard]] T front() {

    if (empty()) {
      ExecEnv::log().critical("Attempted to access empty ringbuffer");
    }

    return std::move(data_[tail_]);

  }

  void push(T item) {

    if (is_full_) {
      ExecEnv::log().critical("Attempted to push full ringbuffer");
    }

    data_[head_] = std::move(item);
    head_ = (head_ + 1) % buffer_size_;
    is_full_ = head_ == tail_;

  }

  void pop() {

    if (empty()) {
      ExecEnv::log().critical("Attempted to pop empty ringbuffer");
    }
    tail_ = (tail_ + 1) % buffer_size_;
    is_full_ = false;

  }

  [[nodiscard]] bool empty() const noexcept { return !is_full_ && tail_ == head_; }

  [[nodiscard]] size_t capacity() { return buffer_size_; }

  [[nodiscard]] size_t size() const {

    if (is_full_) return buffer_size_;

    if (head_ >= tail_) return head_ - tail_;

    return buffer_size_ + head_ - tail_;

  }

private:

  const size_t buffer_size_;
  std::vector<T> data_;
  size_t head_{0};
  size_t tail_{0};
  bool is_full_{false};

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


constexpr static const size_t TIDAL_QUEUE_DEFAULT_HIGH_TIDE{10000};
constexpr static const size_t TIDAL_QUEUE_DEFAULT_LOW_TIDE{2000};

template<typename T> requires std::move_constructible<T>
class TidalQueue {

public:

  explicit TidalQueue( size_t high_tide = TIDAL_QUEUE_DEFAULT_HIGH_TIDE
                     , size_t low_tide = TIDAL_QUEUE_DEFAULT_LOW_TIDE) : high_tide_(high_tide), low_tide_(low_tide), queue_(high_tide) {}

  // This constructor attaches a queue monitor for a 'stalled' queue condition and generates tidal statistics.
  TidalQueue( size_t high_tide
      , size_t low_tide
      , std::string queue_name
      , size_t sample_frequency): high_tide_(high_tide), low_tide_(low_tide), queue_(high_tide) {

    monitor_ptr_ = std::make_unique<BoundedQueueMonitor<TidalQueue<T>>>();
    monitor_ptr_->launchStats(this, sample_frequency, queue_name);

  }
  ~TidalQueue() { monitor_ptr_ = nullptr; }

  // Enqueue function can be called by multiple threads.
  // These threads will block if the queue has reached high-tide size until the queue size reaches low-tide (ebb-tide)..
  // Once the queue has reached low-tide through consumer activity the producer threads are once again unblocked (flood-tide).
  void push(T new_value) {

    { // Mutex
      std::unique_lock<std::mutex> lock(queue_mutex_);
      tide_cond_.wait(lock, [this]()->bool{ return queue_state_; });

      queue_.push(std::move(new_value));
      ++queue_size_;
      ++queue_activity_;
      queue_state_ = queue_size_ < high_tide_;

    } // ~Mutex

    empty_cond_.notify_one();

  }

  // Dequeue function can be called by multiple threads.
  // These threads will only block if the queue is empty.
  [[nodiscard]] T waitAndPop() {

    std::unique_lock<std::mutex> lock(queue_mutex_); // Mutex
    empty_cond_.wait(lock, [this]()->bool{ return queue_size_ > 0; });
    T value(std::move(queue_.front()));
    queue_.pop();
    --queue_size_;
    ++queue_activity_;
    if (not queue_state_) {

      queue_state_ = queue_size_ < low_tide_;

    }
    lock.unlock();  // ~Mutex

    tide_cond_.notify_one();

    return value;

  }

  // All of these functions are thread safe but not guaranteed to show the current state of an active queue.
  [[nodiscard]] bool empty() const { return queue_size_ == 0; }
  [[nodiscard]] size_t size() const { return queue_size_; }
  [[nodiscard]] size_t activity() const { return queue_activity_; }
  [[nodiscard]] bool queueState() const { return queue_state_; }

  [[nodiscard]] size_t highTide() const { return high_tide_; }
  [[nodiscard]] size_t lowTide() const { return low_tide_; }

private:

  // Tidal limits,
  const size_t high_tide_;
  const size_t low_tide_;

  // Actual queue implementation.
//  std::queue<T> queue_;
  RingBuffer<T> queue_;

  // If the queue state is 'true' then producer threads can push() onto the queue ('flood tide').
  // If the queue state is 'false' the producer threads are blocked and are waiting to push() onto the queue ('ebb tide').
  std::atomic<bool> queue_state_{true};
  std::atomic<size_t> queue_size_{0};
  std::atomic<size_t> queue_activity_{0};

  // Condition variable blocks queue producers on 'high tide' and subsequent 'ebb tide' conditions.
  std::condition_variable tide_cond_;
  // Condition variable blocks queue consumers on queue empty.
  std::condition_variable empty_cond_;
  std::mutex queue_mutex_;

  // Held in a pointer for explicit object lifetime.
  std::unique_ptr<BoundedQueueMonitor<TidalQueue<T>>> monitor_ptr_;


};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The bounded queue implemented using an MtQueue.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


constexpr static const size_t BOUNDED_QUEUE_DEFAULT_HIGH_TIDE{10000};
constexpr static const size_t BOUNDED_QUEUE_DEFAULT_LOW_TIDE{2000};

template<typename T> requires std::move_constructible<T> class BoundedMtQueue {

public:

  explicit BoundedMtQueue( size_t high_tide = BOUNDED_QUEUE_DEFAULT_HIGH_TIDE
                         , size_t low_tide = BOUNDED_QUEUE_DEFAULT_LOW_TIDE) : high_tide_(high_tide), low_tide_(low_tide) {}

  // This constructor attaches a queue monitor for a 'stalled' queue condition and generates tidal statistics.
  BoundedMtQueue( size_t high_tide
                , size_t low_tide
                , std::string queue_name
                , size_t sample_frequency): high_tide_(high_tide), low_tide_(low_tide) {

    monitor_ptr_ = std::make_unique<BoundedQueueMonitor<BoundedMtQueue<T>>>();
    monitor_ptr_->launchStats(this, sample_frequency, queue_name);

  }
  ~BoundedMtQueue() { monitor_ptr_ = nullptr; }

  // Enqueue function can be called by multiple threads.
  // These threads will block if the queue has reached high-tide size until the queue size reaches low-tide (ebb-tide)..
  // Once the queue has reached low-tide through consumer activity the producer threads are once again unblocked (flood-tide).
  void push(T new_value) {

    if (queue_state_) {

      if (size() < high_tide_) {  // Possible race condition with size() >= high_tide is considered unimportant.

        mt_queue_.push(std::move(new_value));
        return;

      }
      else {

        queue_state_ = false;

      }

    }

    {
      std::unique_lock<std::mutex> lock(tide_mutex_);
      tide_cond_.wait(lock, [this] ()->bool { return queue_state_; });
    }

    mt_queue_.push(std::move(new_value));
    tide_cond_.notify_one();

  }

  // Dequeue function can be called by multiple threads.
  // These threads will only block if the queue is empty.
  [[nodiscard]] T waitAndPop() {

    T value = mt_queue_.waitAndPop();
    { // Locked block to modify condition variable.
      std::scoped_lock lock(tide_mutex_);
      if (not queue_state_) {

        queue_state_ = size() <= low_tide_;

      }
    } // ~Locked.
    tide_cond_.notify_one();
    return value;

  }

  // All of these functions are thread safe.
  [[nodiscard]] bool empty() const { return mt_queue_.empty(); }
  [[nodiscard]] size_t size() const { return mt_queue_.size(); }
  [[nodiscard]] size_t activity() const { return mt_queue_.activity(); }

  [[nodiscard]] bool queueState() const { return queue_state_; }
  [[nodiscard]] size_t highTide() const { return high_tide_; }
  [[nodiscard]] size_t lowTide() const { return low_tide_; }


private:

  const size_t high_tide_;
  const size_t low_tide_;

  // Multi-thread queue is safe for multiple consumer and producer threads.
  MtQueue<T> mt_queue_;

  // If the queue state is 'true' then producer threads can push() onto the queue ('flood tide').
  // If the queue state is 'false' the producer threads are blocked and are waiting to push() onto the queue ('ebb tide').
  std::atomic<bool> queue_state_{true};

  // Condition variable blocks queue producers on 'high tide' and subsequent 'ebb tide' conditions.
  std::condition_variable tide_cond_;
  mutable std::mutex tide_mutex_;

  // Held in a pointer for explicit object lifetime.
  std::unique_ptr<BoundedQueueMonitor<BoundedMtQueue<T>>> monitor_ptr_;


};



} // namespace

#endif //KEL_BOUND_QUEUE_H
