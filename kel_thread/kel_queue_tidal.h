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


#include "kel_queue_mt_safe.h"
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
// Objects on the queue can be std::move_constructible<T> (std::unique_ptr<...>).
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The current state of the tidal queue.
// If the queue is in the 'FLOOD_TIDE' state, producers can push objects onto the queue until 'high tide' is reached
// then the queue reverts to 'EBB_TIDE' state.
// If the queue is in the 'EBB_TIDE' state then producers are blocked and the queue is drained by consumers until 'low tide'
// is reached. Then the queue reverts to the 'FLOOD_TIDE' state and producers are re-enabled.
enum class QueueTidalState { FLOOD_TIDE, EBB_TIDE };

constexpr static const size_t TIDAL_QUEUE_DEFAULT_HIGH_TIDE{10000};
constexpr static const size_t TIDAL_QUEUE_DEFAULT_LOW_TIDE{2000};


template<typename T>
requires std::move_constructible<T>
class QueueTidal {

public:

  explicit QueueTidal( size_t high_tide = TIDAL_QUEUE_DEFAULT_HIGH_TIDE
                       , size_t low_tide = TIDAL_QUEUE_DEFAULT_LOW_TIDE) : high_tide_(high_tide), low_tide_(low_tide) {}

  // This constructor attaches a queue monitor to check for a 'stalled' condition and generates tidal statistics.
  QueueTidal( size_t high_tide
              , size_t low_tide
              , std::string queue_name
              , size_t sample_frequency): high_tide_(high_tide), low_tide_(low_tide) {

    monitor_ptr_ = std::make_unique<MonitorTidal<QueueTidal<T>>>();
    monitor_ptr_->launchStats(this, sample_frequency, queue_name);

  }
  ~QueueTidal() { monitor_ptr_ = nullptr; }

  // Enqueue function can be called by multiple threads.
  // These threads will block if the queue has reached high-tide size until the queue size reaches low-tide (EBB_TIDE)
  // Once the queue has reached low-tide through consumer activity the producer threads are once again unblocked (FLOOD_TIDE).
  void push(T new_value) {

    { // Mutex
      std::unique_lock<std::mutex> lock(queue_mutex_);
      tide_cond_.wait(lock, [this]()->bool{ return queue_tidal_state_ == QueueTidalState::FLOOD_TIDE; });

      queue_.push(std::move(new_value));

      ++queue_size_;
      ++queue_activity_;

      if (queue_size_ >= high_tide_) {

        queue_tidal_state_ = QueueTidalState::EBB_TIDE;

      }

    } // ~Mutex

    empty_cond_.notify_one();

  }

  // Dequeue function can be called by multiple threads.
  // This function returns the null value if the queue is empty.
  [[nodiscard]] std::optional<T> pop() {

    std::unique_lock<std::mutex> lock(queue_mutex_);

    if (empty()) return std::nullopt;

    T value(std::move(queue_.front()));
    queue_.pop();

    --queue_size_;
    ++queue_activity_;

    if (queue_tidal_state_ == QueueTidalState::EBB_TIDE and queue_size_ <= low_tide_) {

      queue_tidal_state_ = QueueTidalState::FLOOD_TIDE;

    }

    lock.unlock();  // ~Mutex

    tide_cond_.notify_one();

    return value;

  }

  // Dequeue function can be called by multiple threads.
  // These threads will block if the queue is empty.
  [[nodiscard]] T waitAndPop() {

    std::unique_lock<std::mutex> lock(queue_mutex_);
    empty_cond_.wait(lock, [this]()->bool{ return not empty(); });

    T value(std::move(queue_.front()));
    queue_.pop();

    --queue_size_;
    ++queue_activity_;

    if (queue_tidal_state_ == QueueTidalState::EBB_TIDE and queue_size_ <= low_tide_) {

      queue_tidal_state_ = QueueTidalState::FLOOD_TIDE;

    }

    lock.unlock();  // ~Mutex

    tide_cond_.notify_one();

    return value;

  }

  void clear() {

    { // Mutex

      std::scoped_lock<std::mutex> lock(queue_mutex_);
      queue_ = {};
      queue_size_ = 0;
      queue_tidal_state_ = QueueTidalState::FLOOD_TIDE;

    } // ~Mutex

    tide_cond_.notify_one();

  }


  // All of these functions are thread safe.
  [[nodiscard]] bool empty() const { return queue_size_ == 0; }
  [[nodiscard]] size_t size() const { return queue_size_; }
  [[nodiscard]] size_t activity() const { return queue_activity_; }
  [[nodiscard]] QueueTidalState queueTidalState() const { return queue_tidal_state_; }
  [[nodiscard]] bool queueState() const { return queue_tidal_state_ == QueueTidalState::FLOOD_TIDE; }

  [[nodiscard]] size_t highTide() const { return high_tide_; }
  [[nodiscard]] size_t lowTide() const { return low_tide_; }

private:

  // Tidal limits.
  const size_t high_tide_;
  const size_t low_tide_;

  // Actual queue implementation.
  std::queue<T> queue_;

  // If the queue state is FLOOD_TIDE then producer threads can push() onto the queue.
  // If the queue state is EBB_TIDE the producer threads are blocked and are waiting to push() onto the queue.
  std::atomic<QueueTidalState> queue_tidal_state_{QueueTidalState::FLOOD_TIDE};
  std::atomic<size_t> queue_size_{0};
  std::atomic<size_t> queue_activity_{0};

  // Condition variable blocks queue producers on 'high tide' and subsequent 'ebb tide' conditions.
  std::condition_variable tide_cond_;
  // Condition variable blocks queue consumers on queue empty.
  std::condition_variable empty_cond_;
  std::mutex queue_mutex_;

  // Queue monitor asynchronously gathers dynamic producer/consumer statistics.
  // it also (optionally) monitors for a 'stalled' queue where there is a possible deadlock condition and the queue is inactive.
  // Held in a pointer for explicit object lifetime.
  std::unique_ptr<MonitorTidal<QueueTidal<T>>> monitor_ptr_;


};


} // namespace

#endif //KEL_BOUND_QUEUE_H
