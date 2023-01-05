//
// Created by kellerberrin on 19/12/20.
//

#ifndef KEL_BOUND_QUEUE_H
#define KEL_BOUND_QUEUE_H


#include "kel_mt_queue.h"
#include "kel_queue_monitor.h"


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

template<typename T> requires std::move_constructible<T> class BoundedMtQueue {

public:

  explicit BoundedMtQueue( size_t high_tide = DEFAULT_HIGH_TIDE
                         , size_t low_tide = DEFAULT_LOW_TIDE) : high_tide_(high_tide), low_tide_(low_tide) {}

  // This constructor attaches a queue monitor for a 'stalled' queue condition and generates tidal statistics.
  BoundedMtQueue( size_t high_tide
                , size_t low_tide
                , std::string queue_name
                , size_t sample_frequency): high_tide_(high_tide), low_tide_(low_tide) {

    monitor_ptr_ = std::make_unique<BoundedQueueMonitor<T>>();
    monitor_ptr_->launchStats(this, sample_frequency, queue_name);

  }
  ~BoundedMtQueue() { monitor_ptr_ = nullptr; }

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

  T waitAndPop() {

    T value = mt_queue_.waitAndPop();
    { // Locked block to modify condition variable.
      std::scoped_lock lock(tide_mutex_);
      if (not queue_state_) {

        queue_state_ = size() <= low_tide_;

      }
    }
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

  constexpr static const size_t DEFAULT_HIGH_TIDE{10000};
  constexpr static const size_t DEFAULT_LOW_TIDE{2000};

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
  std::unique_ptr<BoundedQueueMonitor<T>> monitor_ptr_;


};



} // namespace

#endif //KEL_BOUND_QUEUE_H
