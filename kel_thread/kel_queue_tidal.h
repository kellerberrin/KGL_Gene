// Copyright 2023 Kellerberrin
//

#ifndef KEL_QUEUE_TIDAL_H
#define KEL_QUEUE_TIDAL_H

#include "kel_exec_env.h"
#include "kel_queue_monitor.h"

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <utility>

namespace kellerberrin {   //  organization level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The tidal multi thread queue has a maximum of high_tide elements and a low tide
// when the producer(s) can again start enqueueing elements after a high tide event.
// Producer threads can enqueue elements onto the queue until high tide is reached, the producer threads are then
// blocked by a condition variable. The queue is then drained by consumer threads (ebbing tide) until the low tide
// level is reached when the producer threads are unblocked and can once again push elements onto the queue (flood tide).
// A low tide of zero is valid: it simply means producers are only re-enabled once the queue has been completely
// drained. When the ebb reaches low tide every blocked producer is woken (tide_cond_.notify_all()) so that all
// producers can flood the queue up to high tide together rather than being released one at a time.
// This automatically balances CPU usage between producer threads and consumer threads.
// It also stops excessive memory usage if producers can queue records faster than consumers can remove them.
//
// Note, if there are no active consumers then this queue will block forever on high tide, a condition known as a 'stalled'
// queue. This stalled condition can be optionally monitored. The monitor also generates statistics on 'high-tide',
// 'low-tide', 'flood-tide', 'ebb-tide' and average queue size. These statistics can be used to efficiently allocate
// CPU resources (threads) between producers and consumers.
//
// Exception safety: push() and waitAndPop() move-construct their payload under the mutex and only update the
// counters after the container operation succeeds, so a throwing T move constructor leaves the queue contents,
// the counters and the tidal state unchanged (the element remains queued / the queue remains non-empty).
// clear() swaps the underlying container into a temporary (noexcept for std::deque) rather than assigning {},
// which would additionally require T to be move-assignable.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The bounded (tidal) queue implemented using a std::queue.
// Objects on the queue must be std::move_constructible<T> (std::unique_ptr<...>).
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
                       , size_t low_tide = TIDAL_QUEUE_DEFAULT_LOW_TIDE) : high_tide_(high_tide), low_tide_(low_tide) {

    validateTide();
    monitor_ptr_ = std::make_unique<MonitorTidal<T>>(this);

  }

  // This constructor attaches a queue monitor to check for a 'stalled' condition and generates tidal statistics.
  // Delegates to the primary constructor so the tide parameters are validated exactly once.
  QueueTidal( size_t high_tide
              , size_t low_tide
              , std::string queue_name
              , size_t sample_frequency) : QueueTidal(high_tide, low_tide) {

    monitor_ptr_->launchStats(sample_frequency, std::move(queue_name));

  }
  ~QueueTidal() noexcept { monitor_ptr_ = nullptr; }

  // Non-copyable, non-movable (std::mutex and std::condition_variable members make the implicit
  // operations deleted anyway; declared explicitly for clarity).
  QueueTidal(const QueueTidal&) = delete;
  QueueTidal& operator=(const QueueTidal&) = delete;
  QueueTidal(QueueTidal&&) = delete;
  QueueTidal& operator=(QueueTidal&&) = delete;

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
  // These threads will block if the queue is empty.
  [[nodiscard]] T waitAndPop() {

    std::unique_lock<std::mutex> lock(queue_mutex_);
    empty_cond_.wait(lock, [this]()->bool{ return not empty(); });

    T value(std::move(queue_.front()));
    queue_.pop();

    --queue_size_;
    ++queue_activity_;

    // Track whether the tide transitioned so we can notify all blocked producers.
    bool transitioned_to_flood = false;

    if (queue_tidal_state_ == QueueTidalState::EBB_TIDE and queue_size_ <= low_tide_) {

      queue_tidal_state_ = QueueTidalState::FLOOD_TIDE;
      transitioned_to_flood = true;

    }

    lock.unlock();  // ~Mutex

    // When the tide flips from EBB to FLOOD, release ALL blocked producers so
    // they can refill the queue in parallel. notify_one() here would serialize
    // producer wakeups and drastically reduce throughput under high producer concurrency.
    // The gate avoids calling notify_all() on every pop when no transition occurred,
    // eliminating spurious wakeups (thundering-herd churn).
    if (transitioned_to_flood) {

      tide_cond_.notify_all();

    }

    return value;

  }

  void clear() {

    {
      std::scoped_lock<std::mutex> lock(queue_mutex_);
      // Swap with a temporary empty queue (noexcept for std::deque) and destroy the old elements when the
      // temporary goes out of scope. The assignment form (queue_ = {}) requires T to be move-assignable,
      // which the class does not guarantee.
      std::queue<T>().swap(queue_);
      queue_size_ = 0;
      queue_tidal_state_ = QueueTidalState::FLOOD_TIDE;
    }

    tide_cond_.notify_all();
    empty_cond_.notify_all();

  }

  [[nodiscard]] MonitorTidal<T>& monitor() const noexcept { return *monitor_ptr_; }
  // All of these functions are thread safe.
  [[nodiscard]] bool empty() const noexcept { return queue_size_ == 0; }
  [[nodiscard]] size_t size() const noexcept { return queue_size_; }
  [[nodiscard]] size_t activity() const noexcept { return queue_activity_; }
  [[nodiscard]] QueueTidalState queueTidalState() const noexcept { return queue_tidal_state_; }
  [[nodiscard]] bool queueState() const noexcept { return queue_tidal_state_ == QueueTidalState::FLOOD_TIDE; }

  [[nodiscard]] size_t highTide() const noexcept { return high_tide_; }
  [[nodiscard]] size_t lowTide() const noexcept { return low_tide_; }

private:

  // Tidal limits.
  size_t high_tide_;
  size_t low_tide_;

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
  std::unique_ptr<MonitorTidal<T>> monitor_ptr_;

  /// Ensures that the tidal parameters are valid. Invalid combinations (high tide <= low tide
  /// or zero high tide) would cause producers to deadlock or the queue to oscillate.
  /// A zero low tide is valid: producers are only re-enabled once the queue is fully drained.
  void validateTide() {

    if (high_tide_ == 0 or high_tide_ <= low_tide_) {

      // Invalid parameters would make the tidal queue unusable. Log and fall back to safe
      // defaults so that the queue continues to function.
      ExecEnv::log().warn( "QueueTidal received invalid tide parameters (high: {}, low: {}); using defaults (high: {}, low: {})."
                         , high_tide_, low_tide_, TIDAL_QUEUE_DEFAULT_HIGH_TIDE, TIDAL_QUEUE_DEFAULT_LOW_TIDE);

      // Since the members are not const, we can safely apply the validated defaults here.
      high_tide_ = TIDAL_QUEUE_DEFAULT_HIGH_TIDE;
      low_tide_ = TIDAL_QUEUE_DEFAULT_LOW_TIDE;

    }

  }


};


} // namespace

#endif // KEL_QUEUE_TIDAL_H
