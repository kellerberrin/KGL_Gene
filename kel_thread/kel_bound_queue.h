//
// Created by kellerberrin on 19/12/20.
//

#ifndef KEL_BOUND_QUEUE_H
#define KEL_BOUND_QUEUE_H


#include "kel_mt_queue.h"
#include "kel_thread_pool.h"
#include "kel_bound_queue_monitor.h"


namespace kellerberrin {   //  organization level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The bounded multithread queue has a maximum of high_tide elements and a low tide
// when the producer(s) can again start pushing elements after a high tide event.
// This stops excessive memory usage (and swapping) if the producer(s) can queue records
// faster than consumer(s) can remove them.
// Note, if there are no active consumers then this queue will block forever on high tide.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T> class BoundedMtQueue {

public:

  BoundedMtQueue(size_t high_tide = DEFAULT_HIGH_TIDE,
                 size_t low_tide = DEFAULT_LOW_TIDE,
                 std::string queue_name = DEFAULT_QUEUE_NAME,
                 size_t sample_frequency = BoundedQueueMonitor<T>::DISABLE_QUEUE_MONITOR):
                  high_tide_(high_tide),
                  low_tide_(low_tide),
                  queue_name_(std::move(queue_name)),
                  sample_frequency_(sample_frequency),
                  empty_size_(high_tide_/EMPTY_PROPORTION_) {

    monitor_ptr_ = std::make_unique<BoundedQueueMonitor<T>>(this, sample_frequency_);
  }
  ~BoundedMtQueue() { monitor_ptr_ = nullptr; }
  BoundedMtQueue(const BoundedMtQueue&) = delete;
  BoundedMtQueue(BoundedMtQueue&&) = delete;
  BoundedMtQueue& operator=(const BoundedMtQueue&) = delete;

  void push(T&& new_value) {

    if (queue_state_) {

      if (size() < high_tide_) {  // Possible race condition with size() >= high_tide is considered unimportant.

        mt_queue_.push(std::forward<T>(new_value));
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


  [[nodiscard]] bool empty() const { return mt_queue_.empty(); }
  [[nodiscard]] size_t size() const { return mt_queue_.size(); }
  [[nodiscard]] size_t activity() const { return mt_queue_.activity(); }

  // Bool 'true' the queue is active and accepting producer threads.
  // Bool 'false' and producer threads are blocked.
  [[nodiscard]] bool queueState() const { return queue_state_; }

  [[nodiscard]] size_t highTide() const { return high_tide_; }

  [[nodiscard]] size_t lowTide() const { return low_tide_; }

  [[nodiscard]] size_t emptySize() const { return empty_size_; }

  [[nodiscard]] size_t sampleFrequency() const { return sample_frequency_; }

  [[nodiscard]] const std::string& queueName() const { return queue_name_; }

  constexpr static const char* DEFAULT_QUEUE_NAME{"BoundedMtQueue"};
  constexpr static const size_t DEFAULT_HIGH_TIDE{10000};
  constexpr static const size_t DEFAULT_LOW_TIDE{2000};

private:

  const size_t high_tide_;
  const size_t low_tide_;
  const std::string queue_name_;
  const size_t sample_frequency_;
  constexpr static const size_t EMPTY_PROPORTION_ = 10; // Queue at 10% is considered empty.
  const size_t empty_size_;

  MtQueue<T> mt_queue_;
  std::atomic<bool> queue_state_{true};
  std::condition_variable tide_cond_;
  mutable std::mutex tide_mutex_;

  // Held in a pointer for explicit object lifetime.
  std::unique_ptr<BoundedQueueMonitor<T>> monitor_ptr_;


};



} // namespace

#endif //KEL_BOUND_QUEUE_H
