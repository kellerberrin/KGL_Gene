//
// Created by kellerberrin on 19/12/20.
//

#ifndef KEL_BOUND_QUEUE_H
#define KEL_BOUND_QUEUE_H

#include "kel_mt_queue.h"
#include "kel_thread_pool.h"

#include <chrono>


namespace kellerberrin {   //  organization level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The bounded queue monitor collects queue statistics to facilitate optimal producer-consumer thread utilization.
// The monitor also detects stalled queue conditions - blocked consumer threads.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Forward queue declaration.
template<typename T> class BoundedMtQueue;

// Realtime queue monitor
template<typename T> class BoundedQueueMonitor {

public:

  BoundedQueueMonitor(BoundedMtQueue<T>* queue_ptr, size_t sample_milliseconds = DISABLE_QUEUE_MONITOR) :
    queue_ptr_(queue_ptr),
    sample_milliseconds_(sample_milliseconds) {

    launchStats();

  }
  ~BoundedQueueMonitor() {

    stopStats();
    if (queue_samples_ > MIN_SAMPLES_) {
      displayQueueStats();
    }

  }

  constexpr static const size_t DISABLE_QUEUE_MONITOR{0};

  [[nodiscard]] size_t sampleFrequency() const { return sample_milliseconds_; }
  [[nodiscard]] size_t cumulativeQueueSize() const { return cumulative_queue_size_; }
  [[nodiscard]] size_t queueSamples() const { return queue_samples_; }


private:

  BoundedMtQueue<T>* queue_ptr_;
  size_t sample_milliseconds_;
  ThreadPool stats_thread_{1};

  std::mutex stats_mutex_;
  std::condition_variable stats_condition_;
  std::atomic<bool> terminate_flag_{false};
  std::atomic<size_t> low_tide_count_{0};
  std::atomic<size_t> high_tide_count_{0};
  std::atomic<size_t> inter_tidal_count_{0};
  std::atomic<size_t> empty_count_{0};
  std::atomic<size_t> previous_count_{0};
  std::atomic<size_t> previous_activity_{0};
  std::atomic<size_t> cumulative_queue_size_{0};
  std::atomic<size_t> queue_samples_{0};

  // Somewhat arbitrary but should work in most cases.
  constexpr static const size_t MIN_SAMPLES_ {100};
  constexpr static const size_t WARN_INACTIVE_COUNT_ {100};

  [[nodiscard]] double averageHighTide() const {

    if (queue_samples_ > 0) {

      return static_cast<double>(high_tide_count_) / static_cast<double>(queue_samples_);

    }

    return 0.0;

  }
  [[nodiscard]] double averageLowTide() const {

    if (queue_samples_ > 0) {

      return static_cast<double>(low_tide_count_) / static_cast<double>(queue_samples_);

    }

    return 0.0;

  }
  [[nodiscard]] double ebbingTide() const {

    if (queue_samples_ > 0) {

      return static_cast<double>(inter_tidal_count_) / static_cast<double>(queue_samples_);

    }

    return 0.0;

  }
  [[nodiscard]] double floodTide() const { return 1.0 - ebbingTide(); }
  [[nodiscard]] double averageEmpty() const {

    if (queue_samples_ > 0) {

      return static_cast<double>(empty_count_) / static_cast<double>(queue_samples_);

    }

    return 0.0;

  }
  [[nodiscard]] double averageSize() const {

    if (queue_samples_ > 0) {

      return static_cast<double>(cumulative_queue_size_) / static_cast<double>(queue_samples_);

    }

    return 0.0;

  }

  [[nodiscard]] double avUtilization() const { return (averageSize() * 100.0) / static_cast<double>(queue_ptr_->highTide()); }

  void displayQueueStats() const {

    ExecEnv::log().info("Queue: {}, Samples: {}; High Tide ({}): {:.2f}%, Flood Tide: {:.2f}%, Ebbing Tide: {:.2f}%, Low Tide ({}): {:.2f}%, Empty (<={}): {:.2f}%, Av. Level: {:.2f}% ({:.2f})",
                        queue_ptr_->queueName(), queueSamples(), queue_ptr_->highTide(), (averageHighTide() * 100.0),
                        (floodTide() * 100.0), (ebbingTide() * 100.0),  queue_ptr_->lowTide(), (averageLowTide() * 100.0),
                        queue_ptr_->emptySize(), (averageEmpty() * 100.0), avUtilization(), averageSize());

  }

  void launchStats() {

    if (sample_milliseconds_ != DISABLE_QUEUE_MONITOR) {

      ExecEnv::log().info("Sampling queue: {} every: {} milliseconds", queue_ptr_->queueName(), sample_milliseconds_);
      stats_thread_.enqueueWork(&BoundedQueueMonitor::SampleQueue, this);

    }

  }

  void stopStats() {

    terminate_flag_ = true;
    stats_condition_.notify_one();

  }

  void SampleQueue() {

    while(not terminate_flag_) {

      std::unique_lock<std::mutex> lock(stats_mutex_);
      stats_condition_.wait_for(lock, std::chrono::milliseconds(sample_milliseconds_));

      ++queue_samples_;
      size_t sample_size = queue_ptr_->size();
      size_t sample_activity = queue_ptr_->activity();
      cumulative_queue_size_ += sample_size;
      if (sample_size <= queue_ptr_->lowTide()) {

        ++low_tide_count_;

      }
      if (sample_size <= queue_ptr_->emptySize()) {

        ++empty_count_;

      }
      if (sample_size >= queue_ptr_->highTide()) {

        ++high_tide_count_;

      }
      if (not queue_ptr_->queueState()) {

        ++inter_tidal_count_;

      }

      // Check for stalled queues (deadlock).
      if (previous_activity_ == sample_activity) {

        ++previous_count_;
        if (previous_count_ >= WARN_INACTIVE_COUNT_) {

          ExecEnv::log().warn("Monitor Queue: {}, size: {}, FSM_State: {}, stalled for samples: {}",
                              queue_ptr_->queueName(), sample_size, (queue_ptr_->queueState() ? "'active'" : "'inactive'"), static_cast<size_t>(previous_count_));

        }

      } else {

        previous_count_ = 0;

      }

      previous_activity_ = sample_activity;

    }

  }


};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The bounded multithreaded queue has a maximum of high_tide elements and a low tide
// when the producer(s) can again start pushing elements after a high tide event.
// This stops excessive memory usage (and swapping) if the producer(s) can queue records
// faster than consumer(s) can remove them.
// Note, if there are no active consumers then this queue will block forever on high tide.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T> class BoundedMtQueue {

public:

  BoundedMtQueue(size_t high_tide,
                 size_t low_tide,
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
