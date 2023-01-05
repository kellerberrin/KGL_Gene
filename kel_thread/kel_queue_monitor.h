//
// Created by kellerberrin on 30/12/22.
//

#ifndef KEL_QUEUE_MONITOR_H
#define KEL_QUEUE_MONITOR_H


#include "kel_exec_env.h"
#include "kel_mt_queue.h"

#include <chrono>


namespace kellerberrin {   //  organization level namespace



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The MtQueue monitor collects queue statistics to facilitate optimal producer-consumer thread utilization.
// The monitor also detects stalled queue conditions - blocked consumer threads.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Forward queue declaration.
template<typename T> requires std::move_constructible<T> class MtQueue;

// Realtime queue monitor
template<typename T> class MtQueueMonitor {

public:

  MtQueueMonitor() = default;
  ~MtQueueMonitor() {

    stopStats();
    if (queue_samples_ > MIN_SAMPLES_) {

      displayQueueStats();

    }

  }


  [[nodiscard]] size_t sampleFrequency() const { return sample_milliseconds_; }

  [[nodiscard]] size_t cumulativeQueueSize() const { return cumulative_queue_size_; }

  [[nodiscard]] size_t queueSamples() const { return queue_samples_; }

  void launchStats(MtQueue<T> *queue_ptr, size_t sample_milliseconds, std::string queue_name) {

    queue_ptr_ = queue_ptr;
    sample_milliseconds_ = sample_milliseconds;
    queue_name_ = std::move(queue_name);

    if (stats_thread_ptr_) {

      stopStats();

    }

    if (sample_milliseconds_ != DISABLE_QUEUE_MONITOR) {

      ExecEnv::log().info("Sampling queue: {} every: {} milliseconds", queue_name_, sample_milliseconds_);
      stats_thread_ptr_ = std::move(std::make_unique<std::thread>(&MtQueueMonitor::SampleQueue, this));

    }

  }

  void stopStats() {

    terminate_flag_ = true;
    stats_condition_.notify_one();
    if (stats_thread_ptr_) {

      stats_thread_ptr_->join();
      stats_thread_ptr_ = nullptr;

    }

  }

  constexpr static const size_t DISABLE_QUEUE_MONITOR{0};
  constexpr static const char* DEFAULT_QUEUE_NAME{"MtQueue"};

private:

  MtQueue<T> *queue_ptr_;
  size_t sample_milliseconds_;
  std::string queue_name_;
  std::unique_ptr<std::thread> stats_thread_ptr_;

  std::mutex stats_mutex_;
  std::condition_variable stats_condition_;
  std::atomic<bool> terminate_flag_{false};
  size_t previous_activity_{0};
  size_t cumulative_queue_size_{0};
  size_t queue_samples_{0};
  size_t previous_count_{0};

  // Somewhat arbitrary but should work in most cases.
  constexpr static const size_t MIN_SAMPLES_{100};
  constexpr static const size_t WARN_INACTIVE_COUNT_{100};


  [[nodiscard]] double averageSize() const {

    if (queue_samples_ > 0) {

      return static_cast<double>(cumulative_queue_size_) / static_cast<double>(queue_samples_);

    }

    return 0.0;

  }


  void displayQueueStats() const {

    ExecEnv::log().info("Queue: {}, Samples: {}; Av. Level: {:.2f}", queue_name_, queueSamples(),  averageSize());

  }


  void SampleQueue() {

    while (not terminate_flag_) {

      std::unique_lock<std::mutex> lock(stats_mutex_);
      stats_condition_.wait_for(lock, std::chrono::milliseconds(sample_milliseconds_));

      ++queue_samples_;
      size_t sample_size = queue_ptr_->size();
      size_t sample_activity = queue_ptr_->activity();
      cumulative_queue_size_ += sample_size;

      // Check for stalled queues (deadlock).
      if (previous_activity_ == sample_activity) {

        if (not queue_ptr_->empty()) {

          ++previous_count_;

        }
        if (previous_count_ >= WARN_INACTIVE_COUNT_) {

          ExecEnv::log().warn("Monitor Queue: {}, size: {}, FSM_State: {}, stalled for samples: {}",
                              queue_name_, sample_size, static_cast<size_t>(previous_count_));

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
// The bounded queue monitor collects queue statistics to facilitate optimal producer-consumer thread utilization.
// The monitor also detects stalled queue conditions - blocked consumer threads.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Forward queue declaration.
template<typename T> requires std::move_constructible<T> class BoundedMtQueue;

// Realtime queue monitor
template<typename T> class BoundedQueueMonitor {

public:

  BoundedQueueMonitor() = default;
  ~BoundedQueueMonitor() {

    stopStats();
    if (queue_samples_ > MIN_SAMPLES_) {

      displayQueueStats();

    }

  }


  [[nodiscard]] size_t sampleFrequency() const { return sample_milliseconds_; }

  [[nodiscard]] size_t cumulativeQueueSize() const { return cumulative_queue_size_; }

  [[nodiscard]] size_t queueSamples() const { return queue_samples_; }

  void launchStats(BoundedMtQueue<T> *queue_ptr, size_t sample_milliseconds, std::string queue_name) {

    queue_ptr_ = queue_ptr;
    sample_milliseconds_ = sample_milliseconds;
    queue_name_ = std::move(queue_name);
    empty_size_ = static_cast<size_t>(queue_ptr_->highTide() * EMPTY_PROPORTION_);

    if (stats_thread_ptr_) {

      stopStats();

    }

    if (sample_milliseconds_ != DISABLE_QUEUE_MONITOR) {

      ExecEnv::log().info("Sampling queue: {} every: {} milliseconds", queue_name_, sample_milliseconds_);
      stats_thread_ptr_ = std::move(std::make_unique<std::thread>(&BoundedQueueMonitor::SampleQueue, this));

    }

  }

  void stopStats() {

    terminate_flag_ = true;
    stats_condition_.notify_one();
    if (stats_thread_ptr_) {

      stats_thread_ptr_->join();
      stats_thread_ptr_ = nullptr;

    }

  }

  constexpr static const size_t DISABLE_QUEUE_MONITOR{0};
  constexpr static const char* DEFAULT_QUEUE_NAME{"BoundedMtQueue"};

private:

  BoundedMtQueue<T> *queue_ptr_;
  size_t sample_milliseconds_;
  std::string queue_name_;
  constexpr static const double EMPTY_PROPORTION_ = 0.1; // Queue at 10% of high tide is considered empty.
  size_t empty_size_;
  std::unique_ptr<std::thread> stats_thread_ptr_;

  std::mutex stats_mutex_;
  std::condition_variable stats_condition_;
  std::atomic<bool> terminate_flag_{false};
  size_t low_tide_count_{0};
  size_t high_tide_count_{0};
  size_t inter_tidal_count_{0};
  size_t empty_count_{0};
  size_t previous_count_{0};
  size_t previous_activity_{0};
  size_t cumulative_queue_size_{0};
  size_t queue_samples_{0};

  // Somewhat arbitrary but should work in most cases.
  constexpr static const size_t MIN_SAMPLES_{100};
  constexpr static const size_t WARN_INACTIVE_COUNT_{100};

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
                        queue_name_, queueSamples(), queue_ptr_->highTide(), (averageHighTide() * 100.0),
                        (floodTide() * 100.0), (ebbingTide() * 100.0), queue_ptr_->lowTide(), (averageLowTide() * 100.0),
                        empty_size_, (averageEmpty() * 100.0), avUtilization(), averageSize());

  }


  void SampleQueue() {

    while (not terminate_flag_) {

      std::unique_lock<std::mutex> lock(stats_mutex_);
      stats_condition_.wait_for(lock, std::chrono::milliseconds(sample_milliseconds_));

      ++queue_samples_;
      size_t sample_size = queue_ptr_->size();
      size_t sample_activity = queue_ptr_->activity();
      cumulative_queue_size_ += sample_size;
      if (sample_size <= queue_ptr_->lowTide()) {

        ++low_tide_count_;

      }
      if (sample_size <= empty_size_) {

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

        if (not queue_ptr_->empty()) {

          ++previous_count_;

        }
        if (previous_count_ >= WARN_INACTIVE_COUNT_) {

          ExecEnv::log().warn("Monitor Queue: {}, size: {}, FSM_State: {}, stalled for samples: {}",
                              queue_name_, sample_size, (queue_ptr_->queueState() ? "'active'" : "'inactive'"), static_cast<size_t>(previous_count_));

        }

      } else {

        previous_count_ = 0;

      }

      previous_activity_ = sample_activity;

    }

  }


};


} // namespace



#endif //KEL_QUEUE_MONITOR_H
