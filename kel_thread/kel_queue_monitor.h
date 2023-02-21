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

#ifndef KEL_QUEUE_MONITOR_H
#define KEL_QUEUE_MONITOR_H

#include "kel_exec_env.h"

#include <mutex>
#include <condition_variable>
#include <thread>
#include <chrono>


namespace kellerberrin {   //  organization level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The QueueMtSafe monitor collects queue statistics to facilitate optimal producer-consumer thread utilization.
// The monitor also detects stalled queue conditions - inactive consumer threads.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Forward queue declaration.
template<typename T> requires std::move_constructible<T> class QueueMtSafe;

// Realtime queue monitor
template<typename T> class MonitorMtSafe {

public:

  explicit MonitorMtSafe(QueueMtSafe<T> *queue_ptr) : queue_ptr_(queue_ptr) {}
  ~MonitorMtSafe() {

    stopStats();
    if (queue_samples_ > MIN_SAMPLES_) {

      displayStats();

    }

  }


  [[nodiscard]] size_t sampleFrequency() const { return sample_milliseconds_; }
  [[nodiscard]] size_t cumulativeQueueSize() const { return cumulative_queue_size_; }
  [[nodiscard]] size_t queueSamples() const { return queue_samples_; }

  void launchStats( size_t sample_milliseconds
                  , std::string queue_name = DEFAULT_QUEUE_NAME
                  , bool monitor_stalled = true) {

    queue_name_ = (std::move(queue_name));
    sample_milliseconds_ = sample_milliseconds;
    monitor_stalled_ = monitor_stalled;

    if (stats_thread_ptr_) {

      stopStats();

    }

    if (sample_milliseconds_ != DISABLE_QUEUE_MONITOR) {

      stats_thread_ptr_ = std::move(std::make_unique<std::thread>(&MonitorMtSafe::SampleQueue, this));
      ExecEnv::log().info("Sampling queue: {}, every milliseconds: {}", queue_name_, sample_milliseconds_);

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

  void displayQueueStats() const {

    if (terminate_flag_) {

      ExecEnv::log().info("Monitored Queue: {}, monitor is not active, no queue statistics available.");

    } else if (queue_samples_ < MIN_SAMPLES_) {

      ExecEnv::log().info( "Monitored Queue: {}, Sample Interval (ms): {}, Samples: {} (min {}); Insufficient for reliable statistics."
                         , queue_name_, sample_milliseconds_, queue_samples_, MIN_SAMPLES_);

    } else {

      displayStats();

    }

  }


  constexpr static const size_t DISABLE_QUEUE_MONITOR{0};
  constexpr static const char* DEFAULT_QUEUE_NAME{"QueueMtSafe"};

private:

  QueueMtSafe<T> *queue_ptr_;
  size_t sample_milliseconds_{0};
  bool monitor_stalled_{true};
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

  void displayStats() const {

    ExecEnv::log().info("Monitored Queue: {}, Sample Interval (ms): {}, Samples: {}; Average Queue Size: {}",
                        queue_name_, sample_milliseconds_, queue_samples_, static_cast<size_t>(averageSize()));

  }

  void SampleQueue() {

    while (true) {

      { // Lock block.
        std::unique_lock<std::mutex> lock(stats_mutex_);
        stats_condition_.wait_for(lock, std::chrono::milliseconds(sample_milliseconds_));
      }
      if (terminate_flag_) break;

      ++queue_samples_;
      size_t sample_size = queue_ptr_->size();
      size_t sample_activity = queue_ptr_->activity();
      cumulative_queue_size_ += sample_size;

      // Check for stalled queues (deadlock).
      if (monitor_stalled_ and previous_activity_ == sample_activity) {

        if (not queue_ptr_->empty()) {

          ++previous_count_;

        }
        if (previous_count_ >= WARN_INACTIVE_COUNT_) {

          ExecEnv::log().warn( "Monitor Queue: {} Size: {} Stalled (no consumer activity) for milliseconds: {}"
                             , queue_name_, sample_size, (queue_samples_ * sample_milliseconds_));

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

constexpr static const size_t TIDAL_QUEUE_MONITOR_DISABLE{0};
constexpr static const bool TIDAL_QUEUE_MONITOR_STALL{true};
constexpr static const char* TIDAL_QUEUE_DEFAULT_NAME{"QueueTidal"};

// Forward queue declaration.
template<typename T> requires std::move_constructible<T> class QueueTidal;

// Realtime queue monitor
template<typename T> class MonitorTidal {

public:

  explicit MonitorTidal(QueueTidal<T> *queue_ptr) : queue_ptr_(queue_ptr) {}
  ~MonitorTidal() {

    stopStats();
    if (queue_samples_ > MIN_SAMPLES_) {

      displayStats();

    }

  }


  [[nodiscard]] size_t sampleFrequency() const { return sample_milliseconds_; }

  [[nodiscard]] size_t cumulativeQueueSize() const { return cumulative_queue_size_; }

  [[nodiscard]] size_t queueSamples() const { return queue_samples_; }

  void launchStats(  size_t sample_milliseconds
                   , std::string queue_name = TIDAL_QUEUE_DEFAULT_NAME
                   , bool monitor_stalled = TIDAL_QUEUE_MONITOR_STALL) {

    sample_milliseconds_ = sample_milliseconds;
    monitor_stalled_ = monitor_stalled;
    queue_name_ = std::move(queue_name);
    empty_size_ = static_cast<size_t>(queue_ptr_->highTide() * EMPTY_PROPORTION_);

    if (stats_thread_ptr_) {

      stopStats();

    }

    if (sample_milliseconds_ != TIDAL_QUEUE_MONITOR_DISABLE) {

      stats_thread_ptr_ = std::move(std::make_unique<std::thread>(&MonitorTidal::SampleQueue, this));
      ExecEnv::log().info("Sampling queue: {}; every milliseconds: {}", queue_name_, sample_milliseconds_);

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

  void displayQueueStats() const {

    if (terminate_flag_) {

      ExecEnv::log().info("Queue monitor is not active, no queue statistics available.");

    } else if (queue_samples_ < MIN_SAMPLES_) {

      ExecEnv::log().info( "Monitored Queue: {}, Sample Interval (ms): {}, Samples: {} (min {}); Insufficient for reliable statistics."
          , queue_name_, sample_milliseconds_, queue_samples_, MIN_SAMPLES_);

    } else {

      displayStats();

    }

  }

private:

  QueueTidal<T> *queue_ptr_;
  size_t sample_milliseconds_{0};
  bool monitor_stalled_{true};
  std::string queue_name_;
  constexpr static const double EMPTY_PROPORTION_ = 0.1; // Queue at 10% of high tide is considered empty.
  size_t empty_size_{0};
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

  void displayStats() const {

    ExecEnv::log().info( "Queue Name: {},  High Tide: {}, Low Tide: {};  Samples: {}  Flood Tide: {:.2f}%, Ebb Tide: {:.2f}%; Empty (<={}): {:.2f}%, Av. Util.: ({}) {:.2f}%"
        , queue_name_,   queue_ptr_->highTide(), queue_ptr_->lowTide()
        , queueSamples(), (floodTide() * 100.0), (ebbingTide() * 100.0)
        , empty_size_, (averageEmpty() * 100.0), static_cast<size_t>(averageSize()), avUtilization());

  }

  void SampleQueue() {

    while (true) {

      { // Lock block.
        std::unique_lock<std::mutex> lock(stats_mutex_);
        stats_condition_.wait_for(lock, std::chrono::milliseconds(sample_milliseconds_));
      }
      if (terminate_flag_) break;

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
        if (monitor_stalled_ and previous_count_ >= WARN_INACTIVE_COUNT_) {

          ExecEnv::log().info( "Stalled Queue: {}, Size: {}, Queue State: {}, Stalled (no consumer activity) for milliseconds: {}"
                             , queue_name_, sample_size, (queue_ptr_->queueState() ? "Flood Tide" : "Ebb Tide")
                             , (queue_samples_ * sample_milliseconds_));

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
