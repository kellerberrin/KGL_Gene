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

#include "kel_mt_queue.h"

#include <mutex>
#include <condition_variable>
#include <thread>
#include <chrono>
#include <iostream>
#include <sstream>
#include <iomanip>

// Replace the spdlog based messaging system.
enum class MessageType { INFO, WARNING, ERROR};
void workflowStreamOut(MessageType type, const std::string& message);

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

  void launchStats( MtQueue<T> *queue_ptr
                   , size_t sample_milliseconds
                   , std::string queue_name = DEFAULT_QUEUE_NAME
                   , bool monitor_stalled = true) {

    queue_ptr_ = queue_ptr;
    sample_milliseconds_ = sample_milliseconds;
    queue_name_ = std::move(queue_name);
    monitor_stalled_ = monitor_stalled;

    if (stats_thread_ptr_) {

      stopStats();

    }

    if (sample_milliseconds_ != DISABLE_QUEUE_MONITOR) {

      stats_thread_ptr_ = std::move(std::make_unique<std::thread>(&MtQueueMonitor::SampleQueue, this));

      std::stringstream begin_message;
      begin_message << "Sampling queue: " << queue_name_ << " every milliseconds: " << sample_milliseconds_;
      workflowStreamOut(MessageType::INFO, begin_message.str());

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


  void displayQueueStats() const {

    std::stringstream stats_message;
    stats_message << "Queue Name: " <<  queue_name_
                  << ", Sample Interval (ms): " << sample_milliseconds_
                  << ", Samples: " << queue_samples_
                  << "; Average Queue Size:" << std::setprecision(2) << averageSize();
    workflowStreamOut(MessageType::INFO, stats_message.str());

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
      if (monitor_stalled_ and previous_activity_ == sample_activity) {

        if (not queue_ptr_->empty()) {

          ++previous_count_;

        }
        if (previous_count_ >= WARN_INACTIVE_COUNT_) {

          std::stringstream stalled_message;
          stalled_message << "Monitor Queue: " <<  queue_name_
                          << " Size: " << sample_size
                          << " Stalled (no consumer activity) for milliseconds: "
                          << (queue_samples_ * sample_milliseconds_);
          workflowStreamOut(MessageType::INFO, stalled_message.str());

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

  void launchStats(BoundedMtQueue<T> *queue_ptr
                   , size_t sample_milliseconds
                   , std::string queue_name = DEFAULT_QUEUE_NAME
                   , bool monitor_stalled = true) {

    queue_ptr_ = queue_ptr;
    sample_milliseconds_ = sample_milliseconds;
    monitor_stalled_ = monitor_stalled;
    queue_name_ = std::move(queue_name);
    empty_size_ = static_cast<size_t>(queue_ptr_->highTide() * EMPTY_PROPORTION_);

    if (stats_thread_ptr_) {

      stopStats();

    }

    if (sample_milliseconds_ != DISABLE_QUEUE_MONITOR) {

      stats_thread_ptr_ = std::move(std::make_unique<std::thread>(&BoundedQueueMonitor::SampleQueue, this));
      std::stringstream begin_message;
      begin_message << "Sampling queue: " << queue_name_ << " every milliseconds: " << sample_milliseconds_;
      workflowStreamOut(MessageType::INFO, begin_message.str());

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
  bool monitor_stalled_{true};
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

    std::stringstream stats_message;
    stats_message << "Queue Name: " << queue_name_
                  << ", Samples :" <<  queueSamples()
                  << "; High Tide (" << queue_ptr_->highTide() <<"): "
                  << std::setprecision(2) << (averageHighTide() * 100.0) << "%"
                  << ", Flood Tide: " << (floodTide() * 100.0) << "%"
                  << ", Ebbing Tide: " << (ebbingTide() * 100.0) << "%"
                  << ", Low Tide (" << std::setprecision(0) << queue_ptr_->lowTide() << "): "
                  << std::setprecision(2) << (averageLowTide() * 100.0) << "%"
                  << ", Empty (<=" << std::setprecision(0) << empty_size_ << "): "
                  << std::setprecision(2) << (averageEmpty() * 100.0) << "%"
                  << "Av. Size: (" << std::setprecision(0) <<  averageSize() << ") "
                  << std::setprecision(0) << avUtilization() << "%";
    workflowStreamOut(MessageType::INFO, stats_message.str());

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
        if (monitor_stalled_ and previous_count_ >= WARN_INACTIVE_COUNT_) {

          std::stringstream stalled_message;
          stalled_message << "Monitor Queue: " <<  queue_name_
                          << " Size: " << sample_size
                          << " Stalled (no consumer activity) for milliseconds: "
                          << (queue_samples_ * sample_milliseconds_);
          workflowStreamOut(MessageType::INFO, stalled_message.str());

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
