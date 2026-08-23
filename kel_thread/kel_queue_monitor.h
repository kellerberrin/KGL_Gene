// Copyright 2023 Kellerberrin
//

#ifndef KEL_QUEUE_MONITOR_H
#define KEL_QUEUE_MONITOR_H

#include "kel_exec_env.h"

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <memory>
#include <mutex>
#include <string>
#include <thread>


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
  ~MonitorMtSafe() noexcept {

    // stopStats() joins a thread and displayStats() calls into the logging framework; both
    // can throw. Contain the failure so a noexcept destructor never calls std::terminate
    // for an ordinary exception (a throwing logger inside the catch still terminates, but
    // that is a far narrower window than an unguarded path).
    try {

      stopStats();
      if (queue_samples_ > MIN_SAMPLES_) {

        displayStats();

      }

    } catch (const std::exception& e) {

      ExecEnv::log().error("~MonitorMtSafe() caught exception: {}", e.what());

    } catch (...) {

      ExecEnv::log().error("~MonitorMtSafe() caught unknown exception");

    }

  }


  [[nodiscard]] size_t sampleFrequency() const noexcept { return sample_milliseconds_; }
  [[nodiscard]] size_t cumulativeQueueSize() const noexcept { return cumulative_queue_size_; }
  [[nodiscard]] size_t queueSamples() const noexcept { return queue_samples_; }

  void launchStats( size_t sample_milliseconds
                  , std::string queue_name = DEFAULT_QUEUE_NAME
                  , bool monitor_stalled = true) {

    // Serialize against concurrent stopStats()/launchStats() to prevent double-join.
    std::scoped_lock lock(lifecycle_mutex_);

    // Stop any existing monitor thread before rewriting the shared members, so a still
    // running SampleQueue cannot observe them mid-update (data race).
    if (stats_thread_ptr_) {

      stopStatsLocked();

    }

    queue_name_ = std::move(queue_name);
    sample_milliseconds_ = sample_milliseconds;
    monitor_stalled_ = monitor_stalled;

    if (sample_milliseconds_ != DISABLE_QUEUE_MONITOR) {

      // stopStats() leaves terminate_flag_ set; clear it or the new thread exits immediately.
      // Also reset the stall-detection accumulators so a relaunch does not carry over stale
      // state from the previous run (a near-threshold stall would otherwise trip a false
      // warning within the first few samples).
      terminate_flag_ = false;
      previous_activity_ = 0;
      stall_start_sample_ = 0;
      stall_warned_ = false;
      stats_thread_ptr_ = std::make_unique<std::thread>(&MonitorMtSafe::SampleQueue, this);
      ExecEnv::log().info("Sampling queue: {}, every milliseconds: {}", queue_name_, sample_milliseconds_);

    }

  }

  void stopStats() {

    std::scoped_lock lock(lifecycle_mutex_);
    stopStatsLocked();

  }

private:

  // Precondition: caller holds lifecycle_mutex_.
  void stopStatsLocked() {

    terminate_flag_ = true;
    stats_condition_.notify_one();
    if (stats_thread_ptr_) {

      stats_thread_ptr_->join();
      stats_thread_ptr_ = nullptr;

    }

  }

public:

  void displayQueueStats() const {

    if (terminate_flag_) {

      ExecEnv::log().info("Monitored Queue: {}, monitor is not active, no queue statistics available.", queue_name_);

    } else if (queue_samples_ < MIN_SAMPLES_) {

      ExecEnv::log().info( "Monitored Queue: {}, Sample Interval (ms): {}, Samples: {} (min {}); Insufficient for reliable statistics."
                         , queue_name_, sample_milliseconds_, queue_samples_.load(), MIN_SAMPLES_);

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

  // Serializes launchStats()/stopStats() to prevent concurrent double-join of the
  // sampling thread.
  std::mutex lifecycle_mutex_;
  std::mutex stats_mutex_;
  std::condition_variable stats_condition_;
  std::atomic<bool> terminate_flag_{false};
  size_t previous_activity_{0};
  std::atomic<size_t> cumulative_queue_size_{0};
  std::atomic<size_t> queue_samples_{0};
  // Stall tracking: start sample of the current stall episode (0 = not stalling).
  size_t stall_start_sample_{0};
  // Rate-limit: warn once per stall episode.
  bool stall_warned_{false};

  // Somewhat arbitrary but should work in most cases.
  constexpr static const size_t MIN_SAMPLES_{100};
  constexpr static const size_t WARN_INACTIVE_COUNT_{100};


  [[nodiscard]] double averageSize() const {

    if (queue_samples_ > 0) {

      return static_cast<double>(cumulative_queue_size_) / static_cast<double>(queue_samples_.load());

    }

    return 0.0;

  }

  void displayStats() const {

    ExecEnv::log().info("Monitored Queue: {}, Sample Interval (ms): {}, Samples: {}; Average Queue Size: {}",
                        queue_name_, sample_milliseconds_, queue_samples_.load(), static_cast<size_t>(averageSize()));

  }

  void SampleQueue() {

    while (true) {

      { // Lock block.
        std::unique_lock<std::mutex> lock(stats_mutex_);
        // Predicate returns true if notified (terminate); false on timeout.
        // This absorbs spurious wakeups: we only sample on a real timeout.
        stats_condition_.wait_for(lock, std::chrono::milliseconds(sample_milliseconds_),
                                  [this]->bool { return terminate_flag_.load(); });
      }
      if (terminate_flag_) break;

      ++queue_samples_;
      size_t sample_size = queue_ptr_->size();
      size_t sample_activity = queue_ptr_->activity();
      cumulative_queue_size_ += sample_size;

      // Check for stalled queues (deadlock). Warn at most once per stall episode to
      // avoid flooding the log with repeated warnings for the same stall.
      if (monitor_stalled_ and previous_activity_ == sample_activity) {

        if (not queue_ptr_->empty()) {

          if (stall_start_sample_ == 0) {

            stall_start_sample_ = queue_samples_;

          }

        }
        if (stall_start_sample_ > 0
            and (queue_samples_ - stall_start_sample_) >= WARN_INACTIVE_COUNT_
            and not stall_warned_) {

          ExecEnv::log().warn( "Monitor Queue: {} Size: {} Stalled (no consumer activity) for milliseconds: {}"
                             , queue_name_, sample_size
                             , ((queue_samples_ - stall_start_sample_) * sample_milliseconds_));
          stall_warned_ = true;

        }

      } else {

        stall_start_sample_ = 0;
        stall_warned_ = false;

      }

      previous_activity_ = sample_activity;

    }

  }

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The tidal queue monitor collects queue statistics to facilitate optimal producer-consumer thread utilization.
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
  ~MonitorTidal() noexcept {

    try {

      stopStats();
      if (queue_samples_ > MIN_SAMPLES_) {

        displayStats();

      }

    } catch (const std::exception& e) {

      ExecEnv::log().error("~MonitorTidal() caught exception: {}", e.what());

    } catch (...) {

      ExecEnv::log().error("~MonitorTidal() caught unknown exception");

    }

  }


  [[nodiscard]] size_t sampleFrequency() const noexcept { return sample_milliseconds_; }

  [[nodiscard]] size_t cumulativeQueueSize() const noexcept { return cumulative_queue_size_; }

  [[nodiscard]] size_t queueSamples() const noexcept { return queue_samples_; }

  void launchStats(  size_t sample_milliseconds
                   , std::string queue_name = TIDAL_QUEUE_DEFAULT_NAME
                   , bool monitor_stalled = TIDAL_QUEUE_MONITOR_STALL) {

    // Serialize against concurrent stopStats()/launchStats() to prevent double-join.
    std::scoped_lock lock(lifecycle_mutex_);

    // Stop any existing monitor thread before rewriting the shared members (data race).
    if (stats_thread_ptr_) {

      stopStatsLocked();

    }

    sample_milliseconds_ = sample_milliseconds;
    monitor_stalled_ = monitor_stalled;
    queue_name_ = std::move(queue_name);
    empty_size_ = static_cast<size_t>(queue_ptr_->highTide() * EMPTY_PROPORTION_);

    if (sample_milliseconds_ != TIDAL_QUEUE_MONITOR_DISABLE) {

      // stopStats() leaves terminate_flag_ set; clear it or the new thread exits immediately.
      // Reset the stall-detection accumulators so a relaunch does not carry over stale state.
      terminate_flag_ = false;
      previous_activity_ = 0;
      stall_start_sample_ = 0;
      stall_warned_ = false;
      stats_thread_ptr_ = std::make_unique<std::thread>(&MonitorTidal::SampleQueue, this);
      ExecEnv::log().info("Sampling queue: {}; every milliseconds: {}", queue_name_, sample_milliseconds_);

    }

  }

  void stopStats() {

    std::scoped_lock lock(lifecycle_mutex_);
    stopStatsLocked();

  }

private:

  // Precondition: caller holds lifecycle_mutex_.
  void stopStatsLocked() {

    terminate_flag_ = true;
    stats_condition_.notify_one();
    if (stats_thread_ptr_) {

      stats_thread_ptr_->join();
      stats_thread_ptr_ = nullptr;

    }

  }

public:

  void displayQueueStats() const {

    if (terminate_flag_) {

      ExecEnv::log().info("Queue monitor is not active, no queue statistics available.");

    } else if (queue_samples_ < MIN_SAMPLES_) {

      ExecEnv::log().info( "Monitored Queue: {}, Sample Interval (ms): {}, Samples: {} (min {}); Insufficient for reliable statistics."
          , queue_name_, sample_milliseconds_, queue_samples_.load(), MIN_SAMPLES_);

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

  // Serializes launchStats()/stopStats() to prevent concurrent double-join of the
  // sampling thread.
  std::mutex lifecycle_mutex_;
  std::mutex stats_mutex_;
  std::condition_variable stats_condition_;
  std::atomic<bool> terminate_flag_{false};
  std::atomic<size_t> low_tide_count_{0};
  std::atomic<size_t> high_tide_count_{0};
  std::atomic<size_t> inter_tidal_count_{0};
  std::atomic<size_t> empty_count_{0};
  size_t previous_activity_{0};
  // Stall tracking: start sample of the current stall episode (0 = not stalling).
  size_t stall_start_sample_{0};
  // Rate-limit: warn once per stall episode.
  bool stall_warned_{false};
  std::atomic<size_t> cumulative_queue_size_{0};
  std::atomic<size_t> queue_samples_{0};

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
        // Predicate returns true if notified (terminate); false on timeout.
        // This absorbs spurious wakeups: we only sample on a real timeout.
        stats_condition_.wait_for(lock, std::chrono::milliseconds(sample_milliseconds_),
                                  [this]->bool { return terminate_flag_.load(); });
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

      // Check for stalled queues (deadlock). Warn at most once per stall episode to
      // avoid flooding the log with repeated warnings for the same stall.
      if (previous_activity_ == sample_activity) {

        if (not queue_ptr_->empty()) {

          if (stall_start_sample_ == 0) {

            stall_start_sample_ = queue_samples_;

          }

        }
        if (monitor_stalled_ and stall_start_sample_ > 0
            and (queue_samples_ - stall_start_sample_) >= WARN_INACTIVE_COUNT_
            and not stall_warned_) {

          ExecEnv::log().info( "Stalled Queue: {}, Size: {}, Queue State: {}, Stalled (no consumer activity) for milliseconds: {}"
                             , queue_name_, sample_size, (queue_ptr_->queueState() ? "Flood Tide" : "Ebb Tide")
                             , ((queue_samples_ - stall_start_sample_) * sample_milliseconds_));
          stall_warned_ = true;

        }

      } else {

        stall_start_sample_ = 0;
        stall_warned_ = false;

      }

      previous_activity_ = sample_activity;

    }

  }


};


} // namespace



#endif //KEL_QUEUE_MONITOR_H
