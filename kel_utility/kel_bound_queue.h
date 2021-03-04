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
                 size_t sample_frequency = DEFAULT_SAMPLE_FREQUENCY): high_tide_(high_tide),
                                                                      low_tide_(low_tide),
                                                                      queue_name_(std::move(queue_name)),
                                                                      sample_frequency_(sample_frequency),
                                                                      empty_size_(high_tide_/EMPTY_PROPORTION_) { launchStats(); }
  ~BoundedMtQueue() { stopStats(); if (queue_samples_ > MIN_SAMPLES_) displayQueueStats(); }
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
      mt_queue_.push(std::move(new_value));
    }

    tide_cond_.notify_one();

  }

  T waitAndPop() {

    T value = mt_queue_.waitAndPop();
    if (not queue_state_) {

      queue_state_ = size() <= low_tide_;

    }
    tide_cond_.notify_one();
    return value;

  }


  [[nodiscard]] bool empty() const { return mt_queue_.empty(); }

  [[nodiscard]] size_t size() const { return mt_queue_.size(); }

  [[nodiscard]] size_t highTide() const { return high_tide_; }

  [[nodiscard]] size_t lowTide() const { return low_tide_; }

  [[nodiscard]] size_t emptySize() const { return empty_size_; }

  [[nodiscard]] size_t sampleFrequency() const { return sample_frequency_; }

  [[nodiscard]] size_t cumulativeQueueSize() const { return cumulative_queue_size_; }
  [[nodiscard]] size_t queueSamples() const { return queue_samples_; }

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
  [[nodiscard]] double averageInterTidal() const {

    if (queue_samples_ > 0) {

      return static_cast<double>(inter_tidal_count_) / static_cast<double>(queue_samples_);

    }

    return 0.0;

  }
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

  [[nodiscard]] double avUtilization() const { return (averageSize() * 100.0) / static_cast<double>(high_tide_); }

  [[nodiscard]] const std::string& queueName() const { return queue_name_; }

  void displayQueueStats() const {

    ExecEnv::log().info("Queue: {}, Samples: {}; High Tide ({}): {:.2f}%, Inter Tidal: {:.2f}%, Low Tide ({}): {:.2f}%, Empty ({<={}}): {:.2f}%, Av. Level: {:.2f}% ({:.2f})",
                        queueName(), queueSamples(), highTide(), (averageHighTide() * 100.0),
                        (averageInterTidal() * 100.0), lowTide(), (averageLowTide() * 100.0),
                        emptySize(), (averageEmpty() * 100.0), avUtilization(), averageSize());

  }

  void launchStats() {

    if (sample_frequency_ > 0) {

      ExecEnv::log().info("Sampling queue: {} every: {} milliseconds", queue_name_, sample_frequency_);
      queue_stats_.enqueueWork(&BoundedMtQueue::SampleQueue, this);

    }

  }

  void stopStats() {

    terminate_ = true;
    stats_cond_.notify_one();

  }


  constexpr static const char* DEFAULT_QUEUE_NAME{"BoundedMtQueue"};
  constexpr static const size_t DEFAULT_SAMPLE_FREQUENCY{0};

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

  ThreadPool queue_stats_{1};
  std::mutex stats_mutex_;
  std::condition_variable stats_cond_;
  std::atomic<bool> terminate_{false};
  std::atomic<size_t> low_tide_count_{0};
  std::atomic<size_t> high_tide_count_{0};
  std::atomic<size_t> inter_tidal_count_{0};
  std::atomic<size_t> empty_count_{0};
  std::atomic<size_t> previous_count_{0};
  std::atomic<size_t> previous_activity_{0};
  std::atomic<size_t> cumulative_queue_size_{0};
  // Somewhat arbitrary.
  constexpr static const size_t MIN_SAMPLES_ {100};
  constexpr static const size_t WARN_INACTIVE_COUNT_ {100};
  std::atomic<size_t> queue_samples_{0};


  void SampleQueue() {

    while(!terminate_) {

      std::unique_lock<std::mutex> lock(stats_mutex_);
      stats_cond_.wait_for( lock, std::chrono::milliseconds(sample_frequency_));

      ++queue_samples_;
      size_t sample_size = mt_queue_.size();
      size_t activity = mt_queue_.activity();
      cumulative_queue_size_ += sample_size;
      if (sample_size <= low_tide_) {

        ++low_tide_count_;

      }
      if (sample_size <= empty_size_) {

        ++empty_count_;

      }
      if (sample_size >= high_tide_) {

        ++high_tide_count_;

      }
      if (not queue_state_) {

        ++inter_tidal_count_;

      }

      // Check for stalled queues (deadlock).
      if (previous_activity_ == activity) {

        ++previous_count_;
        if (previous_count_ >= WARN_INACTIVE_COUNT_) {

          ExecEnv::log().warn("Monitor Queue: {}, size: {} stalled for samples: {}", queueName(), sample_size, static_cast<size_t>(previous_count_));

        }

      } else {

        previous_count_ = 0;

      }

      previous_activity_ = activity;

    }

  }

};



} // namespace

#endif //KEL_BOUND_QUEUE_H
