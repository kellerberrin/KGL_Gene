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
                                                                      queue_state_(QueueState::Normal) { launchStats(); }
  ~BoundedMtQueue() { stopStats(); if (queue_samples_ > MIN_SAMPLES_) displayQueueStats(); }
  BoundedMtQueue(const BoundedMtQueue&) = delete;
  BoundedMtQueue(BoundedMtQueue&&) = delete;
  BoundedMtQueue& operator=(const BoundedMtQueue&) = delete;

  void push(T new_value) {


    if (queue_state_ == QueueState::Normal) {

      if (size() < high_tide_) {  // Possible race condition with size() >= high_tide is considered unimportant.

        mt_queue_.push(std::move(new_value));
        return;

      }
      else {

        queue_state_ = QueueState::HighTide;

      }
    }

    std::unique_lock<std::mutex> lock(tide_mutex_);
    tide_cond_.wait(lock, [this]{return size() <= low_tide_ or queue_state_ == QueueState::Normal;});

    mt_queue_.push(std::move(new_value));
    queue_state_ = QueueState::Normal;
    tide_cond_.notify_one();

  }

  inline void waitAndPop(T& value) {

    mt_queue_.waitAndPop(value);
    tide_cond_.notify_one();

  }

  inline T waitAndPop() {

    T value = std::move(mt_queue_.waitAndPop());
    tide_cond_.notify_one();
    return value;

  }


  [[nodiscard]] bool empty() const { return mt_queue_.empty(); }

  [[nodiscard]] size_t size() const { return mt_queue_.size(); }

  [[nodiscard]] size_t highTide() const { return high_tide_; }

  [[nodiscard]] size_t lowTide() const { return low_tide_; }

  [[nodiscard]] size_t sampleFrequency() const { return sample_frequency_; }

  [[nodiscard]] size_t cumulativeQueueSize() const { return cumulative_queue_size_; }
  [[nodiscard]] size_t queueSamples() const { return queue_samples_; }
  [[nodiscard]] double averageSize() const {

    if (queue_samples_ > 0) {

      return static_cast<double>(cumulative_queue_size_) / static_cast<double>(queue_samples_);

    }

    return 0.0;

  }

  [[nodiscard]] double avUtilization() const { return (averageSize() * 100.0) / static_cast<double>(high_tide_); }

  [[nodiscard]] const std::string& queueName() const { return queue_name_; }

  void displayQueueStats() const {

    ExecEnv::log().info("Monitor Queue: {}, High Tide: {}, Low Tide: {}, Samples: {}, Average Util: {}% ({})",
                        queueName(), highTide(), lowTide(), queueSamples(), avUtilization(), averageSize());

  }

  void launchStats() {

    if (sample_frequency_ > 0) {

      ExecEnv::log().info("Sampling queue: {} size every: {} seconds", queue_name_, sample_frequency_);
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

  MtQueue<T> mt_queue_;
  enum class QueueState { Normal, HighTide } queue_state_;
  std::condition_variable tide_cond_;
  mutable std::mutex tide_mutex_;

  ThreadPool queue_stats_{1};
  std::mutex stats_mutex_;
  std::condition_variable stats_cond_;
  std::atomic<bool> terminate_{false};
  std::atomic<size_t> cumulative_queue_size_{0};
  // Somewhat arbitrary.
  constexpr static const size_t MIN_SAMPLES_ {100};
  std::atomic<size_t> queue_samples_{0};


  void SampleQueue() {

    while(!terminate_) {

      std::unique_lock<std::mutex> lock(stats_mutex_);
      stats_cond_.wait_for( lock, std::chrono::seconds(sample_frequency_));

      cumulative_queue_size_ += mt_queue_.size();
      ++queue_samples_;

    }

  }

};



} // namespace

#endif //KEL_BOUND_QUEUE_H
