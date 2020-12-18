//
// Created by kellerberrin on 6/09/17.
//

#ifndef KEL_MT_QUEUE_H
#define KEL_MT_QUEUE_H

#include "kel_exec_env.h"

#include <queue>
#include <mutex>
#include <condition_variable>
#include <memory>

namespace kellerberrin {   //  organization level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// threadsafe queue for multiple consumers and producers
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename T> class MtQueue {

public:

  MtQueue() = default;
  ~MtQueue() = default;
  MtQueue(const MtQueue&) = delete;
  MtQueue(MtQueue&&) = delete;
  MtQueue& operator=(const MtQueue&) = delete;

  void push(T new_value) {

    std::lock_guard<std::mutex> lock(mutex_);

    data_queue_.push(std::move(new_value));

    data_cond_.notify_one();
  }

  void waitAndPop(T& value) {

    std::unique_lock<std::mutex> lock(mutex_);

    data_cond_.wait(lock,[this]{return not data_queue_.empty();});

    value = std::move(data_queue_.front());

    data_queue_.pop();

  }

  T waitAndPop() {

    std::unique_lock<std::mutex> lock(mutex_);

    data_cond_.wait(lock,[this]{return not data_queue_.empty();});

    T value = std::move(data_queue_.front());

    data_queue_.pop();

    return value;

  }

  bool tryPop(T& value) {

    std::lock_guard<std::mutex> lock(mutex_) ;

    if (data_queue_.empty()) return false;

    value = std::move(data_queue_.front());

    data_queue_.pop();

    return true;
  }

  std::shared_ptr<T> tryPop() {

    std::lock_guard<std::mutex> lock(mutex_) ;

    if (data_queue_.empty()) return std::shared_ptr<T>();

    std::shared_ptr<T> value_ptr(std::make_shared<T>(std::move(data_queue_.front())));

    data_queue_.pop();

    return value_ptr;
  }

  bool empty() const {

    std::lock_guard<std::mutex> lock(mutex_) ;

    return data_queue_.empty();
  }

  size_t size() const {

    std::lock_guard<std::mutex> lock(mutex_) ;

    return data_queue_.size();
  }

private:

  mutable std::mutex mutex_;
  std::queue<T> data_queue_;
  std::condition_variable data_cond_;

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

  BoundedMtQueue(size_t high_tide, size_t low_tide, std::string queue_name = DEFAULT_QUEUE_NAME,
                 size_t sample_frequency = DEFAULT_SAMPLE_FREQUENCY): high_tide_(high_tide),
                                                                      low_tide_(low_tide),
                                                                      queue_name_(std::move(queue_name)),
                                                                      sample_frequency_(sample_frequency),
                                                                      queue_state_(QueueState::Normal) {}
  ~BoundedMtQueue() { if (sample_frequency_ > 0) displayQueueStats(); }
  BoundedMtQueue(const BoundedMtQueue&) = delete;
  BoundedMtQueue(BoundedMtQueue&&) = delete;
  BoundedMtQueue& operator=(const BoundedMtQueue&) = delete;

  void push(T new_value) {

    if (sample_frequency_ > 0) {

      std::unique_lock<std::mutex> lock(mutex_);
      ++queue_count_;

      if ((queue_count_ % sample_frequency_) == 0) {

        cumulative_queue_size_ += mt_queue_.size();
        ++queue_samples_;

      }

    }

    if (queue_state_ == QueueState::Normal) {

      if (size() < high_tide_) {  // Possible race condition with size() >= high_tide is considered unimportant.

        mt_queue_.push(std::move(new_value));
        return;

      }
      else {

        queue_state_ = QueueState::HighTide;

      }
    }

    std::unique_lock<std::mutex> lock(mutex_);
    data_cond_.wait(lock,[this]{return size() <= low_tide_ or queue_state_ == QueueState::Normal;});

    mt_queue_.push(std::move(new_value));
    queue_state_ = QueueState::Normal;
    data_cond_.notify_one();

  }

  inline void waitAndPop(T& value) {

    mt_queue_.waitAndPop(value);
    data_cond_.notify_one();

  }

  inline T waitAndPop() {

    T value = std::move(mt_queue_.waitAndPop());
    data_cond_.notify_one();
    return value;

  }


  [[nodiscard]] bool empty() const { return mt_queue_.empty(); }

  [[nodiscard]] size_t size() const { return mt_queue_.size(); }

  [[nodiscard]] size_t highTide() const { return high_tide_; }

  [[nodiscard]] size_t lowTide() const { return low_tide_; }

  [[nodiscard]] size_t sampleFrequency() const { return sample_frequency_; }

  [[nodiscard]] size_t queueCount() const { return queue_count_; }

  [[nodiscard]] size_t cumulativeQueueSize() const { return cumulative_queue_size_; }
  [[nodiscard]] size_t queueSamples() const { return queue_samples_; }
  [[nodiscard]] double averageSize() const {

    if (queue_samples_ > 0) {

      return static_cast<double>(cumulative_queue_size_) / static_cast<double>(queue_samples_);

    }

    return 0.0;

  }


  [[nodiscard]] const std::string& queueName() const { return queue_name_; }

  void displayQueueStats() const {

    ExecEnv::log().info("Monitor queue: {}, high tide: {}, low tide: {}, count: {}, samples: {}, average size: {}",
                        queueName(), highTide(), lowTide(), queueCount(), queueSamples(), averageSize());

  }

  constexpr static const char* DEFAULT_QUEUE_NAME{"BoundedMtQueue"};
  constexpr static const size_t DEFAULT_SAMPLE_FREQUENCY{0};

private:

  const size_t high_tide_;
  const size_t low_tide_;
  const std::string queue_name_;
  const size_t sample_frequency_;

  std::atomic<size_t> queue_count_{0};
  size_t cumulative_queue_size_{0};
  size_t queue_samples_{0};

  MtQueue<T> mt_queue_;
  enum class QueueState { Normal, HighTide } queue_state_;
  std::condition_variable data_cond_;
  mutable std::mutex mutex_;

};

}   // end namespace


#endif // KGL_MtQueue_H
