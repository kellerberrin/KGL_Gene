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

    std::scoped_lock lock(mutex_);

    data_queue_.push(std::move(new_value));

    data_cond_.notify_one();

    ++size_;
    ++activity_;

  }

  void waitAndPop(T& value) {

    std::unique_lock<std::mutex> lock(mutex_);

    data_cond_.wait(lock,[this]{return not size_ != 0;});

    value = std::move(data_queue_.front());

    data_queue_.pop();

    --size_;
    ++activity_;

  }

  [[nodiscard]] T waitAndPop() {

    std::unique_lock<std::mutex> lock(mutex_);

    data_cond_.wait(lock,[this]{return size_ != 0;});

    T value = std::move(data_queue_.front());

    data_queue_.pop();

    --size_;
    ++activity_;

    return value;

  }

  bool tryPop(T& value) {

    std::scoped_lock lock(mutex_) ;

    if (data_queue_.empty()) return false;

    value = std::move(data_queue_.front());

    data_queue_.pop();

    --size_;
    ++activity_;

    return true;
  }

  [[nodiscard]] std::shared_ptr<T> tryPop() {

    std::scoped_lock lock(mutex_) ;

    if (data_queue_.empty()) return std::shared_ptr<T>();

    std::shared_ptr<T> value_ptr(std::make_shared<T>(std::move(data_queue_.front())));

    data_queue_.pop();

    --size_;
    ++activity_;

    return value_ptr;

  }

  [[nodiscard]] bool empty() const { return size_ == 0; }

  [[nodiscard]] size_t size() const { return size_; }

  [[nodiscard]] size_t activity() const { return activity_; }

private:

  std::mutex mutex_;
  std::queue<T> data_queue_;
  std::condition_variable data_cond_;
  std::atomic<size_t> size_{0};
  std::atomic<size_t> activity_{0};

};


}   // end namespace


#endif // KGL_MtQueue_H
