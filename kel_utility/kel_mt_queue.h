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


}   // end namespace


#endif // KGL_MtQueue_H
