// MIT License
//
// Copyright (c) 2017
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
//
//
// Created by kellerberrin on 6/09/17.
//

#ifndef KGL_MT_QUEUE_H
#define KGL_MT_QUEUE_H

#include <queue>
#include <mutex>
#include <condition_variable>
#include <memory>

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

// threadsafe queue for multiple consumers and producers
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

  std::shared_ptr<T> waitAndPop() {

    std::unique_lock<std::mutex> lock(mutex_);

    data_cond_.wait(lock,[this]{return not data_queue_.empty();});

    std::shared_ptr<T> value_ptr(std::make_shared<T>(std::move(data_queue_.front())));

    data_queue_.pop();

    return value_ptr;
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

// The bounded multithreaded queue has a maximum of high_tide elements and a low tide
// when the producer(s) can again start pushing elements after a high tide event.
// This stops excessive memory usage (and swapping) if the producer(s) can queue records
// faster than consumer(s) can remove them.
// Note, if there are no active consumers then this queue will block forever on high tide.

template<typename T> class BoundedMtQueue {

public:

  BoundedMtQueue(size_t high_tide, size_t low_tide): high_tide_(high_tide)
                                                   , low_tide_(low_tide)
                                                   , queue_state_(QueueState::Normal) {}
  ~BoundedMtQueue() = default;
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

  inline bool empty() const { return mt_queue_.empty(); }

  inline size_t size() const { return mt_queue_.size(); }

private:

  size_t high_tide_;
  size_t low_tide_;

  MtQueue<T> mt_queue_;
  enum class QueueState { Normal, HighTide } queue_state_;
  std::condition_variable data_cond_;
  mutable std::mutex mutex_;

};

}   // namespace genome
}   // namespace kellerberrin


#endif // KGL_MtQueue_H
