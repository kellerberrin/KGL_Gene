///
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
//

#ifndef KEL_WORKFLOW_PIPELINE_H
#define KEL_WORKFLOW_PIPELINE_H


#include "kel_queue_tidal.h"

#include <future>
#include <vector>
#include <optional>
#include <thread>
#include <fstream>


namespace kellerberrin {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A pipeline is similar to a thread pool, the difference being is that homogenous InputObjects are enqueued
// these are then transformed by a single supplied function 'auto f(args..., InputObject)->OutputObject)'
// using multiple threads and the resultant OutputObjects can then be dequeued.
// The Enqueue and Dequeue operations are thread safe, however the sequential ordering of input/output objects
// is not guaranteed if multiple threads are used to Enqueue and Dequeue.
// Conversely, if single threads are used to Enqueue and Dequeue objects, then sequential ordering of input-output
// objects is guaranteed.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename InputObject, typename OutputObject>
requires std::move_constructible<InputObject> && std::move_constructible<OutputObject>
class WorkflowPipeline {

  using WorkflowFunc = std::move_only_function<OutputObject(InputObject)>;
  using WorkflowFuncPtr = std::shared_ptr<WorkflowFunc>;

  // Simple functor object is queued and consumed by the active threads.
  class QueuedFunctor {

  public:

    QueuedFunctor(WorkflowFuncPtr fn_pointer, InputObject&& input_object) : fn_pointer_(fn_pointer), input_object_(std::move(input_object)) {}
    QueuedFunctor(QueuedFunctor&& queued_functor) = delete;
    QueuedFunctor(const QueuedFunctor& queued_functor) = delete;
    ~QueuedFunctor() = default;

    void operator()() { result_promise_.set_value(fn_pointer_->operator()(std::move(input_object_))); }
    [[nodiscard]] std::future<OutputObject> getFuture() { return result_promise_.get_future(); }

  private:

    WorkflowFuncPtr fn_pointer_; // This functional is held by all threads via a shared_ptr.
    InputObject input_object_;
    std::promise<OutputObject> result_promise_;

  };

public:

  explicit WorkflowPipeline(size_t high_tide = HIGH_TIDE_, size_t low_tide = LOW_TIDE_) : high_tide_(high_tide), low_tide_(low_tide) {}
  ~WorkflowPipeline() { joinThreads(); }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe (or made so).
  // If the work function is a non-static class member then the first of the ...args should be a
  // pointer (MyClass* this) to the class instance.
  // The supplied function should be of the form 'auto f(args..., InputObject)->OutputObject'.
  template<typename F, typename... Args>
  bool activatePipeline(size_t threads, F&& f, Args&&... args)
  {

    // Clear any active threads.
    joinThreads();

    auto callback_fn = std::bind_front(std::forward<F>(f), std::forward<Args>(args)...);
    // This callable object is shared by all threads.
    fn_pointer_ = std::make_shared<WorkflowFunc>(std::move(callback_fn));
    // Re-populate the thread pool.
    return queueThreads(threads);

  }

  [[nodiscard]] OutputObject waitAndPop() {

    auto future_output = output_queue_.waitAndPop();
    return future_output.get();

  }

  void push(InputObject input_object) {

    auto func_ptr = std::make_unique<QueuedFunctor>(fn_pointer_, std::move(input_object));
    auto future = func_ptr->getFuture();
    input_queue_.push(std::move(func_ptr));
    output_queue_.push(std::move(future));

  }

  void clear() {

    joinThreads();
    input_queue_.clear();
    output_queue_.clear();

  }
  // Access queue stats.
  [[nodiscard]] const QueueTidal<std::future<OutputObject>>& outputQueue() const { return output_queue_; }
  [[nodiscard]] const QueueTidal<std::unique_ptr<QueuedFunctor>>& inputQueue() const { return input_queue_; }

private:

  // The default tidal IO queue parameters.
  static constexpr const size_t HIGH_TIDE_{10000};          // Maximum QueueTidal size
  static constexpr const size_t LOW_TIDE_{2000};            // Low water mark to begin queueing data records
  size_t high_tide_{HIGH_TIDE_};
  size_t low_tide_{LOW_TIDE_};
  // Tidal queue holds buffered output objects.
  QueueTidal<std::future<OutputObject>> output_queue_{high_tide_, low_tide_};
  // Tidal queue holds buffered input objects.
  QueueTidal<std::unique_ptr<QueuedFunctor>> input_queue_{high_tide_, low_tide_};
  // Thread Pool.
  std::vector<std::thread> threads_;
  // The supplied processing function held in the std::MoveFunction functional object.
  WorkflowFuncPtr fn_pointer_;


  void threadProlog() {

    while(true)
    {

      auto functor_ptr = input_queue_.waitAndPop();

      if (not functor_ptr) {

        input_queue_.push(nullptr);
        break;

      }

      (*functor_ptr)();

    }

  }

  void joinThreads() {

    input_queue_.push(nullptr);

    for(auto& thread : threads_) {

      thread.join();

    }

    threads_.clear();
    input_queue_.clear();

  }

  bool queueThreads(size_t threads)
  {

    // Always have at least one worker thread queued.
    threads = std::max<size_t>(threads, 1);

    // Queue the worker threads,
    for(size_t i = 0; i < threads; ++i) {

      threads_.emplace_back(&WorkflowPipeline::threadProlog, this);

    }

    return true;

  }

};


} // namespace




#endif //KEL_WORKFLOW_PIPELINE_H
