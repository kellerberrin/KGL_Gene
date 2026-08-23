// Copyright 2023 Kellerberrin
//

#ifndef KEL_WORKFLOW_PIPELINE_H
#define KEL_WORKFLOW_PIPELINE_H

#include "kel_exec_env.h"
#include "kel_queue_tidal.h"

#include <algorithm>
#include <concepts>
#include <exception>
#include <functional>
#include <future>
#include <memory>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

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
// Stop tokens: there are two distinct null conventions.
//  - A null InputObject (e.g. nullptr when InputObject is std::unique_ptr<T>) pushed by the user via push() is a
//    *data* stop marker: it is wrapped in an ordinary QueuedFunctor, processed by the worker function (which should
//    forward it, e.g. return a null OutputObject) and appears on the output queue so consumers can detect the end.
//  - A null QueuedFunctor pointer is the *shutdown* token used only by joinThreads()/clear(); it circulates through
//    the worker threads (each thread re-queues it before exiting) and is never passed to the worker function.
//
// push() before activatePipeline() logs an error and drops the input: with no worker threads the functor could
// never execute. clear() joins (drains) the workers - processing every already-queued input - then discards the
// accumulated output queue and deactivates the pipeline; push() afterwards logs an error and drops the input until
// the pipeline is re-activated.
//
// waitAndPop() propagates worker exceptions: any exception thrown by the worker function is captured in the
// per-item promise and rethrown to the caller of waitAndPop() for that item.
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

    void operator()() {
      try {
        result_promise_.set_value(fn_pointer_->operator()(std::move(input_object_)));
      } catch (...) {
        result_promise_.set_exception(std::current_exception());
      }
    }
    [[nodiscard]] std::future<OutputObject> getFuture() noexcept { return result_promise_.get_future(); }

  private:

    WorkflowFuncPtr fn_pointer_; // This functional is held by all threads via a shared_ptr.
    InputObject input_object_;
    std::promise<OutputObject> result_promise_;

  };

public:

  explicit WorkflowPipeline(size_t high_tide = HIGH_TIDE_, size_t low_tide = LOW_TIDE_) : high_tide_(high_tide), low_tide_(low_tide) {}
  ~WorkflowPipeline() noexcept {

    // joinThreads() can throw (queue push, thread::join). Contain the failure so a noexcept
    // destructor never calls std::terminate for an ordinary exception.
    try {

      joinThreads();

    } catch (const std::exception& e) {

      ExecEnv::log().error("~WorkflowPipeline() caught exception: {}", e.what());

    } catch (...) {

      ExecEnv::log().error("~WorkflowPipeline() caught unknown exception");

    }

  }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe (or made so).
  // If the work function is a non-static class member then the first of the ...args should be a
  // pointer (MyClass* this) to the class instance.
  // The supplied function should be of the form 'auto f(args..., InputObject)->OutputObject'.
  template<typename F, typename... Args>
  bool activatePipeline(size_t threads, F&& f, Args&&... args)
  {

    // Stop any active threads and drain all queued input before installing the new function.
    joinThreads();

    auto callback_fn = std::bind_front(std::forward<F>(f), std::forward<Args>(args)...);
    // This callable object is shared by all threads. Stored atomically so that concurrent
    // push() calls and worker threads see a consistent pointer without data races.
    fn_pointer_.store(std::make_shared<WorkflowFunc>(std::move(callback_fn)),
                      std::memory_order_release);
    // Re-populate the thread pool.
    return queueThreads(threads);

  }

  [[nodiscard]] OutputObject waitAndPop() {

    auto future_output = output_queue_.waitAndPop();
    return future_output.get();

  }

  void push(InputObject input_object) {

    // Load the function pointer atomically to avoid data races with activatePipeline().
    auto fn_ptr = fn_pointer_.load(std::memory_order_acquire);

    if (not fn_ptr) {

      // A functor queued without worker threads would never execute and waitAndPop() would
      // block forever, so log an error and drop the input rather than enqueueing it.
      ExecEnv::log().error("WorkflowPipeline::push() called on an inactive pipeline; call activatePipeline() first");
      return;

    }

    auto func_ptr = std::make_unique<QueuedFunctor>(std::move(fn_ptr), std::move(input_object));
    auto future = func_ptr->getFuture();

    // Enqueue the future first so that consumers can always make progress. Then enqueue the
    // input functor so that workers can process it. Reversing this order can deadlock when
    // the input queue is at high tide and the output queue is also full.
    output_queue_.push(std::move(future));
    input_queue_.push(std::move(func_ptr));

  }

  void clear() {

    joinThreads();
    input_queue_.clear();
    output_queue_.clear();
    // Deactivate the pipeline: without worker threads any later push() could never execute.
    fn_pointer_.store(nullptr, std::memory_order_release);

  }

  // Access queue stats.
  [[nodiscard]] const QueueTidal<std::future<OutputObject>>& outputQueue() const noexcept { return output_queue_; }
  [[nodiscard]] const QueueTidal<std::unique_ptr<QueuedFunctor>>& inputQueue() const noexcept { return input_queue_; }

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
  // The supplied processing function, shared by all threads. Accessed via atomic
  // load/store to avoid data races during re-activation (concurrent push() vs activatePipeline()).
  std::atomic<WorkflowFuncPtr> fn_pointer_;


  void threadProlog() {

    while(true)
    {

      auto functor_ptr = input_queue_.waitAndPop();

      if (not functor_ptr) {

        // Shutdown token: re-queue it so the remaining worker threads also stop, then exit.
        // Exactly one token circulates; each thread re-queues it before terminating, so the
        // count is preserved.
        input_queue_.push(nullptr);
        break;

      }

      (*functor_ptr)();

    }

  }

  void joinThreads() {

    // If there are no workers there is nothing to drain: pushing a stop token would simply
    // block if the input queue is at high tide with nobody to consume it.
    if (not threads_.empty()) {

      input_queue_.push(nullptr);

    }

    for(auto& thread : threads_) {

      thread.join();

    }

    threads_.clear();
    input_queue_.clear();

  }

  bool queueThreads(size_t threads)
  {

    // Always have at least one worker thread.
    threads = std::max<size_t>(threads, 1);

    // Queue the worker threads. If thread construction fails partway through (resource
    // exhaustion), roll back the partially started pool so the vector never holds joinable
    // threads that would call std::terminate in its destructor.
    size_t started = 0;
    try {

      for (; started < threads; ++started) {

        threads_.emplace_back(&WorkflowPipeline::threadProlog, this);

      }

    } catch (...) {

      // The freshly started workers have not processed anything (activatePipeline drains
      // the input queue before queueThreads), so pushing one stop token per worker shuts
      // them all down; the tokens are cleared after the joins.
      for (size_t i = 0; i < started; ++i) {

        input_queue_.push(nullptr);

      }
      for (auto& thread : threads_) {

        thread.join();

      }
      threads_.clear();
      input_queue_.clear();
      throw;

    }

    return true;

  }

};


} // namespace




#endif // KEL_WORKFLOW_PIPELINE_H
