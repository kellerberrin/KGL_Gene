//
// Created by kellerberrin on 7/02/23.
//

#ifndef KEL_WORKFLOW_FUTURE_H
#define KEL_WORKFLOW_FUTURE_H

#include "kel_queue_tidal.h"
#include "kel_workflow_threads.h"

#include <functional>
#include <vector>
#include <map>
#include <set>
#include <thread>
#include <optional>


namespace kellerberrin {  //  organization level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is a work pipeline, where input objects are pushed, transformed by multiple work threads and then the resultant
// output objects are (optionally) popped off the pipeline. The order of input objects to output objects is guaranteed.
// Input objects must be std::copy_constructable (such as std::shared_ptr<T>), this is a limitation of C++ 20.
// Output objects can be std::move_constructable (such as std::unique_ptr<T>).
// This workflow guarantees that the output objects are removed from the output queue in exactly the same order in which
// the matching input object was presented to the workflow.
//
//
// In many use cases the InputObject and OutputObject will be the same type (and possibly the same object).
// Note that the processing function std::optional(ly) queues output objects for each input object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// The output queue is a QueueTidal (bounded tidal) queue.
template<typename InputObject, typename OutputObject>
// The input object must be std::copy_constructible, this is a limitation of C++20 (that will be removed in C++ 23).
requires std::copy_constructible<InputObject> && std::move_constructible<OutputObject>
class WorkflowFuture
{

  // User supplied function must std::optional(ly) return an output object (possibly the same object).
  // The output object is std::optional, if missing (std::nullopt) then no object is forwarded to the output.
  using WorkProc = std::function<std::optional<OutputObject>(InputObject)>;
  // Convenience alias
  template<typename Output> using SyncOutputQueue = QueueTidal<std::future<std::optional<Output>>>;

public:

  // A bounded (tidal) queue is used to buffer ouput.
  explicit WorkflowFuture( size_t high_tide = TIDAL_QUEUE_DEFAULT_HIGH_TIDE
                         , size_t low_tide = TIDAL_QUEUE_DEFAULT_LOW_TIDE
                         , std::string workflow_name = BOUNDED_QUEUE_DEFAULT_NAME
                         , size_t sample_frequency = BOUNDED_QUEUE_MONITOR_DISABLE)
                         : output_queue_(high_tide, low_tide, workflow_name + "_output_queue", sample_frequency) {}

  ~WorkflowFuture() { workflow_threads_.joinThreads(); }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe (or made so).
  // If the work function is a non-static class member then the first of the ...args should be a
  // pointer (MyClass* this) to the class instance.
  template<typename F, typename... Args>
  bool activateWorkflow(size_t threads, F&& f, Args&&... args) noexcept
  {

    // Clear any active threads.
    workflow_threads_.joinThreads();
    // Clear the queue.
    output_queue_.clear();
    // Set up the work function and requeue the work threads.
    workflow_callback_ = std::bind_front(std::forward<F>(f), std::forward<Args>(args)...);
    // Re-populate the thread pool.
    workflow_threads_.queueThreads(threads);

    return true;

  }


  void push(InputObject input_obj) {

      auto future_output = workflow_threads_.enqueueFuture(workflow_callback_, std::move(input_obj));
      output_queue_.push(std::move(future_output));

  }


  // If the output queue is empty the calling thread will block.
  [[nodiscard]] OutputObject waitAndPop() {

    while (true) {

      auto future_output = output_queue_.waitAndPop();
      auto opt_output = future_output.get();
      if (opt_output.has_value()) {

        auto out_obj = std::move(opt_output.value());
        return out_obj;

      }

    }

  }

  // Queue const access. All const public functions on the queue are thread safe.
  [[nodiscard]] const SyncOutputQueue<OutputObject>& outputQueue() const { return output_queue_; }

private:

  // Thread pool.
  WorkflowThreads workflow_threads_;
  WorkProc workflow_callback_;

   // Output queue, this queue is tidal.
  SyncOutputQueue<OutputObject> output_queue_;

};


}   // end namespace

#endif //KGL_KEL_WORKFLOW_FUTURE_H
