//
// Created by kellerberrin on 7/02/23.
//

#ifndef KEL_WORKFLOW_FUTURE_H
#define KEL_WORKFLOW_FUTURE_H


#include "kel_queue_mt_safe.h"
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
// A threaded workflow for std::move constructable objects (such as std::unique_ptr<T>).
// This workflow guarantees that the output objects are removed from the output queue in exactly the same order in which
// the matching input object was presented to the workflow.
//
// At the end of processing the WorkflowSynch is STOPPED by a stop_token being pushed onto the workflow queue.
// The stop_token is defined in the WorkflowSynch constructor. This will be a nullptr in the usual case of the
// InputObject being a pointer (std::unique_ptr<T>).
//
// In many use cases the InputObject and OutputObject will be the same type (and possibly the same object).
// Note that the processing function std::optional(ly) queues output objects for each input object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// The input and output objects must be std::move constructable. In addition, the input object should be comparable
// to enable the detection of a stop token placed on the input queue.
// The input queue is a QueueTidal (bounded tidal) queue.
// The output queue is default unbounded but can be constrained to a specified (approximate) size.
template<typename InputObject, typename OutputObject>
requires (std::move_constructible<InputObject> && std::equality_comparable<InputObject>) && std::move_constructible<OutputObject>
class WorkflowFuture
{

  // User supplied function must std::optional(ly) return an output object (possibly the same object).
  // The output object is std::optional, if missing (std::nullopt) then no object is forwarded to the output queue.
  using WorkProc = std::function<std::optional<OutputObject>(InputObject)>;
  // Convenience alias
  template<typename Output> using SyncOutputQueue = QueueTidal<std::future<std::optional<Output>>>;

public:

  // The constructor requires that an input stop token is specified.
  // If the InputObject is a pointer (a typical case is InputObject = std::unique_ptr<T>) then this will be nullptr.
  // A bounded (tidal) queue is used to buffer input.
  // The output queue can be optionally limited in size.
  explicit WorkflowFuture( size_t high_tide = TIDAL_QUEUE_DEFAULT_HIGH_TIDE
                         , size_t low_tide = TIDAL_QUEUE_DEFAULT_LOW_TIDE
                         , std::string workflow_name = BOUNDED_QUEUE_DEFAULT_NAME
                         , size_t sample_frequency = BOUNDED_QUEUE_MONITOR_DISABLE)
                         : output_queue_(high_tide, low_tide, workflow_name + "_output_queue", sample_frequency) {}

  ~WorkflowFuture() { workflow_threads_.joinThreads(); }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe (or made so).
  // If the work function is a non-static class member then the first of the ...args should be a
  // pointer (MyClass* this) to the class instance.
  // If the queue has been STOPPED with a stop token, this function can be called with different workflow functions
  // and thread counts and the workflow will re-activate.
  // Calling this function on a workflow that is not STOPPED will return false and fail.
  template<typename F, typename... Args>
  bool activateWorkflow(size_t threads, F&& f, Args&&... args) noexcept
  {

    // Clear any active threads.
    workflow_threads_.joinThreads();
    // Clear the queues
    output_queue_.clear();
    // Set up the work function and requeue the work threads.
    workflow_callback_ = std::bind_front(f, std::forward<Args>(args)...);
    workflow_threads_.queueThreads(threads);

    return true;

  }


  static std::optional<OutputObject> func(InputObject input_obj) {

    return std::optional<OutputObject>();

  }

/*
  // Returns a std::future holding the work function return value.
  template<typename F, typename... Args>
  [[nodiscard]] auto enqueueFuture(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type>
  {

    using return_type = typename std::result_of<F(Args...)>::type;

    auto callable = std::move(std::bind_front(std::forward<F>(f), std::forward<Args>(args)...));
//    auto d = std::move(callable)();
    auto promise_ptr = std::unique_ptr<std::promise<return_type>>();
    auto return_future = promise_ptr->get_future();
    auto task_lambda = [pr=std::move(promise_ptr), fn=std::move(callable)]() mutable { pr->set_value(std::move(fn)()); };

//    auto task_ptr = std::make_shared<std::packaged_task<return_type()>>(std::move(callable));
//    std::future<return_type> future = task_ptr->get_future();
//    work_queue_.push([task_ptr]()->void{ (*task_ptr)(); });
//    return future;
    return return_future;

  }
*/

  void push(InputObject input_obj) {

      auto future_output = workflow_threads_.enqueueFuture(workflow_callback_, std::move(input_obj));
      output_queue_.push(std::move(future_output));

  }


  // If the output queue is empty the calling thread will block.
  // A suggestion is to push a user defined output stop token onto the output queue when the user workflow function
  // receives an input stop token. The dequeue thread can then stop requesting further output objects.
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

  // Input and output queue const access. All const public functions on these queues are thread safe.
  [[nodiscard]] const SyncOutputQueue<OutputObject>& outputQueue() const { return output_queue_; }

private:

  WorkflowThreads workflow_threads_;
  WorkProc workflow_callback_;

   // Output queue, this queue is tidal.
  SyncOutputQueue<OutputObject> output_queue_;

};


}   // end namespace

#endif //KGL_KEL_WORKFLOW_FUTURE_H
