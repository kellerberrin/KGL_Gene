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


// The workflow is STOPPED if there are no active work threads.
// The workflow is STOPPED on object creation and before activating the workflow with a task function.
// The workflow is in a SHUTDOWN state if a stop token has been received.
// SHUTDOWN is a temporary state that will transition to STOPPED once the input queue has been emptied and all threads are stopped.
// The workflow is ACTIVE if a task function has been supplied and all threads are ready for processing.

enum class FutureWorkflowState { ACTIVE, SHUTDOWN, STOPPED};


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
  template<typename Input> using SyncInputQueue = QueueTidal<std::future<Input>>;

public:

  // The constructor requires that an input stop token is specified.
  // If the InputObject is a pointer (a typical case is InputObject = std::unique_ptr<T>) then this will be nullptr.
  // A bounded (tidal) queue is used to buffer input.
  // The output queue can be optionally limited in size.
  explicit WorkflowFuture(InputObject input_stop_token
      , size_t high_tide = TIDAL_QUEUE_DEFAULT_HIGH_TIDE
      , size_t low_tide = TIDAL_QUEUE_DEFAULT_LOW_TIDE
      , std::string workflow_name = BOUNDED_QUEUE_DEFAULT_NAME
      , size_t sample_frequency = BOUNDED_QUEUE_MONITOR_DISABLE
      , size_t max_output_size = WORKFLOW_OUT_QUEUE_UNBOUNDED_SIZE)
      : input_stop_token_(std::move(input_stop_token))
      , input_queue_(high_tide, low_tide, workflow_name + "_input_queue", sample_frequency)
      , output_queue_(workflow_name + "_output_queue", sample_frequency)
      , max_output_queue_size_(max_output_size) {}

  ~WorkflowFuture() {

    if (workflow_state_ != FutureWorkflowState::STOPPED) {

      // Push a stop token to shut down the active threads so that they can be joined.
      input_queue_.push(std::move(input_stop_token_));

    }

    joinAndDeleteThreads();

  }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe (or made so).
  // If the work function is a non-static class member then the first of the ...args should be a
  // pointer (MyClass* this) to the class instance.
  // If the queue has been STOPPED with a stop token, this function can be called with different workflow functions
  // and thread counts and the workflow will re-activate.
  // Calling this function on a workflow that is not STOPPED will return false and fail.
  template<typename F, typename... Args>
  bool activateWorkflow(size_t threads, F&& f, Args&&... args) noexcept
  {

    { // mutex scope
      std::scoped_lock<std::mutex> lock1(state_mutex_);

      if (workflow_state_ != FutureWorkflowState::STOPPED) {

        return false;

      }

    } // ~mutex scope

    // Clear the queues
    input_queue_.clear();
    output_queue_.clear();

    // Set up the work function and requeue the work threads.
    workflow_callback_ = std::bind_front(f, args...);
    workflow_threads_.queueThreads(threads);

    { // mutex scope
      std::scoped_lock<std::mutex> lock2(state_mutex_);
      workflow_state_ = FutureWorkflowState::ACTIVE;
    } // ~mutex scope
    active_condition_.notify_one();

    return true;

  }


  // This will block if the workflow is not active.
  void push(InputObject input_obj) {

    auto in_obj_opt = pushNoBlock(std::move(input_obj));

    while (in_obj_opt.has_value()) {

      waitUntilActive();
      in_obj_opt = pushNoBlock(std::move(in_obj_opt.value()));

    }

  }


  // The input object is returned to the caller if the workflow is not active.
  // This function can still block (temporarily) if the input queue is at high-tide.
  [[nodiscard]] std::optional<InputObject> pushNoBlock(InputObject input_obj) {

    // Input objects cannot be pushed onto a stopped/shutdown workflow.
    { // mutex scope
      std::scoped_lock<std::mutex> lock(state_mutex_);

      if ( workflow_state_ != FutureWorkflowState::ACTIVE) {

        return std::make_optional<InputObject>(std::move(input_obj));

      }

      if (input_obj == input_stop_token_) {

        workflow_state_ = FutureWorkflowState::SHUTDOWN;

      }

    } // ~mutex.

    input_queue_.push(std::move(input_obj));

    // The input object has been consumed.
    return std::nullopt;

  }

  // If the output queue is empty the calling thread will block.
  // A suggestion is to push a user defined output stop token onto the output queue when the user workflow function
  // receives an input stop token. The dequeue thread can then stop requesting further output objects.
  [[nodiscard]] OutputObject waitAndPop() {

    if (max_output_queue_size_.load() != WORKFLOW_OUT_QUEUE_UNBOUNDED_SIZE) {

      output_condition_.notify_one();

    }
    return output_queue_.waitAndPop();

  }

  // This function is not thread safe! A race condition potentially exists.
  // If the producer thread has pushed a stop token, another thread calling this function may still return ACTIVE.
  // Consider using waitOnStopped() or waitOnActive() instead.
  [[nodiscard]] FutureWorkflowState workflowState() const { return workflow_state_.load(); }

  // Calling thread(s) wait on a std::condition_variable until the workflow queue is STOPPED.
  void waitUntilStopped() const {

    {
      std::unique_lock<std::mutex> lock(state_mutex_);
      stopped_condition_.wait(lock, [this]{ return workflow_state_ == FutureWorkflowState::STOPPED; });
    }
    stopped_condition_.notify_one(); // In case multiple threads are blocked.

  }

  // Calling thread(s) wait on a std::condition_variable until the workflow queue is ACTIVE.
  void waitUntilActive() const {

    {
      std::unique_lock<std::mutex> lock(state_mutex_);
      active_condition_.wait(lock, [this]{ return workflow_state_ == FutureWorkflowState::ACTIVE; });
    }
    active_condition_.notify_one(); // In case multiple threads are blocked.

  }


  // Input and output queue const access. All const public functions on these queues are thread safe.
  [[nodiscard]] const SyncInputQueue<InputObject>& inputQueue() const { return input_queue_; }
  [[nodiscard]] const QueueMtSafe<OutputObject>& outputQueue() const { return output_queue_; }

  // A value of zero means the output queue can grow without bound.
  constexpr static const size_t WORKFLOW_OUT_QUEUE_UNBOUNDED_SIZE{0};

private:

  InputObject input_stop_token_;

  WorkflowThreads workflow_threads_;
  WorkProc workflow_callback_;
  std::mutex process_mutex_;
  std::condition_variable output_condition_;  // Blocks when then output queue (optionally) exceeds a specified size.

  // Input queue. This is a bounded (tidal) queue to control queue size and balance CPU load.
  SyncInputQueue<InputObject> input_queue_;
    // Output queue, note that the default is that this queue can grow without bound if not explicitly constrained.
  QueueMtSafe<OutputObject> output_queue_;
  // Optional limit on the approximate maximum size of the output queue.
  // Placing an upper limit on the output queue size may result in the input queue blocking.
  // Defaults to an unbounded (0 size) output queue size.
  const std::atomic<size_t> max_output_queue_size_{WORKFLOW_OUT_QUEUE_UNBOUNDED_SIZE};

  // Threads can wait on workflow queue state.
  std::atomic<FutureWorkflowState> workflow_state_{FutureWorkflowState::STOPPED };
  mutable std::mutex state_mutex_;
  mutable std::condition_variable stopped_condition_;
  mutable std::condition_variable active_condition_;

  // Workflow threads do not block if they process input objects out of order.
  // The processed results (output objects) are placed on a std::priority queue to facilitate re-ordering so that
  // the ordering of output objects matches that of input objects without blocking.
  // Shuts down the all the active threads if a stop token is found on the input queue.
  void threadProlog() {

    // Loop until a stop token is encountered and then recursively shutdown the active threads.
    while(true) {

      // We control the approximate max output queue size by blocking reading the input queue.
      // The output queue cannot be blocked directly as it is accessed in a critical code section.
      if (max_output_queue_size_.load() != WORKFLOW_OUT_QUEUE_UNBOUNDED_SIZE) {

        {
          std::unique_lock lock(process_mutex_);
          output_condition_.wait(lock, [this]{ return output_queue_.size() < max_output_queue_size_.load(); });
        }

      }

      // Get the next input object. If it is a stop token then recursively shutdown all the work threads.
      // If not a stop token then call the user supplied workflow function and manage the input/output synchronization logic.
      InputObject work_item = input_queue_.waitAndPop();

      // If stop token then shutdown the thread.
      if (work_item.second == input_stop_token_) {

        // Process the stop token and optionally push it onto the output queue.
        std::optional<OutputObject> out_opt = workflow_callback_(std::move(work_item.second));
        if (out_opt.has_value()) {

          output_queue_.push(std::move(out_opt.value()));

        }

         // Workflow is now stopped
        {
            std::scoped_lock<std::mutex> lock(state_mutex_);
            workflow_state_ = FutureWorkflowState::STOPPED;
        }
       stopped_condition_.notify_one();

        return;  // Thread exits.

      }

      // Call the user supplied workflow function. Note that CPU work is done outside the critical code section.
      auto out_opt = workflow_callback_(std::move(work_item.second));

      // Mutex protected critical code.
      {
        std::scoped_lock lock(process_mutex_);


      } // Mutex protected critical code ends.

    } // ~thread while loop

  }

  void joinAndDeleteThreads() noexcept {

    // Join all the threads
    for(auto& thread : threads_) {

      thread.join();

    }

    threads_.clear();

  }

};


}   // end namespace

#endif //KGL_KEL_WORKFLOW_FUTURE_H
