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

#ifndef KEL_WORKFLOW_SYNC_H
#define KEL_WORKFLOW_SYNC_H

#include "kel_mt_queue.h"
#include "kel_bound_queue.h"

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
// These queues guarantee that the output objects are removed from the output queue in exactly the same order in which
// the matching input object was presented to the input workflow queue.
//
// At the end of processing the WorkflowSynchQueue is STOPPED by a stop_token being pushed onto the workflow queue.
// The stop_token is defined in the WorkflowSynchQueue constructor. This will be a nullptr in the usual case of the
// InputObject being a pointer (std::unique_ptr<T>).
//
// In many use cases the InputObject and OutputObject will be the same type (and the same object).
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// The workflow is STOPPED if there are no active work threads.
// The workflow is STOPPED on object creation and before activating the workflow with a task function.
// The workflow is in a SHUTDOWN state if a stop token has been received.
// SHUTDOWN is a temporary state that will transition to STOPPED.
// The workflow is ACTIVE if a task function has been supplied and all threads are ready for processing.

enum class SyncWorkflowState { ACTIVE, SHUTDOWN, STOPPED};

// Some very long-lived (heat death of the universe) applications may need to use a 128 bit unsigned int
// to assign an ordering to the input and output objects.
// If the 128 bit unsigned type is not supported on the target architecture then use uint64_t instead.
using WorkFlowObjectCounter = __uint128_t;
//using WorkFlowObjectCounter = uint64_t;

// The input and output objects must be std::move constructable. In addition, the input object should be comparable
// to enable the detection of a stop token placed on the input queue.
// The input queue can be specified as MtQueue (unbounded) or BoundedMtQueue (bounded tidal).
template<typename InputObject, typename OutputObject, template <typename> typename InputQueue>
requires (std::move_constructible<InputObject> && std::equality_comparable<InputObject>)
         && std::move_constructible<OutputObject>
class WorkflowSyncQueue
{

  // A custom ordering used by std::priority_queue.
  struct CompareProcessed {
    bool operator()(const std::pair<WorkFlowObjectCounter, OutputObject> &lhs, const std::pair<WorkFlowObjectCounter, OutputObject> &rhs) const {
      return lhs.first > rhs.first;
    }
  };

  // User supplied function must return an output object (possibly the same object).
  using WorkProc = std::function<OutputObject(InputObject)>;

public:

  // The constructor requires that an input stop token is specified.
  // If the InputObject is a pointer (a typical case is InputObject = std::unique_ptr<T>) then this will be nullptr.
  // The input queue will be either a MtQueue (unbounded) or BoundedMtQueue (a bounded tidal queue).
  explicit WorkflowSyncQueue(InputObject input_stop_token
      , std::unique_ptr<InputQueue<InputObject>> input_queue_ptr = std::make_unique<InputQueue<InputObject>>())
      : input_stop_token_(std::move(input_stop_token))
      , input_queue_ptr_(std::move(input_queue_ptr)) {}

  ~WorkflowSyncQueue() {

    if (workflow_state_ != SyncWorkflowState::STOPPED) {

      // This stop token is not processed by the workflow function.
      input_queue_ptr_->push({0, std::move(input_stop_token_)});

    }

    joinAndDeleteThreads();

  }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe (or made so).
  // If the work function is a non-static class member then the first of the ...args should be a
  // pointer (MyClass* this) to the class instance.
  // If the queue has been STOPPED, this function can be called with different workflow functions and thread counts.
  // Calling this function on a workflow that is not STOPPED will return false and fail.
  template<typename F, typename... Args>
  bool activateWorkflow(size_t threads, F&& f, Args&&... args) noexcept
  {

    { // mutex scope
      std::scoped_lock<std::mutex> lock1(state_mutex_);

      if (workflow_state_ != SyncWorkflowState::STOPPED) {

        return false;

      }

    } // ~mutex scope

    joinAndDeleteThreads();
    workflow_callback_ = std::bind_front(f, args...);
    queueThreads(threads);

    { // mutex scope
      std::scoped_lock<std::mutex> lock2(state_mutex_);
      workflow_state_ = SyncWorkflowState::ACTIVE;
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


  // Does not block, the input object is returned to the caller if the workflow is not active.
  [[nodiscard]] std::optional<InputObject> pushNoBlock(InputObject input_obj) {

    // Input objects cannot be pushed onto a stopped/shutdown workflow.
    { // mutex scope
      std::scoped_lock<std::mutex> lock(state_mutex_);

      if ( workflow_state_ != SyncWorkflowState::ACTIVE) {

        return std::make_optional<InputObject>(std::move(input_obj));

      }

      if (input_obj == input_stop_token_) {

        workflow_state_ = SyncWorkflowState::SHUTDOWN;

      }

    } // ~mutex.

    // Note, this variable is only used to carry the object ordering tag outside the mutex scope.
    WorkFlowObjectCounter input_tag;
    // Mutex protected critical code.
    {
      std::scoped_lock lock(process_mutex_);

      if (input_obj != input_stop_token_) {

        ++object_counter_;
        input_tag = object_counter_;
        ordered_requests_.push(object_counter_);

      }

    } // End of critical code.

    // The order of objects pushed onto the input queue does not matter as the input objects have already been tagged.
    input_queue_ptr_->push({input_tag, std::move(input_obj)});

    // The input object has been consumed.
    return std::nullopt;

  }

  // If the output queue is empty the calling thread will block.
  // A suggestion is to push a user defined output stop token onto the output queue when the user workflow function
  // receives an input stop token. The dequeue thread can then stop requesting further output objects.
  [[nodiscard]] OutputObject waitAndPop() {

    return output_queue_.waitAndPop();

  }

  // This function is not thread safe! A race condition potentially exists.
  // If a producer thread has pushed a stop token, another thread calling this function may still return ACTIVE.
  // Consider using waitOnStopped() or waitOnActive() instead.
   [[nodiscard]] SyncWorkflowState workflowState() const { return workflow_state_.load(); }

  // Calling thread(s) wait on a std::condition_variable until the workflow queue is STOPPED.
  void waitUntilStopped() const {

    {
      std::unique_lock<std::mutex> lock(state_mutex_);
      stopped_condition_.wait(lock, [this]{ return workflow_state_ == SyncWorkflowState::STOPPED; });
    }
    stopped_condition_.notify_one(); // In case multiple threads are blocked.

  }

  // Calling thread(s) wait on a std::condition_variable until the workflow queue is ACTIVE.
  void waitUntilActive() const {

    {
      std::unique_lock<std::mutex> lock(state_mutex_);
      active_condition_.wait(lock, [this]{ return workflow_state_ == SyncWorkflowState::ACTIVE; });
    }
    active_condition_.notify_one(); // In case multiple threads are blocked.

  }


  // Input and output queue const access. All const public functions on these queues are thread safe.
  [[nodiscard]] const InputQueue<InputObject>& inputQueue() const { return *input_queue_ptr_; }
  [[nodiscard]] const MtQueue<OutputObject>& outputQueue() const { return output_queue_; }

private:

  InputObject input_stop_token_;

  std::atomic<size_t> active_threads_{0};
  std::vector<std::jthread> threads_;

  WorkFlowObjectCounter object_counter_{0};
  WorkProc workflow_callback_;
  std::mutex process_mutex_;

  // Input queue. This can be a bounded (tidal) queue to control queue size and balance CPU load.
  std::unique_ptr<InputQueue<InputObject>> input_queue_ptr_;
  // Custom comparators ensure lower counter values are at the top of the priority queues.
  std::priority_queue< WorkFlowObjectCounter
      , std::vector<WorkFlowObjectCounter>
      , std::greater<>> ordered_requests_;
  std::priority_queue< std::pair<WorkFlowObjectCounter, OutputObject>
      , std::vector<std::pair<WorkFlowObjectCounter, OutputObject>>
      , CompareProcessed > processed_objects_;
  // Output queue, note that this queue can grow without bound.
  MtQueue<OutputObject> output_queue_;

  // Threads can wait on workflow queue state.
  std::atomic<SyncWorkflowState> workflow_state_{SyncWorkflowState::STOPPED };
  mutable std::mutex state_mutex_;
  mutable std::condition_variable stopped_condition_;
  mutable std::condition_variable active_condition_;

  // Can only be called on a STOPPED workflow.
  void queueThreads(size_t threads) noexcept
  {

    // Always have at least one worker thread.
    threads = threads < 1 ? 1 : threads;

    // Queue the worker threads,
    for(size_t i = 0; i < threads; ++i) {

      threads_.emplace_back(&WorkflowSyncQueue::threadProlog,  this);

    }
    active_threads_.store(threads);

  }

  // Workflow threads do not block if they process input objects out of order.
  // The processed results (output objects) are placed on a std::priority queue to facilitate re-ordering so that
  // the ordering of output objects matches that of input objects without blocking.
  // Shuts down the all the active threads if a stop token is found on the input queue.
  void threadProlog() {

  // Loop until a stop token is encountered and then recursively shutdown the active threads.
    while(true) {

      // Get the next input object. If it is a stop token then recursively shutdown all the work threads.
      // If not a stop token then call the user supplied workflow function and manage the input/output synchronization logic.
      std::pair<WorkFlowObjectCounter, InputObject> work_item = input_queue_ptr_->waitAndPop();

      // If stop token then shutdown the thread.
      if (work_item.second == input_stop_token_) {

        // Check if the last thread active.
        if (active_threads_.fetch_sub(1) != 1) {

          input_queue_ptr_->push(std::move(work_item)); // Re-queue the stop token.

        } else {

          output_queue_.push(std::move(workflow_callback_(std::move(work_item.second))));

          // Workflow is now stopped
          {
            std::scoped_lock<std::mutex> lock(state_mutex_);
            workflow_state_ = SyncWorkflowState::STOPPED;
          }
          stopped_condition_.notify_one();
        }

        return;  // Thread exits.

      }

      // Call the user supplied workflow function. Note that CPU work is done outside the critical code section.
      std::pair<WorkFlowObjectCounter, OutputObject> output(work_item.first, std::move(workflow_callback_(std::move(work_item.second))));

      // Mutex protected critical code.
      {
        std::scoped_lock lock(process_mutex_);

        // Is the processed object the next ordered (earliest) request?
        if (output.first == ordered_requests_.top()) {

          ordered_requests_.pop();
          output_queue_.push(std::move(output.second));

          // Dequeue any processed requests that match the request priority queue.
          while(not ordered_requests_.empty() and not processed_objects_.empty()) {

            if (ordered_requests_.top() == processed_objects_.top().first) {

              ordered_requests_.pop();
              // Should be able to std::move a std::unique_ptr from a std::priority_queue without resorting to this unpleasantness.
              std::pair<WorkFlowObjectCounter, OutputObject> processed = std::move(const_cast<std::pair<WorkFlowObjectCounter, OutputObject> &>(processed_objects_.top()));
              processed_objects_.pop();
              output_queue_.push(std::move(processed.second));

            } else {

              break;

            }

          } // ~while

        } else { // Place the processed object onto the internal re-ordering std::priority_queue.

          processed_objects_.emplace(output.first, std::move(output.second));

        }

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

// Convenience synchronous workflow typedefs.

// Implemented with a bounded tidal input queue.
template<typename WorkObj> using BoundedSyncInput = BoundedMtQueue<std::pair<WorkFlowObjectCounter, WorkObj>>;
template<typename WorkInput, typename WorkOutput> using WorkflowSyncBounded = WorkflowSyncQueue<WorkInput, WorkOutput, BoundedSyncInput>;

// Implemented with an unbounded input queue.
template<typename WorkObj> using SyncInputQueue = MtQueue<std::pair<WorkFlowObjectCounter, WorkObj>>;
template<typename WorkInput, typename WorkOutput> using WorkflowSync = WorkflowSyncQueue<WorkInput, WorkOutput, SyncInputQueue>;


}   // end namespace


#endif //KEL_WORKFLOW_SYNC_H
