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


enum class SyncWorkflowState { ACTIVE, STOPPED};

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

  ~WorkflowSyncQueue() noexcept {

    // Shutdown any active threads and join() them.
    // If we have active threads then push an input stop token so that the workflow is STOPPED.
    if (active_threads_ != 0) {

      input_queue_ptr_->push({0, std::move(input_stop_token_)});

    }
    joinAndDeleteThreads();
    // Explicitly remove the workflow lambda to prevent circular references from any captured pointer arguments.
    workflow_callback_ = nullptr;

  }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe (or made so).
  // If the work function is a non-static class member then the first of the ...args should be a
  // pointer (MyClass* this) to the class instance.
  // This function is not multi-threaded and must be called after workflow queue creation for workflow processing to be ACTIVE.
  // If the queue has been STOPPED, this function can be called with different workflow functions and thread counts.
  // Calling this function on an ACTIVE workflow queue will return false.
  template<typename F, typename... Args>
  bool activateWorkflow(size_t threads, F&& f, Args&&... args)
  {

    if (workflow_state_ == SyncWorkflowState::ACTIVE) {

      return false;

    }
    joinAndDeleteThreads();
    workflow_callback_ = [f, args...](InputObject t)->OutputObject{ return std::invoke(f, args..., std::move(t)); };
    queueThreads(threads);
    workflow_state_ = SyncWorkflowState::ACTIVE;
    return true;

  }


  // Although this function is thread-safe, it should really only be called with a single thread.
  // Since calling with multiple threads will lead to undefined input object ordering.
  // This function can be called on a STOPPED workflow queue.
  // However, if the underlying input queue is a BoundedMtQueue, the function will block once high-tide is reached.
  // The function will not block if the underlying input queue is an MtQueue.
  void push(InputObject input_obj) {

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
    // In particular, this queue can block and be implemented as a tidal (automatic load balancing) queue.
    input_queue_ptr_->push({input_tag, std::move(input_obj)});

  }

  // Although this function is thread-safe, it should really only be called with a single thread.
  // Since calling with multiple threads will lead to undefined output object ordering.
  // This function can be called on a STOPPED workflow queue.
  // If the output queue is empty the calling thread will block.
  // A suggestion is to push a user defined output stop token onto the output queue when the user workflow function
  // receives an input stop token. The dequeue thread can then stop requesting further output objects.
  [[nodiscard]] OutputObject waitAndPop() {

    return output_queue_.waitAndPop();

  }

  // Workflow is STOPPED if there are no active work threads.
  // The Workflow queue is STOPPED on object creation and before activating the workflow function and creating active threads.
  // The Workflow queue is also STOPPED after a stop token is placed on the input object queue and there are no longer any active threads.
  [[nodiscard]] SyncWorkflowState workflowState() const { return workflow_state_; }

  // Calling thread(s) wait on a std::condition_variable until the workflow queue is STOPPED.
  void waitUntilStopped() const {

    {
      std::unique_lock<std::mutex> lock(state_mutex_);
      stopped_condition_.wait(lock, [this]{ return workflow_state_ == SyncWorkflowState::STOPPED; });
    }
    stopped_condition_.notify_one(); // In case multiple threads are blocked.

  }

  // Input and output queue const access. All const public functions on these queues are thread safe.
  [[nodiscard]] const InputQueue<InputObject>& inputQueue() const { return *input_queue_ptr_; }
  [[nodiscard]] const MtQueue<OutputObject>& outputQueue() const { return output_queue_; }

private:

  InputObject input_stop_token_;

  std::atomic<size_t> active_threads_{0};
  std::vector<std::thread> threads_;

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

  // Can only be called on a STOPPED workflow.
  void queueThreads(size_t threads)
  {

    // Always have at least one worker thread.
    threads = threads < 1 ? 1 : threads;

    // Queue the worker threads,
    for(size_t i = 0; i < threads; ++i) {

      threads_.emplace_back(&WorkflowSyncQueue::threadProlog, this);
      ++active_threads_;

    }

  }

  // Workflow threads do not block if they process input objects out of order.
  // The processed results (output objects) are placed on a std::priority queue to facilitate re-ordering so that
  // the ordering of output objects matches that of input objects without blocking.
  // Shuts down the all the active threads if a stop token is found on the input queue.
  void threadProlog() {

    // Loop until a stop token is encountered.
    while(true) {

      // Get the next input object. If it is a stop token then recursively shutdown all the work threads.
      // If not a stop token then call the user supplied workflow function and manage the input/output synchronization logic.
      std::pair<WorkFlowObjectCounter, InputObject> work_item = input_queue_ptr_->waitAndPop();

      if (work_item.second == input_stop_token_) {

        if (--active_threads_ != 0) {

          // If multiple threads are active then re-queue the input stop token.
          input_queue_ptr_->push(std::move(work_item));

        } else {

          // Last active thread.
          // Call the workflow function with the stop token
          // to notify the processing function logic that the workflow queue will be STOPPED.
          // Typically, this will push an output stop token onto the output queue.
          output_queue_.push(workflow_callback_(std::move(work_item.second)));
          // Explicitly remove the workflow lambda to prevent circular references from any captured pointer arguments.
          workflow_callback_ = nullptr;
          // Notify any threads waiting on the workflow STOPPED condition.
          {
            std::lock_guard<std::mutex> lock(state_mutex_);
            workflow_state_ = SyncWorkflowState::STOPPED;
          }
          stopped_condition_.notify_one();

        }
        break; // Thread terminates,

      } else {

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

            }

          } else { // Place the processed object onto the internal re-ordering std::priority_queue.

            processed_objects_.emplace(output.first, std::move(output.second));

          }

        } // Mutex protected critical code ends.

      }

    }

  }

  void joinAndDeleteThreads() {

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
