//
// Created by kellerberrin on 30/12/22.
//

#ifndef KEL_WORKFLOW_SYNCH_H
#define KEL_WORKFLOW_SYNCH_H


#include <functional>
#include <vector>
#include <map>
#include <set>
#include <thread>

#include "kel_mt_queue.h"
#include "kel_bound_queue.h"


namespace kellerberrin {  //  organization level namespace



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A threaded workflow for std::move constructable objects (std::unique_ptr).
// These queues guarantee that the output objects are removed from the output queue in exactly the same order in which
// the matching input object was presented to the input workflow queue.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class SynchWorkflowState { ACTIVE, STOPPED};

// Some very long lived (heat death of the universe) applications may need to use a 128 bit unsigned int
// to assign an ordering to the input and output objects.
// If the 128 bit unsigned type is not supported on the target architecture then use uint64_t.

using WorkFlowObjectCounter = __uint128_t;
//using WorkFlowObjectCounter = uint64_t;

template<typename InputObject, typename OutputObject, template <typename> typename InputQueue>
requires (std::move_constructible<InputObject> && std::equality_comparable<InputObject>)
         && (std::move_constructible<OutputObject> && std::equality_comparable<OutputObject>)
class WorkflowSynchQueue
{

  struct CompareProcessed {
    bool operator()(const std::pair<WorkFlowObjectCounter, OutputObject> &lhs, const std::pair<WorkFlowObjectCounter, OutputObject> &rhs) const {
      return lhs.first > rhs.first;
    }
  };

  using WorkProc = std::function<OutputObject(InputObject)>;

public:

  WorkflowSynchQueue(InputObject input_stop_token
      , std::unique_ptr<InputQueue<InputObject>> input_queue_ptr = std::make_unique<InputQueue<InputObject>>())
      : input_stop_token_(std::move(input_stop_token))
      , input_queue_ptr_(std::move(input_queue_ptr)) {}

  ~WorkflowSynchQueue() {

    stopProcessing();
    // Explicitly remove the workflow lambda to prevent circular references from any captured pointer arguments.
    workflow_callback_ = nullptr;

  }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe (or made so).
  // If the work function is a non-static class member then the first of the ...args should be a
  // pointer (MyClass* this) to the class instance.
  template<typename F, typename... Args>
  void activateWorkflow(size_t threads, F&& f, Args&&... args) noexcept
  {

    workflow_callback_ = [f, args...](InputObject t)->OutputObject{ return std::invoke(f, args..., std::move(t)); };
    queueThreads(threads);
    workflow_state_ = SynchWorkflowState::ACTIVE;

  }


  // Although this function is thread-safe, it should really only be called with a single thread.
  // Since calling with multiple threads will lead to undefined input object ordering.
  void push(InputObject input_obj) {

    WorkFlowObjectCounter input_tag{0};
    // RAII mutex protected critical code.
    {
      std::scoped_lock lock(process_mutex_);

      if (input_obj != input_stop_token_) {

        ++object_counter_;
        input_tag = object_counter_;
        ordered_requests_.push(object_counter_);

      }

    } // ~mutex

    // The order of objects pushed onto the input queue does not matter as the input objects have already been tagged.
    // In particular, this queue can block and be implemented as a tidal (automatic load balancing) queue.
    input_queue_ptr_->push({input_tag, std::move(input_obj)});

  }

  // Although this function is thread-safe, it should really only be called with a single thread.
  // Since calling with multiple threads will lead to undefined output object ordering.
  [[nodiscard]] OutputObject waitAndPop() {

    return output_queue_.waitAndPop();

  }

  // Workflow is STOPPED if there are no active work threads.
  // The Workflow queue is STOPPED on object creation and before activating the workflow function and creating active threads.
  // The Workflow queue is also STOPPED after a stop token is placed on the input object queue and there are no longer any active threads.
  [[nodiscard]] SynchWorkflowState workflowState() const { return workflow_state_; }

  // Calling thread(s) wait on a std::condition_variable until the queue is STOPPED.
  void waitUntilStopped() const {

    {
      std::unique_lock<std::mutex> lock(state_mutex_);
      stopped_condition_.wait(lock, [this]{ return workflow_state_ == SynchWorkflowState::STOPPED; });
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

  std::unique_ptr<InputQueue<InputObject>> input_queue_ptr_;
  // Custom comparators ensure lower counter values are at the top of the priority queues.
  std::priority_queue< WorkFlowObjectCounter
      , std::vector<WorkFlowObjectCounter>
      , std::greater<>> ordered_requests_;
  std::priority_queue< std::pair<WorkFlowObjectCounter, OutputObject>
      , std::vector<std::pair<WorkFlowObjectCounter, OutputObject>>
      , CompareProcessed > processed_objects_;
  MtQueue<OutputObject> output_queue_;

  // Threads can wait on workflow queue state.
  SynchWorkflowState workflow_state_{ SynchWorkflowState::STOPPED };
  mutable std::mutex state_mutex_;
  mutable std::condition_variable stopped_condition_;

  void queueThreads(size_t threads)
  {

    // Remove any existing threads.
    stopProcessing();

    // Always have at least one worker thread queued.
    threads = threads < 1 ? 1 : threads;

    // Queue the worker threads,
    for(size_t i = 0; i < threads; ++i) {

      threads_.emplace_back(&WorkflowSynchQueue::threadProlog, this);
      ++active_threads_;

    }

  }

  void threadProlog() {

    // Loop until a stop token is encountered.
    while(true) {

      // Can block.
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
            workflow_state_ = SynchWorkflowState::STOPPED;
          }
          stopped_condition_.notify_one();

        }
        break; // Thread terminates,

      } else {

        // Call the worker function.
        std::pair<WorkFlowObjectCounter, OutputObject> output(work_item.first, std::move(workflow_callback_(std::move(work_item.second))));

        // RAII mutex protected critical code.
        {
          std::scoped_lock lock(process_mutex_);

          // Is the processed object the next ordered (earliest) request?

          if (output.first == ordered_requests_.top()) {

            ordered_requests_.pop();
            // Output queue cannot block since this code is mutex protected.
            output_queue_.push(std::move(output.second));

            // Unqueue any processed requests that match the request priority queue.
            while(not ordered_requests_.empty() and not processed_objects_.empty()) {

              if (ordered_requests_.top() == processed_objects_.top().first) {

                ordered_requests_.pop();
                // Very nasty. Should by able to std::move a std::unique_ptr from a std::priority_queue without resorting to this kind of unpleasantness.
                std::pair<WorkFlowObjectCounter, OutputObject> processed = std::move(const_cast<std::pair<WorkFlowObjectCounter, OutputObject> &>(processed_objects_.top()));
                processed_objects_.pop();
                // Output queue cannot block since this code is mutex protected.
                output_queue_.push(std::move(processed.second));

              } else {

                break;

              }

            }

          } else {

            processed_objects_.emplace(output.first, std::move(output.second));

          }

        } // Mutex protected critical code ends.

      }

    }

  }

  void stopProcessing() {

    // If any active threads then push the stop token onto the input queue.
    if (active_threads_ != 0) {

      input_queue_ptr_->push({0, std::move(input_stop_token_)});

    }

    // Join all the threads
    for(auto& thread : threads_) {

      thread.join();

    }

    threads_.clear();

  }

};


}   // end namespace


#endif //KEL_WORKFLOW_SYNCH_H
