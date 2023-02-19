//
// Created by kellerberrin on 15/02/23.
//

#ifndef KEL_WORKFLOW_SYNC_H
#define KEL_WORKFLOW_SYNC_H


#include "kel_queue_tidal.h"

#include <future>
#include <vector>
#include <optional>
#include <thread>
#include <fstream>


namespace kellerberrin {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T> class Function; // Not defined.

template<typename ReturnType, typename... Args> class Function<ReturnType(Args...)> {

  struct InvokeBase {

    InvokeBase() = default;
    virtual ~InvokeBase() = default;

    virtual ReturnType invoke(Args&&... args) = 0;

  }; // InvokeBase

  template <typename Executable> struct InvokeImpl : InvokeBase {

    explicit InvokeImpl(Executable executable) : executable_(std::move(executable)) {}
    ~InvokeImpl() override = default;

    [[nodiscard]] ReturnType invoke(Args&&... args) override { return std::invoke(executable_, std::forward<Args>(args)...); }

    Executable executable_;

  }; // InvokeImpl

public:

  template<typename FunctionObject> explicit Function(FunctionObject function_object)
  : executable_ptr_(std::make_unique<InvokeImpl<FunctionObject>>(std::move(function_object))) {}
  Function() = default;
  Function(const Function&) = delete;
  Function(Function&& function) noexcept : executable_ptr_(std::move(function.executable_ptr_)) {}
  ~Function() = default;

  [[nodiscard]] ReturnType operator()(Args&&... args) { return executable_ptr_->invoke(std::forward<Args>(args)...); }
  
private:
  
  std::unique_ptr<InvokeBase> executable_ptr_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename InputObject, typename OutputObject>
requires std::move_constructible<InputObject> && std::move_constructible<OutputObject>
class WorkflowSync {

  using WorkflowFunc = Function<std::optional<OutputObject>(InputObject)>;
  using WorkflowFuncPtr = std::shared_ptr<WorkflowFunc>;
  class QueuedFunctor {

  public:

    QueuedFunctor(WorkflowFuncPtr fn_pointer, InputObject&& input_object) : fn_pointer_(fn_pointer), input_object_(std::move(input_object)) {}
    QueuedFunctor(QueuedFunctor&& queued_functor) = delete;
    QueuedFunctor(const QueuedFunctor& queued_functor) = delete;
    ~QueuedFunctor() = default;

    void operator()() { result_promise_.set_value(fn_pointer_->operator()(std::move(input_object_))); }
    [[nodiscard]] std::future<std::optional<OutputObject>> getFuture() { return result_promise_.get_future(); }

  private:

    WorkflowFuncPtr fn_pointer_;
    InputObject input_object_;
    std::promise<std::optional<OutputObject>> result_promise_;

  };

public:

  WorkflowSync() = default;
  ~WorkflowSync() { joinThreads(); }

  // Note that the variadic args... are presented to ALL active threads and must be thread safe (or made so).
  // If the work function is a non-static class member then the first of the ...args should be a
  // pointer (MyClass* this) to the class instance.
  // The supplied function should be of the form 'std::optional<OutputObject> f(args..., OutputObject)'.
  template<typename F, typename... Args>
  bool activateWorkflow(size_t threads, F&& f, Args&&... args)
  {

    // Clear any active threads.
    joinThreads();

    auto callback_fn = std::bind_front(std::forward<F>(f), std::forward<Args>(args)...);
    fn_pointer_ = std::make_shared<WorkflowFunc>(std::move(callback_fn));

    // Re-populate the thread pool.
    return queueThreads(threads);

  }

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

  void push(InputObject input_object) {

    auto func_ptr = std::make_unique<QueuedFunctor>(fn_pointer_, std::move(input_object));
    auto future = func_ptr->getFuture();
    input_queue_.push(std::move(func_ptr));
    output_queue_.push(std::move(future));

  }


  // Access queue stats.
  [[nodiscard]] const QueueTidal<std::future<OutputObject>>& outputQueue() const { return output_queue_; }
  [[nodiscard]] const QueueTidal<std::unique_ptr<QueuedFunctor>>& inputQueue() const { return input_queue_; }

private:

  // The tidal IO queue parameters.
  static constexpr const size_t HIGH_TIDE_{10000};          // Maximum QueueTidal size
  static constexpr const size_t LOW_TIDE_{2000};            // Low water mark to begin queueing data records
  static constexpr const char* QUEUE_NAME_{"Workflow Sync Queue"};      // The queue name
  static constexpr const size_t SAMPLE_RATE_{100};            // The queue monitor sampling rate (ms), zero (0) disables the monitor.
  // Tidal queue holds buffered output objects.
  QueueTidal<std::future<std::optional<OutputObject>>> output_queue_{HIGH_TIDE_, LOW_TIDE_, std::string(QUEUE_NAME_) + "_Output", SAMPLE_RATE_};
  // Tidal queue holds buffered input objects.
  QueueTidal<std::unique_ptr<QueuedFunctor>> input_queue_{HIGH_TIDE_, LOW_TIDE_, std::string(QUEUE_NAME_) + "_Input", SAMPLE_RATE_};
  // Thread Pool.
  std::vector<std::thread> threads_;
  // The supplied processing function.
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

      threads_.emplace_back(&WorkflowSync::threadProlog, this);

    }

    return true;

  }

};


} // namespace




#endif //KEL_WORKFLOW_SYNC_H
