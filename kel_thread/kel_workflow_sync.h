//
// Created by kellerberrin on 15/02/23.
//

#ifndef KEL_WORKFLOW_SYNC_H
#define KEL_WORKFLOW_SYNC_H




#include "kel_basic_io.h"
#include "kel_queue_tidal.h"
#include "kel_workflow_threads.h"

#include <string>
#include <string_view>
#include <vector>
#include <optional>
#include <fstream>


namespace kellerberrin {   //  organization::project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
class WorkflowBase {

public:

  WorkflowBase() = default;
  ~WorkflowBase() = default;

private:


};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
class WorkflowSync : public WorkflowBase<T> {

public:

  // The threads argument is only valid for '.bgz' file types. The argument is ignored for other stream types.
  WorkflowSync() = default;
  ~WorkflowSync() = default;

  // Uses the filename extension heuristic documented in kel_basic_io.h to open the underlying file stream.
  template<typename F, typename... Args>
  bool activateWorkflow(size_t threads, F&& f, Args&&... args)
  {

    workflow_callback_ = std::bind_front(std::forward<F>(f), std::forward<Args>(args)...);

    return true;

  }


  [[nodiscard]] T waitAndPOP() override;


  // Access queue stats.
  [[nodiscard]] const QueueTidal<IOLineRecord>& lineQueue() const { return line_io_queue_; }

private:

  // The tidal IO queue parameters.
  static constexpr const size_t IO_HIGH_TIDE_{10000};          // Maximum QueueTidal size
  static constexpr const size_t IO_LOW_TIDE_{2000};            // Low water mark to begin queueing data records
  static constexpr const char* IO_QUEUE_NAME_{"Workflow Sync Queue"};      // The queue name
  static constexpr const size_t IO_SAMPLE_RATE_{100};            // The queue monitor sampling rate (ms), zero (0) disables the monitor.
  // Tidal queue holds buffered IOLineRecord objects.
  QueueTidal<IOLineRecord> line_io_queue_{IO_HIGH_TIDE_, IO_LOW_TIDE_, IO_QUEUE_NAME_, IO_SAMPLE_RATE_};
  // Thread pool asynchronously queues IOLineRecord objects to the queue.
  static constexpr const size_t WORKER_THREAD_COUNT{1};
  std::vector<std::thread> threads_;
  // The StreamIO object.
  std::unique_ptr<BaseStreamIO> stream_ptr_;
  // Set if the stream is in an EOF condition.
  std::atomic<bool> EOF_received_{false};
  // Mutex for readLine()
  std::mutex mutex_;
  // The thread worker function.
  std::function<std::future<T>()> workflow_callback_;

};




} // namespace




#endif //KEL_WORKFLOW_SYNC_H
