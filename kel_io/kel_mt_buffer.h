//
// Created by kellerberrin on 14/02/23.
//

#ifndef KEL_MT_BUFFER_H
#define KEL_MT_BUFFER_H

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
// A queue buffered multi-thread adapter for text Stream IO.
// A soon as the underlying stream is successfully opened, records are read and stored in a tidal queue and can be
// retrieved using readLine().
// This object cannot be copied.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class StreamMTBuffer {

public:

  StreamMTBuffer() = default;
  ~StreamMTBuffer() = default;

  // Uses the filename extension heuristic documented in kel_basic_io.h to open the underlying file stream.
  // The threads argument is only valid for '.bgz' file types. The argument is ignored for other stream types.
  bool open(const std::string &file_name, size_t decompression_threads = BaseStreamIO::BGZ_DEFAULT_THREADS);
  // After closing the object can be re-opened on another file.
  void close();

  // This is thread safe and buffered.
  // If read by 1 thread then records are guaranteed to presented in file order.
  // If read by multiple threads then there is no line order guarantee.
  // This function will block on an empty queue and no EOF.
  // After an EOF has been received, then subsequent calls will not block and will return EOF objects.
  [[nodiscard]] IOLineRecord readLine();
  // Access queue stats.
  [[nodiscard]] const QueueTidal<IOLineRecord>& lineQueue() const { return line_io_queue_; }

private:

  // The tidal IO queue parameters.
  static constexpr const size_t IO_HIGH_TIDE_{10000};          // Maximum QueueTidal size
  static constexpr const size_t IO_LOW_TIDE_{2000};            // Low water mark to begin queueing data records
  static constexpr const char* IO_QUEUE_NAME_{"StreamMTBuffer Queue"};      // The queue name
  static constexpr const size_t IO_SAMPLE_RATE_{100};            // The queue monitor sampling rate (ms), zero (0) disables the monitor.
  // Tidal queue holds buffered IOLineRecord objects.
  QueueTidal<IOLineRecord> line_io_queue_{IO_HIGH_TIDE_, IO_LOW_TIDE_, IO_QUEUE_NAME_, IO_SAMPLE_RATE_};
  // Thread pool asynchronously queues IOLineRecord objects to the queue.
  static constexpr const size_t WORKER_THREAD_COUNT{1};
  WorkflowThreads line_io_thread_;
  // The StreamIO object.
  std::unique_ptr<BaseStreamIO> stream_ptr_;
  // Set if the stream is in an EOF condition.
  std::atomic<bool> EOF_received_{false};
  // Mutex for readLine()
  std::mutex mutex_;
  // The thread worker function.
  void enqueueIOLineRecord();

};




} // namespace


#endif //KGL_KEL_MT_BUFFER_H
