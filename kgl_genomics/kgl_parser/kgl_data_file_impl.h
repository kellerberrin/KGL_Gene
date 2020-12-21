//
// Created by kellerberrin on 19/4/20.
//

#ifndef KGL_DATA_FILE_IMPL_H
#define KGL_DATA_FILE_IMPL_H



#include "kel_exec_env.h"
#include "kel_bound_queue.h"
#include "kel_basic_io.h"
#include "kgl_genome_types.h"

#include <string>
#include <vector>
#include <fstream>
#include <thread>

namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// IO object (multi-threaded) - enqueues file line reads for further processing
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class FileDataIO  {

public:

  FileDataIO() = default;
  ~FileDataIO() noexcept;

  // Begin reading records, spawns thread(s).
  bool commenceIO(std::string read_file_name);

  [[nodiscard]] const std::string& fileName() const { return read_file_name_; }

  // Called by each reader thread
  [[nodiscard]] IOLineRecord readIORecord() { return raw_io_queue_.waitAndPop(); }

  // Push an eof marker onto the queue
  void enqueueEOF() { raw_io_queue_.push(QUEUED_EOF_MARKER); }

private:

  std::string read_file_name_{"IO Stream Inactive"};

  // The raw tidal IO queue parameters
  static constexpr const size_t IO_HIGH_TIDE_{2000};          // Maximum BoundedMtQueue size
  static constexpr const size_t IO_LOW_TIDE_{1000};            // Low water mark to begin queueing data records
  static constexpr const char* IO_QUEUE_NAME_{"IO Queue"};      // The queue name
  static constexpr const size_t IO_SAMPLE_RATE_{500};            // The queue stats sampling rate.
  BoundedMtQueue<IOLineRecord> raw_io_queue_{ IO_HIGH_TIDE_, IO_LOW_TIDE_, IO_QUEUE_NAME_, IO_SAMPLE_RATE_};

  // VCF queue worker threads
  // Unless there is a good reason, the number of worker threads should be 1.
  constexpr static const size_t IO_THREAD_COUNT_{1};
  // The detached main thread.
  ThreadPool detached_launch_{1};
  // Synchronize shutdown
  std::future<void> launch_token_;

  ThreadPool io_threads_{IO_THREAD_COUNT_};

  // The data input stream.
  std::optional<std::unique_ptr<BaseStreamIO>> file_stream_opt_;


  void launchThreads();
  void rawDataIO(); // Read/decompress from disk and place on the tidal queue.

};



} // end namespace

#endif //KGL_VARIANT_FILE_IMPL_H
