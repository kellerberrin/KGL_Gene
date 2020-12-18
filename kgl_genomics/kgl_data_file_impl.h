//
// Created by kellerberrin on 19/4/20.
//

#ifndef KGL_DATA_FILE_IMPL_H
#define KGL_DATA_FILE_IMPL_H



#include "kel_exec_env.h"
#include "kel_mt_queue.h"
#include "kel_basic_io.h"
#include "kel_thread_pool.h"
#include "kgl_genome_types.h"

#include <string>
#include <vector>
#include <fstream>
#include <thread>

namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IO object (multi-threaded) - enqueues file line reads for further processing

class FileDataIO  {

public:

  explicit FileDataIO(std::string read_file_name, size_t high_tide = HIGH_TIDE_, size_t low_tide = LOW_TIDE_)
    : read_file_name_(std::move(read_file_name)), raw_io_queue_(high_tide, low_tide), io_threads_(IO_THREAD_COUNT) {}
  virtual ~FileDataIO() = default;

  // Begin reading records, spawns thread(s).
  // Specifies how many consumer threads. This is required to enqueue eof tokens.
  bool commenceIO(size_t consumer_threads);

  [[nodiscard]] const std::string& fileName() const { return read_file_name_; }

  // Called by each reader thread
  [[nodiscard]] IOLineRecord readIORecord() { return raw_io_queue_.waitAndPop(); }

  constexpr static const size_t IO_THREAD_COUNT = 1;

private:

  std::string read_file_name_;
  size_t consumer_threads_{0};
  BoundedMtQueue<IOLineRecord> raw_io_queue_; // The raw tidal IO queue
  ThreadPool io_threads_;
  std::optional<std::unique_ptr<BaseStreamIO>> file_stream_opt_;
  bool eof_flag_{false};
  std::mutex eof_mutex_;

  static constexpr const long HIGH_TIDE_{100000};          // Maximum BoundedMtQueue size
  static constexpr const long LOW_TIDE_{10000};            // Low water mark to begin queueing data records

  void rawDataIO(); // The read/decompress from disk and place on the tidal queue.
  void enqueueEOF(); // Shutdown the consumer threads.

};



} // end namespace

#endif //KGL_VARIANT_FILE_IMPL_H
