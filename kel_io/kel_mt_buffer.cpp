//
// Created by kellerberrin on 14/02/23.
//

#include "kel_mt_buffer.h"


namespace kel = kellerberrin;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

// Uses the filename extension heuristic documented above to open the underlying StreamIO
// The threads argument is only valid for '.bgz' file types. The argument is ignored for other stream types.
bool kel::StreamMTBuffer::open(const std::string &file_name, size_t decompression_threads) {

  auto stream_opt = BaseStreamIO::getStreamIO(file_name, decompression_threads);
  if (not stream_opt) {

    ExecEnv::log().error("StreamMTBuffer::open; could not open stream for file: {}", file_name);
    return false;

  }

  EOF_received_ = false;
  stream_ptr_ = std::move(stream_opt.value());
  line_io_thread_.queueThreads(WORKER_THREAD_COUNT);
  line_io_thread_.enqueueVoid(&StreamMTBuffer::enqueueIOLineRecord, this);

  return true;

}

void kel::StreamMTBuffer::close() {

  EOF_received_ = true;
  line_io_thread_.joinThreads();
  stream_ptr_ = nullptr;
  line_io_queue_.clear();

}

void kel::StreamMTBuffer::enqueueIOLineRecord() {

  while(true) {

    auto line_record = stream_ptr_->readLine();

    if (line_record.EOFRecord()) {

      EOF_received_ = true;
      line_io_queue_.push(std::move(line_record));
      break;

    }

    line_io_queue_.push(std::move(line_record));

  }

}

kel::IOLineRecord kel::StreamMTBuffer::readLine() {

  // Dont block on EOF
  std::lock_guard<std::mutex> lock(mutex_);
  if (EOF_received_) return IOLineRecord::createEOFMarker();
  return line_io_queue_.waitAndPop();

}

