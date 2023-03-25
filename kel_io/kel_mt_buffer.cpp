//
// Created by kellerberrin on 14/02/23.
//

#include "kel_mt_buffer.h"


namespace kel = kellerberrin;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

// Uses the filename extension heuristic documented above to open the underlying StreamIO
// The threads argument is only valid for '.bgz' file types. The argument is ignored for other stream types.
bool kel::StreamMTBuffer::open(const std::string &file_name) {

  auto stream_opt = BaseStreamIO::getStreamIO(file_name);
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

std::optional<std::unique_ptr<kel::BaseStreamIO>> kel::StreamMTBuffer::getStreamIO( const std::string& file_name) {

  auto stream_ptr = std::make_unique<StreamMTBuffer>();
  if (stream_ptr->open(file_name)) {

    return stream_ptr;

  }

  return std::nullopt;

}


std::optional<std::unique_ptr<kel::BaseStreamIO>> kel::StreamMTBuffer::getStreamIO(std::unique_ptr<BaseStreamIO> open_stream_ptr) {

  auto stream_ptr = std::make_unique<StreamMTBuffer>();
  stream_ptr->stream_ptr_ = std::move(open_stream_ptr);
  stream_ptr->line_io_thread_.queueThreads(WORKER_THREAD_COUNT);
  stream_ptr->line_io_thread_.enqueueVoid(&StreamMTBuffer::enqueueIOLineRecord, stream_ptr.get());

  return stream_ptr;

}


void kel::StreamMTBuffer::close() {

  {

    std::lock_guard<std::mutex> lock(mutex_);
    EOF_received_ = true;

  }
  line_io_thread_.joinThreads();
  stream_ptr_ = nullptr;
  line_io_queue_.clear();

}

void kel::StreamMTBuffer::enqueueIOLineRecord() {

  while(true) {

    auto line_record = stream_ptr_->readLine();

    if (line_record.EOFRecord()) {

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
  auto line_record = line_io_queue_.waitAndPop();
  if (line_record.EOFRecord()) {

    EOF_received_ = true;

  }
  return line_record;

}

