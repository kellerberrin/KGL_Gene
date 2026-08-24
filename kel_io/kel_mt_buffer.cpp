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

    ExecEnv::log().error("Could not open stream for file: {}", file_name);
    return false;

  }

  return openStream(std::move(stream_opt.value()));

}

bool kel::StreamMTBuffer::open(const std::string &file_name, size_t decompression_threads) {

  auto stream_opt = BaseStreamIO::getStreamIO(file_name, decompression_threads);
  if (not stream_opt) {

    ExecEnv::log().error("Could not open stream for file: {}", file_name);
    return false;

  }

  return openStream(std::move(stream_opt.value()));

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
  if (not stream_ptr->openStream(std::move(open_stream_ptr))) {

    return std::nullopt;

  }

  return stream_ptr;

}


bool kel::StreamMTBuffer::openStream(std::unique_ptr<BaseStreamIO> stream_ptr) {

  // The supplied stream must be a valid, non-null object.
  if (not stream_ptr) {

    ExecEnv::log().error("StreamMTBuffer::openStream; null underlying stream object");
    return false;

  }

  try {

    close_stream_ = false;
    EOF_received_ = false;
    stream_ptr_ = std::move(stream_ptr);
    line_io_thread_.queueThreads(WORKER_THREAD_COUNT);
    line_io_thread_.enqueueVoid(&StreamMTBuffer::enqueueIOLineRecord, this);

  }
  catch (std::exception const &e) {

    ExecEnv::log().error("StreamMTBuffer::openStream; exception starting worker thread: {}", e.what());
    // The worker pool may already be running (queueThreads succeeded) if enqueueVoid() threw.
    // Signal shutdown and join before releasing the underlying stream so the worker thread is
    // never left executing enqueueIOLineRecord() against a destroyed object (leaked thread + UAF).
    close_stream_ = true;
    line_io_thread_.joinThreads();
    stream_ptr_ = nullptr;
    return false;

  }

  return true;

}


void kel::StreamMTBuffer::close() {

  // Signal the worker thread to stop reading from the underlying stream.
  // Note: close() still joins the worker thread, which blocks until the underlying stream
  // returns EOF or the current readLine() returns. It is intended to be called once the
  // stream has been consumed, or concurrently with a consumer that is draining the queue.
  close_stream_ = true;
  EOF_received_ = true;
  line_io_thread_.joinThreads();
  stream_ptr_ = nullptr;
  line_io_queue_.clear();

}

void kel::StreamMTBuffer::enqueueIOLineRecord() {

  while (not close_stream_.load()) {

    try {

      auto line_record = stream_ptr_->readLine();

      if (line_record.EOFRecord()) {

        line_io_queue_.push(std::move(line_record));
        break;

      }

      line_io_queue_.push(std::move(line_record));

    }
    catch (std::exception const &e) {

      ExecEnv::log().error("StreamMTBuffer::enqueueIOLineRecord; exception reading line: {}", e.what());
      break;

    }

  }

}

kel::IOLineRecord kel::StreamMTBuffer::readLine() {

  // Don't block on EOF, return additional EOF markers.
  // Note: if the worker thread is blocked on the underlying stream's readLine() and
  // close() has not been called, then waitAndPop() will block until a record is available.
  std::lock_guard<std::mutex> lock(mutex_);
  if (EOF_received_) return IOLineRecord::createEOFMarker();
  auto line_record = line_io_queue_.waitAndPop();
  if (line_record.EOFRecord()) {

    EOF_received_ = true;

  }
  return line_record;

}
