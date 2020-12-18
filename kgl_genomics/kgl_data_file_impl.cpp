//
// Created by kellerberrin on 18/4/20.
//

#include "kgl_data_file_impl.h"
#include <vector>
#include <memory>
#include <thread>

namespace kgl = kellerberrin::genome;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::FileDataIO::commenceIO(size_t reader_threads) {

  // The number of consumer threads reading the raw record queue; used for shutdown to enqueue eof markers.
  consumer_threads_ = reader_threads;

  file_stream_opt_ = BaseStreamIO::getReaderStream(fileName());

  if (not file_stream_opt_) {

    ExecEnv::log().critical("FileDataIO::rawDataIO; I/O error; could not open Data file: {}", fileName());
    enqueueEOF();
    return false;

  }

  for (size_t thread = 0; thread < io_threads_.threads().size(); ++thread) {

    io_threads_.enqueueWork(&FileDataIO::rawDataIO, this);

  }

  return true;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Pass in the text or gzip stream as an rvalue to a pointer.
void kgl::FileDataIO::rawDataIO() {

  try {

    while (not eof_flag_) {

      // Race condition; readLine() should not block on eof.
      IOLineRecord line_record = file_stream_opt_.value()->readLine();
      if (not line_record) {

        // Enqueue the null eof indicator for each consumer thread.
        enqueueEOF();
        // Terminate the read line loop.
        break;

      }

      // Check we have a valid line record.
      if (line_record.value().second->empty()) {

        ExecEnv::log().error("FileDataIO::rawDataIO; Empty Data record on line: {} - ignored", line_record.value().first);
        continue;

      }

      raw_io_queue_.push(std::move(line_record));

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("FileDataIO::rawDataIO; Data file: {} unexpected I/O exception: {}", fileName(), e.what());

  }

}

void kgl::FileDataIO::enqueueEOF() {

  std::lock_guard lock(eof_mutex_);

  if (eof_flag_) return;

  eof_flag_ = true;
  for (size_t i = 0; i < consumer_threads_; ++i) {

    raw_io_queue_.push(std::nullopt);

  }

}